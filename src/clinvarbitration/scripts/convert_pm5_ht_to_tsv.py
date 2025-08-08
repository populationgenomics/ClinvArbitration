"""
Reads the secondary PM5 ClinVar table from the Zenodo release, and converts it to a TSV format.
All ClinVar entries referenced in this file are Pathogenic missense variants.

The headered columns in the TSV will be:
- transcript: the ENSTnnnnn transcript ID (e.g. "ENST00000367770")
- codon: int, the codon position in the transcript (e.g. 334)
- clinvar_alleles: String

There can be multiple records for the same codon, so the encoding is:

- "ALLELE_IC::GOLD_STARS", e.g. "1341586::1"
- there can be multiple records, so these are separated by "+", e.g. "1341586::1+1234567::2"
- the number of entries is arbitrary, and they should all be linked to a single transcript/codon, so in the TSV export
    these are being retained as a single string, and a user can interact with that as they need

The allele IDs in this compound String are usable directly in URLs:
    http://www.ncbi.nlm.nih.gov/clinvar?term=123456[alleleid]
"""

from argparse import ArgumentParser

import hail as hl


def reformat_table(ht: hl.Table) -> hl.Table:
    """
    Reformat the PM5 ClinVar decisions table to a simpler structure for export.
    """

    # this code splits up the locus field (chr:pos) and the alleles field ([ref, alt]) into separate columns for export
    transcript_codon = ht.newkey.split('::')

    ht = ht.annotate(
        transcript=transcript_codon[0],
        codon=transcript_codon[1],
    )

    # re-key the table on these fields, so everything else can be dropped
    ht = ht.key_by('transcript', 'codon')

    # use select to drop all replaced fields
    return ht.select(ht.clinvar_alleles)


if __name__ == '__main__':
    parser = ArgumentParser(description='Convert ClinVar PM5 Hail Table to TSV format.')
    parser.add_argument('--input', type=str, required=True, help='Input Hail Table path (clinvar_decisions.ht).')
    parser.add_argument('--output', type=str, required=True, help='Output TSV file path.')
    args = parser.parse_args()

    print(f'Converting ClinVar decisions from {args.input} to {args.output}')

    hl.init(quiet=True)
    hl.default_reference('GRCh38')

    ht = hl.read_table(args.input)

    print(f'Loaded ClinVar decisions table with {ht.count()} rows.')

    # reformat the table for export as a TSV
    ht = reformat_table(ht)

    # export the data as a TSV
    ht.export(
        output=args.output,
        header=True,
    )

    print(f'Exported ClinVar decisions to {args.output}')
