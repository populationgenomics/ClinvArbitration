"""
Reads the main ClinVar table from the Zenodo release, and converts it to a TSV format.
This is used for exporting the ClinVar decisions data in a more accessible format, without duplicating the data.

The headered columns in the TSV will be:
- chrom: the chromosome (e.g. "chr1")
- pos: the position on the chromosome (e.g. 123456)
- ref: the reference allele (e.g. "A")
- alt: the alternate allele (e.g. "T")
- clinical_significance: the clinical significance of the variant (e.g. "Pathogenic/Likely Pathogenic", "VUS", "Benign")
- gold_stars: the number of gold stars (0-4)
- allele_id: the ClinVar allele ID (e.g. "123456")

The allele ID is usable directly in URLs:
    http://www.ncbi.nlm.nih.gov/clinvar?term=123456[alleleid]
"""

from argparse import ArgumentParser

import hail as hl


def reformat_table(ht: hl.Table) -> hl.Table:
    """
    Reformat the ClinVar decisions table to a simpler structure for export.

    Some fields are stored in the way Hail requires them (a locus is `chr:pos`, alleles are a list of [ref, alt])
    Which doesn't work well for TSV export.

    Here we split these compound fields into separate columns,
    and re-key the table on these fields to simplify the structure.
    """

    # this code splits up the locus field (chr:pos) and the alleles field ([ref, alt]) into separate columns for export
    ht = ht.annotate(
        chrom=ht.locus.contig,
        pos=ht.locus.position,
        ref=ht.alleles[0],
        alt=ht.alleles[1],
    )

    # re-key the table on these fields, so everything else can be dropped
    ht = ht.key_by('chrom', 'pos', 'ref', 'alt')

    # use select to drop all replaced fields
    ht = ht.select(
        ht.clinical_significance,
        ht.gold_stars,
        ht.allele_id,
    )

    return ht


if __name__ == '__main__':
    parser = ArgumentParser(description='Convert ClinVar decisions Hail Table to TSV format.')
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
