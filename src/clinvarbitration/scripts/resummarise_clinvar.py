"""
read clinvar submissions; identify consensus and disagreement

Requires two files from
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

submission_summary.txt
 - all individual submissions to ClinVar
relevant fields:
1.  VariationID: the identifier assigned by ClinVar
2.  ClinicalSignificance:
7.  ReviewStatus: the level of review for this submission
10. Submitter

variant_summary.txt
 - links clinvar AlleleID, Variant ID, position and alleles

These need to be localised prior to running this script.
"""

import gzip
import re
import zoneinfo
from argparse import ArgumentParser
from collections import defaultdict
from collections.abc import Generator
from dataclasses import dataclass
from datetime import datetime
from enum import Enum

import pysam
from loguru import logger


ASSEMBLY = 'Assembly'
GRCH38 = 'GRCh38'
BENIGN_SIGS = {'Benign', 'Likely benign', 'Benign/Likely benign', 'protective'}
CONFLICTING = 'conflicting data from submitters'
PATH_SIGS = {
    'Pathogenic',
    'Likely pathogenic',
    'Pathogenic, low penetrance',
    'Likely pathogenic, low penetrance',
    'Pathogenic/Likely pathogenic',
}
UNCERTAIN_SIGS = {'Uncertain significance', 'Uncertain risk allele'}

NO_STAR_RATINGS: set[str] = {'no assertion criteria provided'}
USELESS_RATINGS: set[str] = set()

MAJORITY_RATIO: float = 0.6
MINORITY_RATIO: float = 0.2
STRONG_REVIEWS: list[str] = ['practice guideline', 'reviewed by expert panel']
ORDERED_CONTIGS: list[str] = [f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY', 'chrM', 'chrMT']
TSV_KEYS = ['contig', 'position', 'reference', 'alternate', 'clinical_significance', 'gold_stars', 'allele_id']

# contig lengths, used to write a valid VCF header
CONTIG_LENGTHS: dict[str, int] = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chrX': 156040895,
    'chrY': 57227415,
    'chrM': 16569,
    'chrMT': 16569,
}

VCF_INFO_HEADERS = [
    '##INFO=<ID=allele_id,Number=1,Type=Integer,Description="ClinVar Allele ID">',
    '##INFO=<ID=gold_stars,Number=1,Type=Integer,Description="Review confidence, re-calculated by ClinvArbitration">',
    (
        '##INFO=<ID=clinical_significance,Number=1,Type=String,'
        'Description="Clinical significance, re-calculated by ClinvArbitration">'
    ),
]
VCF_COLUMNS = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'

# I really want the linter to just tolerate naive datetimes, but it won't
TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')

# published Nov 2015, available pre-print since March 2015
# assumed to be influential since 2016
ACMG_THRESHOLD = datetime(year=2016, month=1, day=1, tzinfo=TIMEZONE)

# a default date assigned to un-dated entries
VERY_OLD = datetime(year=1970, month=1, day=1, tzinfo=TIMEZONE)

LARGEST_COMPLEX_INDELS = 40
BASES = re.compile(r'[ACGTN]+')

# add the exact name of any submitters whose evidence is not trusted
BLACKLIST: set[str] = set()


class Consequence(Enum):
    """
    csq enumeration
    """

    BENIGN = 'Benign'
    CONFLICTING = 'Conflicting'
    PATHOGENIC = 'Pathogenic'
    UNCERTAIN = 'VUS'
    UNKNOWN = 'Unknown'


# an example of a qualified blacklist - entries of this type and site will be ignored
QUALIFIED_BLACKLIST = [(Consequence.BENIGN, ['illumina laboratory services; illumina'])]


@dataclass
class Submission:
    """
    POPO to store details on each Submission
    """

    date: datetime
    submitter: str
    classification: Consequence
    review_status: str


def get_allele_locus_map(summary_file: str) -> dict:
    """
    Process variant_summary.txt
     - links the allele ID, Locus/Alleles, and variant ID
    relevant fields:
    0 AlleleID
    20 Chromosome
    30 VariationID
    31 Start
    32 ReferenceAllele
    33 AlternateAllele

    Args:
        summary_file (str): path to the gzipped text file

    Returns:
        dictionary of each variant ID to the positional details
    """

    allele_dict = {}

    for line in dicts_from_gzip(summary_file):
        if line[ASSEMBLY] != GRCH38:
            continue

        chromosome = f'chr{line["Chromosome"]}'

        # normalise mitochondrial contig naming
        if chromosome == 'chrMT':
            chromosome = 'chrM'

        ref = line['ReferenceAlleleVCF']
        alt = line['AlternateAlleleVCF']

        # skip over cytogenetic locations
        if any(x == 'na' for x in [ref, alt]) or ref == alt:
            continue

        # skip non-standard chromosomes
        if chromosome not in ORDERED_CONTIGS:
            continue

        # skip chromosomal deletions and insertions, or massive indels
        if len(ref) + len(alt) > LARGEST_COMPLEX_INDELS:
            continue

        # pull values from the line
        allele_id = int(line['AlleleID'])
        var_id = int(line['VariationID'])
        uniq_var_id = f'{chromosome}_{var_id}'
        pos = int(line['PositionVCF'])

        # don't include any of the trash bases in ClinVar
        if BASES.match(ref) and BASES.match(alt):
            allele_dict[uniq_var_id] = {
                'var_id': var_id,
                'allele': allele_id,
                'chrom': chromosome,
                'pos': pos,
                'ref': ref,
                'alt': alt,
            }

    return allele_dict


def dicts_from_gzip(filename: str) -> Generator[dict[str, str], None, None]:
    """
    generator for gzip reading

    Args:
        filename (str): the gzipped input file

    Returns:
        generator; yields each line as a dictionary
    """

    # start with an empty list to please the linter
    header: list[str] = []

    with gzip.open(filename, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                header = line[1:].rstrip().split('\t')
                continue

            yield dict(zip(header, line.rstrip().split('\t'), strict=True))


def consequence_decision(subs: list[Submission]) -> Consequence:
    """
    determine overall consequence assignment based on submissions

    Args:
        subs (): a list of submission objects for this allele

    Returns:
        a single Consequence object
    """

    # start with a default consequence
    decision = Consequence.UNCERTAIN

    # establish counts for this allele
    counts = {
        Consequence.BENIGN: 0,
        Consequence.PATHOGENIC: 0,
        Consequence.UNCERTAIN: 0,
        Consequence.UNKNOWN: 0,
        'total': 0,
    }

    for each_sub in subs:
        # for 3/4-star ratings, don't look any further
        if each_sub.review_status in STRONG_REVIEWS:
            return each_sub.classification

        counts['total'] += 1
        if each_sub.classification in [
            Consequence.PATHOGENIC,
            Consequence.BENIGN,
            Consequence.UNCERTAIN,
            Consequence.UNKNOWN,
        ]:
            counts[each_sub.classification] += 1

    if counts[Consequence.PATHOGENIC] and counts[Consequence.BENIGN]:
        if (max(counts[Consequence.PATHOGENIC], counts[Consequence.BENIGN]) >= (counts['total'] * MAJORITY_RATIO)) and (
            min(counts[Consequence.PATHOGENIC], counts[Consequence.BENIGN]) <= (counts['total'] * MINORITY_RATIO)
        ):
            decision = (
                Consequence.BENIGN
                if counts[Consequence.BENIGN] > counts[Consequence.PATHOGENIC]
                else Consequence.PATHOGENIC
            )

        # both path and benign, but no clear majority - conflicting
        else:
            decision = Consequence.CONFLICTING

    # more than MAJORITY_RATIO are uncertain or unknown, call it that
    elif counts[Consequence.UNKNOWN] > (counts['total'] * MAJORITY_RATIO):
        decision = Consequence.UNKNOWN

    elif counts[Consequence.UNCERTAIN] > (counts['total'] * MAJORITY_RATIO):
        decision = Consequence.UNCERTAIN

    # any pathogenic - call it pathogenic
    elif counts[Consequence.PATHOGENIC]:
        decision = Consequence.PATHOGENIC

    # any benign - call it benign
    elif counts[Consequence.BENIGN]:
        decision = Consequence.BENIGN

    return decision


def check_stars(subs: list[Submission]) -> int:
    """
    processes the submissions, and assigns a 'gold star' rating
    this is a subset of the full ClinVar star system

    The NO_STAR_RATINGS set is ratings which we don't ascribe any
    star rating to, otherwise everything has a floor of 1, with
    an exit for 3 or 4 stars, for those superior review statuses

    Args:
        subs (): list of all submissions at this allele

    Returns:
        integer, summarising the rating
    """
    minimum = 0
    for sub in subs:
        if sub.classification in (Consequence.UNCERTAIN, Consequence.UNKNOWN):
            continue
        if sub.review_status == 'practice guideline':
            minimum = 4
        if sub.review_status == 'reviewed by expert panel':
            minimum = max(minimum, 3)
        if sub.review_status not in NO_STAR_RATINGS:
            minimum = max(minimum, 1)

    return minimum


def process_submission_line(data: dict[str, str]) -> tuple[int, Submission]:
    """
    takes a line, strips out useful content as a 'Submission'. Relevant fields:
    #VariationID
    ClinicalSignificance
    DateLastEvaluated
    ReviewStatus
    Submitter

    Args:
        data (): the array of line content

    Returns:
        the allele ID and corresponding Submission details
    """
    var_id = int(data['VariationID'])
    if data['ClinicalSignificance'] in PATH_SIGS:
        classification = Consequence.PATHOGENIC
    elif data['ClinicalSignificance'] in BENIGN_SIGS:
        classification = Consequence.BENIGN
    elif data['ClinicalSignificance'] in UNCERTAIN_SIGS:
        classification = Consequence.UNCERTAIN
    else:
        classification = Consequence.UNKNOWN
    date = (
        datetime.strptime(data['DateLastEvaluated'], '%b %d, %Y').replace(tzinfo=TIMEZONE)
        if data['DateLastEvaluated'] != '-'
        else VERY_OLD
    )
    sub = data['Submitter'].lower()
    rev_status = data['ReviewStatus'].lower()

    return var_id, Submission(date, sub, classification, rev_status)


def get_all_decisions(submission_file: str, var_ids: set[int]) -> dict[int, list[Submission]]:
    """
    obtains all submissions per-allele which pass basic criteria
        - not a blacklisted submitter
        - not a csq-specific blacklisted submitter

    Args:
        submission_file (): file containing submission-per-line
        var_ids (): only process Var IDs we have pos data for

    Returns:
        dictionary of var IDs and their corresponding submissions
    """

    submission_dict = defaultdict(list)

    for line in dicts_from_gzip(submission_file):
        var_id, line_sub = process_submission_line(line)

        # skip rows where the variantID isn't in this mapping
        # this saves a little effort on haplotypes, CNVs, and SVs
        if (
            (var_id not in var_ids)
            or (line_sub.submitter in BLACKLIST)
            or (line_sub.classification == Consequence.UNKNOWN)
        ):
            continue

        # screen out some submitters per-consequence
        for consequence, submitters in QUALIFIED_BLACKLIST:
            if line_sub.classification == consequence and line_sub.submitter in submitters:
                continue

        submission_dict[var_id].append(line_sub)

    return submission_dict


def acmg_filter_submissions(subs: list[Submission]) -> list[Submission]:
    """
    filter submissions by dates
    if any submissions for this variant occur after the ACMG introduction
        - only return those
    if not
        - return all submissions

    If the submission is an expert panel review or practice guideline, it is always retained.
    """

    # apply the date threshold to all submissions
    date_filt_subs = [sub for sub in subs if sub.date >= ACMG_THRESHOLD or sub.review_status in STRONG_REVIEWS]

    # if this contains results, return only those
    # default to returning everything
    return date_filt_subs or subs


def sort_decisions(all_subs: list[dict]) -> list[dict]:
    """Applies dual-layer sorting to the list of all decisions, on chr & pos."""

    return sorted(all_subs, key=lambda x: (ORDERED_CONTIGS.index(x['contig']), x['position']))


def generate_vcf_header(contigs: list[str]) -> str:
    """Handwrite a VCF header for the given contigs, mirroring the header Hail used to export."""

    header_lines = [
        '##fileformat=VCFv4.2',
        f'##fileDate={datetime.now(tz=TIMEZONE).strftime("%Y-%m-%d")}',
        '##source=ClinvArbitration',
        *VCF_INFO_HEADERS,
        *[
            f'##contig=<ID={contig},length={CONTIG_LENGTHS[contig]},assembly=GRCH38>'
            for contig in contigs
        ],
        VCF_COLUMNS,
    ]
    return '\n'.join(header_lines) + '\n'


def write_vcf(decisions: list[dict], output_vcf: str, pm5_filter: bool = True):
    """Takes the sorted clinvar decisions, optionally filter to SNV & Pathogenic. Writes results to an indexed VCF."""

    if pm5_filter:
        # filter to Pathogenic SNVs
        # there is at least one ClinVar submission which is Pathogenic without being a changed base?
        # https://www.ncbi.nlm.nih.gov/clinvar/variation/1705890/
        # new behaviour - we're not annotating chrM sites, as the default GTF file doesn't have Mito genes, so no csq
        decisions = [
            entry
            for entry in decisions
            if len(entry['alleles'][0]) == 1
            and len(entry['alleles'][1]) == 1
            and entry['clinical_significance'] == Consequence.PATHOGENIC.value
            and entry['contig'] != 'chrM'
        ]

    # only include contig header lines for contigs present in the data, in order of appearance
    contigs = list(dict.fromkeys(entry['contig'] for entry in decisions))

    # write header and rows block-gzipped, so the result can be tabix-indexed
    with pysam.BGZFile(output_vcf, 'wb') as handle:
        handle.write(generate_vcf_header(contigs).encode())
        for entry in decisions:
            ref, alt = entry['alleles']
            info = ';'.join(
                [
                    f'allele_id={entry["allele_id"]}',
                    f'gold_stars={entry["gold_stars"]}',
                    f'clinical_significance={entry["clinical_significance"]}',
                ],
            )
            handle.write(f'{entry["contig"]}\t{entry["position"]}\t.\t{ref}\t{alt}\t.\t.\t{info}\n'.encode())

    pysam.tabix_index(output_vcf, preset='vcf', force=True)
    logger.info(f'Wrote VCF to {output_vcf}')


def write_dicts_as_tsv(dicts: list[dict], output_path: str):
    """
    Writes a list of dictionaries to a TSV file, with headers.
    Args:
        dicts (list[dict]): List of dictionaries to write.
        output_path (str): Path to write the TSV file.
    """

    if not dicts:
        logger.warning('No data to write to TSV.')
        raise ValueError('No ClinVar decisions present.')

    logger.info(f'Writing {len(dicts)} entries to TSV at {output_path}')
    with open(output_path, 'w', encoding='utf-8') as tsv_file:
        # Write header
        tsv_file.write('\t'.join(TSV_KEYS) + '\n')

        # Write each dictionary as a row
        for each_dict in dicts:
            ref, alt = each_dict['alleles']
            each_dict['reference'] = ref
            each_dict['alternate'] = alt
            tsv_file.write('\t'.join(str(each_dict[key]) for key in TSV_KEYS) + '\n')

    logger.info(f'Wrote TSV to {output_path}')


def cli_main():
    parser = ArgumentParser(description='Generates a new clinVar summary from raw submission data')
    parser.add_argument(
        '-s',
        help='submission_summary.txt.gz from NCBI',
        required=True,
    )
    parser.add_argument(
        '-v',
        help='variant_summary.txt.gz from NCBI',
        required=True,
    )
    parser.add_argument(
        '-o',
        help='output root, for tsv and pathogenic-only VCF',
        required=True,
    )
    parser.add_argument(
        '-b',
        help='sites to blacklist',
        nargs='+',
        default=[],
    )
    parser.add_argument(
        '--all_vcf',
        help='if provided, write a VCF containing all entries',
        default=None,
    )

    args = parser.parse_args()

    # if sites are blacklisted on the CLI, update the global BLACKLIST value
    # temporary solution while we continue to validate Talos
    if args.b:
        BLACKLIST.update(args.b)

    main(subs=args.s, variants=args.v, output_root=args.o, all_vcf=args.all_vcf)


def main(subs: str, variants: str, output_root: str, all_vcf: str | None = None):
    """Parse all ClinVar submissions, and re-summarise with new algorithm."""
    logger.info('Getting alleleID-VariantID-Loci from variant summary')
    allele_map = get_allele_locus_map(variants)

    logger.info('Getting all decisions, indexed on clinvar Var ID')

    # the raw IDs - some have ambiguous X/Y mappings
    all_uniq_ids = {x['var_id'] for x in allele_map.values()}
    decision_dict = get_all_decisions(submission_file=subs, var_ids=all_uniq_ids)

    # placeholder to fill wth per-allele decisions
    all_decisions = {}

    # now filter each set of decisions per allele
    for var_id, submissions in decision_dict.items():
        # filter against ACMG date, if appropriate
        filtered_submissions = acmg_filter_submissions(submissions)

        # obtain an aggregate rating
        rating = Consequence.UNCERTAIN if not filtered_submissions else consequence_decision(filtered_submissions)

        # assess stars in remaining entries
        stars = check_stars(filtered_submissions)

        all_decisions[var_id] = (rating, stars)

    # now match those up with the variant coordinates
    logger.info('Matching decisions to variant coordinates')
    complete_decisions = []
    for var_details in allele_map.values():
        var_id = var_details['var_id']

        # we may have found no relevant submissions for this variant
        if var_id not in all_decisions:
            continue

        # add the decision to the list, inc. variant details
        complete_decisions.append(
            {
                'contig': var_details['chrom'],
                'position': var_details['pos'],
                'alleles': [var_details['ref'], var_details['alt']],
                'clinical_significance': all_decisions[var_id][0].value,
                'gold_stars': all_decisions[var_details['var_id']][1],
                'allele_id': var_details['allele'],
            },
        )

    logger.info(f'{len(complete_decisions)} ClinVar entries remain')

    # sort all collected decisions by contig & position, required for tabix indexing
    complete_decisions_sorted = sort_decisions(complete_decisions)

    tsv_path = f'{output_root}.tsv'
    write_dicts_as_tsv(complete_decisions_sorted, output_path=tsv_path)

    # write a VCF containing all variants, not just pathogenic SNV (Echtvar use case)
    if all_vcf:
        write_vcf(complete_decisions_sorted, all_vcf, pm5_filter=False)

    # export the pathogenic SNVs as a tabix-indexed VCF
    vcf_output = f'{output_root}.vcf.bgz'
    logger.info(f'Writing out Pathogenic SNV VCF to {vcf_output}')
    write_vcf(complete_decisions_sorted, vcf_output)


if __name__ == '__main__':
    cli_main()
