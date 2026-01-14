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
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from functools import partial
from itertools import islice

import hail as hl
import pandas as pd
from loguru import logger

ASSEMBLY = 'Assembly'
GRCH37 = 'GRCh37'
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
ORDERED_CONTIGS: dict[str, list[str]] = {
    GRCH38: [f'chr{x}' for x in list(range(1, 23))]
    + ['chrX', 'chrY', 'chrM', 'chrMT'],
    GRCH37: [*list(map(str, range(1, 23))), 'X', 'Y', 'M', 'MT'],
}
TSV_KEYS = [
    'contig',
    'position',
    'reference',
    'alternate',
    'clinical_significance',
    'gold_stars',
    'allele_id',
]

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
    PATHOGENIC = 'Pathogenic/Likely Pathogenic'
    UNCERTAIN = 'VUS'
    UNKNOWN = 'Unknown'


# an example of a qualified blacklist - entries of this type and site will be ignored
QUALIFIED_BLACKLIST = [
    (Consequence.BENIGN, ['illumina laboratory services; illumina'])
]


@dataclass
class Submission:
    """
    POPO to store details on each Submission
    """

    date: datetime
    submitter: str
    classification: Consequence
    review_status: str


def chunk_file(filename: str, chunk_size: int = 10000) -> Generator[list[str], None, None]:
    """
    Generator that yields chunks of lines from a file.
    Args:
        filename: path to file
        chunk_size: number of lines per chunk
    """
    with gzip.open(filename, 'rt') as f:
        # skip header
        header = f.readline()
        chunk = []
        for line in f:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk


def get_header(filename: str) -> list[str]:
    """
    Reads just the header from a gzipped file.
    """
    with gzip.open(filename, 'rt') as handle:
        for line in handle:
            # the submission file has multiple header lines, including two that start with "#VariationID".
            # skip short #-prefixed lines to get to the true header
            if len(line.split('\t')) < 3:
                continue
            if line.startswith('#'):
                return line[1:].rstrip().split('\t')
            return line.rstrip().split('\t')
    return []


def process_variant_chunk(chunk: list[str], header: list[str], assembly: str) -> dict:
    """
    Worker function to process a chunk of variant_summary.txt lines
    """
    allele_dict = {}
    for line_str in chunk:
        line = dict(zip(header, line_str.rstrip().split('\t'), strict=True))
        
        if line[ASSEMBLY] != assembly:
            continue

        chromosome = (
            f'chr{line["Chromosome"]}' if assembly == GRCH38 else line['Chromosome']
        )

        # swap chrM to something Hail will tolerate
        if chromosome == 'chrMT':
            chromosome = 'chrM'

        ref = line['ReferenceAlleleVCF']
        alt = line['AlternateAlleleVCF']

        # skip over cytogenetic locations
        if any(x == 'na' for x in [ref, alt]) or ref == alt:
            continue

        # skip non-standard chromosomes
        if chromosome not in ORDERED_CONTIGS[assembly]:
            continue

        # skip chromosomal deletions and insertions, or massive indels
        if len(ref) + len(alt) > LARGEST_COMPLEX_INDELS:
            continue

        # pull values from the line
        allele_id = int(line['AlleleID'])
        var_id = line['VariationID']
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


def get_allele_locus_map(summary_file: str, assembly: str) -> dict:
    """
    Process variant_summary.txt in parallel
    Args:
        summary_file (str): path to the gzipped text file
        assembly (str): genome build to use

    Returns:
        dictionary of each variant ID to the positional details
    """
    
    header = get_header(summary_file)

    allele_dict = {}
    
    # Use ProcessPoolExecutor to parallelize
    with ProcessPoolExecutor() as executor:
        # Create a partial function with fixed arguments
        worker = partial(process_variant_chunk, header=header, assembly=assembly)
        
        # Submit chunks to the executor
        # Using a larger chunk size for efficiency
        results = executor.map(worker, chunk_file(summary_file, chunk_size=50000))
        
        # Merge results
        for result in results:
            allele_dict.update(result)
            
    return allele_dict


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
        if (
            max(counts[Consequence.PATHOGENIC], counts[Consequence.BENIGN])
            >= (counts['total'] * MAJORITY_RATIO)
        ) and (
            min(counts[Consequence.PATHOGENIC], counts[Consequence.BENIGN])
            <= (counts['total'] * MINORITY_RATIO)
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


def process_submission_line(data: dict[str, str]) -> tuple[str, Submission]:
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
    var_id = data['VariationID']
    if data['ClinicalSignificance'] in PATH_SIGS:
        classification = Consequence.PATHOGENIC
    elif data['ClinicalSignificance'] in BENIGN_SIGS:
        classification = Consequence.BENIGN
    elif data['ClinicalSignificance'] in UNCERTAIN_SIGS:
        classification = Consequence.UNCERTAIN
    else:
        classification = Consequence.UNKNOWN
    date = (
        datetime.strptime(data['DateLastEvaluated'], '%b %d, %Y').replace(
            tzinfo=TIMEZONE
        )
        if data['DateLastEvaluated'] != '-'
        else VERY_OLD
    )
    sub = data['Submitter'].lower()
    rev_status = data['ReviewStatus'].lower()

    return var_id, Submission(date, sub, classification, rev_status)


def process_submission_chunk(
    chunk: list[str], header: list[str], var_ids: set[int]
) -> dict[str, list[Submission]]:
    """
    Worker function to process a chunk of submission_summary.txt lines
    """
    submission_dict = defaultdict(list)
    for line_str in chunk:
        if line_str.startswith('#'):
            continue
        line = dict(zip(header, line_str.rstrip().split('\t'), strict=True))
        var_id, line_sub = process_submission_line(line)

        # skip rows where the variantID isn't in this mapping
        if (
                (var_id not in var_ids)
                or (line_sub.submitter in BLACKLIST)
                or (line_sub.classification == Consequence.UNKNOWN)
        ):
            continue

        # screen out some submitters per-consequence
        skip = False
        for consequence, submitters in QUALIFIED_BLACKLIST:
            if (
                    line_sub.classification == consequence
                    and line_sub.submitter in submitters
            ):
                skip = True
                break
        if skip:
            continue

        submission_dict[var_id].append(line_sub)
    return submission_dict


def get_all_decisions(
    submission_file: str, var_ids: set[int]
) -> dict[int, list[Submission]]:
    """
    obtains all submissions per-allele which pass basic criteria in parallel
        - not a blacklisted submitter
        - not a csq-specific blacklisted submitter

    Args:
        submission_file (): file containing submission-per-line
        var_ids (): only process Var IDs we have pos data for

    Returns:
        dictionary of var IDs and their corresponding submissions
    """
    header = get_header(submission_file)
    submission_dict = defaultdict(list)

    with ProcessPoolExecutor() as executor:
        worker = partial(process_submission_chunk, header=header, var_ids=var_ids)
        
        # larger chunk size for submission entries
        results = executor.map(worker, chunk_file(submission_file, chunk_size=50000))
        
        for result in results:
            for var_id, subs in result.items():
                submission_dict[var_id].extend(subs)

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
    date_filt_subs = [
        sub
        for sub
        in subs
        if sub.date >= ACMG_THRESHOLD or sub.review_status in STRONG_REVIEWS
    ]

    # if this contains results, return only those
    # default to returning everything
    return date_filt_subs or subs


def process_decision_batch(
    batch: list[tuple[int, list[Submission]]]
) -> list[tuple[int, Consequence, int]]:
    """
    Worker function to process a batch of decisions
    """
    results = []
    for var_id, submissions in batch:
        # filter against ACMG date, if appropriate
        filtered_submissions = acmg_filter_submissions(submissions)

        # obtain an aggregate rating
        rating = (
            Consequence.UNCERTAIN
            if not filtered_submissions
            else consequence_decision(filtered_submissions)
        )

        # assess stars in remaining entries
        stars = check_stars(filtered_submissions)
        results.append((var_id, rating, stars))
    return results


def sort_decisions(all_subs: list[dict], assembly: str) -> list[dict]:
    """Applies dual-layer sorting to the list of all decisions, on chr & pos."""

    return sorted(
        all_subs,
        key=lambda x: (ORDERED_CONTIGS[assembly].index(x['contig']), x['position']),
    )


def parse_into_table(tsv_path: str, out_path: str) -> hl.Table:
    """Takes the file of one clinvar variant per line, processes that line into a table."""
    logger.info(f'Importing table from {tsv_path}')
    ht = hl.import_table(
        tsv_path,
        types={'position': hl.tint32, 'gold_stars': hl.tint32, 'allele_id': hl.tint32},
    )

    # create a locus value, and key the table by this. Combine [ref, alt] alleles into a list
    ht = ht.transmute(
        locus=hl.locus(ht.contig, ht.position),
        alleles=[ht.reference, ht.alternate],
    )

    ht = ht.key_by(ht.locus, ht.alleles)

    ht = ht.annotate_globals(
        creation_date=datetime.now(tz=TIMEZONE).strftime('%Y-%m-%d'),
        blacklist=sorted(BLACKLIST) or ['no blacklisted sites'],
    )

    # write out to the specified location
    ht.write(out_path, overwrite=True)

    logger.info(f'Wrote ClinVar decisions HailTable to {out_path}')

    # read the localised version
    return hl.read_table(out_path)


def write_pm5_vcf(clinvar_table: hl.Table, output_vcf: str):
    """Takes a clinvar decisions HailTable, filters to SNV & Pathogenic. Writes results to a VCF file."""

    # filter to Pathogenic SNVs
    # there is at least one ClinVar submission which is Pathogenic without being a changed base?
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/1705890/
    # new behaviour - we're not annotating chrM sites, as the default GTF file doesn't have Mito genes, so no csq
    clinvar_table = clinvar_table.filter(
        (hl.len(clinvar_table.alleles[0]) == 1)
        & (hl.len(clinvar_table.alleles[1]) == 1)
        & (clinvar_table.clinical_significance == Consequence.PATHOGENIC.value)
        & (clinvar_table.locus.contig != 'chrM'),
    )

    # persist the relevant clinvar annotations in INFO (for vcf export)
    clinvar_table = clinvar_table.transmute(
        info=hl.struct(
            allele_id=clinvar_table.allele_id,
            gold_stars=clinvar_table.gold_stars,
            clinical_significance=clinvar_table.clinical_significance,
        ),
    )

    # export this data in VCF format
    hl.export_vcf(clinvar_table, output_vcf, tabix=True)
    logger.info(f'Wrote SNV VCF to {output_vcf}')


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
    parser = ArgumentParser(
        description='Generates a new clinVar summary from raw submission data'
    )
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
        help='output root, for table, tsv, and pathogenic-only VCF',
        required=True,
    )
    parser.add_argument(
        '-b',
        help='sites to blacklist',
        nargs='+',
        default=[],
    )
    parser.add_argument(
        '--assembly',
        help='genome build to use',
        default='GRCh38',
        choices=[GRCH37, GRCH38],
    )

    args = parser.parse_args()

    # if sites are blacklisted on the CLI, update the global BLACKLIST value
    # temporary solution while we continue to validate Talos
    if args.b:
        BLACKLIST.update(args.b)

    main(subs=args.s, variants=args.v, output_root=args.o, assembly=args.assembly)


def batch_iterable(iterable, batch_size=1000):
    """Yields batches of items from an iterable."""
    iterator = iter(iterable)
    while batch := list(islice(iterator, batch_size)):
        yield batch


def main(subs: str, variants: str, output_root: str, assembly: str):
    """Parse all ClinVar submissions, and re-summarise with new algorithm."""

    logger.info('Getting alleleID-VariantID-Loci from variant summary')
    allele_map = get_allele_locus_map(variants, assembly)

    logger.info('Getting all decisions, indexed on clinvar Var ID')

    # the raw IDs - some have ambiguous X/Y mappings
    all_uniq_ids = {x['var_id'] for x in allele_map.values()}
    decision_dict = get_all_decisions(submission_file=subs, var_ids=all_uniq_ids)

    # placeholder to fill wth per-allele decisions
    all_decisions = {}

    # now filter each set of decisions per allele
    logger.info('Processing decisions in parallel')
    
    # Batch the decision processing
    decision_items = list(decision_dict.items())
    
    with ProcessPoolExecutor() as executor:
        # submit batches to workers
        # Using a significantly large batch size to reduce overhead
        results = executor.map(process_decision_batch, batch_iterable(decision_items, batch_size=1000))
        
        for batch_result in results:
            for var_id, rating, stars in batch_result:
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

    # sort all collected decisions, trying to reduce overhead in HT later
    complete_decisions_sorted = sort_decisions(complete_decisions, assembly=assembly)

    tsv_path = f'{output_root}.tsv'
    write_dicts_as_tsv(complete_decisions_sorted, output_path=tsv_path)

    # Initialize Spark here, just before it makes sense to use it
    hl.context.init_spark(master='local[*]')
    hl.default_reference(assembly)
    
    ht_output = f'{output_root}.ht'
    ht = parse_into_table(tsv_path=tsv_path, out_path=ht_output)

    # export the pathogenic SNVs as a tabix-indexed VCF
    vcf_output = f'{output_root}.vcf.bgz'
    logger.info(f'Writing out Pathogenic SNV VCF to {vcf_output}')
    write_pm5_vcf(ht, vcf_output)


if __name__ == '__main__':
    cli_main()
