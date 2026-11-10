import zoneinfo
from copy import deepcopy
from datetime import datetime
from gzip import open as gzip_open
from pathlib import Path

import pytest

from clinvarbitration.scripts.resummarise_clinvar import (
    Consequence,
    Submission,
    check_stars,
    consequence_decision,
    write_vcf,
)

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
BASIC_SUB = Submission(datetime.now(tz=TIMEZONE), 'foo', Consequence.UNKNOWN, 'review')


def test_check_stars_none():
    """
    check that we always get the right number of stars!
    """
    # Unknown, should be skipped over
    dud_1 = deepcopy(BASIC_SUB)
    dud_2 = deepcopy(BASIC_SUB)
    # Pathogenic, but has a no-star rating
    dud_2.classification = Consequence.PATHOGENIC
    dud_2.review_status = 'no assertion criteria provided'
    assert check_stars([dud_1, dud_2]) == 0


def test_check_stars_1():
    """
    check that we always get the right number of stars!
    """
    expect_1_path = deepcopy(BASIC_SUB)
    expect_1_path.classification = Consequence.PATHOGENIC
    assert check_stars([expect_1_path]) == 1
    expect_1_benign = deepcopy(BASIC_SUB)
    expect_1_benign.classification = Consequence.BENIGN
    assert check_stars([expect_1_benign]) == 1
    assert check_stars([expect_1_path, expect_1_benign]) == 1


def test_check_stars_3():
    """
    check that we always get the right number of stars!
    """
    expect_1_path = deepcopy(BASIC_SUB)
    expect_1_path.classification = Consequence.PATHOGENIC
    expect_1_neutral = deepcopy(BASIC_SUB)
    assert check_stars([expect_1_path, expect_1_neutral]) == 1
    expert_panel = deepcopy(BASIC_SUB)
    expert_panel.review_status = 'reviewed by expert panel'
    expert_panel.classification = Consequence.PATHOGENIC
    assert check_stars([expect_1_path, expect_1_neutral, expert_panel]) == 3  # noqa: PLR2004


def test_check_stars_4():
    """
    check that we always get the right number of stars!
    """
    expect_1_path = deepcopy(BASIC_SUB)
    expect_1_path.classification = Consequence.PATHOGENIC
    expect_1_neutral = deepcopy(BASIC_SUB)
    assert check_stars([expect_1_path, expect_1_neutral]) == 1
    expert_panel = deepcopy(BASIC_SUB)
    expert_panel.review_status = 'reviewed by expert panel'
    expert_panel.classification = Consequence.PATHOGENIC
    assert check_stars([expect_1_path, expect_1_neutral, expert_panel]) == 3  # noqa: PLR2004
    practice_guideline = deepcopy(BASIC_SUB)
    practice_guideline.review_status = 'practice guideline'
    practice_guideline.classification = Consequence.PATHOGENIC
    assert check_stars([expect_1_path, expect_1_neutral, expert_panel, practice_guideline]) == 4  # noqa: PLR2004


@pytest.mark.parametrize(
    ('consequences', 'expected'),
    [
        ([Consequence.PATHOGENIC], Consequence.PATHOGENIC),
        ([Consequence.PATHOGENIC, Consequence.PATHOGENIC], Consequence.PATHOGENIC),
        ([Consequence.PATHOGENIC, Consequence.UNCERTAIN], Consequence.PATHOGENIC),
        ([Consequence.BENIGN], Consequence.BENIGN),
        ([Consequence.BENIGN, Consequence.PATHOGENIC], Consequence.CONFLICTING),
    ],
)
def test_consequence_decision_path_single(
    consequences: list[Consequence],
    expected: Consequence,
):
    all_subs = []
    for con in consequences:
        sub = deepcopy(BASIC_SUB)
        sub.classification = con
        all_subs.append(sub)
    assert consequence_decision(all_subs) == expected


def test_write_vcf(tmp_path: Path):
    """
    write a small set of decisions, check the handwritten VCF content and index
    """
    decisions = [
        {
            'contig': 'chr1',
            'position': 100,
            'alleles': ['A', 'T'],
            'clinical_significance': Consequence.PATHOGENIC.value,
            'gold_stars': 1,
            'allele_id': 111,
        },
        # non-SNV, dropped by the pm5 filter
        {
            'contig': 'chr1',
            'position': 200,
            'alleles': ['AC', 'T'],
            'clinical_significance': Consequence.PATHOGENIC.value,
            'gold_stars': 2,
            'allele_id': 222,
        },
        # benign, dropped by the pm5 filter
        {
            'contig': 'chr2',
            'position': 300,
            'alleles': ['G', 'C'],
            'clinical_significance': Consequence.BENIGN.value,
            'gold_stars': 3,
            'allele_id': 333,
        },
        # chrM, dropped by the pm5 filter
        {
            'contig': 'chrM',
            'position': 400,
            'alleles': ['G', 'C'],
            'clinical_significance': Consequence.PATHOGENIC.value,
            'gold_stars': 0,
            'allele_id': 444,
        },
    ]

    vcf_path = str(tmp_path / 'decisions.vcf.bgz')
    write_vcf(decisions, vcf_path)

    with gzip_open(vcf_path, 'rt') as handle:
        lines = [line.rstrip() for line in handle]

    header = [line for line in lines if line.startswith('##')]
    records = [line for line in lines if not line.startswith('#')]

    assert '##fileformat=VCFv4.2' in header
    assert '##contig=<ID=chr1,length=248956422,assembly=GRCh38>' in header
    # filtered contigs don't appear in the header
    assert not any(line.startswith('##contig=<ID=chrM') for line in header)

    assert records == [
        'chr1\t100\t.\tA\tT\t.\t.\tallele_id=111;gold_stars=1;clinical_significance=Pathogenic',
    ]

    # a tabix index was created alongside
    assert (tmp_path / 'decisions.vcf.bgz.tbi').exists()


def test_write_vcf_unfiltered(tmp_path: Path):
    """
    without the pm5 filter, everything is retained
    """
    decisions = [
        {
            'contig': 'chr1',
            'position': 200,
            'alleles': ['AC', 'T'],
            'clinical_significance': Consequence.BENIGN.value,
            'gold_stars': 2,
            'allele_id': 222,
        },
        {
            'contig': 'chrM',
            'position': 400,
            'alleles': ['G', 'C'],
            'clinical_significance': Consequence.PATHOGENIC.value,
            'gold_stars': 0,
            'allele_id': 444,
        },
    ]

    vcf_path = str(tmp_path / 'all.vcf.bgz')
    write_vcf(decisions, vcf_path, pm5_filter=False)

    with gzip_open(vcf_path, 'rt') as handle:
        records = [line.rstrip() for line in handle if not line.startswith('#')]

    assert len(records) == 2  # noqa: PLR2004
    assert records[0].startswith('chr1\t200\t.\tAC\tT')
    assert records[1].startswith('chrM\t400\t.\tG\tC')
