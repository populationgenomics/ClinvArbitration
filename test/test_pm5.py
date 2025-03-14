from pathlib import Path

from clinvarbitration.scripts.clinvar_by_codon import parse_tsv_into_dict

input_path = Path(__file__).parent / 'input'
test_tsv_file = input_path / 'post_annotation.tsv'


def test_pm5_table():
    """
    test that the PM5 table is generated correctly
    """

    result = parse_tsv_into_dict(str(test_tsv_file))
    assert result == {
        'ENST00000338591::561': {'904889::0'},
        'ENST00000341290::663': {'1310278::0'},
        'ENST00000694917::777': {'1310278::0'},
        'ENST00000433179::777': {'1310278::0'},
        'ENST00000649529::28': {'1341586::1'},
        'ENST00000624697::20': {'1341586::1'},
        'ENST00000624652::20': {'1341586::1'},
        'ENST00000379370::76': {'244110::1'},
        'ENST00000620552::310': {'822393::0'},
        'ENST00000651234::343': {'822393::0'},
    }
