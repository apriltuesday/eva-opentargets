from cmat.output_generation.evaluation.ols_utils import fetch_eval_data


def test_fetch_eval_data():
    expected = ('MONDO:0004975', False, {'MONDO:0004975'})
    assert fetch_eval_data(db_iden=('MONDO', 'MONDO:0004975')) == expected
    assert fetch_eval_data(uri='http://purl.obolibrary.org/obo/MONDO_0004975') == expected


def test_fetch_eval_data_include_neighbors():
    expected = ('MONDO:0004975', False, {'MONDO:0004975'},
                {'EFO:0005815', 'MONDO:0001627'}, {'MONDO:0100087', 'MONDO:0014265', 'EFO:1001870'})
    assert fetch_eval_data(db_iden=('MONDO', 'MONDO:0004975'), include_neighbors=True) == expected
