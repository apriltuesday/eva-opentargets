import os
from collections import Counter

from cmat.output_generation.report import Report


resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


def get_test_report():
    report = Report()
    report.evidence_string_count = 42
    report.used_trait_mappings = {'one', 'two', 'three'}
    report.unmapped_trait_names = Counter('abracadabra')
    return report


def test_dump_and_load():
    counts_file = os.path.join(resources_dir, 'counts.yml')
    original_report = get_test_report()
    original_report.dump_to_file(resources_dir)

    loaded_report = Report()
    loaded_report.load_from_file(counts_file)
    assert original_report == loaded_report

    if os.path.exists(counts_file):
        os.remove(counts_file)


def test_add_operator():
    report_1 = get_test_report()

    report_2 = Report()
    report_2.evidence_string_count = 13
    report_2.total_trait_mappings = 100
    report_2.used_trait_mappings = {'three', 'four'}
    report_2.unmapped_trait_names = Counter('alakazam')

    combined = report_1 + report_2
    assert combined.clinvar_total == 0
    assert combined.evidence_string_count == 55
    assert combined.total_trait_mappings == 100
    assert combined.used_trait_mappings == {'one', 'two', 'three', 'four'}
    assert combined.unmapped_trait_names == Counter(
        {'a': 9, 'b': 2, 'r': 2, 'c': 1, 'd': 1, 'l': 1, 'k': 1, 'z': 1, 'm': 1}
    )


def test_check_counts():
    report = get_test_report()
    assert report.check_counts()
    report.clinvar_total = 5
    assert not report.check_counts()
