#!/usr/bin/env python3
"""Tests for the repeat expansion pipeline. Test resources are compressed XML files which contain one or a few records
manually extracted from the main ClinVar XML to check specific cases."""

import os
import tempfile

from consequence_prediction.repeat_expansion_variants import pipeline


def get_test_resource(resource_name):
    """Gets full path to the test resource located in the same directory as the test module."""

    # Full path to this module.
    this_module = os.path.abspath(__file__)

    # Full path to the directory where it is contained.
    module_directory = os.path.dirname(this_module)

    # Full path to the requested resource.
    return os.path.join(module_directory, 'resources', resource_name)


def run_pipeline(resource_name):
    """Runs the pipeline on a given test resource and returns the output consequences as a list of lists."""
    input_filename = get_test_resource(resource_name)
    output_consequences, output_dataframe = [tempfile.NamedTemporaryFile(delete=False) for _ in range(2)]
    pipeline.main(input_filename, output_consequences.name, output_dataframe.name)
    consequences = [line.rstrip().split('\t') for line in open(output_consequences.name).read().splitlines()]
    for temp_file in (output_consequences, output_dataframe):
        os.remove(temp_file.name)
    return consequences


def test_not_microsatellite():
    """Records which are not Microsatellite records should not result in any consequences."""
    assert run_pipeline('not_microsatellite.xml.gz') == []


def test_deletion():
    """Microsatellite deletion (contraction) events should not result in any consequences."""
    # Contains two deletions:
    # 1. With explicit coordinates: RCV000000275, delCTC
    # 2. Without explicit coordinates: RCV000481576, NM_000044.4(AR):c.172_174CAG(7_34) (p.Gln66_Gln80del)
    assert run_pipeline('deletion.xml.gz') == []


def test_short_insertion():
    """Short insertion events should not be treated as proper repeat expansion variants."""
    assert run_pipeline('short_insertion.xml.gz') == []


def test_explicit_coordinates():
    """Repeat expansion events with complete coordinates must be processed with the correct type."""
    assert sorted(run_pipeline('explicit_coords.xml.gz')) == sorted([
        # CGC expansion, trinucleotide.
        ['RCV001051772', 'ENSG00000130711', 'PRDM12', 'trinucleotide_repeat_expansion'],
        # CCGGGACCGAGA (12 base unit) expansion, also to be considered a trinucleotide expansion.
        ['RCV000722291', 'ENSG00000142599', 'RERE', 'trinucleotide_repeat_expansion'],
        # CT expansion, non-trinucleotide.
        ['RCV000292700', 'ENSG00000163554', 'SPTA1', 'short_tandem_repeat_expansion'],
        # TCAT expansion, non-trinucleotide but could be mistaken for one (3 units of 4 = 12, divisible by 3).
        ['RCV000122358', 'ENSG00000135100', 'HNF1A', 'short_tandem_repeat_expansion'],
    ])


def test_no_explicit_coordinates():
    """Repeat expansion events without complete coordinates must also be processed and parsed."""
    assert run_pipeline('no_explicit_coords.xml.gz') == [
        # NM_023067.3(FOXL2):c.661GCN[15_24] (p.Ala221[(15_24)]) should be parsed as a trinucleotide expansion
        ['RCV000192035', 'ENSG00000183770', 'FOXL2', 'trinucleotide_repeat_expansion'],
    ]


def test_alternative_identifiers():
    """Protein HGVS and human-readable identifiers should also be processed."""
    assert sorted(run_pipeline('alternatives.xml.gz')) == sorted([
        # NP_003915.2:p.Ala260(5_9) is protein hgvs and assumed to be trinucleotide expansion
        ['RCV000006377', 'ENSG00000109132', 'PHOX2B', 'trinucleotide_repeat_expansion'],
        # ATXN2, (CAG)n REPEAT EXPANSION is inferred to be trinucleotide from the repeat unit length
        ['RCV000008583', 'ENSG00000204842', 'ATXN2', 'trinucleotide_repeat_expansion']
    ])


def test_missing_names():
    """Records with missing or unparseable variant names are still processed using HGVS identifier."""
    assert sorted(run_pipeline('missing_names.xml.gz')) == sorted([
        ['RCV000006519', 'ENSG00000230223', 'ATXN8OS', 'trinucleotide_repeat_expansion'],
        ['RCV000087738', 'ENSG00000141543', 'EIF4A3', 'trinucleotide_repeat_expansion'],
        ['RCV001356600', 'ENSG00000136869', 'TLR4', 'short_tandem_repeat_expansion'],
        ['RCV001357106', 'ENSG00000214253', 'FIS1', 'short_tandem_repeat_expansion'],
        ['RCV001358198', 'ENSG00000171552', 'BCL2L1', 'short_tandem_repeat_expansion']
    ])


def test_missing_names_and_hgvs():
    """Records that are missing variant names and HGVS should use coordinate span from alleles instead."""
    assert sorted(run_pipeline('missing_names_and_hgvs.xml.gz')) == [
        # ref=T, alt=TGAAAGAAAGAAAGAAAGAAA => correctly classified as short tandem repeat
        ['RCV001355211', 'ENSG00000109861', 'CTSC', 'short_tandem_repeat_expansion'],
        # ref=T, alt=TACACACACACAC => classified as trinucleotide repeat without repeating unit inference.
        ['RCV001356600', 'ENSG00000136869', 'TLR4', 'trinucleotide_repeat_expansion']
    ]
