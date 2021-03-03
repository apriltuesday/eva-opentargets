#!/usr/bin/env python3
"""Tests for the repeat expansion pipeline. Test resources are compressed XML files which contain one or a few records
manually extracted from the main ClinVar XML to check specific cases."""


import gzip
from io import BytesIO
import os
import pandas
import pandas.testing
import pytest
import tempfile

from repeat_expansion_variants import pipeline


def get_test_resource(resource_name):
    """Gets full path to the test resource located in the same directory as the test module."""

    # E. g. repeat_expansion_variants/test/repeat_expansion_pipeline_test.py::test_pipeline (call)
    pytest_current_test = os.getenv('PYTEST_CURRENT_TEST')

    # E.g. repeat_expansion_variants/test/repeat_expansion_pipeline_test.py
    test_module_path = pytest_current_test.split(':')[0]

    # E. g. repeat_expansion_variants/test/
    test_module_dir = os.path.split(test_module_path)[0]

    # E. g.  repeat_expansion_variants/test/input_variant_summary.tsv
    return os.path.join(test_module_dir, 'resources', resource_name)


def run_pipeline(resource_name):
    """Runs the pipeline on a given test resource and returns the output consequences as a list of lists."""
    input_filename = get_test_resource(resource_name)
    output_consequences, output_dataframe = [tempfile.NamedTemporaryFile(delete=False) for _ in range(2)]
    main.pipeline(input_filename, output_consequences.name, output_dataframe.name)
    consequences = [line.rstrip().split('\t') for line in open(output_consequences).read().splitlines()]
    output_consequences.cleanup()
    output_dataframe.cleanup()
    return consequences


def test_not_microsatellite():
    """Records which are not Microsatellite records should not result in any consequences."""
    assert run_pipeline('not_microsatellite.xml.gz') == []


def test_deletion():
    """Microsatellite deletion (contraction) events should not result in any consequences."""
    assert run_pipeline('deletion.xml.gz') == []


def test_short_insertion():
    """Short insertion events should not be treated as proper repeat expansion variants."""
    assert run_pipeline('short_insertion.xml.gz') == []


def test_explicit_coordinates():
    """Repeat expansion events with complete coordinates must be processed with the correct type."""
    assert run_pipeline('explicit_coords.xml.gz') == [
        # CT expansion — non-trinucleotide
        ['RCV000292700', '1', 'ENSG00000163554', 'SPTA1', 'short_tandem_repeat_expansion', '0'],
        # CGC expansion — trinucleotide
        ['RCV001051772', '1', 'ENSG00000130711', 'PRDM12', 'trinucleotide_repeat_expansion', '0'],
    ]


def test_no_explicit_coordinates():
    """Repeat expansion events without complete coordinates must also be processed and parsed."""
    assert run_pipeline('no_explicit_coords.xml.gz') == [
        # NM_023067.3(FOXL2):c.661GCN[15_24] (p.Ala221[(15_24)]) should be parsed as a trinucleotide expansion
        ['RCV000192035', '1', 'ENSG00000183770', 'FOXL2', 'trinucleotide_repeat_expansion', '0']
    ]
