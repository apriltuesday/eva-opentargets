from collections import defaultdict

import pandas as pd

from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT
from eva_cttv_pipeline.evidence_string_generation.consequence_type import get_so_accession_dict
from tests.eva_cttv_pipeline.evidence_string_generation import config


def test_process_gene():
    test_consequence_type_dict = defaultdict(list)
    test_rs_id = "rs121912888"
    test_ensembl_gene_id = "ENSG00000139219"
    test_so_name = "missense_variant"

    test_consequence_type = CT.ConsequenceType(test_ensembl_gene_id, CT.SoTerm(test_so_name))
    CT.process_gene(test_consequence_type_dict, test_rs_id, test_ensembl_gene_id, test_so_name)

    assert test_consequence_type_dict["rs121912888"][0] == test_consequence_type


def test_process_consequence_type_file_tsv():
    test_consequence_type = CT.ConsequenceType("ENSG00000139988", CT.SoTerm("synonymous_variant"))
    consequence_type_dict = CT.process_consequence_type_file(config.snp_2_gene_file)
    assert consequence_type_dict["14:67729241:C:T"][0] == test_consequence_type


def test_process_consequence_type_dataframes():
    dataframe_1 = pd.DataFrame(
        [('NC_000011.10:g.5226797_5226798insGCC', 'ENSG00000244734', 'HBB', 'coding_sequence_variant')],
        columns=('VariantID', 'EnsemblGeneID', 'EnsemblGeneName', 'ConsequenceTerm'))
    dataframe_2 = pd.DataFrame(
        [('RCV001051772', 'ENSG00000130711', 'PRDM12', 'trinucleotide_repeat_expansion')],
        columns=('1', '2', '3', '4'))  # column names can be anything
    consequence_type_dict = CT.process_consequence_type_dataframes(dataframe_1, dataframe_2)
    assert consequence_type_dict['NC_000011.10:g.5226797_5226798insGCC'][0].ensembl_gene_id == 'ENSG00000244734'
    assert consequence_type_dict['RCV001051772'][0].ensembl_gene_id == 'ENSG00000130711'


def test_ensembl_so_term():
    so_term = CT.SoTerm('stop_gained')
    assert so_term.accession == 'SO_0001587'
    assert so_term.rank == 4


def test_nonexistent_so_term():
    so_term = CT.SoTerm('not_real_term')
    assert so_term.accession is None
    assert so_term.rank == 40


def test_repeat_expansion_so_term():
    so_term = CT.SoTerm('short_tandem_repeat_expansion')
    assert so_term.accession == 'SO_0002162'
    assert so_term.rank == 40


def test_get_so_accession_dict():
    results = get_so_accession_dict(page_size=100)
    assert len(results) == 217
