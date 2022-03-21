from collections import defaultdict

from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT
from tests.eva_cttv_pipeline.evidence_string_generation import config


class TestProcessGene:
    def test__process_gene(self):
        test_consequence_type_dict = defaultdict(list)
        test_rs_id = "rs121912888"
        test_ensembl_gene_id = "ENSG00000139219"
        test_so_name = "missense_variant"

        test_consequence_type = CT.ConsequenceType(test_ensembl_gene_id, CT.SoTerm(test_so_name))

        CT.process_gene(test_consequence_type_dict, test_rs_id, test_ensembl_gene_id, test_so_name)

        assert test_consequence_type_dict["rs121912888"][0] == test_consequence_type


class TestProcessConsequenceTypeFileTsv:
    def test__process_consequence_type_file_tsv(self):
        test_consequence_type = CT.ConsequenceType("ENSG00000139988", CT.SoTerm("synonymous_variant"))
        consequence_type_dict = CT.process_consequence_type_file(config.snp_2_gene_file)
        assert consequence_type_dict["14:67729241:C:T"][0] == test_consequence_type


class TestSoTerm:
    @classmethod
    def setup_class(cls):
        cls.test_so_term_a = CT.SoTerm("stop_gained")
        cls.test_so_term_b = CT.SoTerm("not_real_term")

    def test_accession(self):
        assert self.test_so_term_a.accession == "SO_0001587"
        assert self.test_so_term_b.accession is None

    def test_rank(self):
        assert self.test_so_term_a.rank == 6
        assert self.test_so_term_b.rank == 37
