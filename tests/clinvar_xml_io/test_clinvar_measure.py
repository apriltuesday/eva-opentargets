import os

from cmat.clinvar_xml_io import ClinVarDataset

resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


class TestClinvarRecordMeasure:
    @classmethod
    def setup_class(cls):
        input_file = os.path.join(resources_dir, 'clinvar_dataset_v2.xml.gz')
        cls.test_crm = next(iter(ClinVarDataset(input_file))).measure

    def test_hgvs(self):
        text_hgvs = [h.text for h in self.test_crm.all_hgvs]
        assert text_hgvs == ['NM_152443.3:c.677A>G',
                             'NG_008321.1:g.32324A>G',
                             'NC_000014.9:g.67729209A>G',
                             'NC_000014.8:g.68195926A>G',
                             'NM_152443.2:c.677A>G',
                             'Q96NR8:p.Tyr226Cys',
                             'NP_689656.2:p.Tyr226Cys']

    def test_preferred_current_hgvs(self):
        assert self.test_crm.preferred_current_hgvs.text == 'NC_000014.9:g.67729209A>G'

    def test_rs(self):
        assert self.test_crm.rs_id == 'rs28940313'

    def test_nsv(self):
        assert self.test_crm.nsv_id is None

    def test_variant_type(self):
        assert self.test_crm.variant_type == 'single nucleotide variant'

    def test_measure_set_pubmed_refs(self):
        assert self.test_crm.pubmed_refs == []

    def test_so_terms(self):
        assert self.test_crm.existing_so_terms == {'SO:0001583'}
