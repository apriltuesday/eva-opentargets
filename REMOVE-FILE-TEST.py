import os
from eva_cttv_pipeline.clinvar_xml_utils.src import clinvar_xml_utils

def get_test_clinvar_record(filename='test_clinvar_record.xml.gz'):
    """The default test file contains an extract of ClinVar XML for the record RCV000002127."""
    test_clinvar_record_file = "eva_cttv_pipeline/clinvar_xml_utils/data/test_clinvar_record.xml.gz"
    return [r for r in clinvar_xml_utils.ClinVarDataset(test_clinvar_record_file)][0]

a = get_test_clinvar_record()
print(a.measure)

print()
print(a.measure.hgvs)
