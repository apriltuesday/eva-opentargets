import os

test_dir = os.path.dirname(__file__)

snp_2_gene_file = os.path.join(test_dir, 'resources/snp2gene_extract.tsv')
test_clinvar_record_file = os.path.join(test_dir, 'resources/test_clinvar_record.json')
open_targets_schema_gz = os.path.join(test_dir, 'resources/opentargets.1.7.5.json.gz')
expected_genetics_evidence_string = os.path.join(test_dir, 'resources/expected_genetics_evidence_string.json')
expected_somatic_evidence_string = os.path.join(test_dir, 'resources/expected_somatic_evidence_string.json')
