import os

import pytest

from cmat.clinvar_xml_io import ClinVarDataset
from cmat.clinvar_xml_io.clinical_classification import MultipleClinicalClassificationsError

resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


def test_clinvar_record_v1_v2_migration():
    # Same RCV in both versions
    v1_input_file = os.path.join(resources_dir, 'clinvar_dataset_v1.xml.gz')
    v2_input_file = os.path.join(resources_dir, 'clinvar_dataset_v2.xml.gz')
    record_v1 = next(iter(ClinVarDataset(v1_input_file)))
    record_v2 = next(iter(ClinVarDataset(v2_input_file)))
    assert record_v1.xsd_version == 1.6
    assert record_v2.xsd_version == 2.0

    assert record_v1.accession == record_v2.accession
    assert record_v1.score == record_v2.score
    assert set(record_v1.valid_clinical_significances) == set(record_v2.valid_clinical_significances)
    assert set(record_v1.valid_allele_origins) == set(record_v2.valid_allele_origins)


def test_multiple_clinical_classifications_record():
    # input dataset with only one record
    input_file = os.path.join(resources_dir, 'multiple_classifications.xml.gz')
    record = next(iter(ClinVarDataset(input_file)))

    assert len(record.clinical_classifications) == 2
    assert set(cc.type for cc in record.clinical_classifications) == {'GermlineClassification', 'SomaticClinicalImpact'}
    with pytest.raises(MultipleClinicalClassificationsError):
        print(record.valid_clinical_significances)
