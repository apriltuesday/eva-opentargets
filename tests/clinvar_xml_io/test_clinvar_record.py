import os

import pytest

from cmat.clinvar_xml_io import ClinVarDataset
from cmat.clinvar_xml_io.clinical_classification import MultipleClinicalClassificationsError

resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


def test_multiple_clinical_classifications_record():
    # input dataset with only one record
    input_file = os.path.join(resources_dir, 'multiple_classifications.xml.gz')
    record = next(iter(ClinVarDataset(input_file)))

    assert len(record.clinical_classifications) == 2
    assert set(cc.type for cc in record.clinical_classifications) == {'GermlineClassification', 'SomaticClinicalImpact'}
    with pytest.raises(MultipleClinicalClassificationsError):
        print(record.valid_clinical_significances)
