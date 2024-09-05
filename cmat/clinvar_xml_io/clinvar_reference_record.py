import logging
from functools import cached_property

from cmat.clinvar_xml_io.clinical_classification import ClinicalClassification
from cmat.clinvar_xml_io.clinvar_record import ClinVarRecord
from cmat.clinvar_xml_io.xml_parsing import find_mandatory_unique_element, find_elements

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class ClinVarReferenceRecord(ClinVarRecord):
    """Reference records (RCVs) summarise information from submitted records (SCVs) and include additional annotations
    and cross-references supplied by ClinVar."""

    def __init__(self, record_xml, xsd_version):
        super().__init__(record_xml, xsd_version)

    def __str__(self):
        return f'ClinVarReferenceRecord object with accession {self.accession}'

    @property
    def last_updated_date(self):
        return self.record_xml.attrib['DateLastUpdated']

    @property
    def created_date(self):
        return self.record_xml.attrib['DateCreated']

    @cached_property
    def clinical_classifications(self):
        clinical_classifications = []
        if self.xsd_version < 2:
            # V1 only ever has a single clinical classification / clinical significance
            clinical_classifications.append(
                ClinicalClassification(find_mandatory_unique_element(self.record_xml, './ClinicalSignificance'), self))
        else:
            for clin_class in find_elements(self.record_xml, './Classifications/*'):
                clinical_classifications.append(ClinicalClassification(clin_class, self))
        return clinical_classifications
