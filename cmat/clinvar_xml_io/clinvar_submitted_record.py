import logging

from cmat.clinvar_xml_io import ClinVarRecord
from cmat.clinvar_xml_io.xml_parsing import find_mandatory_unique_element

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class ClinVarSubmittedRecord(ClinVarRecord):
    """
    Submitted records (SCVs) are structured similarly to reference records (RCVs), though typically with fewer
    annotations - for example, variant coordinates, HGVS expressions or ontology mappings which are added by curators.
    However, these attributes are also technically optional in the RCVs so the code inheritance is possible.

    SCVs also contain additional information about the actual submission, which we model in this class.
    """

    def __init__(self, record_xml, xsd_version, reference_record):
        super().__init__(record_xml, xsd_version)
        # Each SCV is associated with a single RCV
        self.reference_record = reference_record

    @property
    def submission_date(self):
        return find_mandatory_unique_element(self.record_xml, './ClinVarSubmissionID').attrib['submitterDate']

    @property
    def submitter(self):
        return find_mandatory_unique_element(self.record_xml, './ClinVarSubmissionID').attrib['submitter']

    @property
    def submission_name(self):
        # TODO - check whether this is the correct property to filter on
        return self.record_xml.attrib['SubmissionName']
