from cmat.clinvar_xml_io.clinvar_reference_record import ClinVarReferenceRecord
from cmat.clinvar_xml_io.clinvar_submitted_record import ClinVarSubmittedRecord
from cmat.clinvar_xml_io.xml_parsing import find_mandatory_unique_element, find_elements


class ClinVarSet:
    """
    A ClinVarSet groups together a single reference record (RCV) and one or more submitted records (SCVs).
    """

    def __init__(self, cvs_xml, xsd_version):
        self.cvs_xml = cvs_xml

        rcv_elem = find_mandatory_unique_element(self.cvs_xml, 'ReferenceClinVarAssertion')
        self.rcv = ClinVarReferenceRecord(rcv_elem, xsd_version)

        scv_elems = find_elements(self.cvs_xml, 'ClinVarAssertion', allow_zero=False, allow_multiple=True)
        self.scvs = [ClinVarSubmittedRecord(elem, xsd_version, self.rcv) for elem in scv_elems]

    @property
    def id(self):
        return self.cvs_xml.attrib['ID']

    @property
    def title(self):
        return find_mandatory_unique_element(self.cvs_xml, './Title').text

    @property
    def status(self):
        return find_mandatory_unique_element(self.cvs_xml, './RecordStatus').text
