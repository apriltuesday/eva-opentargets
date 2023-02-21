import logging

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.clinvar_record import ClinVarRecord
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.utils import iterate_rcv_from_xml

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class ClinVarDataset:
    """Iterate through records in ClinVar XML dump and convert them into internal ClinVarRecord representation."""
    def __init__(self, clinvar_xml):
        self.clinvar_xml = clinvar_xml

    def __iter__(self):
        for rcv in iterate_rcv_from_xml(self.clinvar_xml):
            yield ClinVarRecord(rcv)
