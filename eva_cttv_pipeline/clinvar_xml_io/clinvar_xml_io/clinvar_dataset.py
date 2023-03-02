import gzip
import logging
from datetime import date

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.clinvar_record import ClinVarRecord
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.xml_parsing import iterate_rcv_from_xml

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class ClinVarDataset:
    """Iterate through records in ClinVar XML dump and convert them into internal ClinVarRecord representation."""
    def __init__(self, clinvar_xml):
        self.clinvar_xml = clinvar_xml
        self.header = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<ReleaseSet Dated="{self.today()}" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Type="full" xsi:noNamespaceSchemaLocation="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/clinvar_public_1.60.xsd">
'''
        self.header = self.header.encode('utf-8')

    def __iter__(self):
        for rcv in iterate_rcv_from_xml(self.clinvar_xml):
            yield ClinVarRecord(rcv)

    def today(self):
        return date.today().strftime('%Y-%m-%d')

    def write(self, output_xml):
        """Writes the entire ClinVarDataset to a gzipped file at output_xml."""
        logger.info(f'Writing ClinVarDataset to: {output_xml}')
        count = 0
        with gzip.open(output_xml, 'wb') as output_file:
            output_file.write(self.header)
            for record in self:
                output_file.write(b'<ClinVarSet>\n  ')
                record.write(output_file)
                output_file.write(b'</ClinVarSet>\n')
                count += 1

            output_file.write(b'\n</ReleaseSet>')
        logger.info(f'Records written: {count}')
