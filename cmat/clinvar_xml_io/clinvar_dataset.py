import gzip
import logging
import re
from datetime import date

from cmat.clinvar_xml_io.clinvar_record import ClinVarRecord
from cmat.clinvar_xml_io.xml_parsing import iterate_rcv_from_xml, parse_header_attributes

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class ClinVarDataset:
    """Iterate through records (RCVS) in ClinVar XML dump and convert them into internal ClinVarRecord representation."""
    def __init__(self, clinvar_xml):
        self.clinvar_xml = clinvar_xml
        self.header_attr = parse_header_attributes(clinvar_xml)
        self.header_attr['LastProcessed'] = date.today().strftime('%Y-%m-%d')
        self.xsd_version = self.get_xsd_version()

    def __iter__(self):
        for rcv in iterate_rcv_from_xml(self.clinvar_xml):
            yield ClinVarRecord(rcv, self.xsd_version)

    def get_xsd_version(self):
        # For format, see https://github.com/ncbi/clinvar/blob/master/FTPSiteXsdChanges.md
        if 'xsi:noNamespaceSchemaLocation' in self.header_attr:
            url = self.header_attr['xsi:noNamespaceSchemaLocation']
            m = re.search(r'(ClinVar_RCV|clinvar_public)_([0-9.]+)\.xsd', url, flags=re.IGNORECASE)
            if m and m.group(2):
                return float(m.group(2))
        # If we cannot parse, assume v2
        return 2.0

    def write_header(self, output_file):
        header = b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<ReleaseSet'
        for attr, val in self.header_attr.items():
            header += f' {attr}="{val}"'.encode('utf-8')
        header += b'>\n'
        output_file.write(header)

    def write(self, output_xml):
        """Writes the entire ClinVarDataset to a gzipped file at output_xml."""
        logger.info(f'Writing ClinVarDataset to: {output_xml}')
        count = 0
        with gzip.open(output_xml, 'wb') as output_file:
            self.write_header(output_file)
            for record in self:
                output_file.write(b'<ClinVarSet>\n  ')
                record.write(output_file)
                output_file.write(b'</ClinVarSet>\n')
                count += 1

            output_file.write(b'\n</ReleaseSet>')
        logger.info(f'Records written: {count}')
