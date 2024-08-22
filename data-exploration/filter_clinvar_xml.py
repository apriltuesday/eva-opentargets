import argparse
import gzip
import logging
import xml.etree.ElementTree as ElementTree

from cmat.clinvar_xml_io import ClinVarRecord
from cmat.clinvar_xml_io.xml_parsing import find_mandatory_unique_element, iterate_cvs_from_xml
from cmat.output_generation.clinvar_to_evidence_strings import get_consequence_types
from cmat.output_generation.consequence_type import process_consequence_type_file

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# pretty print xml
def pprint(x):
    print(ElementTree.tostring(x, encoding='unicode'))


def filter_xml(input_xml, output_xml, filter_fct, max_num=None):
    """ Filter input_xml by boolean condition defined by filter_fct and write to output_xml.
    If max_num is given, will write at most max_num records, otherwise writes all."""
    header = b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
    <ReleaseSet Dated="." xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Type="full" xsi:noNamespaceSchemaLocation="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/clinvar_public_2.0.xsd">
    '''
    count = 0
    with gzip.open(output_xml, 'wb') as output_file:
        output_file.write(header)
        for raw_cvs_xml in iterate_cvs_from_xml(input_xml):
            rcv = find_mandatory_unique_element(raw_cvs_xml, 'ReferenceClinVarAssertion')
            record = ClinVarRecord(rcv, 2.0)
            if filter_fct(record):
                output_file.write(ElementTree.tostring(raw_cvs_xml))
                count += 1
            if max_num and count == max_num:
                break
        output_file.write(b'</ReleaseSet>')
    logger.info(f'Records written: {count}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Filter clinvar xml to get only variants with no functional consequences')
    parser.add_argument('--clinvar-xml', required=True)
    parser.add_argument('--output-xml', required=True)
    parser.add_argument('--gene-mapping', required=True)
    args = parser.parse_args()

    variant_to_gene_mappings = process_consequence_type_file(args.gene_mapping_file)

    def no_consequences(x: ClinVarRecord):
        return x.traits_with_valid_names and x.measure and not get_consequence_types(x.measure, variant_to_gene_mappings)

    filter_xml(
        input_xml=args.clinvar_xml,
        output_xml=args.output_xml,
        filter_fct=no_consequences,
        max_num=None,
    )
