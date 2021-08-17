import argparse
import gzip
import logging
import xml.etree.ElementTree as ElementTree

from eva_cttv_pipeline.clinvar_xml_utils import ClinVarRecord, find_mandatory_unique_element
from eva_cttv_pipeline.evidence_string_generation.clinvar_to_evidence_strings import get_consequence_types
from eva_cttv_pipeline.evidence_string_generation.consequence_type import process_consequence_type_file


logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def iterate_cvs_from_xml(clinvar_xml):
    """Similar to iterate_rcv_from_xml in clinvar_xml_utils, but keeps the entire ClinVarSet XML element.
    This allows us to construct a valid ClinVar XML for easy future processing."""
    with gzip.open(clinvar_xml, 'rt') as fh:
        for event, elem in ElementTree.iterparse(fh):
            # Wait until we have built a complete ClinVarSet element
            if elem.tag != 'ClinVarSet':
                continue
            yield elem
            elem.clear()


def filter_xml(input_xml, output_xml, filter_fct, max_num=None):
    """ Filter input_xml by boolean condition defined by filter_fct and write to output_xml.
    If max_num is given, will write at most max_num records, otherwise writes all."""
    header = b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
    <ReleaseSet Dated="2020-03-02" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Type="full" xsi:noNamespaceSchemaLocation="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/clinvar_public_1.60.xsd">
    '''
    count = 0
    with gzip.open(output_xml, 'wb') as output_file:
        output_file.write(header)
        for raw_cvs_xml in iterate_cvs_from_xml(input_xml):
            rcv = find_mandatory_unique_element(raw_cvs_xml, 'ReferenceClinVarAssertion')
            record = ClinVarRecord(rcv)
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
