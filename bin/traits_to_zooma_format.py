#!/usr/bin/env python3

import argparse
import itertools
from time import gmtime, strftime

from cmat import clinvar_xml_io
from cmat.clinvar_xml_io.ontology_uri import OntologyUri


def write_zooma_record(clinvar_acc, variant_id, trait_name, ontology_uri, date, outfile):
    zooma_output_list = [clinvar_acc,
                         variant_id,
                         'disease',
                         trait_name,
                         str(ontology_uri),
                         'clinvar-xrefs',
                         date]
    outfile.write('\t'.join(zooma_output_list) + '\n')


def process_clinvar_record(clinvar_record, outfile):
    """Extract the variant, trait and ontology from Clinvar record and write them as a ZOOMA feedback record."""
    if clinvar_record.measure is None:
        return
    variant_ids = [variant_id
                   for variant_id in (clinvar_record.measure.rs_id, clinvar_record.measure.nsv_id)
                   if variant_id is not None]
    traits = clinvar_record.traits
    for variant_id, trait in itertools.product(variant_ids, traits):
        if trait.preferred_or_other_valid_name is None:
            continue
        for db, identifier, status in trait.xrefs:
            if status != 'current' or db.lower() not in OntologyUri.db_to_uri_conversion:
                continue
            ontology_uri = OntologyUri(identifier, db)
            write_zooma_record(clinvar_record.accession, variant_id, trait.preferred_or_other_valid_name, ontology_uri,
                               strftime('%d/%m/%y %H:%M', gmtime()), outfile)


def main(clinvar_xml, zooma_feedback):
    with open(zooma_feedback, 'wt') as outfile:
        outfile.write('STUDY\tBIOENTITY\tPROPERTY_TYPE\tPROPERTY_VALUE\tSEMANTIC_TAG\tANNOTATOR\tANNOTATION_DATE\n')
        for clinvar_record in clinvar_xml_io.ClinVarDataset(clinvar_xml):
            process_clinvar_record(clinvar_record, outfile)


parser = argparse.ArgumentParser('Extracts OMIM and MedGen cross-references from ClinVar for submission to ZOOMA')
parser.add_argument('--clinvar-xml',    help='ClinVar XML release',                 required=True)
parser.add_argument('--zooma-feedback', help='Disease string to ontology mappings', required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.clinvar_xml, args.zooma_feedback)
