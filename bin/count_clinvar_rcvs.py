#!/usr/bin/env python3

import argparse

from cmat.clinvar_xml_io import ClinVarDataset

parser = argparse.ArgumentParser('Count number of RCV records in the XML, print to stdout')
parser.add_argument('--clinvar-xml', help='ClinVar XML release', required=True)


if __name__ == '__main__':
    args = parser.parse_args()
    print(sum(1 for _ in ClinVarDataset(args.clinvar_xml)))
