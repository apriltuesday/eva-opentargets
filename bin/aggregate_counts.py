#!/usr/bin/env python3

import argparse
import glob

from cmat.output_generation.report import Report

parser = argparse.ArgumentParser('Aggregate counts reports')
parser.add_argument('--counts-yml', nargs='+', help='YAML files containing intermediate counts', required=True)


if __name__ == '__main__':
    args = parser.parse_args()

    # Load all the reports
    filenames = [f for files in args.counts_yml for f in glob.glob(files)]
    reports = []
    for filename in filenames:
        r = Report()
        r.load_from_file(filename)
        reports.append(r)

    # Sum them up and output the results
    complete_report = sum(reports, start=Report())
    complete_report.print_report()
    complete_report.dump_to_file(dir_out='.')
    complete_report.write_unmapped_terms(dir_out='.')
    if not complete_report.check_counts():
        raise RuntimeError('Aggregate counts not consistent')
