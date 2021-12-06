import argparse
from collections import Counter
import os

import hgvs.parser

from eva_cttv_pipeline.clinvar_xml_io.clinvar_identifier_parsing import parse_variant_identifier
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io import ClinVarDataset, find_elements

hgvs_parser = hgvs.parser.Parser()


def counts(input_xml, output_dir):
    """
    Gather counts for variant descriptors present in input_xml.
    Returns counts and outputs a text file of all unparseable hgvs descriptors to output_dir.
    """
    dataset = ClinVarDataset(input_xml)

    variant_type_hist = Counter()
    total_count = 0
    spdi_count = 0
    hgvs_count = 0
    toplevel_refseq_hgvs_count = 0
    cytogenetic_count = 0
    start_stop_count = 0
    strict_parseable_hgvs_count = 0
    our_parseable_hgvs_count = 0
    any_parseable_hgvs_count = 0
    hgvs_only_count = 0
    cytogenetic_only_count = 0
    start_stop_only_count = 0

    hgvs_output_file = open(os.path.join(output_dir, 'unparseable-hgvs.txt'), 'w+')

    for record in dataset:
        if not record.measure:
            continue
        m = record.measure
        if m.has_complete_coordinates:
            continue
        total_count += 1
        variant_type_hist[m.variant_type] += 1

        has_hgvs = m.hgvs
        has_cytogenetic = find_elements(m.measure_xml, './CytogeneticLocation')
        has_spdi = find_elements(m.measure_xml, './CanonicalSPDI')
        has_start_stop = m.chr and m.sequence_location_helper('start') and m.sequence_location_helper('stop')

        # exclusive counts
        if has_hgvs and not has_cytogenetic and not has_start_stop:
            hgvs_only_count += 1
        elif has_cytogenetic and not has_hgvs and not has_start_stop:
            cytogenetic_only_count += 1
        elif has_start_stop and not has_hgvs and not has_cytogenetic:
            start_stop_only_count += 1

        # non-exclusive counts
        if has_cytogenetic:
            cytogenetic_count += 1
        if has_spdi:
            spdi_count += 1
        if has_start_stop:
            start_stop_count += 1
        if has_hgvs:
            hgvs_count += 1
            if m.toplevel_refseq_hgvs:
                toplevel_refseq_hgvs_count += 1

            # hgvs parseability
            one_strict_parseable = False
            one_our_parseable = False
            for hgvs in m.hgvs:
                try:
                    hgvs_parser.parse_hgvs_variant(hgvs)
                    one_strict_parseable = True
                except:
                    pass
                try:
                    if any(parse_variant_identifier(hgvs)):
                        one_our_parseable = True
                    else:
                        hgvs_output_file.write(hgvs + '\n')
                except:
                    pass  # these are None
            if one_strict_parseable:
                strict_parseable_hgvs_count += 1
            if one_our_parseable:
                our_parseable_hgvs_count += 1
            if one_strict_parseable or one_our_parseable:
                any_parseable_hgvs_count += 1

    hgvs_output_file.close()

    # collect results
    other_counts = {
        'Canonical SPDI': spdi_count,
        'Any HGVS': hgvs_count,
        'Top level refseq HGVS': toplevel_refseq_hgvs_count,
        'Cytogenetic location': cytogenetic_count,
        'Chr, start, stop': start_stop_count,
        'Official HGVS parser': strict_parseable_hgvs_count,
        'Our HGVS parser': our_parseable_hgvs_count,
        'Any HGVS parser': any_parseable_hgvs_count,
    }
    exclusive_counts = {
        'Only HGVS': hgvs_only_count,
        'Only cytogenetic': cytogenetic_only_count,
        'Only chr/start/stop': start_stop_only_count
    }
    return total_count, variant_type_hist, other_counts, exclusive_counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser('')
    parser.add_argument('--input-xml', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    total_count, variant_type_hist, other_counts, exclusive_counts = counts(args.input_xml, args.output_dir)
