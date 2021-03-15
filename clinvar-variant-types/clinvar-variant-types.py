#!/usr/bin/env python3

import argparse
from collections import Counter
import logging
import re

import eva_cttv_pipeline.clinvar_xml_utils as clinvar_xml_utils

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

SIG_STARS = {
    'practice guideline': 4,
    'reviewed by expert panel': 3,
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, conflicting interpretations': 1,
    'criteria provided, single submitter': 1,
}

parser = argparse.ArgumentParser()
parser.add_argument('--clinvar-xml', required=True)
parser.add_argument('--process-items', type=int)
args = parser.parse_args()


class SankeyDiagram(Counter):

    def __init__(self, name, width, height):
        super().__init__()
        self.name = name
        self.width = width
        self.height = height

    def add_transitions(self, *transition_chain):
        """For example, if a transition chain is (RCV, MeasureSet, Variant), two transitions will be added to the Sankey
        diagram: (RCV → MeasureSet) and (MeasureSet → Variant)."""
        for t_from, t_to in zip(transition_chain, transition_chain[1:]):
            self[(t_from, t_to)] += 1

    def __str__(self):
        lines = [f'========== SANKEY DIAGRAM: {self.name} ==========',
                 f'Build using http://sankeymatic.com/build/ with width={self.width}, height={self.height}']
        for (t_from, t_to), t_count in sorted(self.items(), key=lambda x: -x[1]):
            lines.append(f'    {t_from} [{t_count}] {t_to}')
        return '\n'.join(lines)


class SupplementaryTable:

    def __init__(self, name, fields, sort_lambda=None):
        self.name = name
        self.fields = fields
        self.data = []
        self.sort_lambda = sort_lambda

    def add_row(self, row):
        assert len(row) == len(self.fields), 'Incorrect length of the row supplied'
        self.data.append(row)

    def __str__(self):
        lines = [f'>>>>>>>>>> SUPPLEMENTARY TABLE: {self.name} <<<<<<<<<<',
                 '|'.join(self.fields), '|'.join(':--' for _ in range(len(self.fields)))]
        sorted_data = sorted(self.data, key=self.sort_lambda)
        for row in sorted_data:
            lines.append('|'.join([str(v) for v in row]))
        return '\n'.join(lines)


class SupplementaryTableCounter(SupplementaryTable):

    def __init__(self, name, field_name):
        super().__init__(name=name, fields=(field_name, 'Count'), sort_lambda=lambda x: -x[1])
        self.counter = Counter()

    def add_count(self, item):
        self.counter[item] += 1

    def __str__(self):
        # On demand, tally the counter, add the rows, and return the standard string representation
        for item, count in self.counter.items():
            self.add_row([item, count])
        return super().__str__()


def find_attribute(rcv, xpath, attribute_name):
    """Find an attribute in the RCV record which can have either zero or one occurrence. Return a textual representation
    of the attribute, including special representations for the case of zero or multiple, constructed using the
    attribute_name parameter."""

    attributes = rcv.findall(xpath)
    if len(attributes) == 0:
        return '{} missing'.format(attribute_name)
    elif len(attributes) == 1:
        return attributes[0].text
    else:
        return '{} multiple'.format(attribute_name)


def review_status_stars(review_status):
    black_stars = SIG_STARS.get(review_status, 0)
    white_stars = 4 - black_stars
    return '★' * black_stars + '☆' * white_stars


# Sankey diagrams for visualisation
sankey_variant_types = SankeyDiagram('variant-types.png', 1200, 600)
sankey_clinical_significance = SankeyDiagram('clinical-significance.png', 1200, 600)
sankey_star_rating = SankeyDiagram('star-rating.png', 1200, 600)
sankey_mode_of_inheritance = SankeyDiagram('mode-of-inheritance.png', 1200, 500)
sankey_allele_origin = SankeyDiagram('allele-origin.png', 1200, 600)
sankey_inheritance_origin = SankeyDiagram('inheritance-origin.png', 1200, 400)

# Supplementary tables and counters for the report
counter_clin_sig_complex = SupplementaryTableCounter('Complex clinical significance levels', 'Clinical significance')
counter_clin_sig_all = SupplementaryTableCounter('All clinical significance levels', 'Clinical significance')
counter_star_rating = SupplementaryTableCounter('Distribuion of records by star rating', 'Star rating')
table_multiple_mode_of_inheritance = SupplementaryTable('Multiple mode of inheritance', ['RCV', 'Modes of inheritance'])
counter_multiple_allele_origin = SupplementaryTableCounter('Multiple allele origins', 'Allele origins')
table_inconsistent_moi_ao = SupplementaryTable('Inconsistent mode of inheritance and allele origin values',
                                               ['RCV', 'Modes of inheritance', 'Allele origins'])

# ClinVar XML has the following top-level structure:
#   <ReleaseSet>
#     <ClinVarSet>...</ClinVarSet>
#     <ClinVarSet>...</ClinVarSet>
#     ...
#   </ReleaseSet>
# To not load everything into memory, we use iterparse, wait for a complete ClinVarSet element, and then remove it from
# the tree because we don't need it anymore.

elements_processed = 0
for rcv in clinvar_xml_utils.iterate_rcv_from_xml(args.clinvar_xml):
    # TODO: use clinvar_xml_utils abstractions instead of working with XML directly
    # TODO: https://github.com/EBIvariation/eva-opentargets/issues/170
    rcv_id = rcv.find('./ClinVarAccession').attrib['Acc']

    # RCV can contain either a MeasureSet, or a GenotypeSet. It must not contain both.
    measure_sets = rcv.findall('MeasureSet')
    genotype_sets = rcv.findall('GenotypeSet')
    if len(measure_sets) == 1 and len(genotype_sets) == 0:
        # Most common case. RCV directly contains one measure set.
        measure_set = measure_sets[0]
        measure_set_type = measure_set.attrib['Type']
        sankey_variant_types.add_transitions('RCV', 'MeasureSet', measure_set_type)

        if measure_set_type == 'Variant':
            # Most common case, accounting for >99.97% of all ClinVar records. Here, we go into details on various
            # attribute distributions.

            # Variant type
            measures = measure_set.findall('Measure')
            assert len(measures) == 1, 'MeasureSet of type Variant must contain exactly one Measure'
            sankey_variant_types.add_transitions(measure_set_type, measures[0].attrib['Type'])

            # Clinical significance
            clinical_significance = find_attribute(rcv, 'ClinicalSignificance/Description', 'ClinicalSignificance')
            clin_sig_split = re.split(', |/', clinical_significance)
            for clin_sig in clin_sig_split:  # Count all clinical significance levels after splitting
                counter_clin_sig_all.add_count(clin_sig)
            if len(clin_sig_split) == 1:
                sankey_clinical_significance.add_transitions('Variant', 'Simple', clinical_significance)
            else:
                sankey_clinical_significance.add_transitions('Variant', 'Complex')
                counter_clin_sig_complex.add_count(clinical_significance)

            # Review status
            review_status = find_attribute(rcv, 'ClinicalSignificance/ReviewStatus', 'ReviewStatus')
            star_rating = review_status_stars(review_status)
            sankey_star_rating.add_transitions('Variant', star_rating, review_status)
            counter_star_rating.add_count(star_rating)

            # Mode of inheritance
            modes_of_inheritance = {moi.text for moi in rcv.findall(
                'AttributeSet/Attribute[@Type="ModeOfInheritance"]')}
            modes_of_inheritance_text = ', '.join(sorted(modes_of_inheritance))
            if len(modes_of_inheritance) == 0:
                mode_of_inheritance_category = 'Missing'
            elif 'Somatic mutation' in modes_of_inheritance:
                if len(modes_of_inheritance) > 1:
                    mode_of_inheritance_category = 'Germline & somatic'
                else:
                    mode_of_inheritance_category = 'Somatic'
            else:
                mode_of_inheritance_category = 'Germline'
            sankey_mode_of_inheritance.add_transitions('Variant', mode_of_inheritance_category)
            if mode_of_inheritance_category == 'Germline':
                if len(modes_of_inheritance) == 1:
                    sankey_mode_of_inheritance.add_transitions('Germline', 'Single', modes_of_inheritance_text)
                else:
                    sankey_mode_of_inheritance.add_transitions('Germline', 'Multiple')
            # Log multiple ModeOfInheritance cases in a separate table
            if len(modes_of_inheritance) > 1:
                table_multiple_mode_of_inheritance.add_row([rcv_id, modes_of_inheritance_text])

            # Allele origins
            allele_origins = {origin.text for origin in rcv.findall('ObservedIn/Sample/Origin')}
            allele_origin_text = ', '.join(sorted(allele_origins))
            if len(allele_origins) == 0:
                allele_origin_category = 'Missing'
            elif 'somatic' in allele_origins:
                if len(allele_origins) > 1:
                    allele_origin_category = 'Germline & somatic'
                else:
                    allele_origin_category = 'Somatic'
            else:
                allele_origin_category = 'Germline'
            sankey_allele_origin.add_transitions('Variant', allele_origin_category)
            if allele_origin_category == 'Germline':
                if len(allele_origins) == 1:
                    sankey_allele_origin.add_transitions(allele_origin_category, 'Single', allele_origin_text)
                else:
                    sankey_allele_origin.add_transitions(allele_origin_category, 'Multiple')
            # Log multiple allele of origin values in a separate table
            if len(allele_origins) > 1:
                counter_multiple_allele_origin.add_count(allele_origin_text)

            # Mode of inheritance and allele origin mapping
            if mode_of_inheritance_category != 'Missing' and allele_origin_category != 'Missing':
                sankey_inheritance_origin.add_transitions(
                    f'[MoI] {mode_of_inheritance_category}', f'{allele_origin_category} [AO]')
                if mode_of_inheritance_category != allele_origin_category:
                    table_inconsistent_moi_ao.add_row([rcv_id, modes_of_inheritance_text, allele_origin_text])

    elif len(measure_sets) == 0 and len(genotype_sets) == 1:
        # RCV directly contains one genotype set.
        genotype_set = genotype_sets[0]
        sankey_variant_types.add_transitions('RCV', 'GenotypeSet', genotype_set.attrib['Type'])
    else:
        raise AssertionError('RCV must contain either exactly one measure set, or exactly one genotype set')

    # Track the number of already processed elements
    elements_processed += 1
    if elements_processed % 100000 == 0:
        logger.info(f'Processed {elements_processed} elements')
    if args.process_items and elements_processed >= args.process_items:
        break

logger.info(f'Done. Processed {elements_processed} elements')


# Output the code for Sankey diagrams. Transitions are sorted in decreasing number of counts, so that the most frequent
# cases are on top.
for sankey_diagram in (sankey_variant_types, sankey_clinical_significance, sankey_star_rating,
                       sankey_mode_of_inheritance, sankey_allele_origin, sankey_inheritance_origin):
    print('\n')
    print(sankey_diagram)


# Output the supplementary tables for the report.
for supplementary_table in (counter_clin_sig_complex, counter_clin_sig_all, counter_star_rating,
                            table_multiple_mode_of_inheritance, counter_multiple_allele_origin,
                            table_inconsistent_moi_ao):
    print('\n')
    print(supplementary_table)

print('\n')
