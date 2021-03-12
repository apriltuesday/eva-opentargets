#!/usr/bin/env python3

import argparse
from collections import Counter
import logging
import re

import eva_cttv_pipeline.clinvar_xml_utils as clinvar_xml_utils

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

    def __init__(self, name=None, width=None, height=None):
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


# The dicts store transition counts for the Sankey diagrams. Keys are (from, to), values are transition counts.
# Sankey diagrams can be visualised with SankeyMatic (see http://www.sankeymatic.com/build/).
sankey_variant_types = SankeyDiagram('variant-types.png', 1500, 750)
sankey_clinical_significance = SankeyDiagram('clinical-significance.png', 1000, 1000)
sankey_star_rating = SankeyDiagram('star-rating.png', 1000, 600)
sankey_mode_of_inheritance = SankeyDiagram('mode-of-inheritance.png', 1000, 1000)
sankey_allele_origin = SankeyDiagram('allele-origin.png', 400, 1500)

all_clinical_significance_levels = set()


# ClinVar XML have the following top-level structure:
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
            # Most common case, accounting for >99.95% of all ClinVar records.. Here, we go into details on various
            # attribute distributions.

            # Variant type
            measures = measure_set.findall('Measure')
            assert len(measures) == 1, 'MeasureSet of type Variant must contain exactly one Measure'
            sankey_variant_types.add_transitions(measure_set_type, measures[0].attrib['Type'])

            # Clinical significance
            clinical_significance = find_attribute(
                rcv, 'ClinicalSignificance/Description', 'ClinicalSignificance')
            all_clinical_significance_levels.add(clinical_significance)
            significance_type = 'Complex' if re.search('[,/]', clinical_significance) else 'Simple'
            sankey_clinical_significance.add_transitions('Variant', significance_type, clinical_significance)

            # Review status
            review_status = find_attribute(rcv, 'ClinicalSignificance/ReviewStatus', 'ReviewStatus')
            sankey_star_rating.add_transitions('Variant', review_status_stars(review_status), review_status)

            # Mode of inheritance
            mode_of_inheritance_xpath = 'AttributeSet/Attribute[@Type="ModeOfInheritance"]'
            mode_of_inheritance = find_attribute(rcv, mode_of_inheritance_xpath, 'ModeOfInheritance')
            if mode_of_inheritance.endswith('multiple'):
                # Having multiple ModeOfInheritance is rare. Log them for further investigation
                all_modes = '|'.join(sorted(mode.text for mode in rcv.findall(mode_of_inheritance_xpath)))
                print(f'Multiple ModeOfInheritance\t{rcv_id}\t{all_modes}')
            sankey_mode_of_inheritance.add_transitions(
                'Variant',
                mode_of_inheritance if mode_of_inheritance.endswith('missing') else 'ModeOfInheritance present',
            )
            if not mode_of_inheritance.endswith('missing'):
                sankey_mode_of_inheritance.add_transitions(
                    'ModeOfInheritance present', mode_of_inheritance
                )

    elif len(measure_sets) == 0 and len(genotype_sets) == 1:
        # RCV directly contains one genotype set.
        genotype_set = genotype_sets[0]
        sankey_variant_types.add_transitions('RCV', 'GenotypeSet', genotype_set.attrib['Type'])
    else:
        raise AssertionError('RCV must contain either exactly one measure set, or exactly one genotype set')

    allele_origins = {origin.text for origin in rcv.findall('ObservedIn/Sample/Origin')}
    if len(allele_origins) == 0:
        sankey_allele_origin.add_transitions('RCV', 'No allele origin')
    else:
        allele_origins_count = 'Single allele origin' if len(allele_origins) == 1 else 'Multiple allele origins'
        allele_origins_text = ','.join(sorted(allele_origins))
        sankey_allele_origin.add_transitions('RCV', allele_origins_count, allele_origins_text)

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
                       sankey_mode_of_inheritance, sankey_allele_origin):
    print('\n')
    print(sankey_diagram)

print('\n')
print('All clinical significance levels:')
for clin_sig in sorted(all_clinical_significance_levels):
    print(clin_sig)
