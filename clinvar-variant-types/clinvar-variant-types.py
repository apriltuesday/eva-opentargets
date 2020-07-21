#!/usr/bin/env python3

import argparse
from collections import Counter
import gzip
import re
import sys
import xml.etree.ElementTree as ElementTree

PROCESSED_CLIN_SIG = ['Pathogenic', 'Likely pathogenic', 'protective', 'association', 'risk_factor', 'affects',
                      'drug response']

SIG_STARS = {
    'practice guideline': 4,
    'reviewed by expert panel': 3,
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, conflicting interpretations': 1,
    'criteria provided, single submitter': 1,
}

parser = argparse.ArgumentParser()
parser.add_argument('--clinvar-xml', required=True)
args = parser.parse_args()


def add_transitions(transitions_counter, transition_chain):
    """Increments the count of a particular flow in Sankey diagram."""
    for transition_from, transition_to in zip(transition_chain, transition_chain[1:]):
        transitions_counter[(transition_from, transition_to)] += 1


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


def clin_sig_color(clinical_significance):
    return '{} {}'.format(
        clinical_significance,
        '#ff0000' if clinical_significance in PROCESSED_CLIN_SIG else ''
    )


def review_status_stars(review_status):
    black_stars = SIG_STARS.get(review_status, 0)
    white_stars = 4 - black_stars
    return '★' * black_stars + '☆' * white_stars


# The dicts store transition counts for the Sankey diagrams. Keys are (from, to), values are transition counts.
# Sankey diagrams can be visualised with SankeyMatic (see http://www.sankeymatic.com/build/).
variant_type_transitions, clin_sig_transitions, review_status_transitions, inheritance_mode_transitions \
    = Counter(), Counter(), Counter(), Counter()


# ClinVar XML have the following top-level structure:
#   <ReleaseSet>
#     <ClinVarSet>...</ClinVarSet>
#     <ClinVarSet>...</ClinVarSet>
#     ...
#   </ReleaseSet>
# To not load everything into memory, we use iterparse, wait for a complete ClinVarSet element, and then remove it from
# the tree because we don't need it anymore.

elements_processed = 0
for event, elem in ElementTree.iterparse(gzip.open(args.clinvar_xml)):

    # Wait until we have built a complete ClinVarSet element. Skip
    if elem.tag != 'ClinVarSet':
        continue

    # Go to a ReferenceClinVarAssertion element. This corresponds to a single RCV record, the main unit of ClinVar.
    # There should only be one such record per ClinVarSet.
    rcv_records = elem.findall('ReferenceClinVarAssertion')
    assert len(rcv_records) == 1, 'Found multiple RCV records per ClinVarSet'
    rcv = rcv_records[0]

    # RCV can contain either a MeasureSet, or a GenotypeSet. It must not contain both.
    measure_sets = rcv.findall('MeasureSet')
    genotype_sets = rcv.findall('GenotypeSet')

    if len(measure_sets) == 1 and len(genotype_sets) == 0:
        # Most common case. RCV directly contains one measure set.
        measure_set = measure_sets[0]
        measure_set_type = measure_set.attrib['Type']
        add_transitions(variant_type_transitions, ('RCV', 'MeasureSet', measure_set_type))

        if measure_set_type == 'Variant':
            # Most common case, accounting for >99.95% of all ClinVar records.. Here, we go into details on various
            # attribute distributions.

            # Variant type
            measures = measure_set.findall('Measure')
            assert len(measures) == 1, 'MeasureSet of type Variant must contain exactly one Measure'
            add_transitions(variant_type_transitions, (measure_set_type, measures[0].attrib['Type']))

            # Clinical significance
            clinical_significance = find_attribute(
                rcv, 'ClinicalSignificance/Description', 'ClinicalSignificance')
            if clinical_significance in PROCESSED_CLIN_SIG:
                add_transitions(clin_sig_transitions, (
                    'Variant',
                    'Processed',
                    clinical_significance,
                ))
            else:
                significance_type = 'Complex' if re.search('[,/]', clinical_significance) else 'Simple'
                add_transitions(clin_sig_transitions, (
                    'Variant',
                    'Not processed',
                    significance_type,
                    clinical_significance,
                ))

            # Review status
            review_status = find_attribute(
                rcv, 'ClinicalSignificance/ReviewStatus', 'ReviewStatus')
            add_transitions(review_status_transitions, (
                'Variant',
                review_status_stars(review_status),
                review_status,
            ))

            # Mode of inheritance
            mode_of_inheritance = find_attribute(
                rcv, 'AttributeSet/Attribute[@Type="ModeOfInheritance"]', 'ModeOfInheritance')
            add_transitions(inheritance_mode_transitions, (
                'Variant',
                mode_of_inheritance if mode_of_inheritance.endswith('missing') else 'ModeOfInheritance present',
            ))
            if not mode_of_inheritance.endswith('missing'):
                add_transitions(inheritance_mode_transitions, (
                    'ModeOfInheritance present', mode_of_inheritance
                ))

    elif len(measure_sets) == 0 and len(genotype_sets) == 1:
        # RCV directly contains one genotype set.
        genotype_set = genotype_sets[0]
        add_transitions(variant_type_transitions, ('RCV', 'GenotypeSet', genotype_set.attrib['Type']))

    else:
        raise AssertionError('RCV must contain either exactly one measure set, or exactly one genotype set')

    # Remove the processed element from the tree to save memory
    elem.clear()

    # Track the number of already processed elements
    elements_processed += 1
    if elements_processed % 10000 == 0:
        print('Processed {} elements'.format(elements_processed), file=sys.stderr)

# Output the code for Sankey diagram. Transitions are sorted in decreasing number of counts, so that the most frequent
# cases are on top.
for transitions_counter in (variant_type_transitions, clin_sig_transitions, review_status_transitions,
                            inheritance_mode_transitions):
    print()
    for (transition_from, transition_to), count in sorted(transitions_counter.items(), key=lambda x: -x[1]):
        print('{transition_from} [{count}] {transition_to}'.format(**locals()))
