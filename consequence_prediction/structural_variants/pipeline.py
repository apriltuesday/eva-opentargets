import logging
from itertools import zip_longest

import pandas as pd

from consequence_prediction.vep_mapping_pipeline.consequence_mapping import query_vep, VEP_SHORT_QUERY_DISTANCE, \
    extract_consequences, deduplicate_list
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io import ClinVarDataset
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.hgvs_variant import HgvsVariant, VariantType

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def grouper(iterable, n):
    args = [iter(iterable)] * n
    return [x for x in zip_longest(*args, fillvalue=None) if x is not None]


def accepted_by_vep(seq):
    # TODO be a little more specific here
    return seq.startswith('NC') or seq.startswith('LRG') or seq.startswith('NW') or seq.startswith('AC')


def hgvs_to_vep_identifier(hgvs):
    hgvs_variant = HgvsVariant(hgvs)

    seq = hgvs_variant.sequence_identifier
    if not accepted_by_vep(seq):
        return
    start = hgvs_variant.start
    stop = hgvs_variant.stop
    if not start or not stop or start > stop + 1:
        return
    variant_type = hgvs_variant.variant_type
    if variant_type not in {VariantType.DELETION, VariantType.DUPLICATION, VariantType.INSERTION}:
        return

    # TODO check the strand
    return f'{seq} {start} {stop} {variant_type[:3].upper()} + {hgvs}'


def can_process(record):
    """
    To ensure we don't step on other pipeline's toes, this pipeline will not even attempt to map consequences
    for if the measure has complete VCF-style coordinates.
    """
    return record.measure and record.measure.hgvs and not record.measure.has_complete_coordinates


def get_vep_results(clinvar_xml):
    n = 0
    hgvs_list = []
    for record in ClinVarDataset(clinvar_xml):
        if not can_process(record):
            continue
        # We allow possibly querying VEP multiple times for different HGVS expressions in a single measure
        # to maximise coverage, as consequences will be deduplicated later anyway.
        hgvs_list.extend([h for h in record.measure.current_hgvs if h is not None])
        n += 1
    logger.info(f'{n} records processed with {len(hgvs_list)} HGVS expressions')

    variants = [hgvs_to_vep_identifier(hgvs) for hgvs in hgvs_list]
    variants = [v for v in variants if v]  # v is None if it couldn't be parsed
    logger.info(f'{len(variants)} parsed into chrom/start/end/type')

    # VEP only accepts batches of 200
    i = 0
    vep_results = []
    for group in grouper(variants, n=200):
        vep_results.extend(query_vep(variants=group, search_distance=VEP_SHORT_QUERY_DISTANCE))
        i += 1
        logger.info(f'Done with batch {i}')

    return vep_results


def main(clinvar_xml):
    vep_results = get_vep_results(clinvar_xml)
    results_by_variant = {}
    results_by_variant = extract_consequences(vep_results=vep_results, acceptable_biotypes={'protein_coding', 'miRNA'},
                                              only_closest=False, results_by_variant=results_by_variant,
                                              report_distance=False)
    variant_data = []
    for variant_id, variant_consequences in results_by_variant.items():
        for consequence_to_yield in deduplicate_list(variant_consequences):
            # TODO confirm what the variant id should be
            variant_id, gene_id, gene_symbol, consequence_term, distance = consequence_to_yield
            # The second column, set statically to 1, is not used, and is maintained for compatibility purposes
            variant_data.append((variant_id, '1', gene_id, gene_symbol, consequence_term, distance))

    # Return as a dataframe to be compatible with repeat expansion pipeline
    consequences = pd.DataFrame(variant_data, columns=('VariantID', 'PlaceholderOnes', 'EnsemblGeneID',
                                                       'EnsemblGeneName', 'ConsequenceTerm', 'Distance'))
    return consequences
