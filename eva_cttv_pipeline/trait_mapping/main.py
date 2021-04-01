import csv
import logging
import multiprocessing

from eva_cttv_pipeline.clinvar_xml_utils import ClinVarTrait
from eva_cttv_pipeline.trait_mapping.output import output_trait
from eva_cttv_pipeline.trait_mapping.oxo import get_oxo_results
from eva_cttv_pipeline.trait_mapping.oxo import uris_to_oxo_format
from eva_cttv_pipeline.trait_mapping.trait import Trait
from eva_cttv_pipeline.trait_mapping.trait_names_parsing import parse_trait_names
from eva_cttv_pipeline.trait_mapping.zooma import get_zooma_results

logger = logging.getLogger(__package__)


def get_uris_for_oxo(zooma_result_list: list) -> set:
    """
    For a list of Zooma mappings return a list of uris for the mappings in that list with a high
    confidence.

    :param zooma_result_list: List with elements of class ZoomaResult
    :return: set of uris from high confidence Zooma mappings, for which to query OxO
    """
    uri_set = set()
    for mapping in zooma_result_list:
        # Only use high confidence Zooma mappings for querying OxO
        if mapping.confidence.lower() == "high":
            uri_set.update([entry.uri for entry in mapping.mapping_list])
    return uri_set


def process_trait(trait: Trait, filters: dict, zooma_host: str, oxo_target_list: list, oxo_distance: int) -> Trait:
    """
    Process a single trait. Find any mappings in Zooma. If there are no high confidence Zooma
    mappings that are in EFO then query OxO with any high confidence mappings not in EFO.

    :param trait: The trait to be processed.
    :param filters: A dictionary of filters to use for querying Zooma.
    :param zooma_host: A string with the hostname to use for querying Zooma
    :param oxo_target_list: A list of strings, each being an OxO ID for an ontology. Used to specify
                            which ontologies should be queried using OxO.
    :param oxo_distance: int specifying the maximum number of steps to use to query OxO. i.e. OxO's
                         "distance" parameter.
    :return: The original trait after querying Zooma and possibly OxO, with any results found.
    """
    logger.debug('Processing trait {}'.format(trait.name))

    trait.zooma_result_list = get_zooma_results(trait.name, filters, zooma_host)
    trait.process_zooma_results()
    if (trait.is_finished
            or len(trait.zooma_result_list) == 0
            or any([entry.is_current
                    for mapping in trait.zooma_result_list
                    for entry in mapping.mapping_list])):
        return trait

    uris_for_oxo_set = get_uris_for_oxo(trait.zooma_result_list)
    oxo_input_id_list = uris_to_oxo_format(uris_for_oxo_set)
    if len(oxo_input_id_list) == 0:
        return trait
    trait.oxo_result_list = get_oxo_results(oxo_input_id_list, oxo_target_list, oxo_distance)
    if not trait.oxo_result_list:
        logger.debug('No OxO mapping for trait {}'.format(trait.name))
    trait.process_oxo_mappings()

    return trait


def main(input_filepath, output_mappings_filepath, output_curation_filepath, filters, zooma_host, oxo_target_list,
         oxo_distance):
    logger.info('Started parsing trait names')
    trait_list = parse_trait_names(input_filepath)
    logger.info("Loaded {} trait names".format(len(trait_list)))

    with open(output_mappings_filepath, "w", newline='') as mapping_file, \
            open(output_curation_filepath, "wt") as curation_file:
        mapping_writer = csv.writer(mapping_file, delimiter="\t")
        mapping_writer.writerow(["#clinvar_trait_name", "uri", "label"])
        curation_writer = csv.writer(curation_file, delimiter="\t")

        logger.info('Processing trait names in parallel')
        trait_process_pool = multiprocessing.Pool(processes=24)
        processed_trait_list = [
            trait_process_pool.apply(
                process_trait,
                args=(trait, filters, zooma_host, oxo_target_list, oxo_distance)
            )
            for trait in trait_list
        ]

        logger.info('Writing output with the processed traits')
        for trait in processed_trait_list:
            # Remove non-specific trait names which should never be output
            if trait.name.lower() not in ClinVarTrait.NONSPECIFIC_TRAITS:
                output_trait(trait, mapping_writer, curation_writer)

    logger.info('Finished processing trait names')
