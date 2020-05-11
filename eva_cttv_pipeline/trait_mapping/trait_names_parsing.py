import gzip
import json


def parse_trait_names(filepath: str) -> list:
    """
    For a file containing ClinVar records in the TSV format, return a list of the trait names for the records in the
    file.

    :param filepath: Path to a gzipped file containing ClinVar TSV summary.
    :return: A list of traits. It is important that it is a list and not a set, as it is used to calculate the
        frequencies of each trait later on.
    """

    # This is a dict of sets we use to keep track on which alleles are linked to which traits. The reason this is done
    # is because some (most) alleles appear twice in the document with coordinates for GRCh37 and GRCh38, and we don't
    # want to count them twice.
    trait_names_by_allele_id = dict()

    with gzip.open(filepath, "rt") as clinvar_summary:
        header = clinvar_summary.readline().rstrip().split('\t')
        for line in clinvar_summary:
            line = line.rstrip().split('\t')
            data = dict(header, values)
            # Extract allele ID and list of phenotypes
            allele_id = data['#Allele']
            traits = set(data['PhenotypeList'].split(';'))
            trait_names_by_allele_id.setdefault(allele_id, set())
            trait_names_by_allele_id[allele_id] |= traits

    # Now combine all traits into a single list
    trait_name_list = []
    for traits in trait_names_by_allele_id.values():
        trait_name_list.extend(traits)

    return [t.lower() for t in trait_name_list]
