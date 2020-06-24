import gzip
import json


# TODO: unify the acceptable clinical significance levels project-wide
ACCEPTABLE_CLINICAL_SIGNIFICANCE = ['pathogenic', 'likely pathogenic', 'protective', 'association', 'risk_factor',
                                    'affects', 'drug response']


def parse_trait_names(filepath: str) -> list:
    """
    For a file containing ClinVar records in the TSV format, return a list of the trait names for the records in the
    file.

    :param filepath: Path to a gzipped file containing ClinVar TSV summary.
    :return: A list of traits. It is important that it is a list and not a set, as it is used to calculate the
        frequencies of each trait later on.
    """

    # This set stores all unique (AlleleID, RCV, trait name) tuples. The reason for that is each such tuple will,
    # generally speaking, correspond to one output evidence string. So if we want to gauge which trait names are more
    # important to curate, we need to consider how many such tuples it appears in. The reason we need to keep track of
    # only unique tuples is because some (most) alleles will appear twice in the document with coordinates for GRCh37
    # and GRCh38, and we don't want to count them twice.
    unique_association_records = set()

    with gzip.open(filepath, "rt") as clinvar_summary:
        header = clinvar_summary.readline().rstrip().split('\t')
        for line in clinvar_summary:
            values = line.rstrip().split('\t')
            data = dict(zip(header, values))

            # Check if the record should be processed given its level of clinical significance
            acceptable_clinical_significance_present = False
            for clinical_significance in data['ClinicalSignificance'].split(','):
                if clinical_significance.strip().lower() in ACCEPTABLE_CLINICAL_SIGNIFICANCE:
                    acceptable_clinical_significance_present = True
            if not acceptable_clinical_significance_present:
                continue

            # Extract relevant tuple elements
            allele_id = data['#AlleleID']
            traits = set(data['PhenotypeList'].split(';'))
            rcv_ids = set(data['RCVaccession'].split(';'))
            for trait, rcv_id in zip(traits, rcv_ids):
                unique_association_records |= {(allele_id, trait, rcv_id)}

    # Now combine all traits into a single list
    trait_name_list = []
    for _, trait, _ in unique_association_records:
        trait_name_list.extend(traits)

    return [t.lower() for t in trait_name_list]
