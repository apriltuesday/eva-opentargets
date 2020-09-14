from collections import Counter
import gzip
import json

from eva_cttv_pipeline.trait_mapping.trait import Trait


def parse_trait_names(filepath: str) -> list:
    """
    For a file containing ClinVar records in the TSV format, return a list of Traits for the records in the file. Each
    Trait object contains trait name, how many times it occurs in the input file, and whether it is linked to an NT
    expansion variant.

    Trait occurrence count is calculated based on all unique (AlleleID, RCV, trait name) tuples in the input file. This
    is because each such tuple will, generally speaking, correspond to one output evidence string. So if we want to
    gauge which trait names are more important to curate, we need to consider how many such tuples it appears in.

    The reason we need to keep track of only *unique* tuples is because some (most) alleles will appear twice in the
    document with coordinates for GRCh37 and GRCh38, and we don't want to count them twice.

    Traits which are implicated in "NT expansion" variants are marked using a special field, because their curation is
    of highest importance even if the number of records which they are linked to is low.

    :param filepath: Path to a gzipped file containing ClinVar TSV summary.
    :return: A list of Trait objects.
    """

    # Tracks unique (AlleleID, RCV, trait name) tuples
    unique_association_tuples = set()

    # Tracks all traits which are at least once implicated in "NT expansion", or nucleotide repeat expansion, variants.
    # Their curation is of highest importantce regardless of how many records they are actually associated with.
    nt_expansion_traits = set()

    with gzip.open(filepath, "rt") as clinvar_summary:
        header = clinvar_summary.readline().rstrip().split('\t')
        for line in clinvar_summary:
            values = line.rstrip().split('\t')
            data = dict(zip(header, values))

            # Extract relevant fields
            is_nt_expansion_variant = data['Type'] == 'NT expansion'
            allele_id = data['#AlleleID']
            traits = set(data['PhenotypeList'].split(';'))
            rcv_ids = set(data['RCVaccession'].split(';'))

            # Process all (trait, rcv) records
            for trait, rcv_id in zip(traits, rcv_ids):
                unique_association_tuples.add((trait, rcv_id, allele_id))
                if is_nt_expansion_variant:
                    nt_expansion_traits.add(trait)

    # Count trait occurrences
    trait_names = [t[0] for t in unique_association_tuples]
    traits = []
    for trait_name, trait_frequency in Counter(trait_names).items():
        if trait_name == '-':
            print('Skipped {} missing trait names'.format(trait_frequency))
            continue
        associated_with_nt_expansion = trait_name in nt_expansion_traits
        traits.append(Trait(name=trait_name.lower(), frequency=trait_frequency,
                            associated_with_nt_expansion=associated_with_nt_expansion))

    return traits
