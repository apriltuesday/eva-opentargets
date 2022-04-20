from collections import defaultdict
import logging

from consequence_prediction.vep_mapping_pipeline.consequence_mapping import get_severity_ranking, \
    get_so_accessions

logger = logging.getLogger(__package__)


def process_gene(consequence_type_dict, variant_id, ensembl_gene_id, so_term):
    consequence_type_dict[variant_id].append(ConsequenceType(ensembl_gene_id, SoTerm(so_term)))


def process_consequence_type_dataframes(*dataframes):
    """
    Return a dictionary of consequence information extracted from one or more dataframes.
    Assumes all dataframes are in the same format.
    """
    consequence_type_dict = defaultdict(list)
    for consequences_dataframe in dataframes:
        for row in consequences_dataframe.itertuples():
            variant_id = row[1]
            ensembl_gene_id = row[3]
            so_term = row[5]

            process_gene(consequence_type_dict, variant_id, ensembl_gene_id, so_term)

    return consequence_type_dict


def process_consequence_type_file(snp_2_gene_file, consequence_type_dict=None):
    """
    Return a dictionary of consequence information extracted from the given file.
    If consequence_type_dict is provided then the information will be merge into this dictionary.
    """
    logger.info('Loading mapping rs -> ENSG/SOterms')
    if consequence_type_dict is None:
        consequence_type_dict = defaultdict(list)

    with open(snp_2_gene_file, "rt") as snp_2_gene_file:
        for line in snp_2_gene_file:
            line = line.rstrip()
            line_list = line.split("\t")

            if len(line_list) < 6:
                logger.warning('Skip invalid line in snp_2_gene file: {}'.format(line))
                continue

            variant_id = line_list[0]
            ensembl_gene_id = line_list[2]
            so_term = line_list[4]

            if ensembl_gene_id == 'NA':
                logger.warning('Skip line with missing gene ID: {}'.format(line))
                continue

            process_gene(consequence_type_dict, variant_id, ensembl_gene_id, so_term)

    logger.info('{} rs->ENSG/SOterms mappings loaded'.format(len(consequence_type_dict)))
    return consequence_type_dict


class SoTerm(object):
    """
    Represents a sequence ontology term belonging to a consequence type object.
    Holds information on accession and rank.
    """

    so_accession_name_dict = dict(get_so_accessions(), **{
        'trinucleotide_repeat_expansion': 'SO:0002165',
        'short_tandem_repeat_expansion': 'SO:0002162'
    })

    ranked_so_names_list = get_severity_ranking()

    def __init__(self, so_name):
        self.so_name = so_name
        if so_name in SoTerm.so_accession_name_dict:
            self._so_accession = SoTerm.so_accession_name_dict[so_name]
        else:
            self._so_accession = None

    @property
    def accession(self):
        if self._so_accession is not None:
            return self._so_accession.replace(':', '_')
        else:
            return None

    @property
    def rank(self):
        # If So name not in Ensembl's ranked list, return the least severe rank
        if self.so_name not in SoTerm.ranked_so_names_list:
            return len(SoTerm.ranked_so_names_list)
        else:
            return SoTerm.ranked_so_names_list.index(self.so_name)

    def __eq__(self, other):
        return self.accession == other.accession

    def __hash__(self):
        return hash(self._so_accession)


class ConsequenceType:
    """
    Holds information on the type of consequence related to a variation
    with relationship to ensembl gene IDs and SO terms
    """

    def __init__(self, ensembl_gene_id, so_term):
        self.ensembl_gene_id = ensembl_gene_id
        self.so_term = so_term

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)
