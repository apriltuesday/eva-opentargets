import logging
import re

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.hgvs_variant import HgvsVariant, SequenceType

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Human-readable description notation. This pattern is used for about a dozen variants. From it we only
# extract the repeat unit sequence. Example: 'ATXN8, (CAG)n REPEAT EXPANSION'
re_description = re.compile(
    r'\('
    r'(?P<sequence>[{}]+)'.format(IUPACAmbiguousDNA.letters) +
    r'\)n'
    r'(?: REPEAT)? EXPANSION'
)


def parse_repeat_identifier(variant_name):
    """
    Parses variant identifier and extract certain characteristics:
        * TranscriptID: NCBI RefSeq transcript ID, e.g. NM_000044
        * CoordinateSpan: the distance between start and end coordinates described in the HGVS-like notation.
              E. g., for 'NM_000044.4(AR):c.172_174CAG(7_34) (p.Gln66_Gln80del)', it will be 174 - 172 + 1 = 3
        * RepeatUnitLength: length of the sequence being repeated.
              E. g., for 'NC_000004.11:g.3076606GCA[27_35]', it will be len('GCA') = 3
        * IsProteinHGVS: this field simply reflects whether the variant was defined using a protein HGVS-like notation.
              E.g., for 'NP_002964.3:p.Gln166(>=33)' it will be TRUE, and for the two examples above FALSE.

    This mostly leverages HGVS parsing from HgvsVariant, but also does some bespoke matching on variant names.
    """
    transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs = None, None, None, False
    if variant_name is None:
        return transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs

    hgvs_variant = HgvsVariant(variant_name)

    if hgvs_variant.sequence_type in {SequenceType.GENOMIC, SequenceType.CODING}:
        # We are only interested in RefSeq mRNA transcripts for querying
        if hgvs_variant.sequence_identifier.startswith('NM'):
            transcript_id = hgvs_variant.sequence_identifier
        if hgvs_variant.has_valid_precise_span():
            coordinate_span = hgvs_variant.precise_span()
        # If the repeat unit is present, we can calculate its length directly
        if hgvs_variant.repeat_sequence:
            repeat_unit_length = len(hgvs_variant.repeat_sequence)
        return transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs

    # Is this a protein HGVS-like notation?
    if hgvs_variant.sequence_type == SequenceType.PROTEIN:
        is_protein_hgvs = True
        return transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs

    # Is this a human-readable description?
    match = re_description.search(variant_name)
    if match:
        repeat_unit_length = len(match.group('sequence'))
        return transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs

    logger.warning('ClinVar identifier did not match any of the regular expressions: {}'.format(variant_name))
    return transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs


class RepeatExpansionVariant:
    """Contains information specific to parsing measures representing repeat expansion variants."""

    def __init__(self, name, explicit_insertion_length):
        (transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs) = parse_repeat_identifier(name)
        self.transcript_id = transcript_id
        self.coordinate_span = coordinate_span if coordinate_span is not None else explicit_insertion_length
        self.repeat_unit_length = repeat_unit_length
        self.is_protein_hgvs = is_protein_hgvs
        self.name = name

    @property
    def repeat_type(self):
        """Based on all available information about a variant, determine its type. The resulting type can be:
            * trinucleotide_repeat_expansion, corresponding to SO:0002165
            * short_tandem_repeat_expansion, corresponding to SO:0002162
            * None (not able to determine)
        """
        repeat_type = None

        if self.is_protein_hgvs:
            # For protein HGVS notation, assume that repeat is a trinucleotide one, since it affects entire amino acids
            repeat_type = 'trinucleotide_repeat_expansion'
        else:
            # As a priority, use the repeat unit length determined directly from the HGVS-like base sequence
            # If not available, fall back to using and end coordinate difference
            repeat_unit_length = self.repeat_unit_length
            if repeat_unit_length is None:
                repeat_unit_length = self.coordinate_span
            # Determine repeat type based on repeat unit length
            if repeat_unit_length is not None:
                if repeat_unit_length % 3 == 0:
                    repeat_type = 'trinucleotide_repeat_expansion'
                else:
                    repeat_type = 'short_tandem_repeat_expansion'
        # Check if the HGVS-like name of the variant contains a simple deletion. In this case, it should not be
        # processed as a repeat *expansion* variant. The reason such records are present at this stage is that for
        # records without explicit allele sequences we cannot verify whether they definitely represent expansions.
        if self.name and (self.name.endswith('del') or self.name.endswith('del)')):
            repeat_type = None
        return repeat_type
