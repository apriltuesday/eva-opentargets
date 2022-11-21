import logging
import re

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io import ClinVarRecordMeasure
from eva_cttv_pipeline.clinvar_xml_io.clinvar_xml_io.hgvs_variant import HgvsVariant, SequenceType

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Human-readable description notation. This pattern is used for about a dozen variants. From it we only
# extract the repeat unit sequence.
# Examples: 'ATXN8, (CAG)n REPEAT EXPANSION' or 'TNRC6A, 5-BP INS, TTTCA(n) REPEAT EXPANSION'
re_description = re.compile(
    r'\(?'
    r'(?P<sequence>[{}]+)'.format(IUPACAmbiguousDNA.letters) +
    r'\)?\(?n\)?'
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
        if hgvs_variant.reference_sequence.startswith('NM'):
            transcript_id = hgvs_variant.reference_sequence
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


def repeat_type_from_length(length):
    if length:
        return 'trinucleotide_repeat_expansion' if length % 3 == 0 else 'short_tandem_repeat_expansion'
    return None


def infer_repeat_type_and_transcript_id(identifier):
    """
    Determine transcript ID and repeat type from a variant identifier (HGVS-like or other name).
    The resulting type can be:
        * trinucleotide_repeat_expansion, corresponding to SO:0002165
        * short_tandem_repeat_expansion, corresponding to SO:0002162
        * None (not able to determine)
    """
    (transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs) = parse_repeat_identifier(identifier)

    if is_protein_hgvs:
        # For protein HGVS notation, assume that repeat is a trinucleotide one, since it affects entire amino acids
        repeat_type = 'trinucleotide_repeat_expansion'
    else:
        # As a priority, use the repeat unit length determined directly from the HGVS-like base sequence
        # If not available, fall back to using start and end coordinate difference
        repeat_type = repeat_type_from_length(repeat_unit_length)
        if not repeat_type:
            repeat_type = repeat_type_from_length(coordinate_span)

    # Check if the HGVS-like name of the variant contains a simple deletion. In this case, it should not be
    # processed as a repeat *expansion* variant. The reason such records are present at this stage is that for
    # records without explicit allele sequences we cannot verify whether they definitely represent expansions.
    if identifier and (identifier.endswith('del') or identifier.endswith('del)')):
        repeat_type = None

    return repeat_type, transcript_id


def parse_all_identifiers(clinvar_measure: ClinVarRecordMeasure):
    """Attempt to parse all identifiers (HGVS expressions and names) in a measure to get repeat type and transcript id.
    Repeat type must be present if possible, transcript id is useful if present."""
    repeat_type, transcript_id = infer_repeat_type_and_transcript_id(clinvar_measure.get_variant_name_or_hgvs())
    if repeat_type:
        return repeat_type, transcript_id

    # TODO should we prioritise these in some way?
    for hgvs in clinvar_measure.current_hgvs:
        repeat_type, transcript_id = infer_repeat_type_and_transcript_id(hgvs.text)
        if repeat_type:
            return repeat_type, transcript_id
    for name in clinvar_measure.all_names:
        repeat_type, transcript_id = infer_repeat_type_and_transcript_id(name)
        if repeat_type:
            return repeat_type, transcript_id

    return None, None
