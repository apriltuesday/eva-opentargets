import re
from enum import Enum, auto

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


class SequenceType(Enum):
    """All possible sequence types for HGVS identifiers"""
    CODING = auto()
    GENOMIC = auto()
    NONCODING = auto()
    PROTEIN = auto()
    MITOCHONDRIAL = auto()
    CIRCULAR = auto()
    RNA = auto()


class VariantType(Enum):
    """Most possible variant types for HGVS identifiers"""
    SUBSTITUTION = auto()
    DELETION = auto()
    DUPLICATION = auto()
    INSERTION = auto()
    INVERSION = auto()
    DELINS = auto()
    OTHER = auto()  # not intended to be exhaustive


class HgvsVariant:
    """Representation of variant as derived from HGVS identifier"""

    # Re-usable pieces of regular expressions
    sequence_identifier = (
        r'^([a-zA-Z][a-zA-Z0-9_.]+)'  # Sequence accession, e.g. NM_001256054.2
        r'(?:\([a-zA-Z0-9_.]+\))?'    # Optional gene symbol, e.g. (C9orf72)
        r':'                          # Delimiter, transcript/variant info
    )
    any_sequence_type = f'{sequence_identifier}([cgnpmor])\.'
    any_known_range = f'([0-9]+)_([0-9]+)'

    def __init__(self, hgvs):
        self.hgvs = hgvs
        self.reference_sequence = None
        self.sequence_type = None
        self.variant_type = None

        self.start = None
        self.stop = None

        self._match_sequence_info()
        self._match_simple_range()
        self._match_repeat_with_coordinate_pivots()

    def has_valid_precise_span(self):
        return self.start and self.stop and self.stop > self.start

    def precise_span(self):
        if self.has_valid_precise_span():
            return self.stop - self.start + 1
        return None

    def _match_sequence_info(self):
        regex = re.compile(self.any_sequence_type)
        m = regex.match(self.hgvs)
        if m:
            self.reference_sequence = m.group(1)
            seq_type = m.group(2)
            if seq_type == 'c':
                self.sequence_type = SequenceType.CODING
            elif seq_type == 'g':
                self.sequence_type = SequenceType.GENOMIC
            elif seq_type == 'n':
                self.sequence_type = SequenceType.NONCODING
            elif seq_type == 'p':
                self.sequence_type = SequenceType.PROTEIN
            elif seq_type == 'm':
                self.sequence_type = SequenceType.MITOCHONDRIAL
            elif seq_type == 'o':
                self.sequence_type = SequenceType.CIRCULAR
            elif seq_type == 'r':
                self.sequence_type = SequenceType.RNA

    def _match_simple_range(self):
        regex = re.compile(f'{self.any_sequence_type}{self.any_known_range}([a-zA-Z0-9]*)$')
        m = regex.match(self.hgvs)
        if m and m.group(3) and m.group(4):
            self.start = int(m.group(3))
            self.stop = int(m.group(4))
            if m.group(5):
                var_type = m.group(5)
                if 'del' in var_type and 'delins' not in var_type:
                    self.variant_type = VariantType.DELETION
                elif 'dup' in var_type:
                    self.variant_type = VariantType.DUPLICATION
                elif 'ins' in var_type and 'delins' not in var_type:
                    self.variant_type = VariantType.INSERTION

    def _match_repeat_with_coordinate_pivots(self):
        # Common part for start and end coordinate pivots. Pivots for introns and other cases where a coordinate in a
        # noncoding region of mRNA needs to be addressed relative to the coding regions. For example, c.87 means the
        # 87th coding base of mRNA (not counting introns). In comparison, c.87+14 means that base number 87 is the last
        # base of a particular exon, and the base addressed by the coordinate is 14 bases downstream of the pivot base.
        # In this case, an intron repeat expansion variant might be addressed as c.87+14_c.87+17. In this case, we don't
        # want the pivots (87), but only the actual coordinates (+14 and +17). To do that, we have a regular expression
        # which captures the pivot part.
        coordinate_pivot_part = (
            r'(?:'       # Non-capturing group for coordinate pivot part
            r'[-+]?'     # Coordinate pivot may start with a minus or plus
            r'[0-9]+'    # Then it contains a number
            r'(?=[-+])'  # It then must be followed with either a plus or minus
            r')?'        # The coordinate pivot is always optional
        )
        # HGVS-like notation for coding or genomic coordinates. This pattern is used for most variants. From it we
        # can always extract transcript ID and start coord, and sometimes also end coord and repeat unit sequence.
        # Example: 'NM_001256054.2(C9orf72):c.-45+63_-45+80GGGGCC(2_25)'
        re_hgvs_like_transcript_or_genomic = re.compile(
            self.any_sequence_type
            + coordinate_pivot_part +        # Start coordinate pivot, optional                     -45
            r'\*?'                           # Sometimes there is an asterisk in front of the
                                             # coordinate (special case, not important)
            r'(?P<start_coord>[+-]?[0-9]+)'  # Start coordinate                                     +63

            r'(?:'                           # Non-capturing group for end coordinate
            r'_'                             # Delimiter: start coordinate/end coordinate       _
            + coordinate_pivot_part +        # End coordinate pivot, optional                   -45
            r'\*?'                           # Coordinate may start with an asterisk
            r'(?P<end_coord>[+-]?[0-9]+)'    # End coordinate                                   +80
            r')?'                            # The entire end coordinate part is optional

            r'(?P<sequence>[{}]*)'.format(   # Repeat unit sequence, optional                       GGGGCC or
                IUPACAmbiguousDNA.letters)   # IUPAC ambiguity codes are supported                  CGN
        )
        m = re_hgvs_like_transcript_or_genomic.match(self.hgvs)
        if m:
            # If we don't already have a simple span, we'll take the pivot-based ones
            # TODO I think this is wrong actually, in particular if the pivot is different this yields incorrect spans
            if not self.has_valid_precise_span() and m.group('start_coord') and m.group('end_coord'):
                self.start = int(m.group('start_coord'))
                self.stop = int(m.group('end_coord'))
            # TODO other variant types might be interested in this sequence as well...
            #  (e.g. U43746.1:n.9877-68_9877-65delTTAC)
            self.repeat_sequence = m.group('sequence')
