import re
from enum import Enum, auto


class StrEnum(str, Enum):
    """
    Extension of Enum class that uses lower-cased name as auto value and can be compared using string equality
    (e.g. StrEnum.VALUE == 'value')
    """
    def _generate_next_value_(name, start, count, last_values):
        return name.lower()


class SequenceType(StrEnum):
    """All possible sequence types for HGVS identifiers"""
    CODING = auto()
    GENOMIC = auto()
    NONCODING = auto()
    PROTEIN = auto()
    MITOCHONDRIAL = auto()
    CIRCULAR = auto()
    RNA = auto()


class VariantType(StrEnum):
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
    any_sequence_type = f'{sequence_identifier}([cgnpmor]).'
    any_known_range = f'([0-9]+)_([0-9]+)'

    def __init__(self, hgvs):
        self.hgvs = hgvs
        self.reference_sequence = None
        self.sequence_type = None
        self.variant_type = None

        self.start = None
        self.stop = None
        self.inner_start = None
        self.inner_stop = None
        self.outer_start = None
        self.outer_stop = None

        self._match_sequence_info()
        self._match_simple_range()

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
            # For simple ranges, no inners/outers so we set them equal
            self.inner_start = self.outer_start = self.start
            self.inner_stop = self.outer_stop = self.stop
            if m.group(5):
                var_type = m.group(5)
                if 'del' in var_type and 'delins' not in var_type:
                    self.variant_type = VariantType.DELETION
                elif 'dup' in var_type:
                    self.variant_type = VariantType.DUPLICATION
                elif 'ins' in var_type and 'delins' not in var_type:
                    self.variant_type = VariantType.INSERTION
