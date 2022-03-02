"""Contains utilities and classes to parse the ClinVar XML and convert the records into internal representation via
ClinVarDataset, ClinVarRecord, ClinVarTrait, ClinVarRecordMeasure and ClinVarRecordMeasureHGVS classes."""

import gzip
import logging
import re
import xml.etree.ElementTree as ElementTree
from functools import cached_property

from .clinvar_identifier_parsing import parse_variant_identifier

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def iterate_rcv_from_xml(clinvar_xml):
    """Iterates through the gzipped ClinVar XML and yields complete <ReferenceClinVarAssertion> records."""
    with gzip.open(clinvar_xml, 'rt') as fh:
        for event, elem in ElementTree.iterparse(fh):
            # Wait until we have built a complete ClinVarSet element
            if elem.tag != 'ClinVarSet':
                continue

            # Go to a ReferenceClinVarAssertion element. This corresponds to a single RCV record, the main unit of
            # ClinVar. There should only be one such record per ClinVarSet.
            rcv = find_mandatory_unique_element(elem, 'ReferenceClinVarAssertion')

            # Return the complete record and then remove the processed element from the tree to save memory
            yield rcv
            elem.clear()


def find_elements(node, xpath, allow_zero=True, allow_multiple=True):
    """Attempt to find child elements in a node by xpath. Raise exceptions if conditions are violated. Return a
    (possibly empty) list of elements."""
    all_elements = node.findall(xpath)
    if (len(all_elements) == 0 and not allow_zero) or (len(all_elements) > 1 and not allow_multiple):
        raise AssertionError(f'Found {len(all_elements)} instances of {xpath} in {node}, which is not allowed')
    return all_elements


def find_mandatory_unique_element(node, xpath):
    """Attempt to find a child element by xpath which must have exactly one occurrence, otherwise throw an exception."""
    return find_elements(node, xpath, allow_zero=False, allow_multiple=False)[0]


def find_optional_unique_element(node, xpath):
    """Attempt to find a child element by xpath which must have either 0 or 1 occurrence. Return the element, or None if
    not found."""
    elements = find_elements(node, xpath, allow_zero=True, allow_multiple=False)
    return elements[0] if elements else None


def pubmed_refs_to_urls(pubmed_refs):
    """Convert a list of PubMed identifiers to complete URLs."""
    return [f'http://europepmc.org/abstract/MED/{r}' for r in pubmed_refs]


class ClinVarDataset:
    """Iterate through records in ClinVar XML dump and convert them into internal ClinVarRecord representation."""
    def __init__(self, clinvar_xml):
        self.clinvar_xml = clinvar_xml

    def __iter__(self):
        for rcv in iterate_rcv_from_xml(self.clinvar_xml):
            yield ClinVarRecord(rcv)


class ClinVarRecord:
    """Instances of this class hold data on individual ClinVar records. See also:
    * /data-exploration/clinvar-variant-types/README.md for the in-depth explanation of ClinVar data model;
    * Issue https://github.com/EBIvariation/eva-opentargets/issues/127 for the most recent discussions on changing
      support of different ClinVar record types."""

    # A score for the review status of the assigned clinical significance ranges from 0 to 4 and corresponds to the
    # number of "gold stars" displayed on ClinVar website. See details here:
    # https://www.ncbi.nlm.nih.gov/clinvar/docs/details/#review_status
    score_map = {
        "no assertion provided": 0,
        'no assertion criteria provided': 0,
        'criteria provided, single submitter': 1,
        'criteria provided, conflicting interpretations': 1,
        'criteria provided, multiple submitters, no conflicts': 2,
        'reviewed by expert panel': 3,
        'practice guideline': 4,
    }

    # Some allele origin terms in ClinVar are essentially conveying lack of information and are thus not useful.
    NONSPECIFIC_ALLELE_ORIGINS = {'unknown', 'not provided', 'not applicable', 'tested-inconclusive', 'not-reported'}

    def __init__(self, rcv):
        """Initialise a ClinVar record object from an RCV XML record."""
        self.rcv = rcv

        # Add a list of traits
        self.trait_set = []
        for trait in find_elements(self.rcv, './TraitSet/Trait'):
            self.trait_set.append(ClinVarTrait(trait, self))

        # We are currently only processing MeasureSets of type Variant which are included directly in the RCV record.
        # Some other options (currently not supported) are:
        # * MeasureSet of types "Haplotype", "Phase unknown", or "Distinct chromosomes"
        # * GenotypeSet, which contains an assertion about a group of variants from different chromosome copies, with
        #   the type of be either a "CompoundHeterozygote" or a "Diplotype"
        variant_measure = find_optional_unique_element(self.rcv, './MeasureSet[@Type="Variant"]/Measure')
        if not variant_measure:
            self.measure = None
        else:
            self.measure = ClinVarRecordMeasure(variant_measure, self)

    def __str__(self):
        return f'ClinVarRecord object with accession {self.accession}'

    @property
    def accession(self):
        return find_mandatory_unique_element(self.rcv, './ClinVarAccession').attrib['Acc']

    @property
    def date(self):
        """This tracks the latest update date, counting even minor technical updates."""
        return self.rcv.attrib['DateLastUpdated']

    @property
    def last_evaluated_date(self):
        """This tracks the latest (re)evaluation date for the clinical interpretation.
        See https://github.com/opentargets/platform/issues/1161#issuecomment-683938510 for details."""
        # The DateLastEvaluated attribute is not always present. In this case, this property will be None.
        return find_mandatory_unique_element(self.rcv, './ClinicalSignificance').attrib.get('DateLastEvaluated')

    @property
    def review_status(self):
        """Return a review status text for the assigned clinical significance. See score_map above for the list of
        possible values."""
        review_status = find_mandatory_unique_element(self.rcv, './ClinicalSignificance/ReviewStatus').text
        assert review_status in self.score_map, f'Unknown review status {review_status} in RCV {self.accession}'
        return review_status

    @property
    def score(self):
        """Return a score (star rating) for the assigned clinical significance. See score_map above."""
        return self.score_map[self.review_status]

    @property
    def mode_of_inheritance(self):
        """Return a (possibly empty) list of modes of inheritance for a given ClinVar record."""
        return sorted({
            elem.text for elem in find_elements(self.rcv, './AttributeSet/Attribute[@Type="ModeOfInheritance"]')
        })

    @property
    def trait_set_type(self):
        return find_mandatory_unique_element(self.rcv, './TraitSet').attrib['Type']

    @property
    def traits(self):
        """Returns a list of traits associated with the ClinVar record, in the form of Trait objects."""
        return self.trait_set

    @property
    def traits_with_valid_names(self):
        """Returns a list of traits which have at least one valid (potentially resolvable) name."""
        return [trait for trait in self.trait_set if trait.preferred_or_other_valid_name]

    @property
    def evidence_support_pubmed_refs(self):
        """The references of this type represent evidence support for this specific variant being observed in this
        specific disease. These are the references displayed on the ClinVar website in the "Assertion and evidence
        details" section at the bottom of the page."""
        return [int(elem.text)
                for elem in find_elements(self.rcv, './ObservedIn/ObservedData/Citation/ID[@Source="PubMed"]')]

    @property
    def clinical_significance_raw(self):
        """The original clinical significance string as stored in ClinVar. Example: 'Benign/Likely benign'."""
        return find_mandatory_unique_element(self.rcv, './ClinicalSignificance/Description').text

    @property
    def clinical_significance_list(self):
        """The normalised list of all clinical significance values. The original value is (1) split into multiple values
        by two delimiters: ('/', ', '), (2) converted into lowercase and (3) sorted lexicographically. Example:
        'Benign/Likely benign, risk_factor' → ['benign', 'likely benign', 'risk factor']. See /data-exploration/
        clinvar-variant-types/README.md for further explanation."""
        return sorted(re.split('/|, ', self.clinical_significance_raw.lower().replace('_', ' ')))

    @property
    def allele_origins(self):
        return {elem.text for elem in find_elements(self.rcv, './ObservedIn/Sample/Origin')}

    @property
    def valid_allele_origins(self):
        """Returns all valid allele origins, i.e. ones that are not in the list of nonspecific terms."""
        return {origin for origin in self.allele_origins if origin.lower() not in self.NONSPECIFIC_ALLELE_ORIGINS}


class ClinVarTrait:
    """Represents a single ClinVar trait (usually a disease), with the corresponding database and Pubmed
    cross-references."""

    # Some trait records in ClinVar contain names which are non-specific and cannot possibly be resolved to any
    # meaningful EFO term.
    NONSPECIFIC_TRAITS = {
        '', 'allhighlypenetrant', 'disease', 'none provided', 'not provided', 'not specified',
        'reclassified - variant of unknown significance', 'see cases', 'variant of unknown significance'
    }

    def __init__(self, trait_xml, clinvar_record):
        self.trait_xml = trait_xml
        self.clinvar_record = clinvar_record

    def __str__(self):
        return f'ClinVarTrait object with name {self.preferred_or_other_valid_name} from ClinVar record ' \
               f'{self.clinvar_record.accession}'

    @property
    def identifier(self):
        return self.trait_xml.attrib['ID'].strip()

    @property
    def all_names(self):
        """Returns a lexicographically sorted list of all trait names, including the preferred one (if any)."""
        return sorted(name.text for name in find_elements(self.trait_xml, './Name/ElementValue'))

    @property
    def all_valid_names(self):
        """Returns a lexicographically sorted list of all valid trait names. A valid name is defined as something which
        is not contained in the list of nonspecific traits, which cannot possibly be resolved to a valid EFO mapping."""
        return [name for name in self.all_names if name.lower() not in self.NONSPECIFIC_TRAITS]

    @property
    def preferred_name(self):
        """Returns a single preferred name, as indicated in the ClinVar record."""
        name = find_optional_unique_element(self.trait_xml, './Name/ElementValue[@Type="Preferred"]')
        return None if name is None else name.text

    @property
    def preferred_or_other_valid_name(self):
        """Returns a consistent valid name for a trait, if one is present."""
        if self.preferred_name and self.preferred_name.lower() not in self.NONSPECIFIC_TRAITS:
            return self.preferred_name
        elif self.all_valid_names:
            return self.all_valid_names[0]
        else:
            return None

    @property
    def pubmed_refs(self):
        """Trait-specific PubMed references, contained inside a Trait entity. These are usually reviews or practice
        guidelines related to a disease or a group of diseases."""
        return [int(elem.text) for elem in find_elements(self.trait_xml, './Citation/ID[@Source="PubMed"]')]

    @property
    def xrefs(self):
        return [(elem.attrib['DB'], elem.attrib['ID'].strip(), elem.attrib.get('Status', 'current').lower())
                for elem in find_elements(self.trait_xml, './XRef')]

    @property
    def medgen_id(self):
        """Attempts to resolve a single MedGen ID for a trait. If not present, returns None. If multiple are present,
        returns the first one lexicographically."""
        medgen_ids = []
        for db, id_, status in self.xrefs:
            if db == 'MedGen' and status == 'current':
                medgen_ids.append(id_)
        medgen_ids.sort()
        if len(medgen_ids) == 0:
            logger.warning(f'No MedGen ID for {self}')
        elif len(medgen_ids) == 1:
            return medgen_ids[0]
        else:
            logger.warning(f'Multiple MedGen IDs for {self}: {medgen_ids}')
            return medgen_ids[0]


class ClinVarRecordMeasure:
    """This class represents individual ClinVar record "measures". Measures are essentially isolated variants, which can
    be combined into either MeasureSets (include one or more Measures) or GenotypeSets. For a detailed description of
    ClinVar data model, see data-exploration/clinvar-variant-types/."""

    # For ClinVar Microsatellite events with complete coordinates, require the event to be at least this number of bases
    # long in order for it to be considered a repeat expansion event. Smaller events will be processed as regular
    # insertions. The current value was chosen as a reasonable threshold to separate thousands of very small insertions
    # which are technically microsatellite expansion events but are not long enough to be considered clinically
    # significant repeat expansion variants.
    REPEAT_EXPANSION_THRESHOLD = 12

    # Microsatellite variant types.
    MS_DELETION = 'deletion'
    MS_SHORT_EXPANSION = 'short_expansion'
    MS_REPEAT_EXPANSION = 'repeat_expansion'
    MS_NO_COMPLETE_COORDS = 'no_complete_coords'

    def __init__(self, measure_xml, clinvar_record):
        self.measure_xml = measure_xml
        self.clinvar_record = clinvar_record

    @property
    def all_names(self):
        """Returns a lexicographically sorted list of all measure names, including the preferred one (if any)."""
        return sorted(name.text for name in find_elements(self.measure_xml, './Name/ElementValue'))

    @property
    def preferred_name(self):
        """Returns a single preferred measure name, as indicated in the ClinVar record."""
        name = find_optional_unique_element(self.measure_xml, './Name/ElementValue[@Type="Preferred"]')
        return None if name is None else name.text

    @property
    def preferred_or_other_name(self):
        """Returns a consistent name for a measure, if one is present."""
        if self.preferred_name:
            return self.preferred_name
        elif self.all_names:
            return self.all_names[0]
        else:
            return None

    @property
    def preferred_gene_symbols(self):
        return [elem.text for elem in find_elements(
            self.measure_xml, './MeasureRelationship/Symbol/ElementValue[@Type="Preferred"]')]

    @property
    def hgnc_ids(self):
        return [elem.attrib['ID'] for elem in find_elements(self.measure_xml, './MeasureRelationship/XRef[@DB="HGNC"]')]

    @property
    def rs_id(self):
        """Returns dbSNP rsID found for the record measure."""
        rs_ids = ['rs' + elem.attrib['ID'] for elem in find_elements(self.measure_xml, './XRef[@DB="dbSNP"]')]
        if len(rs_ids) == 0:
            return None
        elif len(rs_ids) == 1:
            return rs_ids[0]
        else:
            logger.warning(f'Found multiple RS IDs for {self.clinvar_record}, this is not yet supported')
            return None

    @property
    def nsv_id(self):
        nsv_ids = [elem.attrib['ID']
                   for elem in find_elements(self.measure_xml, './XRef[@DB="dbVar"]')
                   if elem.attrib['ID'].startswith('nsv')]
        if len(nsv_ids) == 0:
            return None
        elif len(nsv_ids) == 1:
            return nsv_ids[0]
        else:
            logger.warning(f'Found multiple NSV IDs for {self.clinvar_record}, this is not yet supported')
            return None

    @property
    def _hgvs_elems(self):
        return [
            elem
            for elem in find_elements(self.measure_xml, './AttributeSet/Attribute')
            if elem.attrib['Type'].startswith('HGVS')
        ]

    @property
    def hgvs(self):
        return [elem.text for elem in self._hgvs_elems]

    @property
    def current_hgvs(self):
        return [elem.text for elem in self._hgvs_elems if 'previous' not in elem.attrib['Type'].lower()]

    @property
    def toplevel_refseq_hgvs(self):
        refseq_elems = [elem for elem in self._hgvs_elems if elem.attrib['Type'].lower() == 'hgvs, genomic, top level']
        return refseq_elems[0].text if refseq_elems else None

    @property
    def variant_type(self):
        return self.measure_xml.attrib['Type']

    @property
    def explicit_insertion_length(self):
        if self.vcf_alt and self.vcf_ref:
            return len(self.vcf_alt) - len(self.vcf_ref)
        return None

    @property
    def microsatellite_category(self):
        if self.variant_type == 'Microsatellite':
            if self.has_complete_coordinates:
                if self.explicit_insertion_length < 0:
                    return self.MS_DELETION
                elif self.explicit_insertion_length < self.REPEAT_EXPANSION_THRESHOLD:
                    return self.MS_SHORT_EXPANSION
                else:
                    return self.MS_REPEAT_EXPANSION
            else:
                return self.MS_NO_COMPLETE_COORDS
        else:
            return None

    @property
    def is_repeat_expansion_variant(self):
        return self.microsatellite_category in (self.MS_REPEAT_EXPANSION, self.MS_NO_COMPLETE_COORDS)

    @property
    def pubmed_refs(self):
        """Variant-specific PubMed references, contained inside a Measure entity. These are usually large reviews which
        focus on genetics of specific types of variants or genomic regions."""
        return [int(elem.text) for elem in find_elements(self.measure_xml, './Citation/ID[@Source="PubMed"]')]

    @property
    def chr(self):
        return self.sequence_location_helper('Chr')

    @property
    def vcf_pos(self):
        return self.sequence_location_helper('positionVCF')

    @property
    def vcf_ref(self):
        return self.sequence_location_helper('referenceAlleleVCF')

    @property
    def vcf_alt(self):
        return self.sequence_location_helper('alternateAlleleVCF')

    @property
    def has_complete_coordinates(self):
        return bool(self.chr and self.vcf_pos and self.vcf_ref and self.vcf_alt)

    @property
    def vcf_full_coords(self):
        """Returns complete variant coordinates in CHROM_POS_REF_ALT format, if present, otherwise None."""
        if self.has_complete_coordinates:
            return '_'.join([self.chr, self.vcf_pos, self.vcf_ref, self.vcf_alt])

    def sequence_location_helper(self, attr):
        if self.variant_type == 'Translocation':
            # TODO: Translocations have multiple locations and are not supported.
            # TODO: https://github.com/EBIvariation/eva-opentargets/issues/171
            return None
        sequence_locations = find_elements(self.measure_xml, './SequenceLocation[@Assembly="GRCh38"]')
        if len(sequence_locations) != 1:
            # TODO: Support variants with multiple locations (for example, chrX/chrY).
            # TODO: https://github.com/EBIvariation/eva-opentargets/issues/172
            return None
        return sequence_locations[0].attrib.get(attr)

    def get_variant_name_or_hgvs(self):
        if self.preferred_or_other_name:
            return self.preferred_or_other_name
        if self.toplevel_refseq_hgvs:
            return self.toplevel_refseq_hgvs
        return None

    @cached_property
    def hgvs_properties(self):
        return ClinVarRecordMeasureHGVS(self.get_variant_name_or_hgvs(), self.explicit_insertion_length)


class ClinVarRecordMeasureHGVS:
    # TODO integrate this with HgvsVariant

    def __init__(self, name, explicit_insertion_length):
        (transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs) = parse_variant_identifier(name)
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
        # Check if the HGVS-like name of the variant contains a simple deletion. In this case, it should not be processed
        # as a repeat *expansion* variant. The reason such records are present at this stage is that for records without
        # explicit allele sequences we cannot verify whether they definitely represent expansions.
        if self.name and (self.name.endswith('del') or self.name.endswith('del)')):
            repeat_type = None
        return repeat_type
