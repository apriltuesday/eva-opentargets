import gzip
import logging
import xml.etree.ElementTree as ElementTree

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def parse_header_attributes(clinvar_xml):
    """Parses out attributes to the root-level ReleaseSet element and returns them as a dict."""
    attrib = None
    with gzip.open(clinvar_xml, 'rt') as fh:
        for event, elem in ElementTree.iterparse(fh, events=['start']):
            if elem.tag == 'ReleaseSet':
                attrib = elem.attrib
                break
            elem.clear()
    if not attrib:
        return {}
    # Resolve xsi:noNamespaceSchemaLocation="...", which is parsed strangely by ElementTree
    updated_attrib = {}
    for attr, val in attrib.items():
        if attr == '{http://www.w3.org/2001/XMLSchema-instance}noNamespaceSchemaLocation':
            updated_attrib['xmlns:xsi'] = 'http://www.w3.org/2001/XMLSchema-instance'
            updated_attrib['xsi:noNamespaceSchemaLocation'] = val
        else:
            updated_attrib[attr] = val
    return updated_attrib


def iterate_rcv_from_xml(clinvar_xml):
    """Iterates through the gzipped ClinVar XML and yields complete <ReferenceClinVarAssertion> records."""
    for cvs in iterate_cvs_from_xml(clinvar_xml):
        # Go to a ReferenceClinVarAssertion element. This corresponds to a single RCV record, the main unit of
        # ClinVar. There should only be one such record per ClinVarSet.
        rcv = find_mandatory_unique_element(cvs, 'ReferenceClinVarAssertion')
        yield rcv


def iterate_cvs_from_xml(clinvar_xml):
    """Iterates through the gzipped ClinVar XML and yields complete <ClinVarSet> elements."""
    with gzip.open(clinvar_xml, 'rt') as fh:
        for event, elem in ElementTree.iterparse(fh):
            # Wait until we have built a complete ClinVarSet element
            if elem.tag != 'ClinVarSet':
                continue
            # Return the complete record and then remove the processed element from the tree to save memory
            yield elem
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
