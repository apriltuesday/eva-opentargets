#!/usr/bin/env python3
"""A script to extract repeat expansion variants from ClinVar TSV dump. For documentation refer to README.md"""

import argparse
from io import StringIO
import logging
import math
import re

import pandas as pd
import requests
from retry import retry

logging.basicConfig()
logger = logging.getLogger('repeat_expansion_variants')
logger.setLevel(logging.INFO)


# ############################################   ENSEMBL BIOMART QUERIES   ############################################

# Ensembl BioMart is an API which can be used to query Ensembl databases. Here we use it to map external identifiers
# (HGNC gene ID, gene symbol, RefSeq transcript ID) into Ensembl gene ID and name. The BioMart API is a bit quirky in
# that it uses XML to specify the request.

# The idea behind this template is that we query BioMart with a list of identifiers (`identifier_list`) from the
# `key_column`. For example, it can be a "hgnc_id" column, which contains a HGNC ID of a given gene. We then ask
# BioMart to return all mappings from that column to the `query_column`. For example, it can be the "ensembl_gene_id"
# column, which contains stable Ensembl gene ID.
biomart_request_template = """http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "{key_column}" value = "{identifier_list}"/>
        <Attribute name = "{key_column}" />
        <Attribute name = "{query_column}" />
    </Dataset>
</Query>""".replace('\n', '')


# Since we are using an external API call here, the @retry decorator will ensure that any sporadic network errors will
# be handled and the request will be retried.
@retry(tries=10, delay=5, backoff=1.2, jitter=(1, 3), logger=logger)
def query_biomart(key_column, query_column, identifier_list):
    """Query Ensembl BioMart with a list of identifiers (`identifier_list`) from one column (`key_column`) and return
    all mappings from those identifiers to another column (`query_column`) in form of a two-column Pandas dataframe.

    Args:
        key_column: A tuple of key column names in Ensembl and in the resulting dataframe, e. g. ('hgnc_id', 'HGNC_ID')
        query_column: A tuple of query column names, similar to `key_column`, e. g. ('ensembl_gene_id', 'EnsemblGeneID')
        identifier_list: List of identifiers to query, e. g. ['HGNC:10548', 'HGNC:10560']


    Returns:
        A Pandas dataframe with two columns. It will contain at most one row per input identifier. The query column will
        always contain a #list# to support the possibility of multiple mappings. In the example above, this will be:
               HGNC_ID      EnsemblGeneID
            0  HGNC:10548   [ENSG00000124788]
            1  HGNC:10560   [ENSG00000285258, ENSG00000163635]
    """
    biomart_key_column, df_key_column = key_column
    biomart_query_column, df_query_column = query_column
    # Construct BioMart query from the template (see explanation above)
    biomart_query = biomart_request_template.format(
        key_column=biomart_key_column,
        query_column=biomart_query_column,
        identifier_list=','.join([i for i in identifier_list])
    )
    result = requests.get(biomart_query)
    # If there was an HTTP error, raise an exception. This will be caught by @retry
    result.raise_for_status()
    resulting_df = pd.read_table(StringIO(result.text), names=(df_key_column, df_query_column))
    # Group all potential mappings into lists.
    resulting_df = resulting_df.groupby(df_key_column)[df_query_column].apply(list).reset_index(name=df_query_column)
    return resulting_df


# ##########################################   PARSING VARIANT IDENTIFIERS   ##########################################

# These regular expressions are used to parse ClinVar's identifiers for describing repeat expansion variants. The pyhgvs
# module cannot be applied here because not all expression used by ClinVar are actually valid HGVS, and that module
# imposes strict validation.

# Common part for all HGVS-like transcript definitions, e. g. 'NM_001256054.2(C9orf72):'
hgvs_like_transcript_part = (
    r'(?P<transcript_id>[A-Za-z0-9_]+)'   # Transcript accession                      NM_001256054
    r'\.'                                 # Delimiter, transcript accession/version   .
    r'[0-9]+'                             # Transcript version                        2
    r'[A-Za-z0-9_.()]#'                   # Gene symbol in parentheses, optional      (C9orf72)
    r':'                                  # Delimiter, transcript/variant info        :
)

# Common part for start and end coordinate pivots. Pivots for introns and other cases where a coordinate in a noncoding
# region of mRNA needs to be addressed relative to the coding regions. For example, c.87 means the 87th coding base of
# mRNA (not counting introns). In comparison, c.87+14 means that base number 87 is the last base of a particular exon,
# and the base addressed by the coordinate is 14 bases downstream of the pivot base. In this case, an intron repeat
# expansion variant might be addressed as c.87+14_c.87+17. In this case, we don't want the pivots (87), but only the
# actual coordinates (+14 and +17). To do that, we have a regular expression which captures the pivot part.
coordinate_pivot_part = (
    r'(?:'       # Non-capturing group for coordinate pivot part
    r'[-+]?'         # Coordinate pivot may start with a minus or plus
    r'[0-9]+'        # Then it contains a number
    r'(?=[-+])'      # It then must be followed with either a plus or minus.  This is a positive lookahead and is #not#
                     # part of the coordinate pivot, hence the (?=...) notation
    r')?'        # The coordinate pivot is always optional
)

# Pattern 1. HGVS-like notation for coding or genomic coordinates. This pattern is used for most variants. From it we
# can always extract transcript ID and start coord, and sometimes also end coord and repeat unit sequence.
# Example: 'NM_001256054.2(C9orf72):c.-45+63_-45+80GGGGCC(2_25)'
re_hgvs_like_transcript_or_genomic = re.compile(
    hgvs_like_transcript_part +       # Transcript definition                                NM_001256054.2(C9orf72):
    r'[gc]'                           # Coordinate type, genomic or coding                   c
    r'\.'                             # Delimiter, coordinate type/coordinate                .

    + coordinate_pivot_part +         # Start coordinate pivot, optional                     -45
    r'\#?'                            # Sometimes there is an asterisk in front of the
                                      # coordinate (special case, not important)
    r'(?P<start_coord>[+-]?[0-9]+)'   # Start coordinate                                     +63

    r'(?:'                            # Non-capruting group for end coordinate
    r'_'                                  # Delimiter: start coordinate/end coordinate       _
    + coordinate_pivot_part +             # End coordinate base, optional                    -45
    r'\#?'                                # Coordinate may start with an asterisk
    r'(?P<end_coord>[+-]?[0-9]#)'         # End coordinate                                   +80
    r')?'                             # The entire end coordinate part is optional
    
    r'(?P<sequence>[ATGC]#)'          # Repeat unit sequence, optional                       GGGGCC
)

# Pattern 2. HGVS-like notation for protein coordinates. This pattern is used for a few variants. From it we do not
# extract any useful information, but always assume that such a notation signified a trinucleotide repeat type.
# Example: NP_002964.3:p.Gln166(>=33)
re_hgvs_like_protein = re.compile(hgvs_like_transcript_part + r'p\.')

# Pattern 3. Human-readable description notation. This pattern is used for about a dozen variants. From it we only
# extract the repeat unit sequence. Example: 'ATXN8, (CAG)n REPEAT EXPANSION'
re_description = re.compile(
    r'\('
    r'(?P<sequence>[ATGC]+)'
    r'\)n'
    r' REPEAT EXPANSION'
)


# #############################################   MAIN PROCESSING LOGIC   #############################################

def parse_variant_identifier(row):
    """Parses variant identifier and extract certain characteristics into separate columns:
        * TranscriptID: NCBI RefSeq transcript ID, e. g. NM_000044
        * CoordinateSpan: the distance between start and end coordinates described in the HGVS-like notation.
              E. g., for 'NM_000044.4(AR):c.172_174CAG(7_34) (p.Gln66_Gln80del)', it will be 174 - 172 + 1 = 3
        * RepeatUnitLength: length of the sequence being repeated.
              E. g., for 'NC_000004.11:g.3076606GCA[27_35]', it will be len('GCA') = 3
        * IsProteinHGVS: this field simply reflects whether the variant was defined using a protein HGVS-like notation.
              E.g., for 'NP_002964.3:p.Gln166(>=33)' it will be TRUE, and for the two examples above FALSE.
    """
    # A variant "name", or identifier, can come from either of three patterns described in the section above.
    variant_name = str(row.Name)
    row['TranscriptID'], row['CoordinateSpan'], row['RepeatUnitLength'], row['IsProteinHGVS'] = None, None, None, False

    # Try to match HGVS-like transcript/genomic ID
    match = re_hgvs_like_transcript_or_genomic.search(variant_name)
    if match:
        transcript_id, start_coord, end_coord, sequence = \
            match.group('transcript_id'), match.group('start_coord'), match.group('end_coord'), match.group('sequence')
        if transcript_id.startswith('NM'):  # We are only interested in RefSeq mRNA transcripts for querying
            row['TranscriptID'] = match.group('transcript_id')
        if start_coord and end_coord:  # If both start #and# end coordinates are present, we can calculate the span
            row['CoordinateSpan'] = int(end_coord) - int(start_coord) + 1
        if sequence:  # If the repeat unit is present, we can calculate its length directly
            row['RepeatUnitLength'] = len(sequence)
        return row

    # Is this a protein HGVS-like notation?
    if re_hgvs_like_protein.search(variant_name):
        row['IsProteinHGVS'] = True
        return row

    # Is this a human-readable description?
    match = re_description.search(variant_name)
    if match:
        row['RepeatUnitLength'] = len(match.group('sequence'))
    return row


def determine_repeat_type(row):
    """Based on all available information about a variant, determine its type. The resulting type can be:
        * trinucleotide_repeat_expansion, corresponding to SO:0002165
        * short_tandem_repeat_expansion, corresponding to SO:0002162
        * None (not able to determine)
    """
    repeat_type = None
    if row['IsProteinHGVS']:
        # For protein HGVS notation, assume that repeat is a trinucleotide one, since it affects entire amino acids
        repeat_type = 'trinucleotide_repeat_expansion'
    else:
        # As a priority, use the repeat unit length determined directly from base sequence
        repeat_unit_length = row['RepeatUnitLength']
        # If not available luck, fall back to using and end coordinate difference
        if math.isnan(repeat_unit_length):
            repeat_unit_length = row['CoordinateSpan']
        # Determine repeat type based on repeat unit length
        if not math.isnan(repeat_unit_length):
            if repeat_unit_length % 3 == 0:
                repeat_type = 'trinucleotide_repeat_expansion'
            else:
                repeat_type = 'short_tandem_repeat_expansion'
    row['RepeatType'] = repeat_type
    return row


def main(clinvar_summary_tsv, output_consequences, output_dataframe):
    """Process data and generate output files.

    Args:
        clinvar_summary_tsv: filepath to the ClinVar variant summary file.
        output_consequences: filepath to the output file with variant consequences. The file uses a 6-column format
            compatible with the VEP mapping pipeline (see /vep-mapping-pipeline/README.md).
        output_dataframe: filepath to the output file with the full dataframe used in the analysis. This will contain
            all relevant columns and can be used for review or debugging purposes.
    """
    # Load variants
    variants = pd.read_table(clinvar_summary_tsv)

    # Filter only NT expansion variants (in case we got anything else in the file)
    variants = variants[variants['Type'] == 'NT expansion']

    # Drop all columns except: name; RCV accessions; HGNC gene ID
    variants = variants[['Name', 'RCVaccession', 'GeneSymbol', 'HGNC_ID']]

    # Records may contain multiple RCVs per row, delimited by semicolon. Here we explode them into separate rows
    variants['RCVaccession'] = variants['RCVaccession'].str.split(';')
    variants = variants.explode('RCVaccession')

    # The same is true for having multiple #gene symbols# per record, they should also be split
    variants['GeneSymbol'] = variants['GeneSymbol'].str.split(';')
    variants = variants.explode('GeneSymbol')

    # Because the same (Name, RCV) pair can have coordinate in multiple builds, it can be repeated. Remove duplicates
    variants = variants.drop_duplicates()

    # Sort values by variant name
    variants = variants.sort_values(by=['Name'])

    # Assign "Repeat Unit" and "Repeat Type" columns
    variants = variants.apply(lambda row: parse_variant_identifier(row), axis=1)

    # OK, now we need to match each (RCV, variant) pair to Ensembl gene ID and name.
    # And it's so much more complicated than you could imagine.

    gene_annotation_sources = (
        # Dataframe column  Biomart column        Filtering function
        ('HGNC_ID',        'hgnc_id',             lambda i: i.startswith('HGNC:')),
        ('GeneSymbol', 'external_gene_name', lambda i: i != '-'),
        ('TranscriptID',   'refseq_mrna',         lambda i: i is not None),
    )
    variants_original = variants.copy(deep=True)
    for column_name_in_dataframe, column_name_in_biomart, filtering_function in gene_annotation_sources:
        # Step 1: get all identifiers we want to query BioMart with
        identifiers_to_query = sorted({
            i for i in variants[column_name_in_dataframe]
            if filtering_function(i)
        })

        # Step 2: query BioMart for Ensembl Gene IDs
        annotation_info = query_biomart(
            key_column=(column_name_in_biomart, column_name_in_dataframe),
            query_column=('ensembl_gene_id', 'EnsemblGeneID'),
            identifier_list=identifiers_to_query,
        )

        # Step 3: make note where the annotations came from
        annotation_info['GeneAnnotationSource'] = column_name_in_dataframe

        # Step 4: combine the information we received with the #original# dataframe
        annotation_df = pd.merge(variants_original, annotation_info, on=column_name_in_dataframe, how='left')

        # Step 5: update main dataframe with the new values
        variants = variants \
            .set_index([column_name_in_dataframe]) \
            .combine_first(annotation_df.set_index([column_name_in_dataframe]))

    # Reset index to default
    variants.reset_index(inplace=True)

    # Some records are being annotated to multiple Ensembl genes. For example, HGNC:10560 is being resolved to
    # ENSG00000285258 and ENSG00000163635. We need to explode dataframe by that column.
    variants = variants.explode('EnsemblGeneID')

    # Fetch Ensembl gene name based on Ensembl gene ID
    annotation_info = query_biomart(
        key_column=('ensembl_gene_id', 'EnsemblGeneID'),
        query_column=('external_gene_name', 'EnsemblGeneName'),
        identifier_list=sorted({str(i) for i in variants['EnsemblGeneID'] if str(i).startswith('ENSG')}),
    )
    variants = pd.merge(variants, annotation_info, on='EnsemblGeneID', how='left')

    # Check that there are no multiple gene name mappings for any given gene ID
    assert variants['EnsemblGeneName'].str.len().dropna().max() == 1, 'Multiple gene ID → gene name mappings found!'

    # Convert the one-item list into a plain column
    variants = variants.explode('EnsemblGeneName')

    # Populate variant type
    variants = variants.apply(lambda row: determine_repeat_type(row), axis=1)

    # Based on all information, mark records as either complete or incomplete
    variants['RecordIsComplete'] = (variants['EnsemblGeneID'].notnull()
                                    & variants['EnsemblGeneName'].notnull()
                                    & variants['RepeatType'].notnull())

    # Rearrange order of dataframe columns
    variants = variants[
        ['Name', 'RCVaccession', 'HGNC_ID', 'GeneSymbol',
         'RepeatUnitLength', 'CoordinateSpan', 'IsProteinHGVS', 'TranscriptID',
         'EnsemblGeneID', 'EnsemblGeneName', 'GeneAnnotationSource',
         'RepeatType', 'RecordIsComplete']
    ]

    variants.to_csv(output_dataframe, sep='\t', index=False)

    # Generate consequences table
    consequences = variants[variants['RecordIsComplete']] \
        .groupby(['RCVaccession', 'EnsemblGeneID', 'EnsemblGeneName'])['RepeatType'] \
        .apply(set).reset_index(name='RepeatType')

    # Check that for every (RCV, gene) pair there is only one consequence type
    assert consequences['RepeatType'].str.len().dropna().max() == 1, 'Multiple (RCV, gene) → variant type mappings!'

    # Get rid of sets
    consequences['RepeatType'] = consequences['RepeatType'].apply(list)
    consequences = consequences.explode('RepeatType')

    # Form six-column file compatible with the consequence mapping pipeline
    # RCV000005966    1       ENSG00000156475 PPP2R2B trinucleotide_repeat_expansion  0
    consequences['PlaceholderOnes'] = 1
    consequences['PlaceholderZeroes'] = 0
    consequences = consequences[['RCVaccession', 'PlaceholderOnes', 'EnsemblGeneID', 'EnsemblGeneName', 'RepeatType',
                                 'PlaceholderZeroes']]
    consequences.sort_values(by=['RepeatType', 'RCVaccession', 'EnsemblGeneID'], inplace=True)

    # Write the consequences table
    consequences.to_csv(output_consequences, sep='\t', index=False, header=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--clinvar-summary-tsv', required=True,
        help='ClinVar summary TSV dump file (variant_summary.txt.gz)'
    )
    parser.add_argument(
        '--output-consequences', required=True,
        help='File to output functional consequences to. Format is compatible with the main VEP mapping pipeline.'
    )
    parser.add_argument(
        '--output-dataframe', required=True,
        help='File to output full dataframe for subsequent analysis and debugging.'
    )
    args = parser.parse_args()
    main(args.clinvar_summary_tsv, args.output_consequences, args.output_dataframe)
