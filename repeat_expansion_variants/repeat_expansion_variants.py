#!/usr/bin/env python3
"""A script to extract repeat expansion variants from ClinVar TSV dump. For documentation, refer to README.md"""

import argparse
from io import StringIO
import itertools
import json
import logging
import math
import os
import pandas as pd
import re
import requests
import sys

from retry import retry

logging.basicConfig()
logger = logging.getLogger('repeat_expansion_variants')
logger.setLevel(logging.INFO)

# Regular expressions to parse ClinVar's HGVS-like notation
# Rest assured I know about the existence of a pyhgvs module... Unfortunately those are HGVS-*like*, not real HGVS.

# Common part for HGVS-like transcripts: transcript ID, version and delimiters
hgvs_like_transcript_part = (
    r'(?P<transcript_id>[A-Za-z0-9_]+)'  # Transcript accession                                   NM_001256054
    r'\.'                                # Delimiter: transcript accession/version                .
    r'[0-9]+'                            # Transcript version                                     2
    r'[A-Za-z0-9_.()]*'                  # Gene symbol in parentheses, optional                   (C9orf72)
    r':'                                 # Delimiter: transcript/variant info                     :
)

# Common part for start and end coordinates: base coordinate which _preceses_ the actual coordinate, and needs to be
# skipped
coordinate_base_part = (
    r'(?:'       # Non-capturing group for coordinate base part
    r'[-+]?'         # Coordinate base may start with a minus or plus
    r'[0-9]+'        # Then it contains a number
    r'(?=[-+])'      # It then must be followed with either a plus or minus (a positive lookahead)
    r')?'        # The coordinate base is always optional
)

# HGVS-like notation for transcript or genomic coordinates.
# Example: NM_001256054.2(C9orf72):c.-45+163_-45+180GGGGCC(2_25)
# This allows us to extract: transcript ID; coordinates (start always, end sometimes) and repeat unit (sometimes)
re_hgvs_like_transcript_or_genomic = re.compile(
    hgvs_like_transcript_part +
    r'[gc]'                              # Coordinate type, genomic or coding                     c
    r'\.'                                # Delimiter: coordinate type/coordinate                  .
    
    + coordinate_base_part +             # Start coordinate base, optional                        -45
    r'\*?'                              # Coordinate may start with an asterisk; ignore that
    r'(?P<start_coord>[+-]?[0-9]+)'      # Start coordinate                                       +163

    r'(?:'                               # The entire end coordinate part is optional:
    r'_'                                     # Delimiter: start coordinate/end coordinate         _
    + coordinate_base_part +                          # End coordinate base, optional                      -45
    r'\*?'                              # Coordinate may start with an asterisk; ignore that
    r'(?P<end_coord>[+-]?[0-9]*)'            # End coordinate                                     +180
    r')?'
    
    r'(?P<sequence>[ATGC]*)'             # Repeat unit, optional                                  GGGGCC
)

# HGVS-like notation for protein coordinates.
# Sometimes we only have a protein HGVS-like notation, not a coding one. (Why, ClinVar, why?)
# In this case we don't parse any coordinates or anything, and just assume that the repeat is a trinucleotide one
re_hgvs_like_protein = re.compile(hgvs_like_transcript_part + 'p\.')  # E. g. NP_002964.3:p.Gln166(>=33)

# Human-readable description notation. Luckily it's consistent.
# Example: ATXN8, (CAG)n REPEAT EXPANSION
re_description = re.compile(
    r'\('
    r'(?P<sequence>[ATGC]+)'
    r'\)n'
    r' REPEAT EXPANSION'
)

biomart_request_template = """http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "{key_column}" value = "{identifier_list}"/>
        <Attribute name = "{key_column}" />
        <Attribute name = "{column_to_query}" />
    </Dataset>
</Query>""".replace('\n', '')


@retry(tries=10, delay=5, backoff=1.2, jitter=(1, 3), logger=logger)
def query_biomart(key_column, column_to_query, identifier_list, df_key_column, df_column_to_query):
    """Query Ensembl BioMart and return information about specific columns."""
    identifier_list_string = ','.join([i for i in identifier_list])
    biomart_query = biomart_request_template.format(
        key_column=key_column, column_to_query=column_to_query, identifier_list=identifier_list_string)
    result = requests.get(biomart_query)
    # If there was an HTTP error, raise an exception. This will be caught by @retry
    result.raise_for_status()
    resulting_df = pd.read_table(StringIO(result.text), names=(df_key_column, df_column_to_query))
    # In case of multiple mappings, group them into the list
    resulting_df = resulting_df.groupby(df_key_column)[df_column_to_query]\
        .apply(list).reset_index(name=df_column_to_query)
    return resulting_df


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


def determine_repeat_unit(row):
    variant_name = str(row.Name)
    row['TranscriptID'], row['CoordinateSpan'], row['RepeatUnitLength'], row['IsProteinHGVS'] = None, None, None, False

    # Try to match HGVS-like transcript/genomic ID
    match = re_hgvs_like_transcript_or_genomic.search(variant_name)
    if match:
        transcript_id, start_coord, end_coord, sequence = \
            match.group('transcript_id'), match.group('start_coord'), match.group('end_coord'), match.group('sequence')
        if transcript_id.startswith('NM'):  # We are only interested in RefSeq mRNA transcripts for querying
            row['TranscriptID'] = match.group('transcript_id')
        if start_coord and end_coord:  # If both start *and* end coordinates are present, we can calculate the span
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
    repeat_type = None
    if row['IsProteinHGVS']:
        # For protein HGVS notation, always assume that repeat is trinucleotide
        repeat_type = 'trinucleotide_repeat_expansion'
    else:
        # Try to determine repeat unit length, prioritising the one obtained directly from HGVS sequence
        repeat_unit_length = row['RepeatUnitLength']
        # If no luck, fall back to start and end coordinate difference
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


def main():
    # Parse command line arguments
    args = parser.parse_args()

    # Load variants
    variants = pd.read_table(args.clinvar_summary_tsv)

    # Filter only NT expansion variants (in case we got anything else in the file)
    variants = variants[variants['Type'] == 'NT expansion']

    # Drop all columns except: name; RCV accessions; HGNC gene ID
    variants = variants[['Name', 'RCVaccession', 'GeneSymbol', 'HGNC_ID']]

    # Records may contain multiple RCVs per row, delimited by semicolon. Here we explode them into separate rows
    variants['RCVaccession'] = variants['RCVaccession'].str.split(';')
    variants = variants.explode('RCVaccession')

    # The same is true for having multiple *gene symbols* per record, they should also be split
    variants['GeneSymbol'] = variants['GeneSymbol'].str.split(';')
    variants = variants.explode('GeneSymbol')

    # Because the same (Name, RCV) pair can have coordinate in multiple builds, it can be repeated. Remove duplicates
    variants = variants.drop_duplicates()

    # Sort values by variant name
    variants = variants.sort_values(by=['Name'])

    # Assign "Repeat Unit" and "Repeat Type" columns
    variants = variants.apply(lambda row: determine_repeat_unit(row), axis=1)

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
            key_column=column_name_in_biomart,
            column_to_query='ensembl_gene_id',
            identifier_list=identifiers_to_query,
            df_key_column=column_name_in_dataframe,
            df_column_to_query='EnsemblGeneID',
        )

        # Step 3: make note where the annotations came from
        annotation_info['GeneAnnotationSource'] = column_name_in_dataframe

        # Step 4: combine the information we received with the *original* dataframe
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
        key_column='ensembl_gene_id',
        column_to_query='external_gene_name',
        identifier_list=sorted({str(i) for i in variants['EnsemblGeneID'] if str(i).startswith('ENSG')}),
        df_key_column='EnsemblGeneID',
        df_column_to_query='EnsemblGeneName'
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

    variants.to_csv(args.output_dataframe, sep='\t', index=False)

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
    consequences.to_csv(args.output_consequences, sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()
