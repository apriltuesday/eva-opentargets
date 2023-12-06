#!/usr/bin/env python3
"""A pipeline to extract repeat expansion variants from ClinVar XML dump. For documentation refer to README.md"""

import logging
from collections import Counter

import numpy as np
import pandas as pd

from cmat import clinvar_xml_io
from cmat.clinvar_xml_io.repeat_variant import parse_all_identifiers, repeat_type_from_length
import cmat.consequence_prediction.common.biomart as biomart

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


STANDARD_CHROMOSOME_NAMES = {str(c) for c in range(1, 23)} | {'X', 'Y', 'M', 'MT'}


def none_to_nan(value):
    """Converts arguments which are None to np.nan, for consistency inside a Pandas dataframe."""
    return np.nan if value is None else value


def load_clinvar_data(clinvar_xml):
    """Load ClinVar data, preprocess, and return it as a Pandas dataframe."""
    # Iterate through ClinVar XML records
    variant_data = []  # To populate the return dataframe (see columns below)
    stats = Counter()
    for i, clinvar_record in enumerate(clinvar_xml_io.ClinVarDataset(clinvar_xml)):
        if i and i % 100000 == 0:
            total_repeat_expansion_variants = stats[clinvar_xml_io.ClinVarRecordMeasure.MS_REPEAT_EXPANSION] + \
                                              stats[clinvar_xml_io.ClinVarRecordMeasure.MS_NO_COMPLETE_COORDS]
            logger.info(f'Processed {i} records, collected {total_repeat_expansion_variants} repeat expansion variant '
                        f'candidates')

        # Skip a record if it does not contain variant information
        if not clinvar_record.measure:
            continue
        measure = clinvar_record.measure

        # Repeat expansion events come in two forms: with explicit coordinates and allele sequences (CHROM/POS/REF/ALT),
        # or without them. In the first case we can compute the explicit variant length as len(ALT) - len(REF). In the
        # second case, which is more rare but still important, we have to resort to parsing HGVS-like variant names.
        if measure.microsatellite_category:
            stats[measure.microsatellite_category] += 1
        # Skip the record if it's a deletion or a short insertion
        if not measure.is_repeat_expansion_variant:
            continue

        # Extract gene symbol(s). Here and below, dashes are sometimes assigned to be compatible with the variant
        # summary format which was used previously.
        gene_symbols = measure.preferred_gene_symbols
        if not gene_symbols:
            gene_symbols = ['-']

        # Extract HGNC ID
        hgnc_ids = measure.hgnc_ids
        hgnc_id = hgnc_ids[0] if len(hgnc_ids) == 1 and len(gene_symbols) == 1 else '-'

        repeat_type, transcript_id = parse_all_identifiers(measure)
        # If no identifier yields a repeat type, try to infer from elsewhere in the measure
        if not repeat_type and measure.explicit_insertion_length:
            repeat_type = repeat_type_from_length(measure.explicit_insertion_length)

        # Append data strings
        for gene_symbol in gene_symbols:
            variant_data.append([
                measure.preferred_or_other_name,
                clinvar_record.accession,
                gene_symbol,
                hgnc_id,
                none_to_nan(transcript_id),
                none_to_nan(repeat_type)
            ])
    total_repeat_expansion_variants = stats[clinvar_xml_io.ClinVarRecordMeasure.MS_REPEAT_EXPANSION] + \
                                      stats[clinvar_xml_io.ClinVarRecordMeasure.MS_NO_COMPLETE_COORDS]
    logger.info(f'Done. A total of {i} records, {total_repeat_expansion_variants} repeat expansion variant candidates')

    variants = pd.DataFrame(variant_data, columns=('Name',
                                                   'RCVaccession',
                                                   'GeneSymbol',
                                                   'HGNC_ID',
                                                   'TranscriptID',
                                                   'RepeatType'))

    # Since the same record can have coordinates in multiple builds, it can be repeated. Remove duplicates
    variants = variants.drop_duplicates()
    # Sort values by variant name
    return variants.sort_values(by=['Name']), stats


def annotate_ensembl_gene_info(variants, include_transcripts):
    """Annotate the `variants` dataframe with information about Ensembl gene ID and name"""

    # Ensembl gene ID can be determined using three ways, listed in the order of decreasing priority. Having multiple
    # ways is necessary because no single method works on all ClinVar variants.
    gene_annotation_sources = (
        # Dataframe column    Biomart column   Filtering function
        ('HGNC_ID',                 'hgnc_id', lambda i: i.startswith('HGNC:')),
        ('GeneSymbol',   'external_gene_name', lambda i: i != '-'),
        ('TranscriptID',        'refseq_mrna', lambda i: pd.notnull(i)),
    )
    query_columns = [('ensembl_gene_id', 'EnsemblGeneID')]
    if include_transcripts:
        query_columns.append(('ensembl_transcript_id', 'EnsemblTranscriptID'))

    # Make a copy of the variants dataframe to query. Variants will be removed from this copy as Ensembl gene IDs are
    # found, ensuring that we don't make unnecessary queries.
    variants_to_query = variants

    # Empty dataframe to collect annotated rows
    annotated_variants = pd.DataFrame()

    for column_name_in_dataframe, column_name_in_biomart, filtering_function in gene_annotation_sources:
        # Get all identifiers we want to query BioMart with
        identifiers_to_query = sorted({
            i for i in variants_to_query[column_name_in_dataframe]
            if filtering_function(i)
        })
        # Query BioMart for Ensembl Gene IDs
        annotation_info = biomart.query_biomart(
            key_column=(column_name_in_biomart, column_name_in_dataframe),
            query_columns=query_columns,
            identifier_list=identifiers_to_query,
        )

        # Make note where the annotations came from
        annotation_info['GeneAnnotationSource'] = column_name_in_dataframe
        # Combine the information we received with the *original* dataframe using an outer join.
        total_merge = pd.merge(variants, annotation_info, on=column_name_in_dataframe, how='outer', indicator=True)

        # Variants annotated in this round are those present in both tables in the join, add these to the
        # cumulative results.
        annotated_variants = pd.concat((annotated_variants, total_merge[total_merge['_merge'] == 'both']))
        annotated_variants.drop('_merge', axis=1, inplace=True)

        # Variants that still need to be queried are those present only in the left table in the join.
        variants_to_query = total_merge[total_merge['_merge'] == 'left_only']
        variants_to_query.drop('_merge', axis=1, inplace=True)
        if len(variants_to_query) == 0:
            break

    # Based on the Ensembl gene ID, annotate (1) gene name and (2) which chromosome it is on
    gene_query_columns = (
        ('external_gene_name', 'EnsemblGeneName'),
        ('chromosome_name', 'EnsemblChromosomeName'),
    )
    annotation_info = biomart.query_biomart(
        key_column=('ensembl_gene_id', 'EnsemblGeneID'),
        query_columns=gene_query_columns,
        identifier_list=sorted({str(i) for i in annotated_variants['EnsemblGeneID'] if str(i).startswith('ENSG')}),
    )
    annotated_variants = pd.merge(annotated_variants, annotation_info, on='EnsemblGeneID', how='left')
    # Check that there are no multiple mappings for any given ID
    for _, column_name_in_dataframe in gene_query_columns:
        assert_uniqueness(annotated_variants, ['EnsemblGeneID'], column_name_in_dataframe,
                          'Found multiple gene ID → gene attribute mappings!')

    return annotated_variants


def determine_complete(row, include_transcripts):
    """Depending on the information, determine whether the record is complete, i.e., whether it has all necessary
        fields to be output for the final "consequences" table."""
    row['RecordIsComplete'] = (
        pd.notnull(row['EnsemblGeneID']) and
        pd.notnull(row['EnsemblGeneName']) and
        pd.notnull(row['RepeatType']) and
        (pd.notnull(row['EnsemblTranscriptID']) if include_transcripts else True) and
        row['EnsemblChromosomeName'] in STANDARD_CHROMOSOME_NAMES
    )
    return row


def generate_consequences_file(consequences, output_consequences):
    """Output final table."""
    if consequences.empty:
        logger.info('There are no records ready for output')
        return
    # Write the consequences table. This is used by the main annotation pipeline.
    consequences.to_csv(output_consequences, sep='\t', index=False, header=False)
    # Output statistics
    logger.info(f'Generated {len(consequences)} consequences in total:')
    logger.info(f'  {sum(consequences.RepeatType == "trinucleotide_repeat_expansion")} trinucleotide repeat expansion')
    logger.info(f'  {sum(consequences.RepeatType == "short_tandem_repeat_expansion")} short tandem repeat expansion')


def assert_uniqueness(df, uniqueness_columns, target_column, error_msg):
    """
    Check uniqueness of values in target_column with respect to tuples in uniqueness_columns.

    Args:
        df: dataframe to check
        uniqueness_columns: iterable of column names
        target_column: name of column that should be unique
        error_msg: message to report if uniqueness check fails
    Returns:
        result_df: input df containing only uniqueness_columns and target_column
    """
    result_df = df.groupby(uniqueness_columns, group_keys=False)[target_column].apply(set).reset_index(name=target_column)
    if result_df.empty:
        return result_df
    assert result_df[target_column].apply(len).dropna().max() == 1, error_msg
    # Get rid of sets
    result_df[target_column] = result_df[target_column].apply(list)
    result_df = result_df.explode(target_column)
    return result_df


def extract_consequences(variants, include_transcripts):
    """Generate consequences table"""
    # Check that for every (RCV, gene) pair or (RCV, gene, transcript) triple there is only one consequence type
    unique_repeat_type_columns = ['RCVaccession', 'EnsemblGeneID', 'EnsemblGeneName']
    if include_transcripts:
        unique_repeat_type_columns.append('EnsemblTranscriptID')
    consequences = assert_uniqueness(variants[variants['RecordIsComplete']], unique_repeat_type_columns, 'RepeatType',
                                     'Multiple (RCV, gene) → variant type mappings!')
    if consequences.empty:
        return consequences
    # Form a four-column file compatible with the consequence mapping pipeline, for example:
    # RCV000005966    ENSG00000156475    PPP2R2B    trinucleotide_repeat_expansion
    output_columns = ['RCVaccession', 'EnsemblGeneID', 'EnsemblGeneName', 'RepeatType']
    if include_transcripts:
        output_columns.append('EnsemblTranscriptID')
    consequences = consequences[output_columns]
    consequences.sort_values(by=['RepeatType', 'RCVaccession', 'EnsemblGeneID'], inplace=True)
    # Check that there are no empty cells in the final consequences table
    assert consequences.isnull().to_numpy().sum() == 0
    return consequences


def generate_all_variants_file(output_dataframe, variants):
    # Rearrange order of dataframe columns
    variants = variants[
        ['Name', 'RCVaccession', 'GeneSymbol', 'HGNC_ID', 'TranscriptID',
         'EnsemblGeneID', 'EnsemblGeneName', 'EnsemblChromosomeName', 'GeneAnnotationSource',
         'RepeatType', 'RecordIsComplete']
    ]
    # Write the full dataframe. This is used for debugging and investigation purposes.
    variants.sort_values(by=['Name', 'RCVaccession', 'GeneSymbol'])
    variants.to_csv(output_dataframe, sep='\t', index=False)


def main(clinvar_xml, include_transcripts, output_consequences=None, output_dataframe=None):
    """Process data and generate output files.

    Args:
        clinvar_xml: filepath to the ClinVar XML file.
        include_transcripts:
        output_consequences: filepath to the output file with variant consequences. The file uses a 6-column format
            compatible with the VEP mapping pipeline (see /consequence_prediction/README.md).
        output_dataframe: filepath to the output file with the full dataframe used in the analysis. This will contain
            all relevant columns and can be used for review or debugging purposes."""

    logger.info('Load and preprocess variant data')
    variants, s = load_clinvar_data(clinvar_xml)

    # Output ClinVar record statistics
    logger.info(f'''
        Microsatellite records: {sum(s.values())}
            With complete coordinates: {s[clinvar_xml_io.ClinVarRecordMeasure.MS_DELETION] +
                                        s[clinvar_xml_io.ClinVarRecordMeasure.MS_SHORT_EXPANSION] +
                                        s[clinvar_xml_io.ClinVarRecordMeasure.MS_REPEAT_EXPANSION]}
                Deletions: {s[clinvar_xml_io.ClinVarRecordMeasure.MS_DELETION]}
                Short insertions: {s[clinvar_xml_io.ClinVarRecordMeasure.MS_SHORT_EXPANSION]}
                Repeat expansions: {s[clinvar_xml_io.ClinVarRecordMeasure.MS_REPEAT_EXPANSION]}
            No complete coordinates: {s[clinvar_xml_io.ClinVarRecordMeasure.MS_NO_COMPLETE_COORDS]}
    '''.replace('\n' + ' ' * 8, '\n'))

    if variants.empty:
        logger.info('No variants to process')
        return None

    logger.info('Match each record to Ensembl gene ID and name')
    variants = annotate_ensembl_gene_info(variants, include_transcripts)

    logger.info('Determine variant type and whether the record is complete')
    variants = variants.apply(lambda row: determine_complete(row, include_transcripts), axis=1)

    logger.info('Postprocess data and output the two final tables')
    if output_dataframe is not None:
        generate_all_variants_file(output_dataframe, variants)
    consequences = extract_consequences(variants, include_transcripts)
    if output_consequences is not None:
        generate_consequences_file(consequences, output_consequences)

    logger.info('Completed successfully')
    return consequences
