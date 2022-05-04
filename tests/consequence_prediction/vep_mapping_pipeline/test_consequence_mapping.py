from consequence_prediction.vep_mapping_pipeline.consequence_mapping import extract_consequences


def test_extract_consequences():
    """Verifies behaviour of extract_consequences."""
    vep_results = [
        {
            'input': '10 27169969 . C A',
            'transcript_consequences': [
                {
                    'gene_symbol': 'MASTL',
                    'gene_id': 'ENSG00000120539',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['missense_variant']
                },
                {  # missing gene symbol is okay
                    'gene_id': 'ENSG00000120538',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['missense_variant']
                },
                {  # less severe consequences should be filtered out
                    'gene_symbol': 'ACBD5',
                    'gene_id': 'ENSG00000107897',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['3_prime_UTR_variant']
                },
                {  # other biotypes should be filtered out
                    'gene_symbol': 'MIR125A',
                    'gene_id': 'ENSG00000208008',
                    'biotype': 'miRNA',
                    'consequence_terms': ['missense_variant']
                }
            ]
        }
    ]
    results = extract_consequences(
        vep_results=vep_results,
        acceptable_biotypes=['protein_coding'],
    )
    assert results == {'10 27169969 . C A': [
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant'),
        ('10 27169969 . C A', 'ENSG00000120538', '', 'missense_variant'),
    ]}
