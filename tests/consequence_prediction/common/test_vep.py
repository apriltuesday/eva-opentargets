from cmat.consequence_prediction.common.vep import extract_consequences, overall_most_severe_consequence, \
    most_severe_consequence_per_gene


def get_vep_results():
    return [
        {
            'input': '10 27169969 . C A',
            'transcript_consequences': [
                {  # most severe consequence
                    'gene_symbol': 'MASTL',
                    'gene_id': 'ENSG00000120539',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['missense_variant'],
                    'transcript_id': 'ENST00000342386'
                },
                {  # same consequence and gene, different transcript
                    'gene_symbol': 'MASTL',
                    'gene_id': 'ENSG00000120539',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['missense_variant'],
                    'transcript_id': 'ENST00000477034'
                },
                {  # missing gene symbol
                    'gene_id': 'ENSG00000120538',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['missense_variant'],
                    'transcript_id': 'ENST00000342386'
                },
                {  # less severe consequence on same gene
                    'gene_id': 'ENSG00000120538',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['stop_retained_variant'],
                    'transcript_id': 'ENST00000342386'
                },
                {  # less severe consequence on distinct gene
                    'gene_id': 'ENSG00000120537',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['stop_retained_variant'],
                    'transcript_id': 'ENST00000342385'
                },
                {  # non-overlapping consequences
                    'gene_symbol': 'ACBD5',
                    'gene_id': 'ENSG00000107897',
                    'biotype': 'protein_coding',
                    'distance': 5,
                    'consequence_terms': ['3_prime_UTR_variant'],
                    'transcript_id': 'ENST00000676731'
                },
                {  # other biotypes
                    'gene_symbol': 'MIR125A',
                    'gene_id': 'ENSG00000208008',
                    'biotype': 'miRNA',
                    'consequence_terms': ['missense_variant'],
                    'transcript_id': 'ENST00000385273'
                }
            ]
        }
    ]


def test_extract_consequences():
    """Verifies behaviour of extract_consequences."""
    vep_results = get_vep_results()
    results = extract_consequences(
        vep_results=vep_results,
        acceptable_biotypes=['protein_coding'],
        include_transcripts=False
    )
    assert results == {'10 27169969 . C A': [
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant'),
        ('10 27169969 . C A', 'ENSG00000120538', '', 'missense_variant'),
        ('10 27169969 . C A', 'ENSG00000120537', '', 'stop_retained_variant'),
    ]}


def test_extract_consequences_with_transcripts():
    """Verifies behaviour of extract_consequences with transcript flag set."""
    vep_results = get_vep_results()
    results = extract_consequences(
        vep_results=vep_results,
        acceptable_biotypes=['protein_coding'],
        include_transcripts=True
    )
    assert results == {'10 27169969 . C A': [
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant', 'ENST00000342386'),
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant', 'ENST00000477034'),
        ('10 27169969 . C A', 'ENSG00000120538', '', 'missense_variant', 'ENST00000342386'),
        ('10 27169969 . C A', 'ENSG00000120537', '', 'stop_retained_variant', 'ENST00000342385'),
    ]}


def test_overall_most_severe_consequence():
    vep_results = get_vep_results()
    variant_identifier = vep_results[0]['input']
    consequences = vep_results[0]['transcript_consequences']
    results = overall_most_severe_consequence(variant_identifier, consequences, include_transcripts=False)
    assert results == [
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant'),
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant'),
        ('10 27169969 . C A', 'ENSG00000120538', '', 'missense_variant'),
        ('10 27169969 . C A', 'ENSG00000208008', 'MIR125A', 'missense_variant')
    ]


def test_overall_most_severe_consequence_with_transcripts():
    vep_results = get_vep_results()
    variant_identifier = vep_results[0]['input']
    consequences = vep_results[0]['transcript_consequences']
    results = overall_most_severe_consequence(variant_identifier, consequences, include_transcripts=True)
    assert results == [
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant', 'ENST00000342386'),
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant', 'ENST00000477034'),
        ('10 27169969 . C A', 'ENSG00000120538', '', 'missense_variant', 'ENST00000342386'),
        ('10 27169969 . C A', 'ENSG00000208008', 'MIR125A', 'missense_variant', 'ENST00000385273')
    ]


def test_most_severe_consequence_per_gene():
    vep_results = get_vep_results()
    variant_identifier = vep_results[0]['input']
    consequences = vep_results[0]['transcript_consequences']
    results = most_severe_consequence_per_gene(variant_identifier, consequences, include_transcripts=False)
    assert results == [
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant'),
        ('10 27169969 . C A', 'ENSG00000120538', '', 'missense_variant'),
        ('10 27169969 . C A', 'ENSG00000120537', '', 'stop_retained_variant'),
        ('10 27169969 . C A', 'ENSG00000107897', 'ACBD5', '3_prime_UTR_variant'),
        ('10 27169969 . C A', 'ENSG00000208008', 'MIR125A', 'missense_variant')
    ]


def test_most_severe_consequence_per_gene_with_transcripts():
    vep_results = get_vep_results()
    variant_identifier = vep_results[0]['input']
    consequences = vep_results[0]['transcript_consequences']
    results = most_severe_consequence_per_gene(variant_identifier, consequences, include_transcripts=True)
    assert results == [
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant', 'ENST00000342386'),
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant', 'ENST00000477034'),
        ('10 27169969 . C A', 'ENSG00000120538', '', 'missense_variant', 'ENST00000342386'),
        ('10 27169969 . C A', 'ENSG00000120537', '', 'stop_retained_variant', 'ENST00000342385'),
        ('10 27169969 . C A', 'ENSG00000107897', 'ACBD5', '3_prime_UTR_variant', 'ENST00000676731'),
        ('10 27169969 . C A', 'ENSG00000208008', 'MIR125A', 'missense_variant', 'ENST00000385273')
    ]
