/*
 * Extract target ontology from mappings file header. Defaults to EFO if missing.
 */
process getTargetOntology {
    label 'short'
    label 'small'

    input:
    val mappingsFile

    output:
    env ONTOLOGY, emit: targetOntology

    script:
    """
    ONTOLOGY=\$(grep '^#ontology=' ${mappingsFile} | sed 's/#ontology=//g')
    ONTOLOGY=\${ONTOLOGY:-EFO}
    """
}
