#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def helpMessage() {
    log.info"""
    Generate ClinVar evidence strings for Open Targets, or annotated ClinVar XML.
    
    Params:
        --output_dir    Directory for output
        --schema        Open Targets JSON schema version (optional, will output XML if omitted)
        --clinvar       ClinVar XML file (optional, will download latest if omitted)
    """
}

params.help = null
params.output_dir = null
params.schema = null
params.clinvar = null

if (params.help) {
    exit 0, helpMessage()
}
if (!params.output_dir) {
    exit 1, helpMessage()
}
batchRoot = params.output_dir


/*
 * Main workflow.
 */
workflow {
    if (params.clinvar != null) {
        clinvarXml = Channel.fromPath(params.clinvar)
    } else {
        clinvarXml = downloadClinvar()
    }

    // Functional consequences
    runSnpIndel(clinvarXml)
    runRepeat(clinvarXml)
    runStructural(clinvarXml)
    combineConsequences(runSnpIndel.out.consequencesSnp,
                        runRepeat.out.consequencesRepeat,
                        runStructural.out.consequencesStructural)

    if (params.schema != null) {
        // Open Targets evidence string output
        downloadJsonSchema()
        generateEvidence(clinvarXml,
                         downloadJsonSchema.out.jsonSchema,
                         combineConsequences.out.consequencesCombined)

        checkDuplicates(generateEvidence.out.evidenceStrings)
        convertXrefs(clinvarXml)

    } else {
        // Annotated ClinVar XML output
        generateAnnotatedXml(clinvarXml, combineConsequences.out.consequencesCombined)
    }
}

/*
 * Download ClinVar data, using the most recent XML dump.
 */
process downloadClinvar {
    output:
    path "clinvar.xml.gz", emit: clinvarXml

    script:
    """
    wget -O clinvar.xml.gz \
        https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz
    """
}

/*
 * Download the Open Targets JSON schema.
 */
process downloadJsonSchema {
    output:
    path "opentargets-${params.schema}.json", emit: jsonSchema

    script:
    """
    wget -O opentargets-${params.schema}.json \
        https://raw.githubusercontent.com/opentargets/json_schema/${params.schema}/schemas/disease_target_evidence.json
    """
}

/*
 * Run simple variants (SNPs and other variants with complete coordinates) through VEP and map them
 * to genes and functional consequences.
 */
process runSnpIndel {
    clusterOptions "-o ${batchRoot}/logs/consequence_snp.out \
                    -e ${batchRoot}/logs/consequence_snp.err"

    publishDir "${batchRoot}/gene_mapping",
        overwrite: true,
        mode: "copy",
        pattern: "*.tsv"

    input:
    path clinvarXml

    output:
    path "consequences_snp.tsv", emit: consequencesSnp

    script:
    """
    \${PYTHON_BIN} "\${CODE_ROOT}/consequence_prediction/extract_variants_for_vep.py" --clinvar-xml ${clinvarXml} \
    | sort -u \
    | parallel \
        --halt now,fail=1    `# If any job fails, kill the remaining ones immediately and report failure` \
        --pipe               `# Input is read from STDIN and split by chunks`                             \
        -j 20                `# Number of concurrent workers`                                             \
        -N 200               `# Number of records (lines) per worker`                                     \
        --tmpdir .           `# Store temporary files in the current directory to avoid /tmp overflow`    \
        \${PYTHON_BIN} "\${CODE_ROOT}/consequence_prediction/vep_mapping_pipeline/consequence_mapping.py" \
    | sort -u > consequences_snp.tsv
    """
}

/*
 * Extract repeat expansion variants from ClinVar and map them to genes.
 */
process runRepeat {
   clusterOptions "-o ${batchRoot}/logs/consequence_repeat.out \
                   -e ${batchRoot}/logs/consequence_repeat.err"

   publishDir "${batchRoot}/gene_mapping",
       overwrite: true,
       mode: "copy",
       pattern: "*.tsv"

   input:
   path clinvarXml

   output:
   path "consequences_repeat.tsv", emit: consequencesRepeat

   script:
   """
   \${PYTHON_BIN} \${CODE_ROOT}/consequence_prediction/run_repeat_expansion_variants.py \
        --clinvar-xml ${clinvarXml} \
        --output-consequences consequences_repeat.tsv

    # create an empty file if nothing generated
    [[ -f consequences_repeat.tsv ]] || touch consequences_repeat.tsv
   """
}

/*
 * Run consequence and gene mapping for structural variants (i.e. no complete coordinates and not
 * known repeat expansions).
 */
process runStructural {
   clusterOptions "-o ${batchRoot}/logs/consequence_structural.out \
                   -e ${batchRoot}/logs/consequence_structural.err"

   publishDir "${batchRoot}/gene_mapping",
       overwrite: true,
       mode: "copy",
       pattern: "*.tsv"

   input:
   path clinvarXml

   output:
   path "consequences_structural.tsv", emit: consequencesStructural

   script:
   """
   \${PYTHON_BIN} \${CODE_ROOT}/consequence_prediction/run_structural_variants.py \
        --clinvar-xml ${clinvarXml} \
        --output-consequences consequences_structural.tsv

    # create an empty file if nothing generated
    [[ -f consequences_structural.tsv ]] || touch consequences_structural.tsv
   """
}

/*
 * Unite results of consequence mapping.
 */
process combineConsequences {
    input:
    path consequencesSnp
    path consequencesRepeat
    path consequencesStructural

    output:
    path "consequences_combined.tsv", emit: consequencesCombined

    script:
    """
    cat ${consequencesRepeat} ${consequencesSnp} ${consequencesStructural} > consequences_combined.tsv
    """
}

/*
 * Generate annotated ClinVar XML
 */
process generateAnnotatedXml {
    clusterOptions "-o ${batchRoot}/logs/annotated_xml_generation.out \
                    -e ${batchRoot}/logs/annotated_xml_generation.err"

    publishDir "${batchRoot}",
        overwrite: true,
        mode: "copy",
        pattern: "*.xml.gz"

    input:
    path clinvarXml
    path consequenceMappings

    output:
    path "annotated_clinvar.xml.gz"

    script:
    """
    \${PYTHON_BIN} \${CODE_ROOT}/bin/generate_annotated_xml.py \
        --clinvar-xml ${clinvarXml} \
        --efo-mapping \${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv \
        --gene-mapping ${consequenceMappings} \
        --output-xml annotated_clinvar.xml.gz
    """
}

/*
 * Generate the evidence strings for submission to Open Targets.
 */
process generateEvidence {
    clusterOptions "-o ${batchRoot}/logs/evidence_string_generation.out \
                    -e ${batchRoot}/logs/evidence_string_generation.err"

    publishDir "${batchRoot}/evidence_strings",
        overwrite: true,
        mode: "copy",
        pattern: "*.json"

    input:
    path clinvarXml
    path jsonSchema
    path consequenceMappings

    output:
    path "evidence_strings.json", emit: evidenceStrings

    script:
    """
    \${PYTHON_BIN} \${CODE_ROOT}/bin/evidence_string_generation.py \
        --clinvar-xml ${clinvarXml} \
        --efo-mapping \${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv \
        --gene-mapping ${consequenceMappings} \
        --ot-schema ${jsonSchema} \
        --out .
    """
}

/*
 * Check that the generated evidence strings do not contain any duplicated evidence strings.
 */
process checkDuplicates {
    input:
    path evidenceStrings

    script:
    """
    jq --arg sep \$'\t' -jr \
        '.datatypeId,\$sep,.studyId,\$sep,.targetFromSourceId,\$sep,.variantId,\$sep,.variantFunctionalConsequenceId,\$sep,.diseaseFromSourceMappedId,\$sep,.diseaseFromSource,"\n"' \
        ${evidenceStrings} | sort | uniq -d > duplicates.tsv
    [[ ! -s duplicates.tsv ]]
    """
}

/*
 * Convert MedGen and OMIM cross-references into ZOOMA format.
 */
process convertXrefs {
    clusterOptions "-o ${batchRoot}/logs/traits_to_zooma_format.out \
                    -e ${batchRoot}/logs/traits_to_zooma_format.err"

    publishDir "${batchRoot}/clinvar",
        overwrite: true,
        mode: "copy",
        pattern: "*.txt"

    input:
    path clinvarXml

    output:
    path "clinvar_xrefs.txt", emit: clinvarXrefs

    """
    \${PYTHON_BIN} \${CODE_ROOT}/bin/traits_to_zooma_format.py \
        --clinvar-xml ${clinvarXml} \
        --zooma-feedback clinvar_xrefs.txt
    """
}
