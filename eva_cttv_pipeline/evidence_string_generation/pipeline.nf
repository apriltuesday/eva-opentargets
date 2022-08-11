#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def helpMessage() {
    log.info"""
    
    Params:
        --batch_root     Directory for current batch
        --schema        Open Targets JSON schema version
        --clinvar       ClinVar XML file (optional, will download latest if omitted)
    """
}

params.help = null
params.batch_root = null
params.schema = null
params.clinvar = null

if (params.help) {
    exit 0, helpMessage()
}
if (!params.batch_root || !params.schema) {
    exit 1, helpMessage()
}
batchRoot = params.batch_root

workflow {
    if (params.clinvar != null) {
        clinvarXml = Channel.fromPath(params.clinvar)
    } else {
        clinvarXml = downloadClinvar()
    }
    downloadJsonSchema()

    runSnp(clinvarXml)
    generateEvidence(clinvarXml,
                     downloadJsonSchema.out.jsonSchema,
                     runSnp.out.consequenceMappingsSnp)
    checkDuplicates(generateEvidence.out.evidenceStrings)

    convertXrefs(clinvarXml)
}

process downloadClinvar {
    output:
    path "clinvar.xml.gz", emit: clinvarXml

    script:
    """
    wget -O clinvar.xml.gz \
        https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz
    """
}

process downloadJsonSchema {
    output:
    path "opentargets-${params.schema}.json", emit: jsonSchema

    script:
    """
    wget -O opentargets-${params.schema}.json \
        https://raw.githubusercontent.com/opentargets/json_schema/${params.schema}/opentargets.json
    """
}

process runSnp {
    input:
    path clinvarXml

    output:
    path "consequences_snp.tsv", emit: consequenceMappingsSnp

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
    // TODO logs for this
}

process generateEvidence {
    publishDir "${batchRoot}/logs",
        overwrite: true,
        mode: "copy",
        pattern: "*.log"

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
        --out . \
        > evidence_string_generation.log
    """
}

process checkDuplicates {
    input:
    path evidenceStrings

    script:
    """
    (jq --arg sep \$'\t' -jr \
        '.datatypeId,\$sep,.studyId,\$sep,.targetFromSourceId,\$sep,.variantId,\$sep,.variantFunctionalConsequenceId,\$sep,.diseaseFromSourceMappedId,\sep,.diseaseFromSource,"\n"' \
        ${evidenceStrings} | sort | uniq -d > duplicates.tsv) \
    || [[ -z duplicates.tsv ]]
    """
}

process convertXrefs {
    publishDir "${batchRoot}/logs",
        overwrite: true,
        mode: "copy",
        pattern: "*.log"

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
        --clinvar-xml    ${clinvarXml} \
        --zooma-feedback clinvar_xrefs.txt \
        > traits_to_zooma_format.log
    """
}