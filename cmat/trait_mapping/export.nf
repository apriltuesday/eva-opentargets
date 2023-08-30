#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def helpMessage() {
    log.info"""
    Export manual curation table.

    Params:
        --curation_root     Directory for current batch
        --input_csv         Input csv file
        --mappings          Current mappings file (optional, will use a default path if omitted)
        --with_feedback     Whether to generate EFO/Zooma feedback and final symlinking
    """
}

params.help = null
params.curation_root = null
params.input_csv = null
params.mappings = "\${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv"
params.with_feedback = false

if (params.help) {
    exit 0, helpMessage()
}
if (!params.curation_root or !params.input_csv) {
    exit 1, helpMessage()
}
curationRoot = params.curation_root


/*
 * Main workflow.
 */
workflow {
    exportTable(params.input_csv)
    combineManualAndAutomated(exportTable.out.finishedMappings)
    mergeWithLatestMappings(combineManualAndAutomated.out.newMappings)
    checkDuplicates(mergeWithLatestMappings.out.newMappings)
    if (params.with_feedback) {
        createEfoTable(exportTable.out.importTerms)
        generateZoomaFeedback(mergeWithLatestMappings.out.newMappings)
        updateLinks(checkDuplicates.out.duplicatesOk, generateZoomaFeedback.out.zoomaFeedback)
    }
}

/*
 * Extract the relevant columns from the input CSV.
 */
process exportTable {
    publishDir "${curationRoot}",
        overwrite: true,
        mode: "copy",
        pattern: "curator_comments.tsv"

    input:
    path inputCsv

    output:
    path "finished_mappings_curation.tsv", emit: finishedMappings
    path "terms_for_efo_import.txt", emit: importTerms
    path "curator_comments.tsv", emit: curatorComments

    script:
    """
    \${PYTHON_BIN} \${CODE_ROOT}/bin/trait_mapping/export_curation_table.py \
        -i ${inputCsv} \
        -d finished_mappings_curation.tsv \
        -m terms_for_efo_import.txt \
        -c curator_comments.tsv
    """
}

/*
 * Concatenate finished automated and manual mappings into a single file.
 */
process combineManualAndAutomated {
    input:
    path finishedMappings

    output:
    path "trait_names_to_ontology_mappings.tsv", emit: newMappings

    script:
    """
    cat ${curationRoot}/automated_trait_mappings.tsv ${finishedMappings} \
        | sort -u > trait_names_to_ontology_mappings.tsv
    """
}

/*
 * Add all mappings from the database which are *not* present in the results of the current curation iteration (automated
 * + manually curated). This is done in order to never lose mappings, even if they are not present in ClinVar during the
 * latest curation iteration.
 */
process mergeWithLatestMappings {
    publishDir "${curationRoot}",
        overwrite: true,
        mode: "copy",
        pattern: "*.tsv"

    input:
    path newMappings

    output:
    path newMappings, emit: newMappings

    script:
    """
    # The first file operand is the list of mappings in the current database; and the second is the list of trait names
    # which are only present in the existing database and not in the new mappings.
    export LC_ALL=C
    join -j 1 -t \$'\t' \
        <(sort -t \$'\t' -k 1,1 ${params.mappings}) \
        <(comm -23 <(cut -d \$'\t' -f 1 ${params.mappings} | sort -u) <(cut -d \$'\t' -f 1 ${newMappings} | sort -u)) \
    >> ${newMappings}
    """
}

/*
 * Prepare the table for EFO import.
 */
process createEfoTable {
    publishDir "${curationRoot}",
        overwrite: true,
        mode: "copy",
        pattern: "*.tsv"

    input:
    path importTerms

    output:
    path "efo_import_table.tsv", emit: efoImportTable

    script:
    """
    \${PYTHON_BIN} \${CODE_ROOT}/bin/trait_mapping/create_efo_table.py \
        -i ${importTerms} \
        -o efo_import_table.tsv
    """
}

/*
 * Generate ZOOMA feedback.
 */
process generateZoomaFeedback {
    publishDir "${curationRoot}",
        overwrite: true,
        mode: "copy",
        pattern: "*.txt"

    input:
    path newMappings

    output:
    path "eva_clinvar.txt", emit: zoomaFeedback

    script:
    """
    echo -e 'STUDY\tBIOENTITY\tPROPERTY_TYPE\tPROPERTY_VALUE\tSEMANTIC_TAG\tANNOTATOR\tANNOTATION_DATE' \
        > eva_clinvar.txt
    tail -n+2 ${newMappings} \
        | cut -f-2 \
        | sort -t\$'\t' -k1,1 \
        | awk -F\$'\t' -vDATE="\$(date +'%y/%m/%d %H:%M')" '{print "\t\tdisease\t" \$1 "\t" \$2 "\teva\t" DATE}' \
    >> eva_clinvar.txt
    """
}

/*
 * Check there are no complete duplicates in the final mappings file.
 */
process checkDuplicates {
    input:
    path newMappings

    output:
    val true, emit: duplicatesOk  // ensure we don't do the final linking if this check fails

    script:
    """
    sort ${newMappings} | uniq -c | awk '\$1 > 1' > duplicates.tsv
    [[ ! -s duplicates.tsv ]]
    """
}

/*
 * Update the symbolic links pointing to the location of the most recent curation result and ZOOMA feedback dataset.
 */
process updateLinks {
    input:
    val duplicatesOk
    path zoomaFeedback

    script:
    """
    ln -s -f ${curationRoot}/trait_names_to_ontology_mappings.tsv \${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv
    ln -s -f ${curationRoot}/eva_clinvar.txt \${BATCH_ROOT_BASE}/manual_curation/eva_clinvar.txt
    ln -s -f ${curationRoot}/curator_comments.tsv \${BATCH_ROOT_BASE}/manual_curation/latest_comments.tsv
    """
}
