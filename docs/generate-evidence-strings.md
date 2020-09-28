# Evidence string generation protocol
_Issue template: https://www.ebi.ac.uk/panda/jira/browse/EVA-1912_

## 1. Set up the environment
First, [set up the common environment.](environment.md)

Next, set up the protocol-specific environment. The variables in it are specific to each Open Targets release. They are either announced the e-mail which they send a few weeks before the data submission deadline, or can be derived from the information in it:
* Year and month of the upcoming Open Targets release (`OT_RELEASE`). For example, if you're processing data for “20.02” release, this variable will be set to `2020-02`.
* Open Targets JSON schema version (`OT_SCHEMA_VERSION`)
* Open Targets validator package version (`OT_VALIDATOR_VERSION`)

You will also need to determine which ClinVar release to use as source (`CLINVAR_RELEASE`). Each Open Targets release is synchronised with a certain Ensembl release version, which is also announced in the aforementioned e-mail. Each Ensembl release is, in turn, synchronised with a certain ClinVar version. Based on Ensembl version, we can find the ClinVar release associated to an Ensembl release in its [sources page](http://www.ensembl.org/info/genome/variation/species/sources_documentation.html). For example, if Ensembl is using ClinVar version “07/2019”, this variable will be set to `2019-07`. Note that this is generally *different* from the Open Targets release year and month.

```bash
# Year and month for the upcoming Open Targets release
export OT_RELEASE=YYYY-MM

# Year and month of ClinVar release used
export CLINVAR_RELEASE_YEAR=YYYY
export CLINVAR_RELEASE_MONTH=MM

# Open Targets JSON schema version
export OT_SCHEMA_VERSION=1.7.3

# Open Targets validator schema version
export OT_VALIDATOR_VERSION=0.6.0

# The root directory for all files for the current batch
export BATCH_ROOT=${BATCH_ROOT_BASE}/batch-${OT_RELEASE}
export CLINVAR_RELEASE=${CLINVAR_RELEASE_YEAR}-${CLINVAR_RELEASE_MONTH}
```

## 2. Set up directories and download the data
```bash
# Create the necessary directories
mkdir -p ${BATCH_ROOT}
cd ${BATCH_ROOT}
mkdir -p clinvar gene_mapping evidence_strings logs

# Download ClinVar data. It is available in several formats. Different formats are convenient for different use cases,
# hence we download three of them: XML dump; VCF summary; and TSV summary.
# ClinVar FTP paths are different depending on whether the current year is the same as the year of the ClinVar release
# which we are trying to download.
if [[ `date +"%Y"` = ${CLINVAR_RELEASE_YEAR} ]]
then
  # Same year
  CLINVAR_XML="xml/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz"
  CLINVAR_TSV="tab_delimited/archive/variant_summary_${CLINVAR_RELEASE}.txt.gz"
else
  # Different (past) year
  CLINVAR_XML="xml/archive/${CLINVAR_RELEASE_YEAR}/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz"
  CLINVAR_TSV="tab_delimited/archive/${CLINVAR_RELEASE_YEAR}/variant_summary_${CLINVAR_RELEASE}.txt.gz"
fi

# ClinVar VCF uses a different path pattern
CLINVAR_VCF_NAME=$(
  curl -sS -l "${CLINVAR_PATH_BASE}/vcf_GRCh38/weekly/" \
  | grep "^clinvar_${CLINVAR_RELEASE_YEAR}${CLINVAR_RELEASE_MONTH}..\.vcf\.gz$" \
  | head -n1
)
CLINVAR_VCF="vcf_GRCh38/weekly/${CLINVAR_VCF_NAME}"

# Download ClinVar XML + TSV + VCF
wget --directory-prefix ${BATCH_ROOT}/clinvar/ \
  "${CLINVAR_PATH_BASE}/${CLINVAR_XML}" \
  "${CLINVAR_PATH_BASE}/${CLINVAR_TSV}" \
  "${CLINVAR_PATH_BASE}/${CLINVAR_VCF}"

# Download the Open Targets JSON schema and update the validator
wget \
  -O ${BATCH_ROOT}/evidence_strings/opentargets-${OT_SCHEMA_VERSION}.json \
  https://raw.githubusercontent.com/opentargets/json_schema/${OT_SCHEMA_VERSION}/opentargets.json
pip3 install --upgrade pip
pip3 install --upgrade opentargets-validator==${OT_VALIDATOR_VERSION}
```

## 3. Manual steps & checks before running the pipeline

### Update ClinVar XML schema version
Schema of ClinVar XML files changes from time to time. The schema version can be obtained by inspecting the XML file header:
```bash
zcat ${BATCH_ROOT}/clinvar/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  | tail -n+2 | head -n1 | sed -e 's|.*clinvar_public_\(.*\)\.xsd.*|\1|'
```

The current supported version is **1.59**. **If the version changes (and only in that case),** we have to regenerate the JAXB binding classes to be able to parse the XML. It can be done using the following commands:
```bash
cd ${CODE_ROOT} && ${BSUB_CMDLINE} \
  -o ${BATCH_ROOT}/logs/update_clinvar_schema.out \
  -e ${BATCH_ROOT}/logs/update_clinvar_schema.err \
  python3 bin/update_clinvar_schema.py \
  -i ${BATCH_ROOT}/clinvar/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  -j clinvar-xml-parser/src/main/java
```

A new Java package will be generated in the directory `clinvar-xml-parser/src/main/java/uk/ac/ebi/eva/clinvar/model`. You will need to rebuild the parser; see [Build instructions](build.md#building-java-clinvar-parser).

Update the test data. All the test does (for the moment) is checking that parsing 10 records from the XML will (1) not crash and (2) provide 10 records after parsing. So to regenerate test data, we just have to extract any 10 records (can just be the first 10 records) from the ClinVar XML file:
```bash
zcat ${BATCH_ROOT}/clinvar/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  | awk 'BEGIN {RS="</ClinVarSet>\n\n"; ORS=RS} {print} NR==10 {exit} END {print "</ReleaseSet>"}' \
  | tee ${CODE_ROOT}/clinvar-xml-parser/src/test/resources/ClinvarExample.xml \
  | gzip -c >${CODE_ROOT}/clinvar-xml-parser/src/test/resources/ClinvarExample.xml.gz
```

Eyeball input & output files to ensure that the ClinVar format has not changed sufficiently enough to render this snippet invalid.

On the cluster, check out a new Git branch, commit all resulting file changes (which will include parser classes and test data), and push them to GitHub (you will need to enter your user name and either a password or an authentication token). Create a pull request to merge this code into the main repository.

Update the currently supported version mentioned in this document.

### Update the code to reflect the Open Targets schema change
Open Targets have a [JSON schema](https://github.com/opentargets/json_schema) used to validate submitted data. Validation of generated evidence strings is carried out during generation.

When Open Targets schema version changes, test files in `tests/evidence_string_generation/resources` also need to be updated, as well as all mentions of old schema version throughout the code (use global search for that). The currently supported schema version is **1.7.3.**

### Check and correct known problematic mappings
There is a [spreadsheet](https://docs.google.com/spreadsheets/d/1m4ld3y3Pfust5JSOJOX9ZmImRCKRGi-fGYj_dExoGj8/edit) which was created to track trait-to-ontology mappings which were especially problematic in the past to users of Open Targets platform. Prior to running subsequent steps, make sure that all traits mentioned in that spreadsheet are mapped to the correct ontology terms in `${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv`.

## 4. Process data
Here, we run several steps in sequence:
* **Transform ClinVar data.** XML file is converted into a JSON file which can be parsed by the downstream tools, using an XML parser which we (if necessary) updated during the previous step. Result is `clinvar.json.gz`.
* **Map variants to gene and functional consequences.** Each evidence string must include the variant, the gene it affects, and the functional effect it has in that gene. We obtain this information separately for repeat expansion variants and for all other types. Detailed information can be found in the [subdirectory where the pipelines are contained](../vep-mapping-pipeline).
  + Run repeat expansion pipeline
  + Run main VEP pipeline
  + Concatenate the results into a single file
* **Generate the evidence strings.** Here, we integrate all of the information produced to generate the evidence strings. This outputs multiple files:
  + `evidence_strings.json` — evidence strings for submitting to Open Targets
  + `eva_clinvar.txt` — ZOOMA feedback with trait names
* Prepare ClinVar xRefs ZOOMA feedback (see more on that below)
* Additionally validate the evidence strings using `opentargets_validator` module

```bash
cd ${CODE_ROOT} && \
${BSUB_CMDLINE} -K -n 8 -M 16G \
  -o ${BATCH_ROOT}/logs/convert_clinvar_files.out \
  -e ${BATCH_ROOT}/logs/convert_clinvar_files.err \
  java -Xmx15G -jar ${CODE_ROOT}/clinvar-xml-parser/target/clinvar-parser-1.0-SNAPSHOT-jar-with-dependencies.jar \
  -i ${BATCH_ROOT}/clinvar/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  -o ${BATCH_ROOT}/clinvar && \
${BSUB_CMDLINE} -K \
  -o ${BATCH_ROOT}/logs/consequence_repeat_expansion.out \
  -e ${BATCH_ROOT}/logs/consequence_repeat_expansion.err \
  python3 ${CODE_ROOT}/vep-mapping-pipeline/run_repeat_expansion_variants.py \
  --clinvar-summary-tsv ${BATCH_ROOT}/clinvar/variant_summary_${CLINVAR_RELEASE}.txt.gz \
  --output-consequences ${BATCH_ROOT}/gene_mapping/consequences_1_repeat.tsv \
  --output-dataframe ${BATCH_ROOT}/gene_mapping/repeat_dataframe.tsv && \
${BSUB_CMDLINE} -K \
  -o ${BATCH_ROOT}/logs/consequence_vep.out \
  -e ${BATCH_ROOT}/logs/consequence_vep.err \
  bash ${CODE_ROOT}/vep-mapping-pipeline/run_consequence_mapping.sh \
  ${BATCH_ROOT}/clinvar/clinvar_${CLINVAR_RELEASE_YEAR}${CLINVAR_RELEASE_MONTH}??.vcf.gz \
  ${BATCH_ROOT}/gene_mapping/consequences_2_vep.tsv && \
cat \
  ${BATCH_ROOT}/gene_mapping/consequences_1_repeat.tsv \
  ${BATCH_ROOT}/gene_mapping/consequences_2_vep.tsv \
  > ${BATCH_ROOT}/gene_mapping/consequences_3_combined.tsv && \
${BSUB_CMDLINE} -K \
  -M 10G \
  -o ${BATCH_ROOT}/logs/evidence_string_generation.out \
  -e ${BATCH_ROOT}/logs/evidence_string_generation.err \
  python3 ${CODE_ROOT}/bin/evidence_string_generation.py \
  -e ${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv \
  -g ${BATCH_ROOT}/gene_mapping/consequences_3_combined.tsv \
  -j ${BATCH_ROOT}/clinvar/clinvar.json.gz \
  --ot-schema ${BATCH_ROOT}/evidence_strings/opentargets-${OT_SCHEMA_VERSION}.json \
  --out ${BATCH_ROOT}/evidence_strings/ && \
${BSUB_CMDLINE} -K \
  -o ${BATCH_ROOT}/logs/traits_to_zooma_format.out \
  -e ${BATCH_ROOT}/logs/traits_to_zooma_format.err \
  python3 ${CODE_ROOT}/bin/clinvar_jsons/traits_to_zooma_format.py \
  -i ${BATCH_ROOT}/clinvar/clinvar.json.gz \
  -o ${BATCH_ROOT}/clinvar/clinvar_xrefs.txt && \
${BSUB_CMDLINE} -K -n 8 -M 16G \
  -o ${BATCH_ROOT}/logs/opentargets_validator.out \
  -e ${BATCH_ROOT}/logs/opentargets_validator.err \
  python3 -m opentargets_validator.cli \
  --schema ${BATCH_ROOT}/evidence_strings/opentargets-${OT_SCHEMA_VERSION}.json \
  ${BATCH_ROOT}/evidence_strings/evidence_strings.json
```

## 5. Manual follow-up actions

### Update summary metrics
After the evidence strings have been generated, summary metrics need to be updated in the Google Sheets [table](https://docs.google.com/spreadsheets/d/1g_4tHNWP4VIikH7Jb0ui5aNr0PiFgvscZYOe69g191k/) on the “Raw statistics” sheet.

### Submit evidence strings
The evidence string file (`evidence_strings.json`) must be uploaded to the [Open Targets Google Cloud Storage](https://console.cloud.google.com/storage/browser/otar012-eva/) and be named in the format `cttv012-[dd]-[mm]-[yyyy].json.gz` (e.g. `cttv012-12-06-2017.json.gz`).

More details can be found on [Open Targets Github wiki](https://github.com/opentargets/data_release/wiki/OT006-Data-Submission#ot009-evidence-string-generation-json-schema-validation--submission).

### Submit feedback to ZOOMA
The idea with [ZOOMA](http://www.ebi.ac.uk/spot/zooma/) is that we not only use it, but also provide feedback to help improve it. The evidence string generation pipeline generates two files with such a feedback:
1. **clinvar_xrefs.** ClinVar data already includes some cross-links from trait names to disease ontologies. Unfortunately, it almost exclusively uses MedGen and OMIM, which are not acceptable for Open Targets (since they're using EFO). However, mappings from trait names to MedGen and OMIM might still be useful to other users. Hence, we extract and submit them to ZOOMA under the evidence handle “ClinVar_xRefs”.
1. **eva_clinvar.** This contains the trait mappings (to EFO) created during the evidence string generation, including automated and manual mappings.

The files are uploaded to the FTP, where ZOOMA will pick it up. At this stage, you only need to upload the **clinvar_xrefs** dataset (the *eva_clinvar* dataset is updated in the process of the manual curation).

To make changes to the FTP, you will need to log in to the cluster using your **personal account** and then run `become <FTP administrative user> /bin/bash`. (Please see [this document](https://www.ebi.ac.uk/seqdb/confluence/display/VAR/Simplified+EVA+FTP+SOP) for details on the FTP administrative user.) After you do this, the environment will be wiped clean, so you will need to set the `BATCH_ROOT` variable again.

```bash
# EXECUTE UNDER FTP ADMINISTRATIVE USER
# DON'T FORGET TO SET THE TWO VARIABLES BELOW AGAIN
export BATCH_ROOT=...
export FTP_PATH_BASE=...

# Create the folder, copy the file to FTP, and update the “latest” folder
FTP_PATH=${FTP_PATH_BASE}/`date +%Y/%m/%d`
mkdir -p ${FTP_PATH}
cp ${BATCH_ROOT}/clinvar/clinvar_xrefs.txt ${FTP_PATH}
cp ${FTP_PATH}/clinvar_xrefs.txt ${FTP_PATH_BASE}/latest/clinvar_xrefs.txt
```

After uploading both files, confirm that the changes have propagated to the FTP:
```bash
md5sum ${BATCH_ROOT}/clinvar/clinvar_xrefs.txt
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/`date +%Y/%m/%d`/clinvar_xrefs.txt | md5sum
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/clinvar_xrefs.txt | md5sum
```

If everything has been done correctly, hash sums will be the same. Note that the FTP may take a few minutes to update after you make changes on the cluster.

## Review checklist
* General
  + All relevant logs are present and contain no unexpected warnings or errors
  + All intermediate and output files
    - Are present
    - In the correct format without obvious errors (compare with the previous OT batch)
    - Line/record count is the same or more compared to the previous OT batch
    - Files are terminated correctly (not ending abruptly after an incomplete record)
* Step 2 “Set up directories and download the data”
  + Working directory exists and its structure is as described in the instructions for the step
* Step 3 “Manual steps & checks before running the pipeline”
  + ClinVar schema version is either unchanged, *or* schema is changed and the PR for the changes has been submitted
    - If the PR is submitted, check that there are no breaking changes in the schema version
    - If the PR is submitted, tests must be updated as well
  + References to the Open Targets schema version are updated throughout the code and in the test files
* Step 4 “Process data”
  + Functional consequences
    - There shouldn't be any warnings such as “Error on last attempt, skipping” in the logs. This might mean that one of the servers was down during processing, and potentially not all information has been gathered.
    - The final file with the genes & consequences contains both “normal variants”, and a few dozen repeat variants, including `short_tandem_repeat_expansion` and `trinucleotide_repeat_expansion` ones.
    - The `${BATCH_ROOT}/logs/consequence_repeat_expansion.err` log will record all repeat expansion variants which could not be parsed using any of the regular expressions. Verify that there are no new such variants compared to the previous batch.
  + Evidence stringts
    - Version of JSON schema is the same as specified in the Open Targets e-mail
    - All traits mentioned in the [spreadsheet](https://docs.google.com/spreadsheets/d/1m4ld3y3Pfust5JSOJOX9ZmImRCKRGi-fGYj_dExoGj8/edit) are mapped to the correct ontology terms in `${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv`.
* Step 5 “Manual follow-up actions”
  + The summary metrics
    - Are present in the [spreadsheet](https://docs.google.com/spreadsheets/d/1g_4tHNWP4VIikH7Jb0ui5aNr0PiFgvscZYOe69g191k/)
    - Are calculated correctly (re-calculate using the commands in the spreadsheet as required)
    - Are not worse than for the previous release
  + Generated evidence strings validate against the schema
    - Have actually been submitted to the Open Targets cloud storage
    - The file is the same as on the cluster (check md5)
    - (After review) E-mail has been sent to Open Targets, with eva-dev in copy
  + ZOOMA feedback (the FTP path is http://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest; to see where files are located on the cluster, see variable `BATCH_ROOT_BASE` on [this page](https://github.com/EBIvariation/configuration/blob/master/open-targets-configuration.md))
    - The changes have been propagated to the FTP, and the files available over FTP are the same as on the cluster
    - The files in the `YYYY/MM/DD` and in the `latest` folders are identical (using either symlinks or copied)
