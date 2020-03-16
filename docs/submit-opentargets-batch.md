# How to submit an Open Targets batch
This protocol describes how to process and submit data for an Open Targets batch. Additional diagrams and background explanations can be found in [this presentation](https://docs.google.com/presentation/d/1nai1dvtfow4RkolyITcymXAsQqEwPJ8pUPcgjLDCntM).

Batch submission process consists of three major parts. For each step, create a single JIRA ticket. This can be done by cloning a corresponding ticket from one of the previous batches and adjusting its fields. List of suitable ticket templates for each of the steps:
1. https://www.ebi.ac.uk/panda/jira/browse/EVA-1910
2. https://www.ebi.ac.uk/panda/jira/browse/EVA-1911
3. https://www.ebi.ac.uk/panda/jira/browse/EVA-1912

At the end of each step there is a list of checks do be done during review of the ticket.

Log in to the LSF cluster, where all data processing must take place. You must use a common EVA production user instead of your personal account. Follow the [Build instructions](build.md). In particular, you'll need to install Python 3.8 (if you don't have it already), build the Java ClinVar parser, build and install the Python pipeline.

## Set up environment
Commands throughout the protocol depend on a number of environment variables. It makes sense to set them all at once before executing any of the steps.

First you will need to set up variables which are specific to your installation. For EVA use case, see “configuration” repository.
```bash
# This variable should point to the directory where the clone of this repository is located on the cluster
export CODE_ROOT=

# Location of Python installation which you configured using build instructions
export PYTHON_INSTALL_PATH=

# Location of bcftools installation path
export BCFTOOLS_INSTALL_PATH=

# The directory where subdirectories for each batch will be created
export BATCH_ROOT_BASE=

# Base path of FTP directory on the cluster
export FTP_PATH_BASE=
```

Next set of variables is specific to Open Targets release. They are announced the e-mail which they send a few weeks before the data submission deadline, or can be derived from it:
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
export OT_SCHEMA_VERSION=1.6.6

# Open Targets validator schema version
export OT_VALIDATOR_VERSION=0.6.0
```

Finally, we define some environment variables which are either constant or based on the above two sets:
```bash
# Base bsub command line for all commands. For example, you can specify your e-mail to receive a notification once
# the job has been completed.
export BSUB_CMDLINE="bsub -u your_email@example.com"

# Setting up Python paths
export PATH=${PYTHON_INSTALL_PATH}:${PYTHON_INSTALL_PATH}/bin:${BCFTOOLS_INSTALL_PATH}:$PATH
export PYTHONPATH=${PYTHON_INSTALL_PATH}

# The root directory for all files for the current batch
export BATCH_ROOT=${BATCH_ROOT_BASE}/batch-${OT_RELEASE}
export CLINVAR_RELEASE=${CLINVAR_RELEASE_YEAR}-${CLINVAR_RELEASE_MONTH}
```

Before proceeding with executing the commands, update code on the cluster and create the necessary directories:
```bash
cd $CODE_ROOT
git fetch
git checkout master
git reset --hard origin/master
git submodule update --init --recursive
python3 setup.py install
mkdir -p ${BATCH_ROOT}
cd ${BATCH_ROOT}
mkdir -p clinvar gene_mapping trait_mapping evidence_strings logs
```

Note that most of the commands below use `bsub` to submit jobs to LSF cluster rather than executing them directly. If you're not using LSF, then omit `bsub` and its arguments.

## Step 1. Process ClinVar data

### 1.1 Download data: XML dump & VCF summary

Two files are required (full release XML and VCF variant summary) and are downloaded from the ClinVar FTP into the `clinvar` subfolder:
```bash
wget --directory-prefix ${BATCH_ROOT}/clinvar/ \
  ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/archive/${CLINVAR_RELEASE_YEAR}/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
```

### 1.2. Update ClinVar XML schema version (if necessary)
Schema of ClinVar XML files changes from time to time. The schema version can be obtained by inspecting the XML file header:
```bash
zcat ${BATCH_ROOT}/clinvar/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  | tail -n+2 | head -n1 | sed -e 's|.*clinvar_public_\(.*\)\.xsd.*|\1|'
```

The current supported version is **1.59**. **If the version changes (and only in that case),** we have to regenerate the JAXB binding classes to be able to parse the XML. It can be done using the following commands:
```bash
# Update parser
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

### 1.3. Convert and filter ClinVar files; run the trait mapping pipeline
Here we transform ClinVar's XML file into a JSON file which can be parsed by the downstream tools, using an XML parser which we (if necessary) updated during the previous step.

Resulting Clinvar JSON (`clinvar.json.gz`) is then filtered, extracting only records with allowed levels of clinical significance (as provided by ClinVar). For example, this step filters out records where the clinical significance is “Benign”, meaning that the variant *does not* contribute to a disease. Result is output to `clinvar.filtered.json.gz`.

```bash
cd ${CODE_ROOT} && \
${BSUB_CMDLINE} -K -n 8 -M 16G \
  -o ${BATCH_ROOT}/logs/convert_clinvar_files.out \
  -e ${BATCH_ROOT}/logs/convert_clinvar_files.err \
  java -Xmx15G -jar ${CODE_ROOT}/clinvar-xml-parser/target/clinvar-parser-1.0-SNAPSHOT-jar-with-dependencies.jar \
  -i ${BATCH_ROOT}/clinvar/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  -o ${BATCH_ROOT}/clinvar && \
${BSUB_CMDLINE} -K \
  -o ${BATCH_ROOT}/logs/filter_clinvar_json.out \
  -e ${BATCH_ROOT}/logs/filter_clinvar_json.err \
  python3 bin/clinvar_jsons/extract_pathogenic_and_likely_pathogenic_variants.py \
  -i ${BATCH_ROOT}/clinvar/clinvar.json.gz \
  -o ${BATCH_ROOT}/clinvar/clinvar.filtered.json.gz && \
${BSUB_CMDLINE} -K -M 4G \
  -o ${BATCH_ROOT}/logs/trait_mapping.out \
  -e ${BATCH_ROOT}/logs/trait_mapping.err \
  python3 bin/trait_mapping.py \
  -i ${BATCH_ROOT}/clinvar/clinvar.filtered.json.gz \
  -o ${BATCH_ROOT}/trait_mapping/automated_trait_mappings.tsv \
  -c ${BATCH_ROOT}/trait_mapping/traits_requiring_curation.tsv
```

See additional information about the trait mapping pipeline [here](trait-mapping-pipeline.md).

### 1.4. Perform gene & functional consequence mapping
This step can be done in parallel with the previous one (1.3), because their input and output files are independent.

To generate the evidence strings, it is necessary to have the gene mapping and functional consequence annotation for each variant. This step was previously performed on the Open Targets side, but is now integrated into our workflow:

```bash
cd ${CODE_ROOT} && \
${BSUB_CMDLINE} -K \
  -o ${BATCH_ROOT}/logs/consequence_mapping.out \
  -e ${BATCH_ROOT}/logs/consequence_mapping.err \
  bash ${CODE_ROOT}/vep-mapping-pipeline/vep_mapping_pipeline/run_consequence_mapping.sh \
  ${BATCH_ROOT}/clinvar/clinvar.vcf.gz \
  ${BATCH_ROOT}/gene_mapping/consequence_mapping_result.tsv
```

### Review checklist
* (See general review checklist at the bottom of the page)
* Working directory exists and its structure is as described in step 1.1 
* ClinVar schema version is either unchanged, *or* schema is changed and the PR for the changes has been submitted
  + If the PR is submitted, check that there are no breaking changes in the schema version
  + If the PR is submitted, tests must be updated as well
* There shouldn't be any warnings such as “Error on last attempt, skipping” in the logs. This might mean that one of the servers was down during processing, and potentially not all information has been gathered.

## Step 2. Manual curation
See separate protocol, [Manual curation](manual-curation.md).

After the manual curation has been completed, traits remaining unmapped or poorly mapped should be submitted to EFO if a suitable parent term is available.

### IMPORT terms
Open the curation spreadsheets and use filters to display only terms with the status of `IMPORT`. Copy just the ontology URLs into the file `${BATCH_ROOT}/trait_mapping/efo_ontology_terms.txt`, one URL per line. Example:
```
http://purl.obolibrary.org/obo/HP_0002647
http://purl.obolibrary.org/obo/MONDO_0000727
http://www.orpha.net/ORDO/Orphanet_199306
...
```

Run the helper script to prepare the table for import:
```bash
python3 ${CODE_ROOT}/bin/trait_mapping/create_efo_table.py \
  -i ${BATCH_ROOT}/trait_mapping/efo_ontology_terms.txt \
  -o ${BATCH_ROOT}/trait_mapping/efo_import_table.tsv
```

The file `${BATCH_ROOT}/trait_mapping/efo_import_table.tsv` will contain a partially ready table for EFO import. Copy its contents into the “Add EFO disease” sheet in the curation spreadsheet.

The table needs to be amended manually:
* Some terms will lack descriptions, because ontologies don't always contain a description field for a particular term. If possible, descriptions should be added for all traits.
* Some terms (or their parent terms) might be marked as obsolete. Although an effort is made to exclude such traits during upstream analysis (inside the trait mapping pipeline), sometimes a trait is not properly marked as obsolete in the ontology but its obsoleteness is indicated in its name or in its parent term. The easy way to detect such issues is to search the table for the term “obsolete”. They must be corrected manually by selecting another term or just removed from the import table.

Open a new git issue with EFO to review and import these novel trait names, e.g. [https://github.com/EBISPOT/efo/issues/223](https://github.com/EBISPOT/efo/issues/223).

### NEW terms
Terms which don't have a suitable mapping cannot be added to the “Add EFO disease“ sheet and must be specified manually in PR description.

### Review checklist
* (See general review checklist at the bottom of the page)
* The mappings selected for each trait are adequate
* Good/bad criteria for curation are observed (see the manual curation protocol, section “Criteria to manually evaluate mapping quality”)
* The number of traits in the `finished_mappings_curation.tsv` file is the same as in the spreadsheet after applying all relevant filters
* _Important:_ spreadhseet does not contain line endings, or extraneous space symbols, in trait names (can be checked by a regexp search)
* For submitting terms to EFO
  + Cross-references has been populated for as many traits as possible
  + GitHub issue has been created and linked in the issue

## Step 3. Generate evidence strings
Here, we integrate all of the information produced to generate the evidence strings. 

### 3.1. Fetch Open Targets JSON schema
Open Targets have a [JSON schema](https://github.com/opentargets/json_schema) used to validate submitted data. Validation of generated evidence strings is carried out during generation. To fetch schema, use the following command. `$VERSION` needs to be filled with the version number recommended by the Open Targets in their announcement e-mail.

```bash
wget \
  -O ${BATCH_ROOT}/evidence_strings/opentargets-${OT_SCHEMA_VERSION}.json \
  https://raw.githubusercontent.com/opentargets/json_schema/${OT_SCHEMA_VERSION}/opentargets.json
```

When Open Targets schema version changes, test files in `tests/evidence_string_generation/resources` also need to be updated, as well as all mentions of old schema version throughout the code (use global search for that). 

### 3.2. Generate evidence strings
In order to generate the evidence strings, run the following command:

```bash
cd ${CODE_ROOT} && ${BSUB_CMDLINE} -K \
  -M 10G \
  -o ${BATCH_ROOT}/logs/evidence_string_generation.out \
  -e ${BATCH_ROOT}/logs/evidence_string_generation.err \
  python3 bin/evidence_string_generation.py \
  -e ${BATCH_ROOT}/trait_mapping/trait_names_to_ontology_mappings.tsv \
  -g ${BATCH_ROOT}/gene_mapping/consequence_mapping_result.tsv \
  -j ${BATCH_ROOT}/clinvar/clinvar.filtered.json.gz \
  --ot-schema ${BATCH_ROOT}/evidence_strings/opentargets-${OT_SCHEMA_VERSION}.json \
  --out ${BATCH_ROOT}/evidence_strings/
```

This outputs multiple files, including the evidence strings (`evidence_strings.json`) for submitting to Open Targets and the file of trait mappings for submitting to ZOOMA (`eva_clinvar.txt`).

Generated evidence strings must be additionally validated using tool provided by Open Targets:

```bash
pip3 install --upgrade pip
pip3 install --upgrade opentargets-validator==${OT_VALIDATOR_VERSION}
python3 -m opentargets_validator.cli \
  --schema https://raw.githubusercontent.com/opentargets/json_schema/${OT_SCHEMA_VERSION}/opentargets.json \
  < ${BATCH_ROOT}/evidence_strings/evidence_strings.json
```

### 3.3. Update summary metrics
After the evidence strings have been generated, summary metrics need to be updated in the Google Sheets [table](https://docs.google.com/spreadsheets/d/1g_4tHNWP4VIikH7Jb0ui5aNr0PiFgvscZYOe69g191k/) on the “Raw statistics” sheet.

### 3.4. Submit evidence strings
The evidence string file (`evidence_strings.json`) must be uploaded to the [Open Targets Google Cloud Storage](https://console.cloud.google.com/storage/browser/otar012-eva/) and be named in the format `cttv012-[dd]-[mm]-[yyyy].json.gz` (e.g. `cttv012-12-06-2017.json.gz`).

More details can be found on [Open Targets Github wiki](https://github.com/opentargets/data_release/wiki/OT006-Data-Submission#ot009-evidence-string-generation-json-schema-validation--submission).

### 3.5. Submit feedback to ZOOMA
The idea with [ZOOMA](http://www.ebi.ac.uk/spot/zooma/) is that we not only use it, but also provide feedback to help improve it. The evidence string generation pipeline generates two files with such a feedback:
1. **clinvar_xrefs.** ClinVar data already includes some cross-links from trait names to disease ontologies. Unfortunately, it almost exclusively uses MedGen and OMIM, which are not acceptable for Open Targets (since they're using EFO). However, mappings from trait names to MedGen and OMIM might still be useful to other users. Hence, we extract and submit them to ZOOMA under the evidence handle “ClinVar_xRefs”.
1. **eva_clinvar.** This contains the trait mappings (to EFO) created during the evidence string generation, including automated and manual mappings.

#### 3.5.1. Prepare ClinVar xRefs file
The mappings are parsed from the ClinVar JSON file into a TSV suitable for submitting to ZOOMA, excluding any traits which already have mappings in ZOOMA from trusted data sources (EVA, Open Targets, GWAS, Uniprot). In order to do so, execute the following command:

```bash
# Convert trait mappings to the ZOOMA format
cd ${CODE_ROOT} && ${BSUB_CMDLINE} \
  -o ${BATCH_ROOT}/logs/traits_to_zooma_format.out \
  -e ${BATCH_ROOT}/logs/traits_to_zooma_format.err \
  python3 bin/clinvar_jsons/traits_to_zooma_format.py \
  -i ${BATCH_ROOT}/clinvar/clinvar.filtered.json.gz \
  -o ${BATCH_ROOT}/clinvar/clinvar_xrefs.txt
```

### 3.5.2. Upload ClinVar xRefs and trait mappings to the FTP
You need to upload two files to the FTP as feedback to ZOOMA: `clinvar_xrefs` (prepared during the previous step) and `eva_clinvar`.

To make changes to the FTP, you will need to log in to the cluster using your **personal account** and then run `become <FTP administrative user> /bin/bash`. (Please see [this document](https://www.ebi.ac.uk/seqdb/confluence/display/VAR/Simplified+EVA+FTP+SOP) for details on the FTP administrative user.) After you do this, the environment will be wiped clean, so you will need to set the `BATCH_ROOT` variable again.

```bash
# EXECUTE UNDER FTP ADMINISTRATIVE USER
# DON'T FORGET TO SET BATCH_ROOT AGAIN
export BATCH_ROOT=...

# Create the folder

FTP_PATH=${FTP_PATH_BASE}/`date +%Y/%m/%d`
mkdir -p ${FTP_PATH}

# Copy both files to FTP
cp ${BATCH_ROOT}/clinvar/clinvar_xrefs.txt ${BATCH_ROOT}/evidence_strings/eva_clinvar.txt ${FTP_PATH}

# Update files in the “latest” folder
for ZOOMA_FILE in clinvar_xrefs eva_clinvar; do
    # ln -f -s ${FTP_PATH}/${ZOOMA_FILE}.txt ${FTP_PATH_BASE}/latest/${ZOOMA_FILE}.txt
    # # Note: as of August 2019, for some reason symbolic links aren't working consistently; using copying for now
    cp ${FTP_PATH}/${ZOOMA_FILE}.txt ${FTP_PATH_BASE}/latest/${ZOOMA_FILE}.txt
done
```

After uploading both files, confirm that the changes have propagated to the FTP:
```bash
md5sum ${BATCH_ROOT}/evidence_strings/eva_clinvar.txt
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/eva_clinvar.txt | md5sum
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/`date +%Y/%m/%d`/eva_clinvar.txt | md5sum

md5sum ${BATCH_ROOT}/clinvar/clinvar_xrefs.txt
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/`date +%Y/%m/%d`/clinvar_xrefs.txt | md5sum
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/clinvar_xrefs.txt | md5sum
```

If everything has been done correctly, hash sums will be the same. Note that the FTP may take a few minutes to update after you make changes on the cluster.

### Review checklist
* (See general review checklist at the bottom of the page)
* Version of JSON schema is the same as specified in the Open Targets e-mail
* The summary metrics
  + Are present in the [spreadsheet](https://docs.google.com/spreadsheets/d/1g_4tHNWP4VIikH7Jb0ui5aNr0PiFgvscZYOe69g191k/)
  + Are calculated correctly (re-calculate using the commands in the spreadsheet as required)
  + Are not worse than for the previous release
* Generated evidence strings validate against the schema
  + Have actually been submitted to the Open Targets cloud storage
  + The file is the same as on the cluster (check md5)
  + (After review) E-mail has been sent to Open Targets, with eva-dev in copy
* ZOOMA feedback
  + The changes have been propagated to the FTP, and the files available over FTP are the same as on the cluster
  + The files in the `YYYY/MM/DD` and in the `latest` folders are identical (using either symlinks or copied)

## General review checklist for all steps
* All relevant logs are present and contain no unexpected warnings or errors
* All intermediate and output files
  + Are present
  + In the correct format without obvious errors (compare with the previous OT batch)
  + Line/record count is the same or more compared to the previous OT batch
  + Files are terminated correctly (not ending abruptly after an incomplete record)
