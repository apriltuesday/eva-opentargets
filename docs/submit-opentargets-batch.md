# How to submit an OpenTargets batch
Please follow the instructions below in order to create and submit an OpenTargets batch using ClinVar data. Additional diagrams and background explanations can be found in [this presentation](https://docs.google.com/presentation/d/1nai1dvtfow4RkolyITcymXAsQqEwPJ8pUPcgjLDCntM).

Before starting the process, follow the [Build instructions](build.md). In particular, you'll need to install Python 3.5 (if you don't have it already), build the Java ClinVar parser, build and install the Python pipeline.

Typical OpenTargets submission process consists of five major parts; for each one, an issue is created in the JIRA tracker. Issue template is linked for each of the steps, as well as the list of checks do be done during review.

## Set up environment
Commands below depend on a number of environment variables. It makes sense to set them all at once before executing any of the steps.

```bash
# Year and month for the upcoming OpenTargets release. This is announced in their e-mail.
export OT_RELEASE=YYYY-MM

# Year and month of ClinVar release used (see step 1.1 Download ClinVar data).
# Note that this is *different* from the OpenTargets release year/month.
export CLINVAR_RELEASE=YYYY-MM

# This variable should point to the directory where this repository clone is located on the cluster.
export CODE_ROOT=/nfs/production3/eva/software/eva-cttv-pipeline

# Setting up Python version
PYTHON_VERSION=3.5.6
INSTALL_PATH=/nfs/production3/eva/software/python-${PYTHON_VERSION}
export PATH=${INSTALL_PATH}:$PATH
export PYTHONPATH=${INSTALL_PATH}

# Base bsub command line for all commands. For example, you can specify your e-mail to receive a notification once
# the job has been completed.
export BSUB_CMDLINE="bsub -u your_email@example.com"

# The following variables do not require modification.
export BATCH_ROOT=/nfs/production3/eva/opentargets/batch-${OT_RELEASE}
```

Before proceeding with executing the commands, make sure to update code on the cluster:
```bash
cd $CODE_ROOT
git pull origin master
python setup.py install
```

## Step 1. Create JSON file from ClinVar [(issue template)](https://www.ebi.ac.uk/panda/jira/browse/EVA-1469) 

### 1.1 Download ClinVar data
The working directory for the processing is `/nfs/production3/eva/opentargets`. Note that most of the commands below use `bsub` to submit jobs to LSF cluster rather than executing them directly. If you're not using LSF, then omit `bsub` and its arguments.

Given the year and month the batch is to be released on, run the following command to create the appropriate folders:

```bash
mkdir ${BATCH_ROOT}
cd ${BATCH_ROOT}
mkdir clinvar gene_mapping trait_mapping evidence_strings logs
```

Each OpenTargets release is synchronised with a certain Ensembl release version. The specific version is announced in the e-mail which they send a few weeks before the data submission deadline. Each Ensembl release is, in turn, synchronised with a certain ClinVar version. Based on Ensembl version, we can find the ClinVar release associated to an Ensembl release in its [sources page](http://www.ensembl.org/info/genome/variation/species/sources_documentation.html).

You need to download two files (full release XML and variant summary) from the ClinVar FTP into the `clinvar` subfolder:
```bash
wget --directory-prefix ${BATCH_ROOT}/clinvar/ \
  ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_${CLINVAR_RELEASE}.txt.gz
```

### 1.2. Update ClinVar schema version (if necessary)
Schema of ClinVar XML files changes from time to time. The schema version can be obtained by inspecting the XML file header:
```bash
zcat ${BATCH_ROOT}/clinvar/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  | tail -n+2 | head -n1 | sed -e 's|.*clinvar_public_\(.*\)\.xsd.*|\1|'
```

The current supported version is **1.57**. If the version changes, we have to regenerate the JAXB binding classes to be able to parse the XML. It can be done using the following command:

```bash
# Set the variable below for year and month of the ClinVar release used
cd ${CODE_ROOT} && ${BSUB_CMDLINE} \
  -o ${BATCH_ROOT}/logs/update_clinvar_schema.out \
  -e ${BATCH_ROOT}/logs/update_clinvar_schema.err \
  python bin/update_clinvar_schema.py \
  -i ${BATCH_ROOT}/clinvar/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  -j clinvar-xml-parser/src/main/java
```

A new Java package will be generated in the directory `clinvar-xml-parser/src/main/java/uk/ac/ebi/eva/clinvar/model`. With each schema version change, test data must be updated as well. See details in [Build instructions](build.md#regenerating-test-data). Create a pull request to merge this code into the main repository. It must contain both the updated schema version as well as the new test data. After a schema update, you'll also need to rebuild Java parser (see [Build instructions](build.md#building-java-clinvar-parser)).

### 1.3. Convert ClinVar files
Here we transform ClinVar's XML file into a JSON file which can be parsed by the downstream tools, using an XML parser which we (if necessary) updated during the previous step.

```bash
cd ${CODE_ROOT} && ${BSUB_CMDLINE} -n 8 -M 4G \
  -o ${BATCH_ROOT}/logs/convert_clinvar_files.out \
  -e ${BATCH_ROOT}/logs/convert_clinvar_files.err \
  java -jar ${CODE_ROOT}/clinvar-xml-parser/target/clinvar-parser-1.0-SNAPSHOT-jar-with-dependencies.jar \
  -i ${BATCH_ROOT}/clinvar/ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  -o ${BATCH_ROOT}/clinvar
```

A file named `clinvar.json.gz` will be created in the output directory.

### 1.4. Filter Clinvar JSON file
Clinvar JSON file obtained on the previous step is then filtered, extracting only records with allowed levels of clinical significance (as provided by ClinVar). For example, this step filters out records where the clinical significance is “Benign”, meaning that the variant *does not* contribute to a disease.

```bash
cd ${CODE_ROOT} && ${BSUB_CMDLINE} \
  -o ${BATCH_ROOT}/logs/filter_clinvar_json.out \
  -e ${BATCH_ROOT}/logs/filter_clinvar_json.err \
  python bin/clinvar_jsons/extract_pathogenic_and_likely_pathogenic_variants.py \
  -i ${BATCH_ROOT}/clinvar/clinvar.json.gz \
  -o ${BATCH_ROOT}/clinvar/clinvar.filtered.json.gz
```

### Review checklist
* (See general review checklist at the bottom of the page)
* Working directory exists and its structure is as described in step 1.1 
* ClinVar schema version is either unchanged, *or* schema is changed and the PR for the changes has been submitted
  + If the PR is submitted, check that there are no breaking changes in the schema version
  + If the PR is submitted, tests must be updated as well

## Step 2. Submit files to OpenTargets for gene & functional consequence mapping [(issue tempalte)](https://www.ebi.ac.uk/panda/jira/browse/EVA-1470)
To generate the evidence strings, it is necessary to have the gene mapping and functional consequence annotation for each variant. Currently, this step is performed (partially) by OpenTargets: we generate the necessary file using the command below and submit it to OpenTargets. This needs to be done well in advance of the final submission deadline, so that they have the time to run the gene & functional consequence mapping pipeline and return the results to us. 

```bash
cd ${CODE_ROOT} && ${BSUB_CMDLINE} \
  -o ${BATCH_ROOT}/logs/gene_mapping.out \
  -e ${BATCH_ROOT}/logs/gene_mapping.err \
  python bin/gene_mapping/gene_map_coords.py \
  -i ${BATCH_ROOT}/clinvar/variant_summary_${CLINVAR_RELEASE}.txt.gz \
  -o ${BATCH_ROOT}/gene_mapping/clinvar_${CLINVAR_RELEASE}_coords.tsv.gz
```

The columns in the TSV output file are:
* Chromosome
* Start position
* Stop position
* Reference allele
* Alternate allele
* Strand
* Structural variant type (calculated in script)
* rs ID
* RCV ID (ID for a ClinVar record)
* NCBI gene ID
* sv ID
* Variant type (as specified by ClinVar)

This output file must be uploaded to the [OpenTargets Google Cloud Storage](https://console.cloud.google.com/storage/browser/otar012-eva/). They will fed it into their mapping pipeline, which consists of VEP and some internal custom software.

After uploading the file to the Cloud Storage bucket, you also need to notify the OpenTargets team via e-mail so that they can start their mapping pipeline (they won't receive an automatic notification). The reply from OpenTargets will take some time (usually around a week); however, it doesn't block most of the subsequent steps, which can be done while waiting for the reply. The only step which _is_ blocked is Step 6, “Evidence string generation”.

The TSV file eventually returned by OpenTargets has these columns:

* Variant (rs ID, sv ID, coordinate and alleles, or RCV)
* na (unused)
* Ensembl gene ID
* HGNC gene symbol
* Functional consequence
* na (unused)

### Review checklist (for the submitted file only)
* (See general review checklist at the bottom of the page)
* Submission
  + File has actually been uploaded to the OpenTargets cloud storage
  + File in the cloud storage is the same as on the cluster (compare md5)
  + (After review) e-mail has been sent to OpenTargets with eva-dev in copy
 

## Step 3. Trait mapping pipeline [(issue template)](https://www.ebi.ac.uk/panda/jira/browse/EVA-1471)

See information about the trait mapping pipeline [here](trait-mapping-pipeline.md). It is run with the following command:

```bash
cd ${CODE_ROOT} && ${BSUB_CMDLINE} -M 4G \
  -o ${BATCH_ROOT}/logs/trait_mapping.out \
  -e ${BATCH_ROOT}/logs/trait_mapping.err \
  python bin/trait_mapping.py \
  -i ${BATCH_ROOT}/clinvar/clinvar.filtered.json.gz \
  -o ${BATCH_ROOT}/trait_mapping/automated_trait_mappings.tsv \
  -c ${BATCH_ROOT}/trait_mapping/traits_requiring_curation.tsv
```

### Review checklist
* (See general review checklist at the bottom of the page)
* There shouldn't be any warnings such as “Error on last attempt, skipping” in the logs. This might mean that one of the servers was down during processing, and potentially not all information has been gathered.

## Step 4. Manual curation [(issue template)](https://www.ebi.ac.uk/panda/jira/browse/EVA-1472)
See separate protocol, [Manual curation](manual-curation.md).

### Review checklist
* (See general review checklist at the bottom of the page)
* The mappings selected for each trait are adequate
* Good/bad criteria for curation are observed (see the manual curation protocol, section “Criteria to manually evaluate mapping quality”)
* The number of traits in the `finished_mappings_curation.tsv` file is the same as in the spreadsheet after applying all relevant filters
* _Important:_ spreadhseet does not contain line endings, or extraneous space symbols, in trait names (can be checked by a regexp search)

## Step 5. Generate evidence strings [(issue template)](https://www.ebi.ac.uk/panda/jira/browse/EVA-1473)

Here, we integrate all of the information produced to generate the evidence strings. 

### 5.1. Fetch OpenTargets JSON schema
OpenTargets have a [JSON schema](https://github.com/opentargets/json_schema) used to validate submitted data. Validation of generated evidence strings is carried out during generation. To fetch schema, use the following command. `$VERSION` needs to be filled with the version number recommended by the OpenTargets in their announcement e-mail.

```bash
VERSION=1.6.2
wget \
  -O ${BATCH_ROOT}/evidence_strings/opentargets-$VERSION.json \
  https://raw.githubusercontent.com/opentargets/json_schema/$VERSION/opentargets.json
```

### 5.2. Generate evidence strings
This step depends on OpenTargets providing the output file after processing the data submitted to them on Step 2, “Gene and consequence type mappings”. This file can be named in any fashion but usually has `.out.gz` suffix. It must be downloaded from Google Cloud storage bucket and saved to `${BATCH_ROOT}/gene_mapping/` directory.

In order to generate the evidence strings, run the following command.

```bash
# Set the variable for the name of output file provided by OpenTargets, without the path
export OT_OUTPUT_FILE=mergeEVA_uniq_clinvar_2019-04_19_09_manually_corrected.out.gz
zcat ${BATCH_ROOT}/gene_mapping/${OT_OUTPUT_FILE} > ${BATCH_ROOT}/gene_mapping/ot_mapping_result.out
cd ${CODE_ROOT} && ${BSUB_CMDLINE} \
  -M 10G \
  -o ${BATCH_ROOT}/logs/evidence_string_generation.out \
  -e ${BATCH_ROOT}/logs/evidence_string_generation.err \
  python bin/evidence_string_generation.py \
  -e ${BATCH_ROOT}/trait_mapping/trait_names_to_ontology_mappings.tsv \
  -g ${BATCH_ROOT}/gene_mapping/ot_mapping_result.out \
  -j ${BATCH_ROOT}/clinvar/clinvar.filtered.json.gz \
  --ot-schema ${BATCH_ROOT}/evidence_strings/opentargets-$VERSION.json \
  --out ${BATCH_ROOT}/evidence_strings/
```

This outputs multiple files, including the evidence strings (`evidence_strings.json`) for submitting to OpenTargets and the file of trait mappings for submitting to ZOOMA (`eva_clinvar.txt`).

Generated evidence strings must be additionally validated using tool provided by OpenTargets. _Note: as of August 2019, there is a problem running `opentargets_validator` module using Python 3; as a workaround you can install and run it locally using Python 2.__

```bash
python -m pip install --upgrade pip
python -m pip install --upgrade opentargets-validator
python -m opentargets_validator.cli \
  --schema https://raw.githubusercontent.com/opentargets/json_schema/$VERSION/opentargets.json \
  < ${BATCH_ROOT}/evidence_strings/evidence_strings.json
```

### 5.3. Update summary metrics
After the evidence strings have been generated, summary metrics need to be updated in the Google Sheets [table](https://docs.google.com/spreadsheets/d/1g_4tHNWP4VIikH7Jb0ui5aNr0PiFgvscZYOe69g191k/) on the “Raw statistics” sheet.

## 5.4. Submit evidence strings
The evidence string file (`evidence_strings.json`) must be uploaded to the [OpenTargets Google Cloud Storage](https://console.cloud.google.com/storage/browser/otar012-eva/) and be named in the format `cttv012-[dd]-[mm]-[yyyy].json.gz` (e.g. `cttv012-12-06-2017.json.gz`).

More details can be found on [OpenTargets Github wiki](https://github.com/opentargets/data_release/wiki/OT006-Data-Submission#ot009-evidence-string-generation-json-schema-validation--submission).

### Review checklist
* (See general review checklist at the bottom of the page)
* Version of JSON schema is the same as specified in the OpenTargets e-mail
* The summary metrics
  + Are present in the spreadsheet
  + Are calculated correctly (re-calculate using the commands in the spreadsheet as required)
  + Are not worse than for the previous release
* Generated evidence strings validate against the schema
  + Have actually been submitted to the OpenTargets cloud storage
  + The file is the same as on the cluster (check md5)
  + (After review) E-mail has been sent to OpenTargets, with eva-dev in copy

## Step 6. Submitting feedback to ZOOMA [(issue template)](https://www.ebi.ac.uk/panda/jira/browse/EVA-1490)

The idea with ZOOMA is that we not only use it, but also provide feedback to help improve it. The evidence string generation pipeline generates two files with such a feedback.

### 6.1. Prepare ClinVar xRefs file
Each ClinVar record is associated with one or more traits. When submitting the data to OpenTargets, the trait needs to be specified as an ontological term present in [the Experimental Factor Ontology (EFO)](http://www.ebi.ac.uk/efo/).

Traits in ClinVar may specify ontological terms for the trait name, but not to a term present in EFO. Cross-references (xrefs) to ontological terms can help in finding a good EFO term for a trait, whether an xref is in EFO itself or an xref can be used as a starting point for searching for a term in EFO later on.

Within the ClinVar records are ontology xRefs that link the trait name to a controlled vocabulary. We submit these to [ZOOMA](http://www.ebi.ac.uk/spot/zooma/) under the evidence handle “ClinVar_xRefs”.

The mappings from the ClinVar trait name to the specified ontology xrefs are parsed from the ClinVar JSON file into a tsv suitable for submitting to ZOOMA, excluding any traits which already have mappings in ZOOMA from trusted data sources (EVA, Open Targets, GWAS, Uniprot). In order to do so, execute the following command:

```bash
# Convert trait mappings to the ZOOMA format
cd ${CODE_ROOT} && ${BSUB_CMDLINE} \
  -o ${BATCH_ROOT}/logs/traits_to_zooma_format.out \
  -e ${BATCH_ROOT}/logs/traits_to_zooma_format.err \
  python bin/clinvar_jsons/traits_to_zooma_format.py \
  -i ${BATCH_ROOT}/clinvar/clinvar.filtered.json.gz \
  -o ${BATCH_ROOT}/clinvar/clinvar_xrefs.txt
```

### 6.2. Upload ClinVar xRefs and trait mappings to the FTP

Two files need to be uploaded to the FTP as feedback to ZOOMA: `clinvar_xrefs` (generated during the previous step) and `eva_clinvar` (containing the trait mappings, created during the evidence string generation). This file needs to be uploaded to the FTP. To do this, you will need to run `become ftpadmin /bin/bash` first from your personal account (*not* from `eva_etl` user). After you do this, the environment will be wiped clean, so you'll need to re-apply it again (see “Set up environment” section at the top).

```bash
# EXECUTE UNDER FTPADMIN

# Create the folder
FTP_PATH=/nfs/ftp/pub/databases/eva/ClinVar/`date +%Y/%m/%d`
mkdir -p ${FTP_PATH}

# Copy both files to FTP
cp ${BATCH_ROOT}/clinvar/clinvar_xrefs.txt ${BATCH_ROOT}/evidence_strings/eva_clinvar.txt ${FTP_PATH}

# Update files in the “latest” folder
for ZOOMA_FILE in clinvar_xrefs eva_clinvar; do
    # ln -f -s ${FTP_PATH}/${ZOOMA_FILE}.txt /nfs/ftp/pub/databases/eva/ClinVar/latest/${ZOOMA_FILE}.txt
    # # Note: as of August 2019, for some reason symbolic links aren't working consistently; using copying for now
    cp ${FTP_PATH}/${ZOOMA_FILE}.txt /nfs/ftp/pub/databases/eva/ClinVar/latest/${ZOOMA_FILE}.txt
done
```

After uploading both files, confirm that the changes have propagated to the FTP:
```bash
md5sum ${BATCH_ROOT}/evidence_strings/eva_clinvar.txt
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/eva_clinvar.txt | md5sum
md5sum ${BATCH_ROOT}/clinvar/clinvar_xrefs.txt
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/clinvar_xrefs.txt | md5sum
```

If everything has been done correctly, hash sums will be the same. Note that the FTP will take a few minutes to update after you make changes on the cluster.

### Review checklist
* (See general review checklist at the bottom of the page)
* The changes have been propagated to the FTP, and the files available over FTP are the same as on the cluster
* The files in the `YYYY/MM/DD` and in the `latest` folders are identical (using either symlinks or copied)


## Step 7. Adding novel trait names to EFO
Traits remaining unmapped or poorly mapped can be submitted to EFO if a suitable parent term is available. Any new terms added will be picked up by ZOOMA in the next iteration of the trait mapping pipeline.

Novel traits can be submitted to EFO using the [Webulous templates](https://www.ebi.ac.uk/efo/webulous/) **Add EFO disease** and **Import HP term**. Open a new Google spreadsheet and connect with the server using the Webulous Add-on.

There is a helper script available for preparing the table. `ontology_mappings` must contain a list of ontology identifiers for EFO import (such as `MONDO_123456`, `Orphanet:123456`, etc.), one entry per line. The output file will contains a partially ready table for EFO import.
```bash
python bin/trait_mapping/create_efo_table.py \ 
  -i ontology_mappings.txt \
  -o efo_table.tsv
```

The table needs to be amended manually:
* Some terms will lack descriptions, because ontologies don't always contain a description field for a particular term. If possible, descriptions should be added for all traits.
* Some terms (or their parent terms) might be marked as obsolete. Although an effort is made to exclude such traits during upstream analysis (inside the trait mapping pipeline), sometimes a trait is not properly marked as obsolete in the ontology but its obsoleteness is indicated in its name or in its parent term. The easy way to detect such issues is to search the table for the term “obsolete”. They must be corrected manually by selecting another term or just removed from the import table.

Open a new git issue with EFO to review and import these novel trait names e.g. [https://github.com/EBISPOT/efo/issues/223](https://github.com/EBISPOT/efo/issues/223)

### Review checklist
* (See general review checklist at the bottom of the page)
* Cross-references has been populated for as many traits as possible
* GitHub issue has been created and linked in the issue

## General review checklist for all issues

Applies to all issues:
* All relevant logs are present and contain no unexpected warnings or errors
* All intermediate and output files
  + Are present
  + In the correct format without obvious errors (compare with the previous OT batch)
  + Line/record count is the same or more compared to the previous OT batch
  + Files are terminated correctly (not ending abruptly after an incomplete record) 
