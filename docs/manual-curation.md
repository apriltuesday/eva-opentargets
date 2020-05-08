# Manual trait name curation protocol
_Issue template: https://www.ebi.ac.uk/panda/jira/browse/EVA-1911_

## Introduction
Data which we submit to Open Targets represents associations between traits and variants. In the case of ClinVar, traits are almost always inherited diseases. In the Open Targets evidence strings, they must be represented using terms from the [Experimental Factor Ontology](https://www.ebi.ac.uk/efo).

Mapping free-text trait names to ontology terms is the only part of the Open Targets submission process which cannot be automated, and this is the reason it is contained as a separate protocol.

The idea is to run this protocol periodically and independently of the main submission protocol. While at submission time Open Targets data is always synchronised to an older, fixed ClinVar release, all new ClinVar data will be *eventually* incorporated into future Open Targets releases. Hence, it makes sense to do the curation work ahead of time and to decouple this process from the main evidence string generation.

This protocol, which is currently performed with TSV files and spreadsheets, will be replaced by a web interface for manual curation before the end of 2020.

## Protocol overview
The protocol consists of four parts which are done in sequence by different people. Some parts are “technical” and are performed by a developer, other are “biological” and are performed by a curator.
1. **Fetch data (technical).** The latest ClinVar data is downloaded and the trait names are extracted. They are attempted to be automatically mapped to ontology terms using ZOOMA. The traits which cannot be mapped automatically are output as a separate file, which is loaded into a Google spreadsheet.
1. **Curate (biological).** The curator goes through the spreadsheet and fills it in, performing manual curation. Other people review the results and comment on them.
1. **Extract results (technical).** Curation results are extracted from the spreadsheet into a TSV file. Some accompanying data is prepared for providing feedback to EFO.
1. **Provide feedback (biological).** The curator, using the data generated on the previous steps, submits feedback to EFO and follows up on this. 

## Setting up environment
First, set up the common environment as explained in the [build instructions](build.md).

Next, set up the protocol-specific environment, **filling in `${CURATION_RELEASE}` first:**
```bash
# Identifier of the current manual curation iteration, which is just the current date.
# Be sure to write this down and to set it to the same value in later parts of this protocol. 
export CURATION_RELEASE=YYYY-MM-DD
export CURATION_RELEASE_ROOT=${BATCH_ROOT_BASE}/manual_curation/${CURATION_RELEASE}
```



# Part I, technical: fetch the latest ClinVar data, attempt automatic mapping, extract unmapped traits

## Run the automated protocol
_Before running, set up environment as described above._

```bash
# Create directories for data processing
mkdir -p ${CURATION_RELEASE_ROOT}

# Download the latest ClinVar variant summary data
wget -q \
  --directory-prefix ${CURATION_RELEASE_ROOT} \
  ${CLINVAR_PATH_BASE}/tab_delimited/variant_summary.txt.gz

# Run the trait mapping pipeline
cd ${CODE_ROOT} && ${BSUB_CMDLINE} -K -M 4G \
  -o ${CURATION_RELEASE_ROOT}/log.trait_mapping.out \
  -e ${CURATION_RELEASE_ROOT}/log.trait_mapping.err \
  python3 bin/trait_mapping.py \
  -i ${CURATION_RELEASE_ROOT}/variant_summary.txt.gz \
  -o ${CURATION_RELEASE_ROOT}/automated_trait_mappings.tsv \
  -c ${CURATION_RELEASE_ROOT}/traits_requiring_curation.tsv

# Download the latest eva_clinvar release from FTP. At this step, mappings produced by the pipeline on the previous
# iteration (including automated and manual) are downloaded to be used to aid the manual curation process.
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/eva_clinvar.txt \
  | cut -f4-5 | sort -u > ${CURATION_RELEASE_ROOT}/previous_mappings.tsv

# Create the final table for manual curation
cd ${CODE_ROOT} && python3 bin/trait_mapping/create_table_for_manual_curation.py \
  --traits-for-curation ${CURATION_RELEASE_ROOT}/traits_requiring_curation.tsv \
  --previous-mappings ${CURATION_RELEASE_ROOT}/previous_mappings.tsv \
  --output ${CURATION_RELEASE_ROOT}/table_for_manual_curation.tsv

# Sort and export to Google Sheets. Note that the number of columns in the output table is limited to 50, because only a
few traits have that many mappings, and in virtually all cases these extra mappings are not meaningful. However, having
a very large table degrades the performance of Google Sheets substantially.
cut -f-50 ${CURATION_RELEASE_ROOT}/table_for_manual_curation.tsv \
  | sort -t$'\t' -k2,2rn > ${CURATION_RELEASE_ROOT}/google_sheets_table.tsv
```

## Create a Google table for curation

Duplicate a [template](https://docs.google.com/spreadsheets/d/1PyDzRs3bO1klvvSv9XuHmx-x7nqZ0UAGeS6aV2SQ2Yg/edit?usp=sharing). Paste the contents of `${CURATION_RELEASE_ROOT}/google_sheets_table.tsv` file into it, starting with column H “ClinVar label”. Example of a table fully populated with data can be found [here](https://docs.google.com/spreadsheets/d/1HQ08UQTpS-0sE9MyzdUPO7EihMxDb2e8N14s1BknjVo/edit?usp=sharing).



# Part II, biological: perform manual curation

The goal is for traits with occurence ≥ 10 to have 100% coverage after the manual curation. For the rest of the traits, curate as many as possible.

Good mappings must be eyeballed to ensure they are actually good. Alternative mappings for medium or low quality mappings can be searched for using OLS. If a mapping can't be found in EFO, look for a mapping to a HP, ORDO, or MONDO trait name. Most HP/ORDO/MONDO terms will also be in EFO but some are not. These can be imported to EFO using the Webulous submission service.

## Criteria to manually evaluate mapping quality
* Exact string for string matches are _good_
* Slight modifications are _good_ e.g. IRAK4 DEFICIENCY → Immunodeficiency due to interleukin-1 receptor-associated kinase-4 deficiency
* Subtype to parent are _good_ e.g ACHROMATOPSIA 3 → Achromatopsia
* Parent to subtype are _bad_ e.g. HEMOCHROMATOSIS → Hemochromatosis type 3
* Familial / congenital represented on only one half are _bad_ e.g. Familial renal glycosuria → Renal glycosuria
* Susceptibility on only one half is _bad_ e.g Alcohol dependence, susceptibility to → alcohol dependence
* Early / late onset on only one half is _bad_ e.g. Alzheimer disease, early-onset → Alzheimer's disease

## Unmapped trait names
Trait names that haven't been automatically mapped against any ontology term can also be searched for using OLS. If a mapping can't be found in EFO, look for a mapping to a HP, ORDO, or MONDO trait name. If these are not already in EFO they should be imported to EFO using the Webulous submission service.

## Curation workflow
Curation should be done by subsequently applying filters to appropriate columns, then making decisions for the traits in the filtered selection.

* 1\. **There is a previously assigned mapping for this trait.** All of these are the decisions that we made in the past, so we trust them (to an extent). Copy and paste previously used mappings into “Mapping to use”. Then review them according to the following steps.
  * 1.1. **The previously assigned mapping is in EFO**
    * 1.1.1. **The previously assigned mapping is in EFO and is exact.** Mark as finished immediately. (It's extremely unlikely that a better mapping could be found).
    * 1.1.2. **The previously assigned mapping is in EFO and IS NOT exact.** Review the mappings to see if a better (more accurate/specific) mapping is available. Then mark as finished.
  * 1.2. **The previously assigned mapping is not contained in EFO.** We need to either find a mapping which is already in EFO, or import these terms into EFO.
    * 1.2.1. **The previously used mapping IS NOT contained in EFO and is exact.** These are good candidates to mark as finished and them import in EFO afterwards. However, quickly check whether there are non-exact matches which are already in EFO are are as good as exact mappings.
      * E. g. if the exact mapping is “erythrocytosis 6, familial” and not in EFO, but there is an inexact mapping “familial erythrocytosis 6” which *is* in EFO, we should use the inexact mapping.
      * If a trait does not have any EFO mappings, it's probably safe to mark it as finished (with subsequent import to EFO).
    * 1.2.2. **The previously assigned mapping IS NOT contained in EFO and IS NOT exact.** Similarly to 1.2.1, attempt to find an acceptable EFO mapping; if not found, use any acceptable mapping (with subsequent import to EFO).
* 2\. **There is no previously assigned mappings for the trait, but exact mappings are available.** Because letter-to-letter matches are extremely likely to be correct, we can use them after eyeballing for correctness.
  * 2.1. **The exact mapping in the EFO.** Mark as finished immediately.
  * 2.2. **The exact mapping IS NOT in the EFO.** Similarly to 1.2.1, attempt to find an acceptable EFO mapping; if not found, use any acceptable mapping (with subsequent import to EFO).
* 3\. **There are no previously used or exact mappings for the trait.** Curate manually as usual.

### Time-saving options
The new manual workflow can be shortened if necessary, while the quality of the results will be _at least as good as for the old workflow_ (because we're reusing the results of previous curations):
* All subsections 1.\* involve review of mappings previously selected by ourselves. Because we trust them (to an extent), this review can be applied not to all mappings, but only to some (selected on a basis of frequency, or just randomly sampled/eyeballed).
* If necessary, section 1 can be skipped completely, i. e. copy-paste previous mappings into “Mapping to use” column, but skip the review.
* Sections 2.2 and 3 can only be applied to some variants (e. g. based on frequency), depending on the time available.

## Entering the curation results

### Adding new mappings
To select a new mapping which does not appear in the list of automatically generated mappings, use the following format: `URL|LABEL|||EFO_STATUS`. Example: `http://www.ebi.ac.uk/efo/EFO_0006329|response to citalopram|||EFO_CURRENT`. The last value can be either `EFO_CURRENT` (trait is present in the latest EFO version available in OLS), or `NOT_CONTAINED` if the term is not contained in the EFO.

### Marking the status of curated terms
The “Status” column has the following acceptable values:
* **DONE** — an acceptable trait contained in EFO has been found for the trait
* **IMPORT** — an acceptable trait has been found from the MONDO/ORDO/HP ontologies which is not contained in EFO and must be imported
* **NEW** — new term must be created in EFO
* **SKIP** — trait is going to be skipped in this iteration, due to being too non-specific, or just having a low frequency
* **UNSURE** — temporary status; traits to be discussed with reviewers/the team

“Comment” field can contain arbitrary additional information.

### Note on spaces and line breaks
Sometimes, especially when copy-pasting information from external sources, a mapping label or URL can contain an additional space symbol (at the beginning or end) or an accidental line break. This causes problems in the downstream processing and must be manually removed. To minimise the occurences of this, Google Sheets template includes a validation formula for the first two columns (“URI of selected mapping” and “Label of selected mapping”). If it detects an extra space symbol or a line break, the cell will be highlighted in red.



# Part III, technical: export curation results
_Before running, set up environment as described above._

## Extract curation results from the spreadsheet

### Finished mappings
Once the manual curation is completed, apply a spreadsheet filter so that only traits with **Status = DONE** are visible. Copy data for all non-empty rows from three columns: “ClinVar label”; “URI of selected mapping”; “Label of selected mapping”, in that order. **Do not include header lines.** Save the data to a file `${CURATION_RELEASE_ROOT}/finished_mappings_curation.tsv`.

### Terms requiring import into EFO
After the manual curation has been completed, traits remaining unmapped or poorly mapped should be submitted to EFO if a suitable parent term is available. Open the curation spreadsheet and use filters to display only terms with **Status = IMPORT.** Copy _just the ontology_ URLs into the file `${CURATION_RELEASE_ROOT}/terms_for_efo_import.txt`, one URL per line, for example Example:
```
http://purl.obolibrary.org/obo/HP_0002647
http://purl.obolibrary.org/obo/MONDO_0000727
http://www.orpha.net/ORDO/Orphanet_199306
```

## Run the automated protocol

```bash
# Concatenate finished automated and manual mappings into a single file
cat \
  ${CURATION_RELEASE_ROOT}/automated_trait_mappings.tsv \
  ${CURATION_RELEASE_ROOT}/finished_mappings_curation.tsv \
> ${CURATION_RELEASE_ROOT}/trait_names_to_ontology_mappings.tsv

# Update the symbolic link pointing to the location of the most recent curation result. This will be used by the main
# evidence string generation protocol.
ln -s -f \
  ${CURATION_RELEASE_ROOT}/trait_names_to_ontology_mappings.tsv \
  ${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv

# Run the helper script to prepare the table for import
python3 ${CODE_ROOT}/bin/trait_mapping/create_efo_table.py \
  -i ${CURATION_RELEASE_ROOT}/terms_for_efo_import.txt \
  -o ${CURATION_RELEASE_ROOT}/efo_import_table.tsv
```

## Copy the table for EFO import
The file `${CURATION_RELEASE_ROOT}/efo_import_table.tsv` will contain a partially ready table for EFO import. Copy its contents into the “Add EFO disease” sheet in the curation spreadsheet.



# Part IV, biological: submit feedback to EFO

## IMPORT terms
The partially ready table for EFO import is available in the "Add EFO disease" sheet in the curation spreadsheet. The table needs to be amended manually:
* Some terms will lack descriptions, because ontologies don't always contain a description field for a particular term. If possible, descriptions should be added for all traits.
* Some terms (or their parent terms) might be marked as obsolete. Although an effort is made to exclude such traits during upstream analysis (inside the trait mapping pipeline), sometimes a trait is not properly marked as obsolete in the ontology but its obsoleteness is indicated in its name or in its parent term. The easy way to detect such issues is to search the table for the term “obsolete”. They must be corrected manually by selecting another term or just removed from the import table.

Open a new git issue with EFO to review and import these novel trait names, e.g. [https://github.com/EBISPOT/efo/issues/223](https://github.com/EBISPOT/efo/issues/223).

## NEW terms
Terms which don't have a suitable mapping cannot be added to the “Add EFO disease“ sheet and must be specified manually in PR description.



# Review checklist
* The mappings selected for each trait are adequate
* Good/bad criteria for curation are observed (see the manual curation protocol, section “Criteria to manually evaluate mapping quality”)
* The number of traits in the `finished_mappings_curation.tsv` file is the same as in the spreadsheet after applying all relevant filters
* _Important:_ spreadhseet does not contain line endings, or extraneous space symbols, in trait names (can be checked by a regexp search)
* For submitting terms to EFO
  + Cross-references has been populated for as many traits as possible
  + GitHub issue has been created and linked in the issue



# Addendum: in-depth description of the automated trait mapping pipeline

When specifying which disease the variant is associated with, ClinVar uses a free-text description, such as “Alzheimer's disease”. However, it is required by OpenTargets that diseases (traits) are specified as EFO terms. EFO is an ontology (basically a controlled hierarchical dictionary) developed at EBI to standartise nomenclature of diseases and other terms.

To convert free-text descriptions to EFO terms, a trait mapping pipeline is used. You can find a diagram of the whole workflow [here](https://docs.google.com/presentation/d/1nai1dvtfow4RkolyITcymXAsQqEwPJ8pUPcgjLDCntM/edit#slide=id.g24b2b34015_0_531). The detailed description is presented below.

## Querying ZOOMA
ZOOMA is first queried using the trait name. ZOOMA is a EBI-developed resource to take free-text trait descriptions as input and provide ontology terms as output. For example, if you query “Alzheimer's disease”, it should return http://www.ebi.ac.uk/efo/EFO_0000249. The pipeline uses ZOOMA results in the following manner:

* If there are any high confidence mappings from EFO then they are output for use and no further processing is taken for that trait name.
* If there are lower confidence mappings from EFO then these are output in a separate curation file, and the process stops here for that trait name.
* If there are high confidence mappings not from EFO then their ontology IDs are used as input for querying OxO (see below)

## Querying OxO
[OxO](http://www.ebi.ac.uk/spot/oxo/) is a tool which attempts to find ontology cross references for provided identifiers. It is used by the pipeline as follows:

* Any EFO mappings found within a distance of 1 step are output for use, no further processing for this trait name.
* Any EFO mappings within a distance greater than 1 step are output for curation.

## Output
The output of the trait mapping pipeline consists of 2 files.

### Automated trait mappings
The file `automated_trait_mappings.tsv` contains the following columns:

* Trait name from ClinVar
* Ontology URI
* Ontology term label

### Manual curation required
The file `traits_requiring_curation.tsv` contains the entries that require manual curation, including any trait name
with no mappings. It contains the following columns:

* Trait name from ClinVar
* Frequency of this trait in ClinVar
* ZOOMA mappings:
  * URI
  * Ontology label
  * Confidence
  * Datasource/annotator
* OxO mappings:
  * URI
  * Ontology label
  * Confidence
  * Distance
