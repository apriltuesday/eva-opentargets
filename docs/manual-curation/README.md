# Manual trait name curation protocol
_Issue template: https://www.ebi.ac.uk/panda/jira/browse/EVA-1911_

## Introduction
Data which we submit to Open Targets represents associations between traits and variants. In the case of ClinVar, traits are almost always inherited diseases. In the Open Targets evidence strings, they must be represented using terms from the [Experimental Factor Ontology](https://www.ebi.ac.uk/efo).

Mapping free-text trait names to ontology terms is the only part of the Open Targets submission process which cannot be automated, and this is the reason it is contained as a separate protocol.

The idea is to run this protocol periodically and independently of the main submission protocol. While at submission time Open Targets data is always synchronised to an older, fixed ClinVar release, all new ClinVar data will be *eventually* incorporated into future Open Targets releases. Hence, it makes sense to do the curation work ahead of time and to decouple this process from the main evidence string generation.

Detailed description of the current approach can be found [in the addendum.](detailed-description.md) This protocol, which is currently performed with TSV files and spreadsheets, will be replaced by a web interface for manual curation before the end of 2020.

## Protocol overview
The protocol consists of four parts which are done in sequence by different people. Some parts are “technical” and are performed by a developer, other are “biological” and are performed by a curator.
1. [**Fetch data**](step1-fetch-clinvar-data.md) (technical). The latest ClinVar data is downloaded and the trait names are extracted. They are attempted to be automatically mapped to ontology terms using ZOOMA. The traits which cannot be mapped automatically are output as a separate file, which is loaded into a Google spreadsheet.
1. [**Curate**](step2-manual-curation.md) (biological). The curator goes through the spreadsheet and fills it in, performing manual curation. Other people review the results and comment on them.
1. [**Extract results**](step3-export-results.md) (technical). Curation results are extracted from the spreadsheet into a TSV file. Some accompanying data is prepared for providing feedback to EFO.
1. [**Provide feedback**](step4-submit-efo-feedback.md) (biological). The curator, using the data generated on the previous steps, submits feedback to EFO and follows up on this. 

## Setting up environment
To follow the technical steps of the protocol, you will need to set up the environment.

First, set up the common environment as explained in the [build instructions](../build.md#setting-up-the-common-environment).

Next, set up the protocol-specific environment, **filling in `${CURATION_RELEASE}`:**
```bash
# Identifier of the current manual curation iteration, which is just the current date.
# Be sure to write this down and to set it to the same value in later parts of this protocol. 
export CURATION_RELEASE=YYYY-MM-DD
export CURATION_RELEASE_ROOT=${BATCH_ROOT_BASE}/manual_curation/${CURATION_RELEASE}
```

## Review checklist
* The mappings selected for each trait are adequate
* Good/bad criteria for curation are observed (see the manual curation protocol, section “Criteria to manually evaluate mapping quality”)
* The number of traits in the `finished_mappings_curation.tsv` file is the same as in the spreadsheet after applying all relevant filters
* _Important:_ spreadhseet does not contain line endings, or extraneous space symbols, in trait names (can be checked by a regexp search)
* For submitting terms to EFO
  + Cross-references has been populated for as many traits as possible
  + GitHub issue has been created and linked in the issue
