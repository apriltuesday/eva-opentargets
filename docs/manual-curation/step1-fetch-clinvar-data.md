# Manual curation, part I, technical: fetch the latest ClinVar data, attempt automatic mapping, extract unmapped traits

## Run the automated protocol
Before running, set up the environment:
* [Common environment](../environment.md)
* [Protocol-specific environment](README.md#setting-up-environment)

```bash
# Create directories for data processing
mkdir -p ${CURATION_RELEASE_ROOT}
cd ${CURATION_RELEASE_ROOT}

# Run the nextflow pipeline, resuming execution of previous attempt if possible.
nextflow run ${CODE_ROOT}/eva_cttv_pipeline/trait_mapping/pipeline.nf \
  --curation_root ${CURATION_RELEASE_ROOT} \
  -resume
```

## Create a Google spreadsheet for curation

Duplicate a [template](https://docs.google.com/spreadsheets/d/1PyDzRs3bO1klvvSv9XuHmx-x7nqZ0UAGeS6aV2SQ2Yg/edit?usp=sharing). Paste the contents of `${CURATION_RELEASE_ROOT}/google_sheets_table.tsv` file into it, starting with column H “ClinVar label”. Example of a table fully populated with data can be found [here](https://docs.google.com/spreadsheets/d/1HQ08UQTpS-0sE9MyzdUPO7EihMxDb2e8N14s1BknjVo/edit?usp=sharing).
