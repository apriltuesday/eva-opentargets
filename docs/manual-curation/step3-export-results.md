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