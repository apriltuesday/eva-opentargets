# Part IV, biological: submit feedback to EFO

## IMPORT terms
The partially ready table for EFO import is available in the "Add EFO disease" sheet in the curation spreadsheet. The table needs to be amended manually:
* Some terms will lack descriptions, because ontologies don't always contain a description field for a particular term. If possible, descriptions should be added for all traits.
* Some terms (or their parent terms) might be marked as obsolete. Although an effort is made to exclude such traits during upstream analysis (inside the trait mapping pipeline), sometimes a trait is not properly marked as obsolete in the ontology but its obsoleteness is indicated in its name or in its parent term. The easy way to detect such issues is to search the table for the term “obsolete”. They must be corrected manually by selecting another term or just removed from the import table.

Open a new git issue with EFO to review and import these novel trait names, e.g. [https://github.com/EBISPOT/efo/issues/223](https://github.com/EBISPOT/efo/issues/223).

## NEW terms
Terms which don't have a suitable mapping cannot be added to the “Add EFO disease“ sheet and must be specified manually in PR description.
