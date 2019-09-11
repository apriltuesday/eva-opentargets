# Overview of trait mapping pipeline

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
