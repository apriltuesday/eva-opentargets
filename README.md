[![Build Status](https://travis-ci.com/EBIvariation/eva-cttv-pipeline.svg?branch=master)](https://travis-ci.com/EBIvariation/eva-cttv-pipeline)
[![Coverage Status](https://coveralls.io/repos/github/EBIvariation/eva-cttv-pipeline/badge.svg?branch=master)](https://coveralls.io/github/EBIvariation/eva-cttv-pipeline?branch=master)

# Introduction

## ClinVar
[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) is a curated database of clinically relevant genetic variation in humans, maintaned by the National Center for Biotechnology Information in the USA. For each variant, it stores a handful of information:
* **Variation location,** e. g.: *NM_007294.3(BRCA1):c.2706delA* (using [HGVS nomenclature](https://varnomen.hgvs.org/) in this example)
* **Gene** which the variant impacts: *BRCA1*
* **Condition** which is associated with this variant: *Hereditary breast and ovarian cancer syndrome*
* **Clinical significance** of the variant. It is most frequently evaluated using the [ACMG guidelines](https://www.acmg.net/docs/standards_guidelines_for_the_interpretation_of_sequence_variants.pdf); however, other approaches also exist. Most common values of this field are:
  * *Pathogenic* generally means that variation is causing the disease or making it worse
  * *Benign*: it has been established that the variant is *not* implicated in this disease
  * *Likely pathogenic* and *Likely benign* have the same meaning as “Pathogenic” and “Benign”, but signify that the connection isn't as certain (usually due to limited information being available for this variant)
  * *Uncertain significance* means that either there is contradictory information on the impact of this variant, or that information is insufficient to draw the conclusion
* **Review status** shows how many submitters provided information for this variant. Example values:
  + no assertion criteria
  + criteria provided, single submitter
  + criteria provided, multiple submitters, no conflicts

ClinVar is continuously updated and holds monthly releases of its database contents.

## OpenTargets
[OpenTargets](https://www.opentargets.org/) is a collaboration between academia and industry. Among other things, it combines associations between genetic variation and human traits (most notably, diseases) into a single integrated resource. This information is then used to provide evidence on the biological validity of therapeutic targets and an initial assessment of the likely effectiveness of pharmacological intervention on these targets.

OpenTargets also holds periodic releases, which happen approximately every two months. Data for every release comes from several data providers. There are several requirements for submitting data to OpenTargets:
* It must be represented in the form of “evidence strings”. These are JSON strings describing:
  + Genes the variant affects
  + Functional consequence of the variant on the gene
  + Traits (usually diseases) associated with the variant, such as “parkinson disease” or “age-related macular degeneration”.
  + Other information about the variant and source, such as related publications
* The variant data must be synchronised with a specific version of external data sources, for example Ensembl.

## Role of EVA
ClinVar data is highly valuable, but in its original form is not suitable for submission to OpenTargets. EVA is registered as one of the submitters for OpenTargets. For every OpenTargets release, the EVA processes ClinVar records (variants), curates the result and submits it to OpenTargets in the form of evidence strings. This allows for the up-to-date ClinVar data to be integrated into the OpenTargets platform.

Approximately one month before the submission deadline, OpenTargets will contact their submitters and specify the requirements for the next release. At this point the EVA can start executing the main submission protocol (see below). Once the data is ready, it is submitted to OpenTargets, and then the same will happen with the next release. Most of the actions in the pipeline are automated.

# Documentation
* [Build instructions](docs/build.md)
* [How to submit an OpenTargets batch](docs/submit-opentargets-batch.md)
  + [Trait mapping pipeline](docs/trait-mapping-pipeline.md)
  + [Manual curation of traits](docs/manual-curation.md)
* [Comparing the evidence strings](compare-evidence-strings/)

# Workflow diagram

![](docs/workflow.png)

Vector version available [here](https://docs.google.com/presentation/d/1cv-_S1A6j6O3PWbQkVhi4IY7cgLH3o0os5hEN9HVmJA/edit?usp=sharing) (for EBI employees only)