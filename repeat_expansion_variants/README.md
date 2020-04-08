# Repeat expansion variants pipeline for ClinVar data

This pipeline is as a companion to the main [functional consequence & gene mapping pipeline](../vep_mapping_pipeline) and is intended to be used in conjunction with it. The main pipeline uses Ensembl's Variant Effect Predictor to assign functional consequences to variants. The problem with VEP is that it does not support repeat expansion variants, which are a special subclass of duplication variants which cause a small repeat unit to increase in number. As a result, the main pipeline does not process these variants.

This pipeline gets data about repeat expansion variants from ClinVar's [TSV summary file](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz). The pipeline can run with the following command:
```bash
bash run_repeat_expansion_variants.sh \
  variant_summary.txt.gz \
  output_consequences.tsv \
  output_dataframe.tsv
```