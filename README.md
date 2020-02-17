# Gene & functional consequence mapping pipeline

## Outline

The pipeline located in the [vep_mapping_pipeline](/vep_mapping_pipeline) directory performs gene & functional consequence mapping. That is, it takes variants as input and maps each of them to a gene (Ensembl gene ID & gene name) and to the most severe functional consequence which a variant causes in that gene, according to Variant Effect Predictor.

## Running the pipeline, input & output formats

### Input file
The pipeline uses succinct, VCF-compatible variant identifiers, in the format of `CHROM:POS:REF:ALT`. Input file must contain variants in exactly this format, one entry per line. Example:
```
10:110821920:G:A
10:110821933:C:T
10:110821934:G:A
```

### Output file
Output is a TSV file consisting of six columns:
1. Variant identifier: the same one as in the input files.
2. The second column is not used and is always set to 1. It is retained for compabitility purposes (see more about that below).
3. Ensembl gene ID.
4. Ensembl gene name.
5. Most severe functional consequence of that variant for that gene.
6. Distance from a variant to a gene. This will be nonzero only for upstream and downstream variants, and zero for all others.

In general, each variant can have more than one record in the output file, because it may affect multiple genes.

Example (note that tabs have been replaced with spaces here for readability):
```
10:110821920:G:A    1    ENSG00000203867    RBM20    missense_variant         0
10:110821933:C:T    1    ENSG00000203867    RBM20    missense_variant         0
10:110821934:G:A    1    ENSG00000203867    RBM20    splice_region_variant    0
```

### Running the script directly (on a small number of entries, for testing/debugging only)
The pipeline core module is [consequence_mapping.py](/bin/consequence_mapping/consequence_mapping.py). It reads data from STDIN and writes to STDOUT in the formats described above. It can be run as follows:
```bash
python3 consequence_mapping.py <input_variants.txt >output_mappings.tsv
```

This should only be done for testing purposes and only for a small number of variants, because when querying VEP API, the script submits of all of the variants it receives in a single query (this is more efficient than submitting them one by one).

### Running the pipeline using a wrapper script
In production environment the pipeline should be run using a wrapper script which would take care of preprocessing and parallelisation. There is a simple wrapper script available, [run_consequence_mapping.sh](/bin/consequence_mapping/run_consequence_mapping.sh). It can be run as follows:
```bash
bash run_consequence_mapping.sh input_variants.vcf output_mappings.tsv
``` 

It performs the following steps:
1. Take regular VCF as input.
1. Extract the necessary fields using `bcftools query`.
1. Sort and remove duplicates.
1. Using `parallel`, split the input file into chunks at most **200** records long and process them using **20** workers in parallel. These values have been selected to process the variants in reasonable time, but not to overload VEP API.
1. Sort the result and remove any duplicate values if present, and save it to a specified output file. 

### Dependencies
The core pipeline module depends only on `retry` module, which can be installed using `pip`: `pip3 install retry`.

The wrapper script depends on `bcftools` and `parallel` (GNU Parallel).

## Mapping process
For each variant:
* Query VEP with default distance to the gene (5,000 bases either way). Consider only consequences which affect **either a protein coding transcript or a miRNA**. Find the most severe consequence for either of those transcripts, according to the list of priorities described [on Ensembl website](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html). Output **all non-redundant** consequences of this (most severe) type.
* If we found nothing during the previous step, repeat VEP search for this variant, now with a distance up to 500,000 bases either way. If there are any consequences which affect a **protein coding** transcript, choose the most severe consequence type (usually this will be either an upstream or a downstream gene variant) and output **a single consequence with the smallest distance**.

## Note on porting & changes from the original pipeline
This pipeline originated as a Python port of the original Open Targets "SNP to gene" pipeline (see https://github.com/opentargets/snp_to_gene). An effort has been taken to retain backwards compatibility where possible; however, many important changes have been introduced.

### Breaking changes
* Original pipeline supported both RS IDs and full variant description for querying. Because RS IDs are in general not allele specific, their support has been dropped. New pipeline only accepts complete, VCF-compatible variant descriptions.
* Input file for the original pipeline was a TSV consisting of 10 columns, most of which were not used, at least in the EVA/ClinVar use case. New pipeline uses a simpler VCF-derived format, which can be used to query VEP directly and can be easily produced from VCF (wrapper script does this already).
* Output format for the new pipeline is a 6 column TSV, mostly the same as for the old pipeline. Changes:
  + Column 1 (variant identifier) is using a different, VCF-compatible notation.
  + Column 5 (functional consequence): special consequence type “nearest_gene_five_prime_end” has been dropped and replaced by conventional “upstream_gene_variant” / “downstream_gene_variant”.
  + Column 6 (distance from variant to gene) is now always non-negative, for both upstream and downstream gene variants. It also always denotes a distance to the gene as reported by VEP, not to the nearest gene 5' end.

### Changes in handling upstream and downstream gene variants
The original mapping process was similar, but not identical, to the new one. The second step did not attempt to queue VEP, but instead searched for the nearest gene 5' end and, if found, output that gene with the special, non-standard consequence type of “nearest_gene_five_prime_end”, along with the computed distance.

However, this consequence type was not used (at least in the EVA/ClinVar use case), and the variants in question are essentially downstream/upstream variants. Hence, the new pipeline handles all downstream/upstream variants in a similar manner, and outputs distance for all of them.

### Changes in handling severity of transcript consequences
The original pipeline contained a serious bug in determining the most severe consequence for a given variant. It worked in the following way:
1. Query VEP with default parameters for a given variant, obtain a list of results.
2. Take note of the `most_severe_consequence` reported by VEP.
3. Filter the list of results based on biotype, leaving only protein coding and miRNA transcripts.
4. Output all consequences in the list from step 3, where type matches the `most_severe_consequence` determined during step 2.

The problem with this approach is that sometimes the `most_severe_consequence` calculated by VEP comes from transcripts of other biotypes (not transcript coding or miRNA), which are filtered out during step 3. This results in the pipeline not outputting any results for such variants.

The new approach does not use VEP's `most_severe_consequence` field. Instead, it first filters the consequences based on the list of acceptable biotypes, and then scans the list of the _remaining_ ones the most severe consequence, based on the severity list described on [Ensembl website](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html). (This is the same list used by VEP internally.) 

### Technical changes (should not affect the results)
* Ensembl REST API (https://rest.ensembl.org/) is used instead of Perl API.
* VEP is queried with multiple variants at once (`vep/:species/region`), rather than querying them one by one (`vep/:species/id/:id` and `vep/:species/region/:region/:allele/`), which greatly speeds up the pipeline and lowers the strain on VEP servers.

## Note on external use
The pipeline is implemented to be as standalone as possible, in hope that it could in future become useful to other people (including, but not limited to, other Open Targets submitters). Feel free to contact the maintainers of this repository if you need any assistance in adopting this pipeline for your use case.
