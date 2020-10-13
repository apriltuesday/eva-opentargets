# ClinVar data model and attribute value distributions

The script in this directory parses ClinVar XML types and calculates statistics on all possible ways the variants can be represented. The results are described below.

## Running the script

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz
python3 \
  clinvar-variant-types.py \
  --clinvar-xml ClinVarFullRelease_00-latest.xml.gz
```

## Results

Graphs can be enlarged by clicking on them. Dates in parentheses specify when the graph was last updated.

### Data model and variant types (2020-07-06)

![](variant-types.png)

**RCV** is the top level of ClinVar data organisation. It is a record which associates one or more traits (usually diseases) with exactly one _VCV record,_ which can be one of two types:
* **MeasureSet** contains one or more _Measures_ (which are basically individual, isolated variants). Its type can be one of four values:
  - **Variant.** This means that the measure “set” has the length of 1 and contains just a single isolated variant. This variant can be one of the following subtypes, listed in the decreasing order of occurrence:
    + single nucleotide variant
    + Deletion
    + copy number loss
    + copy number gain
    + Duplication
    + Microsatellite
    + Indel
    + Insertion
    + Variation
    + Inversion
    + Translocation
    + protein only
    + Complex
    + fusion
    + Tandem duplication
  - Three other complex types, which were not investigated further in this analysis. They may contain multiple Measures (variants), which must all be interpreted together:
    + **Haplotype.** A collection of variants phased on the same chromosome copy and usually inherited together.
    + **Phase unknown**
    + **Distinct chromosomes**
* **GenotypeSet** represents the cases when the variants which are interpreted together are located on different chromosomal copies (paternal/maternal), that is, when they include _trans_ phasing. The GenotypeSet can be one of two types, which were not investigated further in this analysis:
  - **CompoundHeterozygote.** Presumably this should include exactly two variants which are _trans_ phased and interpreted together.
  - **Diplotype.** Similar, but at least one of the _trans_ phased alleles includes a haplotype. An example of this would be three variants located on one copy of the gene, and one variant in the second one, all interpreted together.

As of July 2020, the most common case is the MeasureSet/Variant one, accounting for 1114689 out of 1117817 RCV records, or >99.7%. **Currently, this is the only type being processed by this pipeline.**

### Clinical significance (2020-07-06)

Calculated from `ClinVarFullRelease_2020-0706.xml.gz`.

![](clinical-significance.png)

Under the current criteria, 188,518 out of 1,114,689 (17%) records are being processed.

For the situations where multiple clinical significance levels were reported for a given association, they are converted into a single composite string, e.g. `Benign/Likely benign, other`.

### Star rating and review status (2020-07-06)

![](star-rating.png)

The distribution of records by star rating is:
* ☆☆☆☆ 142,855 (13%)
* ★☆☆☆ 894,109 (80%)
* ★★☆☆ 66,107 (6%)
* ★★★☆ 11,583 (1%)
* ★★★★ 35 (< 0.01%)

### Mode of inheritance (2020-07-21)

![](mode-of-inheritance.png)

Only a small fraction of all records specify their mode of inheritance: 35,009 out of 1,114,689, or about 3%.

### Allele origin (2020-10-13)

![](allele-origin.png)

All records specify an allele origin. It can be either a single value (the majority of them) or multiple ones.