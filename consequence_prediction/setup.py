from setuptools import setup

description = """Contains two packages for processing ClinVar data for Open Targets submission.
* `vep_mapping_pipeline` maps variants (CHROM:POS:REF:ALT) to their most severe functional consequence according to
  Ensembl VEP, as well as their Ensembl gene ID and name.
* `repeat_expansion_variants` parses ClinVar variant_summary file and extracts information about repeat expansion
  variants."""

setup(
    name='vep-mapping-pipeline',
    version='',
    packages=['vep_mapping_pipeline', 'repeat_expansion_variants', 'repeat_expansion_variants.test'],
    url='',
    license='',
    author='European Variation Archive, EMBL-EBI',
    author_email='eva-dev@ebi.ac.uk',
    description=description
)
