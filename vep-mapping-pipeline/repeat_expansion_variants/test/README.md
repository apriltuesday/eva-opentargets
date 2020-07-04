# Tests for the repeat expansion pipeline

ClinVar repeat expansion data includes a number of peculiarities. To check them all in separate unit tests would be
expensive to develop and maintain. Hence, the pipeline uses a hybrid integration test with an annotated dataset.

The dataset includes the input file [`input_variant_summary.tsv`](input_variant_summary.tsv) and two expected output
files: [`output_dataframe.tsv`](output_dataframe.tsv) and [`output_consequences.tsv`](output_consequences.tsv). The
input file is not a sample, but rather a complete selection of “NT expansion” variants from ClinVar data as of
2020-04-08. The expected output files were produced by the pipeline and checked manually for correctness. The idea
behind including the entire dataset is that it will make the tests sensitive to even minor changes.

The test files are annotated using comments, which are removed by the testing function prior to using those files. The records of special interest are listed on top, and their peculiarities are documented. This allows to trace the fate
of each such record from input to full dataframe to the collapsed final output.

In addition to the hybrid integration test, the code of the pipeline itself performs sanity checks whenever possible.