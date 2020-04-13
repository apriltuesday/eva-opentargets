#!/bin/bash
# A script to compare evidence strings. Please see README.md for details on using it.

# A function to sort key in evidence strings. This makes them more readable and helps comparison through word diff.
# Also removes the "validated_aginst_schema_version" field entirely, because it frequently changes between the
# versions, and this change is not important.
# The two arguments are input and output JSON files.
sort_keys () {
  ./jq -S "." --tab <"$1" \
    | tr -d '\t\n' \
    | sed -e 's|}{|}~{|g' \
    | tr '~' '\n' \
    | sed -e 's|,"validated_against_schema_version": "[0-9.]*"||g' \
  > "$2"
}

# A function to extract some fields from the evidence strings. The fields being extracted are:
# * ClinVar RCV accession
# * Variant ID (rsID or, if absent, RCV accession)
# * Ensembl gene ID
# * Functional consequence SO code
extract_fields () {
  ./jq '
    .unique_association_fields.clinvarAccession + "|" +
    .unique_association_fields.variant_id + "|" +
    .unique_association_fields.gene + "|" +
    .evidence.gene2variant.functional_consequence
  ' < "$1" \
  | tr -d '"' | tr '|' '\t' \
  | sed -e 's|http://purl.obolibrary.org/obo/||g' \
        -e 's|http://identifiers.org/clinvar.record/||g' \
        -e 's|http://identifiers.org/dbsnp/||g' \
        -e 's|http://www.ncbi.nlm.nih.gov/clinvar/||g' \
  > "$2"
}

echo "Set up environment & parse parameters"
export -f sort_keys extract_fields
# To ensure that the sort results are consistent, set the sort order locale directly
export LC_COLLATE=C
# The realpath is required to make the paths work after the working directory change
OLD_EVIDENCE_STRINGS=$(realpath "$1")
NEW_EVIDENCE_STRINGS=$(realpath "$2")
mkdir comparison && cd comparison || exit 1

echo "Install JQ — a command line JSON processor"
wget -q -O jq https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64
chmod a+x jq

echo "Install aha — HTML report generator"
wget -q https://github.com/theZiz/aha/archive/0.5.zip
unzip -q 0.5.zip && cd aha-0.5 && make && mv aha ../ && cd .. && rm -rf aha-0.5 0.5.zip

echo "Sort keys and remove schema version from the evidence strings"
sort_keys "${OLD_EVIDENCE_STRINGS}" 01.keys-sorted.old.json \
  & sort_keys "${NEW_EVIDENCE_STRINGS}" 01.keys-sorted.new.json \
  & wait

echo "Extract some fields to pair old and new evidence strings together"
extract_fields 01.keys-sorted.old.json 02.fields.old \
  & extract_fields 01.keys-sorted.new.json 02.fields.new \
  & wait

echo "Paste fields & original strings into the same table and sort"
paste 02.fields.old 01.keys-sorted.old.json | sort -u > 03.fields-and-strings.old \
  & paste 02.fields.new 01.keys-sorted.new.json | sort -u > 03.fields-and-strings.new \
  & wait

# Comparing functional consequences *correctly* will require investigating the nature of relationships (one-to-one,
# one-to-many etc.) for all ClinVar objects. This will be done in subsequent tickets. Until then, functional consequence
# calculations are disabled so as not to be misleading.

#echo "Extract all unique variants"
#cut -f1-2 03.fields-and-strings.old | sort -u > 04.all-variants.old \
#  & cut -f1-2 03.fields-and-strings.new | sort -u > 04.all-variants.new \
#  & wait
#
#echo "Group variants into deleted, new and common"
#comm -23 04.all-variants.old 04.all-variants.new > 05.variants.deleted \
#  & comm -13 04.all-variants.old 04.all-variants.new > 05.variants.added \
#  & comm -12 04.all-variants.old 04.all-variants.new > 05.variants.common \
#  & wait
#
#echo "Find variants which appear in both files but where functional mappings have changed"
#cut -f1-3 03.fields-and-strings.old | sort -u > 06.consequences.old \
#  & cut -f1-3 03.fields-and-strings.new | sort -u > 06.consequences.new \
#  & wait
#join 3_common old.consequences -j 1 | join /dev/stdin new.consequences -j1 | awk '$2 != $3' > 3_common_changed
#
#There are $(wc -l 1_deleted) records present only in the first file:
#$(sed 's|^|  |' 1_deleted)
#
#There are $(wc -l 2_added) records present only in the second file:
#$(sed 's|^|  |' 2_added)
#
#There are also $(wc -l 3_common) records present in both files.
#Of them, $(wc -l 3_common_changed) records have different consequences between the first and the second files:
#$(sed 's|^|  |' 3_common_changed)

echo "Compute difference"
git diff --minimal -U0 --color=always --word-diff=color old.sorted new.sorted &> "diff"

echo "Producing the report"
cat << EOF > report
Compared:
File 1, ${OLD_EVIDENCE_STRINGS}, contained $(wc -l 03.fields-and-strings.old) records
File 2, ${NEW_EVIDENCE_STRINGS}, contained $(wc -l 03.fields-and-strings.new) records

The full diff between two files follows.

$(cat diff)
EOF

./aha --word-wrap < report > report.html
cd ..
