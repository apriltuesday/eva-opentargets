#!/bin/bash
# A script to compare evidence strings. Please see README.md for details on using it.

# A function to sort key in evidence strings. This makes them more readable and helps comparison through word diff.
# The two arguments are input and output JSON files.
sort_keys () {
  ./jq -S "." --tab <"$1" \
    | tr -d '\t\n' \
    | sed -e 's|}{|}~{|g' \
    | tr '~' '\n' \
  > "$2"
}

# A function to extract some fields from the evidence strings. The fields being extracted are:
# * ClinVar RCV accession
# * Variant ID (rsID or, if absent, RCV accession)
# * Functional consequence SO code
extract_fields () {
  ./jq '.evidence.variant2disease.provenance_type.database.dbxref.url + ">" + .variant.id + "|" + .evidence.gene2variant.functional_consequence' <"$1" \
  | tr -d '"' | tr '|' '\t' \
  | sed -e 's|http://purl.obolibrary.org/obo/||g' \
        -e 's|http://identifiers.org/clinvar.record/||g' \
        -e 's|http://identifiers.org/dbsnp/||g' \
        -e 's|http://www.ncbi.nlm.nih.gov/clinvar/||g' \
  > "$1.fields"
}

echo "Set up environment & parse parameters"
export -f sort_keys extract_fields
# To ensure that the sort results are consistent, set the sort order locale directly
export LC_COLLATE=C
export OLD_EVIDENCE_STRINGS="$1"
export NEW_EVIDENCE_STRINGS="$2"
mkdir comparison
ln -s "${OLD_EVIDENCE_STRINGS}" comparison/old.input.json
ln -s "${NEW_EVIDENCE_STRINGS}" comparison/new.input.json
cd comparison || exit 1

echo "Install JQ — a command line JSON processor"
rm -rf jq
wget -q -O jq https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64
chmod a+x jq

echo "Sort keys to make comparison easier"
sort_keys old.input.json old.json \
  & sort_keys new.input.json new.json \
  & wait

echo "Extract some fields to pair old and new evidence strings together"
extract_fields old.json \
  & extract_fields new.json \
  & wait

echo "Paste fields & original strings into the same document"
paste old.json.fields old.json > old \
  & paste new.json.fields new.json > new \
  & wait

echo "Classify records into (disappeared; new; common)"
cut -f1 old.json.fields | sort -u > old.variants \
  & cut -f1 new.json.fields | sort -u > new.variants \
  & wait
comm -23 old.variants new.variants > 1_deleted \
  & comm -13 old.variants new.variants > 2_added \
  & comm -12 old.variants new.variants > 3_common \
  & wait

echo "Find variants which appear in both files but where functional mappings have changed"
cut -f1-2 old.json.fields | sort -u > old.consequences \
  & cut -f1-2 new.json.fields | sort -u > new.consequences \
  & wait
join 3_common old.consequences -j 1 | join /dev/stdin new.consequences -j1 | awk '$2 != $3' > 3_common_changed

echo "Sort the results"
sort -u old > old.sorted \
  & sort -u new > new.sorted \
  & wait

echo "Compute difference"
git diff --minimal -U0 --color=always --word-diff=color old.sorted new.sorted &> "diff"

cd ..