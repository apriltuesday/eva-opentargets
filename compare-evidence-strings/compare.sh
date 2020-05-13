#!/bin/bash
# A script to compare evidence strings. Please see README.md for details on using it.

# A function to sort keys in the evidence strings. This makes them more readable and helps comparison through word diff.
# Also removes several version- and date-related fields which are changed frequently, but do not reflect actual change:
# * validated_aginst_schema_version
# * date_asserted
# The two arguments are input and output JSON files.
sort_keys () {
  ./jq -S "." --tab <"$1" \
    | tr -d '\t\n' \
    | sed -e 's|}{|}~{|g' \
    | tr '~' '\n' \
    | sed -e 's|,"validated_against_schema_version": "[0-9.]*"||g' \
    | sed -e 's|"date_asserted": ".\{19\}"||g' \
  > "$2"
}

# A function to extract all unique identifying fields from the evidence strings. The fields being extracted are:
# * ClinVar RCV accession
# * Phenotype
# * Allele origin
# * Variant ID (rsID or, if absent, RCV accession)
# * Ensembl gene ID
extract_fields () {
  ./jq '
    .unique_association_fields.clinvarAccession + "|" +
    .unique_association_fields.phenotype + "|" +
    .unique_association_fields.alleleOrigin + "|" +
    .unique_association_fields.variant_id + "|" +
    .unique_association_fields.gene
  ' < "$1" | tr -d '"' > "$2"
}

# Computes a word diff between two files using git diff
compute_git_diff () {
  # The --no-index option is important, because otherwise git will refuse to compare the files if you're running this
  # script from right inside the repository (because the files are untracked).
  git diff \
  --minimal \
  -U0 \
  --color=always \
  --word-diff=color \
  --no-index \
  "$1" "$2"
}

echo "Set up environment & parse parameters"
export -f sort_keys extract_fields compute_git_diff
# To ensure that the sort results are consistent, set the sort order locale explicitly
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
unzip -q 0.5.zip && cd aha-0.5 && make &>/dev/null && mv aha ../ && cd .. && rm -rf aha-0.5 0.5.zip

echo "Sort keys and remove schema version from the evidence strings"
sort_keys "${OLD_EVIDENCE_STRINGS}" 01.keys-sorted.old.json \
  & sort_keys "${NEW_EVIDENCE_STRINGS}" 01.keys-sorted.new.json \
  & wait

echo "Extract the unique identifying fields to pair old and new evidence strings together"
extract_fields 01.keys-sorted.old.json 02.fields.old \
  & extract_fields 01.keys-sorted.new.json 02.fields.new \
  & wait

echo "Paste the unique identifying fields & original strings into the same table and sort"
paste 02.fields.old 01.keys-sorted.old.json | sort -k1,1 > 03.fields-and-strings.old \
  & paste 02.fields.new 01.keys-sorted.new.json | sort -k1,1 > 03.fields-and-strings.new \
  & wait

echo "Compute sets of all non-unique identifying fields in each evidence string set"
cut -f1 03.fields-and-strings.old | uniq -c | awk '$1>1 {print $2}' > 04.non-unique-fields.old \
  & cut -f1 03.fields-and-strings.new | uniq -c | awk '$1>1 {print $2}' > 04.non-unique-fields.new \
  & wait
cat 04.non-unique-fields.old 04.non-unique-fields.new | sort -u > 05.all-non-unique-fields

echo "Extract evidence strings with *non-unique* identifying fields into a separate group"
join -j 1 05.all-non-unique-fields 03.fields-and-strings.old > 06.non-unique.old \
  & join -j 1 05.all-non-unique-fields 03.fields-and-strings.new > 06.non-unique.new \
  & wait

echo "Extract evidence strings with *unique* identifying fields into a separate group"
# -v 2 means "print only records from file #2 which cannot be paired"
# If records cannot be paired, it means their identifying fields are *not* in the list of duplicates
join -j 1 -v 2 05.all-non-unique-fields 03.fields-and-strings.old > 07.unique.old \
  & join -j 1 -v 2 05.all-non-unique-fields 03.fields-and-strings.new > 07.unique.new \
  & wait

echo "Separate unique evidence strings into (deleted, common, new)"
join -j 1 07.unique.old 07.unique.new > 08.common \
  & join -j 1 -v 1 07.unique.old 07.unique.new > 08.deleted \
  & join -j 1 -v 2 07.unique.old 07.unique.new > 08.added \
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
DIFF_FILE="diff.$(date +'%Y%m%d%H%M%S')"
export DIFF_FILE
# The --no-index option is important, because otherwise git will refuse to compare the files if you're running this
# script from right inside the repository (because the files are untracked).
git diff \
  --minimal \
  -U0 \
  --color=always \
  --word-diff=color \
  --no-index \
  03.fields-and-strings.old \
  03.fields-and-strings.new &> "${DIFF_FILE}"

echo "Producing the report"
cat << EOF > report
Compared:

File 1
${OLD_EVIDENCE_STRINGS}
Total unique records: $(wc -l <03.fields-and-strings.old)

File 2
${NEW_EVIDENCE_STRINGS}
Total unique records: $(wc -l <03.fields-and-strings.new)

The full diff between two files follows.

$(tail -n+5 "${DIFF_FILE}" | awk '{if ($0 !~ /@@/) {print $0 "\n"}}')
EOF

./aha --word-wrap < report > report.html

cd ..
