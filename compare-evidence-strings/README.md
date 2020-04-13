# Comparison of evidence strings

When major updates to the pipeline are implemented, an important measure of control is comparing evidence strings before and after the update for the same input data. This protocol contains commands which help do this.

Compare the evidence strings with the script:
```bash
bash compare.sh \
  old_evidence_strings.json \
  new_evidence_strings.json
```

The script should take a few minutes to run and will create the `comparison/` subdirectory in the current working directory. It will contain several files:
* `old.json` and `new.json` — evidence strings with their keys sorted lexicographically
* `old.json.fields` and `new.json.fields` — TSV files with fields extracted from the evidence strings
  + Column 1: RCV accession and the associated variant, e. g. `RCV000162096>rs730882195`. In case the record is addressed using its RCV ID (in case of repeat expansion variants), it will look like `RCV000075958>RCV000075958`.
  + Column 2: functional consequence code, e. g. SO_0001583
* `old` and `new` — previous two files pasted together; extracted fields plus the original evidence string with the sorted keys
* `1_deleted`, `2_added` and `3_common` — the variants which were deleted or added compared to the old evidence strings, or which appear in both files
* `3_common_changed` — a subsection of `3_common` where the variant is present in both old and new evidence strings, but its functional consequence has changed
* `old.sorted` and `new.sorted` — same as `old` and `new`, but this time sorted for comparison
* `diff` — the difference between `old.sorted` and `new.sorted`, should be viewed using `less -r`

## Future improvements
There is a [json-diff](https://pypi.org/project/json-diff/) module which allows detailed comparison of JSON objects. If this protocol is going to be updated in the future, this module might be helpful. It provides structured overview of differences; however, it has a few limitations:
 * It can only compare individual evidence strings (so they must be sorted and matched beforehand)
 * When a field's value is updated, `json-diff` only reports the new value of the field, but not the old one, for example:
```json
{
    "_update": {
        "evidence": {
            "_update": {
                "gene2variant": {
                    "_update": {
                        "functional_consequence": "http://purl.obolibrary.org/obo/SO_0001575"
                    }
                }
            }
        }
    }
}
```

Here, the change was from SO_0001589 to SO_0001575, but only the second value is reported.