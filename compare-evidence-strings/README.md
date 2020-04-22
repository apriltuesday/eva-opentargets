# Comparison of evidence strings

When major updates to the pipeline are implemented, an important measure of control is comparing evidence strings before and after the update for the same input data. This protocol contains commands which help do this.

Compare the evidence strings with the script:
```bash
bash compare.sh \
  old_evidence_strings.json \
  new_evidence_strings.json
```

The script should take a few minutes to run and will create the `comparison/` subdirectory in the current working directory. It will contain several files, the most important of which is `report.html`, which can be viewed in a browser.

## Future improvements

### Alternative library for producing diffs
The [diff2html-cli](https://github.com/rtfpessoa/diff2html-cli) is a more advanced library which can be used to replace `aha` in the future.

### json-diff
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