# Filtering functions that can be used in multiple pipelines.

# Identified as problematic submissions, e.g. too many unmappable trait names.
submission_names_to_exclude = ['SUB14299258', 'SUB14767656']


def filter_by_submission_name(clinvar_set):
    """Return False (i.e. filter out) if every submitted record in the set has submission_name in the exclusion list."""
    for submitted_record in clinvar_set.scvs:
        if submitted_record.submission_name not in submission_names_to_exclude:
            return True
    return False
