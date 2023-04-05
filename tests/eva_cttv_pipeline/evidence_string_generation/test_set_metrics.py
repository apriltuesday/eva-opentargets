from eva_cttv_pipeline.evidence_string_generation.set_metrics import SetComparisonMetrics


# Tolerance for float comparisons
EPSILON = 1e-6


def test_f1_score():
    assert abs(SetComparisonMetrics.get_f1_tp_fp_fn({1, 2}, {1})[0] - 2 / 3) < EPSILON
    assert abs(SetComparisonMetrics.get_f1_tp_fp_fn({1, 2}, {2, 1})[0] - 1) < EPSILON
    assert SetComparisonMetrics.get_f1_tp_fp_fn({}, {})[0] == 0


def test_count_and_score():
    metrics = SetComparisonMetrics()
    metrics.count_and_score({}, {})
    metrics.count_and_score({1, 2}, {1})
    metrics.count_and_score({1, 2}, {2, 1})
    metrics.count_and_score({1, 2}, {2, 3})
    metrics.count_and_score({1, 2}, {2, 3, 4})
    metrics.count_and_score({1}, {2})
    metrics.count_and_score([], {2, 3})
    metrics.finalise()

    assert metrics.counts == {
        'exact_match': 1,
        'cmat_superset': 0,
        'cmat_subset': 1,
        'divergent_match': 2,
        'mismatch': 1,
        'cv_missing': 1,
        'cmat_missing': 0,
        'both_missing': 1,
    }

    assert metrics.scores['exact_match'] == 1
    assert metrics.scores['mismatch'] == 0
    # mean(2/(2+1+1), 2/(2+2+1)) = mean(0.5, 0.4) = 0.45
    assert abs(metrics.scores['divergent_match'] - 0.45) < EPSILON
