class SetComparisonMetrics:
    """Records counts and scores for comparing sets of arbitrary values from two sources (ClinVar and CMAT)."""

    match_keys = [
        'exact_match',
        'cmat_superset',
        'cmat_subset',
        'divergent_match'
    ]
    both_present_keys = match_keys + ['mismatch']
    some_missing_keys = [
        'cv_missing',
        'cmat_missing',
        'both_missing'
    ]
    all_keys = both_present_keys + some_missing_keys

    def __init__(self):
        self.counts = {k: 0 for k in self.all_keys}
        self.scores = {k: 0 for k in self.all_keys}
        # To keep counts in self.counts disjoint, add a separate field (calculated last) for both_present and match
        self.match_count = 0
        self.match_score = 0
        self.both_present_count = 0
        self.both_present_score = 0

    def count_and_score(self, cv_set, cmat_set):
        cv_set = set(cv_set)
        cmat_set = set(cmat_set)

        # First check if either set is empty - if so we count these but don't compute f-scores, as we do not assume
        # either ClinVar or CMAT has complete coverage of the data.
        if not cv_set and cmat_set:
            self.counts['cv_missing'] += 1
        elif cv_set and not cmat_set:
            self.counts['cmat_missing'] += 1
        elif not cv_set and not cmat_set:
            self.counts['both_missing'] += 1

        # Both sets present, compute the score and place in correct category.
        else:
            f1_score, true_pos, false_pos, false_neg = self.get_f1_tp_fp_fn(cv_set, cmat_set)
            if false_pos and not false_neg:
                k = 'cmat_superset'
            elif not false_pos and false_neg:
                k = 'cmat_subset'
            elif not false_pos and not false_neg:
                k = 'exact_match'
            elif true_pos:
                k = 'divergent_match'
            else:
                k = 'mismatch'
            self.counts[k] += 1
            self.scores[k] += f1_score

    def finalise(self):
        """Averages scores over counts, and also calculate both_present count & score."""
        self.match_count = sum(self.counts[k] for k in self.match_keys)
        if self.match_count:
            self.match_score = sum(self.scores[k] for k in self.match_keys) / self.match_count
        self.both_present_count = sum(self.counts[k] for k in self.both_present_keys)
        if self.both_present_count:
            self.both_present_score = sum(self.scores[k] for k in self.both_present_keys) / self.both_present_count
        for k in self.all_keys:
            self.scores[k] = self.scores[k] / self.counts[k] if self.counts[k] else 0

    def report(self):
        total = sum(v for v in self.counts.values())
        print(f'Total = {total}')
        self.pretty_print(
            ('Category', 'Count', 'Percent', 'F1 Score'),
            [(k, self.counts[k], f'{self.counts[k] / total:.1%}', f'{self.scores[k]:.2f}') for k in self.both_present_keys]
            + [('-->match', self.match_count, f'{self.match_count / total:.1%}', f'{self.match_score:.2f}')]
            + [('-->both_present', self.both_present_count, f'{self.both_present_count / total:.1%}', f'{self.both_present_score:.2f}')]
            + [(k, self.counts[k], f'{self.counts[k] / total:.1%}', f'{self.scores[k]:.2f}') for k in self.some_missing_keys]
        )

    @staticmethod
    def get_f1_tp_fp_fn(cv_set, cmat_set):
        """Returns f1-score and number of true positives, false positives and false negatives."""
        if len(cv_set) == 0 and len(cmat_set) == 0:
            return 0, 0, 0, 0
        tp = len(cmat_set & cv_set)
        fp = len(cmat_set - cv_set)
        fn = len(cv_set - cmat_set)
        return 2 * tp / (2 * tp + fp + fn), tp, fp, fn

    @staticmethod
    def pretty_print(header, table):
        cell_widths = [len(h) for h in header]
        for row in table:
            for i, cell in enumerate(row):
                cell_widths[i] = max(cell_widths[i], len(str(cell)))
        format_string = '  '.join('{%s:>%s}' % (i, w) for i, w in enumerate(cell_widths))
        print(' ' + format_string.format(*header) + ' ')
        for row in table:
            print(' ' + format_string.format(*row) + ' ')
