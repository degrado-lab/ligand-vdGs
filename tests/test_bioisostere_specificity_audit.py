"""Falsifiers for matched-case assignment and held-out rank artifacts."""
import csv
import json
import os
import tempfile
import unittest

from scripts.bioisostere_specificity_audit import CASES, validate
from tests.vacuity import assert_discriminates

class TestSpecificityAudit(unittest.TestCase):
    def test_validator_refuses_swapped_head_and_bad_percentile(self):
        with tempfile.TemporaryDirectory() as root:
            rows = [dict(case=name, label=label, query_fold=str(fold), query_key=a, target_key=b,
                         query_support='100', target_support='100', pearson_r='.8',
                         rank_query_to_target='2', pool_query_to_target='100', pct_query_to_target='.02',
                         rank_target_to_query='3', pool_target_to_query='100', pct_target_to_query='.03',
                         reciprocal_top5='True')
                    for name, label, a, b, *_ in CASES for fold in (0, 1)]
            with open(os.path.join(root, 'summary.json'), 'w') as handle:
                json.dump(dict(excluded_entries=sorted({pdb for case in CASES for pdb in case[4:]}),
                               entry_counts=[2, 2], profile_counts=[2, 2]), handle)
            def write():
                with open(os.path.join(root, 'rows.tsv'), 'w') as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
                    writer.writeheader()
                    writer.writerows(rows)
            write()
            self.assertEqual(len(validate(root)), 2 * len(CASES))
            kmo = next(row for row in rows if row['case'] == 'KMO' and row['query_fold'] == '0')
            correct = (kmo['query_key'], kmo['target_key'])
            kmo['target_key'] = next(case[3] for case in CASES if case[0] == 'TTR_role_control')
            assert_discriminates(lambda pair: pair == correct, accepts=[correct],
                                 rejects=[(kmo['query_key'], kmo['target_key'])],
                                 label='target head must remain assigned to KMO')
            write()
            self.assertRaisesRegex(ValueError, 'wrong heads', validate, root)
            kmo['target_key'] = correct[1]
            kmo['pct_query_to_target'] = '.01'
            write()
            self.assertRaisesRegex(ValueError, 'rank or percentile', validate, root)

if __name__ == '__main__': unittest.main()
