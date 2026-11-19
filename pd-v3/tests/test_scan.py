import unittest
import numpy as np
import pandas as pd
from pdv3.scan import scan_pairs


def pair(**changes):
    row = dict(gene='A', assay='A1', variant='rs1', phenotype='MRI1', family='volume',
               ea_x='A', oa_x='G', ea_y='G', oa_y='A', beta_x=.5, se_x=.01,
               beta_y=-.05, se_y=.01, build_x='GRCh38', build_y='GRCh38',
               chrom_x=1, chrom_y=1, pos_x=100, pos_y=100)
    return dict(row, **changes)


class ScanTests(unittest.TestCase):
    def test_missing_outcome_remains_in_global_denominator(self):
        full, audit = scan_pairs(pd.DataFrame([pair(), pair(phenotype='MRI2', beta_y=np.nan)]))
        self.assertEqual(len(full), 2)
        self.assertEqual(len(audit), 2)
        self.assertAlmostEqual(full.iloc[0].beta, .1)
        self.assertEqual(full.iloc[1].status, 'no_usable_instrument')
        self.assertEqual(audit.iloc[1].status, 'missing_outcome')
        self.assertAlmostEqual(full.iloc[0].q_bh_global, 2 * full.iloc[0].p)
        self.assertTrue(np.isnan(full.iloc[1].q_bh_global))

    def test_coordinates_and_builds_cannot_be_silently_joined(self):
        for changes, reason in [({'pos_y': 101}, 'coordinate_mismatch'),
                                ({'build_y': 'GRCh37'}, 'build_mismatch_or_unknown')]:
            full, audit = scan_pairs(pd.DataFrame([pair(**changes)]))
            self.assertEqual(audit.iloc[0].status, reason)
            self.assertEqual(full.iloc[0].n_instruments, 0)

    def test_duplicate_variants_fail_before_pooling(self):
        with self.assertRaises(ValueError):
            scan_pairs(pd.DataFrame([pair(), pair()]), independent=True)

    def test_ivw_requires_explicit_independence(self):
        pairs = pd.DataFrame([pair(), pair(variant='rs2', pos_x=200, pos_y=200)])
        full, audit = scan_pairs(pairs)
        self.assertEqual(full.iloc[0].status, 'independence_not_attested')
        self.assertTrue(audit.usable.all())
        full, _ = scan_pairs(pairs, independent=True)
        self.assertEqual(full.iloc[0].method, 'IVW_independent')
        self.assertAlmostEqual(full.iloc[0].beta, .1)

    def test_assay_mixing_and_truthy_strings_fail(self):
        with self.assertRaises(ValueError):
            scan_pairs(pd.DataFrame([pair()]), independent='False')
        with self.assertRaises(ValueError):
            scan_pairs(pd.DataFrame([pair(), pair(variant='rs2', assay='A2')]), independent=True)

    def test_incomplete_family_plan_fails(self):
        with self.assertRaises(ValueError):
            scan_pairs(pd.DataFrame([pair(), pair(gene='B', assay='B1', phenotype='MRI2')]))


if __name__ == '__main__':
    unittest.main()
