"""Boundary cases for streaming precision verification, not a new analysis."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'complete_precision', ROOT/'papers/01_mr_imaging_scan/verify_complete_precision.py')
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


def empty_trait():
    a = pd.DataFrame(dict(exposure_id=['a', 'b'], phenotype=['trait']*2,
                          family=['family']*2, status=['missing_outcome', 'no_cis_instrument']))
    for col in V.SOURCE_COLUMNS[4:]:
        a[col] = np.nan
    b = a[V.KEYS].copy()
    cols = ['conditional_se', 'MDE80_nominal', 'MDE90_nominal',
            'MDE80_full_plan_bonferroni_reference', 'MDE90_full_plan_bonferroni_reference',
            'Fieller95_rho0_lower', 'Fieller95_rho0_upper', 'Fieller95_grid_lower_min',
            'Fieller95_grid_upper_max', 'delta_SE_grid_min', 'delta_SE_grid_max']
    for col in cols:
        b[col] = np.nan
    # Distinct pandas labels must not displace already-verified scientific rows.
    a.index = [101, 102]
    b.index = [901, 902]
    mde = []
    for reference, alpha in V.ALPHAS.items():
        row = dict(policy='drop', phenotype='trait', family='family', planned=2,
                   estimated=0, reference=reference, alpha=alpha)
        row.update({f'MDE{target}_{q}': np.nan for target in (80, 90) for q in ('p05', 'median', 'p95')})
        mde.append(row)
    cov = []
    for rho in V.RHOS:
        row = dict(policy='drop', phenotype='trait', planned=2, estimated=0,
                   rho_error=rho, intervals_including_zero=0)
        row.update({f'{name}_{q}': np.nan for name in ('Fieller_width_ratio_to_rho0', 'delta_SE_ratio_to_rho0')
                    for q in ('p05', 'median', 'p95')})
        cov.append(row)
    return (a, b, pd.DataFrame(mde).set_index(['policy', 'phenotype', 'reference']),
            pd.DataFrame(cov).set_index(['policy', 'phenotype', 'rho_error']))


class CompletePrecisionVerifierTests(unittest.TestCase):
    def test_all_missing_trait_with_different_reader_indices(self):
        a, b, mde, cov = empty_trait()
        V.paired_keys(a, b, dict(phenotype='trait', family='family'), ['a', 'b'])
        report = V.verify_block(a, b, 'drop', mde, cov)
        self.assertEqual(report['estimated'], 0)
        self.assertEqual(report['summary_max_abs_difference'], 0.)

    def test_unestimated_slot_cannot_be_filled_with_zero(self):
        a, b, mde, cov = empty_trait()
        b.loc[901, 'conditional_se'] = 0.
        with self.assertRaises(ValueError):
            V.verify_block(a, b, 'drop', mde, cov)

    def test_duplicate_or_reordered_assays_rejected(self):
        a, b, _, _ = empty_trait()
        b.exposure_id = ['a', 'a']
        with self.assertRaises(ValueError):
            V.paired_keys(a, b, dict(phenotype='trait', family='family'), ['a', 'b'])

    def test_weak_denominator_is_not_given_a_finite_fieller_interval(self):
        with self.assertRaises(ValueError):
            V.quadratic(np.array([1.]), np.array([1.]), np.array([0.]), np.array([1.]), 0.)
        lo, hi, delta, residual = V.quadratic(np.array([1.]), np.array([.1]), np.array([0.]), np.array([1.]), 0.)
        self.assertLess(lo[0], 0.)
        self.assertGreater(hi[0], 0.)
        self.assertEqual(delta[0], 1.)
        self.assertLess(residual, 1e-12)

    def test_incomplete_or_changed_scientific_scope_rejected(self):
        audit = dict(complete=True, mode='full', original_full_plan_pairs=V.TOTAL,
                     analysed_pairs_per_policy=V.TOTAL, alpha_references=V.ALPHAS,
                     rho_error_grid=V.RHOS, policies=dict(drop={}, infer={}),
                     actual_overlap_estimated=False, overlap_bias_corrected=False, new_MR_discovery=False)
        self.assertEqual(V.completed_audit(audit, 'full'), 3963)
        for field, value in [('complete', False), ('analysed_pairs_per_policy', 132345),
                             ('rho_error_grid', [0.]), ('alpha_references', {'nominal': .05}),
                             ('actual_overlap_estimated', True), ('overlap_bias_corrected', True)]:
            with self.subTest(field=field), self.assertRaises(ValueError):
                V.completed_audit({**audit, field: value}, 'full')

    def test_unavailable_full_result_creates_no_output(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(V, 'BASE', Path(temp)):
            output = Path(temp)/'new/report.json'
            with self.assertRaises(ValueError):
                V.verify('full', output)
            self.assertFalse(output.parent.exists())


if __name__ == '__main__':
    unittest.main()
