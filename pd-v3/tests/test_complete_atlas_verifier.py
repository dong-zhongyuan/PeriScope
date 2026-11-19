"""Defects the complete-atlas verifier must detect before result delivery."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'complete_verifier', ROOT/'papers/01_mr_imaging_scan/verify_complete_atlas.py')
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class CompleteAtlasVerifierTests(unittest.TestCase):
    def test_missing_hypotheses_keep_their_planned_denominator(self):
        p = np.array([.01, .04, np.nan])
        np.testing.assert_allclose(V.external_adjust(p, 'fdr_bh'), [.03, .06, np.nan], equal_nan=True)
        np.testing.assert_allclose(V.external_adjust(p, 'fdr_by'), [.055, .11, np.nan], equal_nan=True)
        with self.assertRaises(AssertionError):
            V.close([.02, .04, np.nan], V.external_adjust(p, 'fdr_bh'), 'dropped missing slot', True)

    def test_storey_scales_missing_slots_before_restoring_na(self):
        p = np.array([.001, .002, .003, .004, .005, .006, .007, .008, 1., np.nan])
        # lambda=.5 gives pi0=(1+2)/(10*.5)=.6 for this declared family.
        np.testing.assert_allclose(V.external_adjust(p, 'fdr_bh', .6),
                                   [.006]*8 + [.6, np.nan], equal_nan=True)

    def test_missingness_and_small_probability_errors_are_detected(self):
        with self.assertRaises(ValueError):
            V.close([0.], [np.nan], 'invented estimate', True)
        with self.assertRaises(AssertionError):
            V.close([2e-40], [1e-40], 'wrong small P', True)

    def test_all_four_nonpalindrome_orientations(self):
        rows = []
        for ya, yo, sign in [('A', 'C', 1), ('C', 'A', -1), ('T', 'G', 1), ('G', 'T', -1)]:
            rows.append(dict(ea_x='A', oa_x='C', ea_y=ya, oa_y=yo, beta_y=2.,
                             beta_y_aligned=2.*sign, eaf_x=.2, eaf_y=.2,
                             eaf_y_aligned=.2 if sign == 1 else .8, allele_flipped=sign == -1))
        data = pd.DataFrame(rows)
        self.assertEqual(V.independent_alignment(data, 'drop'), 0)
        data.loc[1, 'beta_y_aligned'] = 2.
        with self.assertRaises(AssertionError):
            V.independent_alignment(data, 'drop')

    def test_palindrome_policy_and_frequency_ambiguity(self):
        data = pd.DataFrame([dict(ea_x='A', oa_x='T', ea_y='A', oa_y='T',
                                  beta_y=2., beta_y_aligned=2., eaf_x=.2, eaf_y=.2,
                                  eaf_y_aligned=.2, allele_flipped=False)])
        self.assertEqual(V.independent_alignment(data, 'infer'), 1)
        with self.assertRaises(ValueError):
            V.independent_alignment(data, 'drop')
        data[['eaf_x', 'eaf_y']] = .49
        with self.assertRaises(ValueError):
            V.independent_alignment(data, 'infer')

    @staticmethod
    def pair():
        a = pd.DataFrame(dict(exposure_id=['a', 'b', 'c'], phenotype=['x']*3,
                              status=['estimated', 'palindromic_dropped', 'missing_outcome']))
        for col in V.NUMERIC:
            a[col] = [.25, np.nan, np.nan]
        b = a.copy()
        b.loc[1, 'status'] = 'estimated'
        b.loc[1, V.NUMERIC] = .25
        return a, b

    def test_policy_changes_only_recover_palindromes(self):
        a, b = self.pair()
        self.assertEqual(sum(V.policy_pair(a, b).values()), 3)
        b.loc[2, 'status'] = 'estimated'
        with self.assertRaises(ValueError):
            V.policy_pair(a, b)

    def test_original_estimates_and_keys_must_be_preserved(self):
        a, b = self.pair()
        b.loc[0, 'beta'] = .5
        with self.assertRaises(AssertionError):
            V.policy_pair(a, b)
        a, b = self.pair()
        b.loc[0, 'exposure_id'] = 'other'
        with self.assertRaises(ValueError):
            V.policy_pair(a, b)

    def test_full_completion_guard_before_output_creation(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            output = folder/'new_outputs/report.json'
            with self.assertRaises(ValueError):
                V.verify(folder, 'full', output)
            self.assertFalse(output.parent.exists())
            run = dict(complete=True, status='complete', mode='full', assays=2941,
                       phenotypes=3963, completed_traits=3963, planned_pairs=11655183,
                       policies=dict(drop={}, infer={}))
            for key, bad in [('complete', False), ('completed_traits', 3962), ('planned_pairs', 11655182),
                             ('assays', 1954), ('mode', 'validation45'), ('policies', dict(drop={}))]:
                (folder/'RUN.json').write_text(json.dumps({**run, key: bad}))
                with self.assertRaises(ValueError):
                    V.completed_run(folder, 'full')

    def test_input_hash_mismatch(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as temp:
            path = Path(temp)/'input'
            path.write_bytes(b'original')
            expected = V.sha(path)
            V.checked(path, expected, {})
            path.write_bytes(b'changed')
            with self.assertRaises(ValueError):
                V.checked(path, expected, {})


if __name__ == '__main__':
    unittest.main()
