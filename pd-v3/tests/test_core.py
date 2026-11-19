"""Arithmetic, alignment and full-family regression tests, independent of real data."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from pdv3.mr import harmonize, wald, ivw, variance_explained
from pdv3.multiplicity import adjust_pvalues, adjust_table, storey_qvalues, two_stage_bh
from pdv3.sumstats import validate_sumstats, read_sumstats
from pdv3.findings import rank_effects, direction_summary
from pdv3.permute import null_calibration


class HarmonizationTests(unittest.TestCase):
    def test_four_nonpalindromic_orientations(self):
        for ea, oa, by, expected in [("A", "G", 0.2, "direct"),
                                     ("G", "A", -0.2, "swap"),
                                     ("T", "C", 0.2, "strand"),
                                     ("C", "T", -0.2, "strand_swap")]:
            r = harmonize("A", "G", ea, oa, by, .2, .8 if by < 0 else .2,
                          exposure_build="hg38", outcome_build="GRCh38")
            self.assertTrue(r["usable"])
            self.assertEqual(r["status"], expected)
            self.assertAlmostEqual(r["beta_y"], .2)
            self.assertAlmostEqual(r["eaf_y"], .2)
            self.assertEqual(r["build_status"], "matched")

    def test_palindromes_require_explicit_frequency_inference(self):
        self.assertFalse(harmonize("A", "T", "A", "T", .2)["usable"])
        inferred = harmonize("A", "T", "A", "T", .2, .1, .9, palindrome="infer")
        self.assertTrue(inferred["usable"])
        self.assertAlmostEqual(inferred["beta_y"], -.2)
        self.assertAlmostEqual(inferred["eaf_y"], .1)
        for ex, ey in [(None, .1), (.5, .5), (.1, .3)]:
            self.assertFalse(harmonize("A", "T", "A", "T", .2, ex, ey,
                                       palindrome="infer")["usable"])

    def test_mismatches_and_invalid_data_are_auditable(self):
        self.assertEqual(harmonize("A", "G", "A", "C", .2)["status"], "allele_mismatch")
        self.assertEqual(harmonize("A", "G", "A", "G", .2, exposure_build=37,
                                   outcome_build=38)["status"], "build_mismatch")
        self.assertEqual(harmonize("A", "G", "A", "G", .2, exposure_build="foo")["status"], "unknown_build")
        self.assertEqual(harmonize("A", "G", "A", "G", np.nan)["status"], "invalid_beta")
        self.assertEqual(harmonize("A", "G", "A", "G", .2, eaf_y=1.5)["status"], "invalid_eaf")
        self.assertEqual(harmonize("AT", "G", "AT", "G", .2)["status"], "unsupported_allele")


class MRTests(unittest.TestCase):
    def test_wald_golden_values_and_sign_invariance(self):
        r = wald(.5, .01, -.05, .01)
        self.assertAlmostEqual(r["beta"], -.1)
        self.assertAlmostEqual(r["se"], .02)
        self.assertAlmostEqual(r["se_delta"], np.sqrt(.0004 + .000004))
        self.assertAlmostEqual(r["F"], 2500)
        self.assertAlmostEqual(r["p"], 5.733031437583866e-7)
        self.assertGreater(r["p_delta"], r["p"])
        flipped = wald(-.5, .01, .05, .01)
        for key in ("beta", "se", "p", "se_delta", "p_delta"):
            self.assertEqual(r[key], flipped[key])

    def test_invalid_wald_inputs_fail(self):
        for vals in [(0, .1, .2, .1), (.1, 0, .2, .1), (.1, .1, np.inf, .1),
                     (.1, .1, .2, -1)]:
            with self.assertRaises(ValueError):
                wald(*vals)

    def test_ivw_closed_form_and_underdispersion_floor(self):
        with self.assertRaises(ValueError):
            ivw([1, 2], [.1, .1], [.3, .6], [.1, .1])
        r = ivw([1, 2], [.1, .1], [.3, .6], [.1, .1], independent=True)
        self.assertAlmostEqual(r["beta"], .3)
        self.assertAlmostEqual(r["se"], 1 / np.sqrt(500))
        self.assertEqual(r["se_random"], r["se"])
        self.assertAlmostEqual(r["Q"], 0)
        self.assertEqual(r["df"], 1)

    def test_ivw_heterogeneity_increases_uncertainty(self):
        r = ivw([1, 1], [.01, .01], [0, 1], [.1, .1], independent=True)
        self.assertAlmostEqual(r["beta"], .5)
        self.assertAlmostEqual(r["Q"], 50)
        self.assertAlmostEqual(r["se_random"], .5)
        self.assertGreater(r["p_random"], r["p"])

    def test_variance_scale_must_be_declared(self):
        with self.assertRaises(ValueError):
            variance_explained(.5, .2)
        self.assertAlmostEqual(variance_explained(.5, .2, standardized=True), .08)
        self.assertAlmostEqual(variance_explained(.5, .2, trait_variance=2), .04)
        with self.assertRaises(ValueError):
            variance_explained(10, .5, standardized=True)


class MultiplicityTests(unittest.TestCase):
    def test_bh_and_by_against_independent_library(self):
        p = np.array([.9, .001, .04, .04, .02, 0, 1])
        for name in ("bh", "by"):
            expected = multipletests(p, method="fdr_" + name)[1]
            np.testing.assert_allclose(adjust_pvalues(p, name), expected)

    def test_missing_p_retains_declared_family_size(self):
        r = adjust_pvalues([.01, np.nan, .03])
        np.testing.assert_allclose(r[[0, 2]], [.03, .045])
        self.assertTrue(np.isnan(r[1]))
        df = pd.DataFrame({"p": [.01, np.nan, .03], "family": ["a", "a", "b"]}, index=[5, 5, 8])
        out = adjust_table(df)
        self.assertEqual(list(out.index), [5, 5, 8])
        self.assertEqual(list(out.family_n), [2, 2, 1])
        np.testing.assert_allclose(out.q_bh_family.iloc[[0, 2]], [.02, .03])

    def test_storey_reports_pi0_and_monotonic_qvalues(self):
        p = np.array([.00001, .0001, .001, .005, .01, .03, .07, .1, .2, .9])
        q, pi0 = storey_qvalues(p, .5)
        self.assertAlmostEqual(pi0, .4)
        self.assertTrue(np.all(np.diff(q) >= 0))
        self.assertTrue(np.all(q <= 1))
        self.assertAlmostEqual(q[-1], .36)

    def test_invalid_probabilities_fail_and_empty_families_work(self):
        for p in ([1.01], [-.001], [np.inf]):
            with self.assertRaises(ValueError):
                adjust_pvalues(p)
        self.assertEqual(len(adjust_pvalues([])), 0)
        self.assertTrue(np.isnan(adjust_pvalues([np.nan])).all())

    def test_two_stage_is_alpha_specific_and_labeled(self):
        r = two_stage_bh([.001, .01, .03, .5, .9], alpha=.05)
        self.assertEqual(r["alpha"], .05)
        self.assertIn("m0", r)
        out = adjust_table(pd.DataFrame({"p": [.01, .8], "family": ["a", "a"]}), two_stage_alpha=.1)
        self.assertTrue(out.two_stage_alpha.eq(.1).all())
        self.assertNotIn("significant", out)


class InputAndPresentationTests(unittest.TestCase):
    def test_all_bad_records_survive_with_reasons(self):
        df = pd.DataFrame({"variant_id": ["rs1", "rs1", "rs2", "rs3"],
                           "effect_allele": ["a", "A", "G", "T"],
                           "other_allele": ["G", "G", "C", "C"],
                           "beta": [.1, .1, "bad", .2], "se": [.01, .01, 0, .03]})
        r = validate_sumstats(df)
        self.assertEqual(len(r), 4)
        self.assertEqual(r.valid.tolist(), [False, False, False, True])
        self.assertIn("duplicate_key", r.validation_errors.iloc[0])
        self.assertIn("invalid_beta", r.validation_errors.iloc[2])
        self.assertIn("nonpositive_se", r.validation_errors.iloc[2])
        self.assertEqual(r.effect_allele.iloc[0], "A")

    def test_source_mapping_and_long_table_keys(self):
        df = pd.DataFrame({"snp": ["rs1", "rs1"], "protein": ["A", "B"],
                           "EA": ["A", "A"], "OA": ["C", "C"], "B": [.1, .2], "SE": [.1, .1]})
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "sumstats.tsv.gz"
            df.to_csv(path, sep="\t", index=False)
            r = read_sumstats(path, column_map={"snp": "variant_id", "EA": "effect_allele", "OA": "other_allele", "B": "beta", "SE": "se"},
                              key_cols=("protein", "variant_id"))
        self.assertTrue(r.valid.all())

    def test_direction_is_descriptive_for_correlated_outcomes(self):
        df = pd.DataFrame({"protein": ["TFRC"] * 10, "beta": [-.1] * 8 + [0, np.nan]})
        r = direction_summary(df).iloc[0]
        self.assertEqual(r.n_negative, 8)
        self.assertEqual(r.n_zero, 1)
        self.assertEqual(r.n_missing, 1)
        self.assertTrue(np.isnan(r.sign_p))
        explicit = direction_summary(df, independent=True).iloc[0]
        self.assertAlmostEqual(explicit.sign_p, 2 / 256)
        with self.assertRaises(ValueError):
            direction_summary(df, independent="False")

    def test_ranking_preserves_null_rows_and_scales(self):
        df = pd.DataFrame({"beta": [-.1, .3, np.nan, 4], "effect_unit": ["SD"] * 3 + ["arbitrary"]})
        r = rank_effects(df)
        self.assertEqual(len(r), 4)
        self.assertEqual(r.loc[1, "effect_rank"], 1)
        self.assertEqual(r.loc[3, "effect_rank"], 1)
        self.assertTrue(np.isnan(r.loc[2, "effect_rank"]))

    def test_synthetic_calibration_is_reproducible_and_labeled(self):
        r = null_calibration(n=10000, seed=123)
        self.assertEqual(r, null_calibration(n=10000, seed=123))
        self.assertFalse(r["real_data_used"])
        self.assertLess(abs(r["mean_p"] - .5), .02)
        self.assertLess(abs(r["fraction_p_below_0p05"] - .05), .015)


if __name__ == "__main__":
    unittest.main(verbosity=2)
