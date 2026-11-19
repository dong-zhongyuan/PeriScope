"""Assay-level full-family, identity matching and provenance regression tests."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from pdv3.expanded import build_pairs, expanded_scan
from pdv3.multiplicity import adjust_pvalues


def fixture():
    coverage = pd.DataFrame({"exposure_id": ["G:P1:OID1:v1", "G:P1:OID2:v1", "H:P2:OID3:v1"],
                             "assay": ["OID1", "OID2", "OID3"], "gene": ["G", "G", "H"],
                             "has_cis_primary": [True, True, False],
                             "description": ["assay one", "assay two", "no cis"]})
    instruments = pd.DataFrame({"exposure_id": coverage.exposure_id[:2], "assay": ["OID1", "OID2"],
                                "gene": ["G", "G"], "variant": ["rs1", "rs1"], "chrom": [1, 1],
                                "pos37": [100, 100], "pos38": [110, 110],
                                "ea_x": ["A", "A"], "oa_x": ["G", "G"], "beta_x": [.5, .25],
                                "se_x": [.01, .01], "eaf_x": [.2, .2], "source_x": ["UKBPPP", "UKBPPP"],
                                "units_x": ["reported protein unit"] * 2, "n_x": [30000, 30000],
                                "original_row": [37, 42]})
    phenotypes = pd.DataFrame({"phenotype": ["left_T2", "right_vol"], "family": ["T2star", "volume"],
                              "trait": ["left structure T2*", "right structure volume"],
                              "n_y": [20043, 22128], "source_y": ["BIG40", "BIG40"],
                              "catalog_id": [1, 2]})
    outcomes = pd.DataFrame({"phenotype": ["left_T2"], "variant": ["rs1"], "chrom_y": [1],
                            "pos37": [100], "build_y": ["GRCh37"], "ea_y": ["A"], "oa_y": ["G"],
                            "beta_y": [-.05], "se_y": [.01], "eaf_y": [.2], "n_y": [20043],
                            "source_y": ["BIG40 original discovery"], "source_file": ["idp_1.tsv.gz"],
                            "n_y_status": ["known_autosomal"]})
    return coverage, instruments, outcomes, phenotypes


class ExpandedTests(unittest.TestCase):
    def test_same_gene_assays_are_separate_and_all_planned_slots_retained(self):
        c, i, o, m = fixture()
        pairs = build_pairs(c, i, o, m)
        full, audit = expanded_scan(pairs)
        self.assertEqual(len(full), 6)
        self.assertEqual(len(audit), 6)
        self.assertEqual(full.status.value_counts().to_dict(), {"estimated": 2, "missing_outcome": 2, "no_cis_instrument": 2})
        estimates = full[full.status == "estimated"].set_index("assay")
        self.assertAlmostEqual(estimates.loc["OID1", "beta"], -.1)
        self.assertAlmostEqual(estimates.loc["OID2", "beta"], -.2)
        self.assertEqual(estimates.gene.tolist(), ["G", "G"])
        self.assertTrue(full.family_n_global.eq(6).all())
        self.assertTrue(full.family_n.eq(3).all())
        self.assertTrue(full.cis_family_n.eq(4).all())
        self.assertTrue(full.cis_family_n_assays.eq(2).all())
        np.testing.assert_allclose(full.q_bh_global, adjust_pvalues(full.p), equal_nan=True)
        mask = full.has_cis_instrument
        np.testing.assert_allclose(full.loc[mask, "q_bh_cis_family"], adjust_pvalues(full.loc[mask, "p"]), equal_nan=True)
        self.assertTrue(full.loc[~mask, "q_bh_cis_family"].isna().all())
        self.assertLess(estimates.loc["OID1", "q_bh_cis_family"], estimates.loc["OID1", "q_bh_global"])

    def test_source_columns_and_trait_survive(self):
        c, i, o, m = fixture()
        full, audit = expanded_scan(build_pairs(c, i, o, m))
        row = full.iloc[0]
        self.assertEqual(row.trait, "left structure T2*")
        self.assertEqual(row["assay_meta__description"], "assay one")
        self.assertEqual(row["instrument__original_row"], 37)
        self.assertEqual(row["phenotype_meta__catalog_id"], 1)
        self.assertEqual(row["outcome__source_file"], "idp_1.tsv.gz")
        self.assertEqual(row.source_y, "BIG40 original discovery")
        self.assertEqual(row["outcome_row_index"], 0)
        self.assertTrue(set(build_pairs(c, i, o, m).columns).issubset(full.columns))

    def test_unknown_variant_sample_n_never_falls_back_to_catalog(self):
        c, i, o, m = fixture()
        o["n_y"] = np.nan
        o["n_y_status"] = "unknown_chrX"
        full, _ = expanded_scan(build_pairs(c, i, o, m))
        rows = full[full.status == "estimated"]
        self.assertTrue(rows.n_y.isna().all())
        self.assertTrue(rows.n_y_status.eq("unknown_chrX").all())
        self.assertTrue(rows["phenotype_meta__n_y"].eq(20043).all())

    def test_rsid_optional_position_and_alleles_are_required(self):
        c, i, o, m = fixture()
        o = o.drop(columns="variant")
        full, _ = expanded_scan(build_pairs(c, i, o, m))
        self.assertEqual((full.status == "estimated").sum(), 2)
        self.assertTrue(full.loc[full.status == "estimated", "variant_match"].eq("same_build_position_and_alleles").all())
        with self.assertRaisesRegex(ValueError, "rsID-only"):
            build_pairs(c, i, o.drop(columns="pos37"), m)

    def test_same_rsid_wrong_position_is_rejected(self):
        c, i, o, m = fixture()
        o["pos37"] = 101
        full, _ = expanded_scan(build_pairs(c, i, o, m))
        rows = full[(full.phenotype == "left_T2") & full.has_cis_instrument]
        self.assertTrue(rows.status.eq("coordinate_mismatch").all())
        self.assertTrue(rows.p.isna().all())
        self.assertTrue(rows["outcome__pos37"].eq(101).all())

    def test_same_position_multiallelic_lookup_uses_alleles_not_row_order(self):
        c, i, o, m = fixture()
        incompatible = o.copy()
        incompatible["oa_y"] = "C"
        incompatible["variant"] = "rsDifferent"
        o = pd.concat([incompatible, o], ignore_index=True)
        full, _ = expanded_scan(build_pairs(c, i, o, m))
        rows = full[full.status == "estimated"]
        self.assertEqual(len(rows), 2)
        self.assertTrue(rows.outcome_row_index.eq(1).all())
        self.assertTrue(rows.outcome_candidate_count.eq(2).all())
        with self.assertRaisesRegex(ValueError, "multiple compatible"):
            build_pairs(c, i, pd.concat([o, o.iloc[[1]]], ignore_index=True), m)

    def test_build38_uses_build38_exposure_position(self):
        c, i, o, m = fixture()
        o = o.drop(columns="pos37")
        o["pos38"] = 110
        o["build_y"] = "GRCh38"
        full, _ = expanded_scan(build_pairs(c, i, o, m))
        rows = full[full.status == "estimated"]
        self.assertEqual(len(rows), 2)
        self.assertTrue(rows.pos_x.eq(110).all())
        self.assertTrue(rows.build_x.eq("GRCh38").all())

    def test_chromosome_x_aliases_and_padded_autosomes(self):
        for exposure_chrom, outcome_chrom, canonical in [(23, "0X", "X"), ("X", "23", "X"),
                                                        ("chrX", "0X", "X"), (1, "01", "1")]:
            c, i, o, m = fixture()
            i["chrom"] = exposure_chrom
            o["chrom_y"] = outcome_chrom
            full, _ = expanded_scan(build_pairs(c, i, o, m))
            rows = full[full.status == "estimated"]
            self.assertEqual(len(rows), 2)
            self.assertTrue(rows.chrom_x.eq(canonical).all())
            self.assertTrue(rows.chrom_y.eq(canonical).all())
            self.assertTrue(rows["instrument__chrom"].eq(exposure_chrom).all())
            self.assertTrue(rows["outcome__chrom_y"].eq(outcome_chrom).all())

    def test_selected_indels_are_audited_without_ambiguity_or_fallback(self):
        for count in (0, 1, 2):
            c, i, o, m = fixture()
            i["oa_x"] = "AG"
            if count == 0:
                o = o.iloc[:0]
            elif count == 2:
                other = o.copy()
                other["oa_y"] = "AT"
                o = pd.concat([o, other], ignore_index=True)
            full, _ = expanded_scan(build_pairs(c, i, o, m))
            cis = full[full.has_cis_instrument]
            self.assertTrue(cis.status.eq("unsupported_allele").all())
            self.assertTrue(cis.p.isna().all())
            self.assertTrue(cis.variant.eq("rs1").all())
            self.assertTrue(cis["instrument__oa_x"].eq("AG").all())
            self.assertTrue(cis.loc[cis.phenotype == "left_T2", "outcome_candidate_count"].eq(count).all())
            self.assertTrue(full.cis_family_n.eq(4).all())
            if count == 1:
                self.assertTrue(cis.loc[cis.phenotype == "left_T2", "outcome__source_file"].eq("idp_1.tsv.gz").all())

    def test_palindrome_frequency_sensitivity_is_explicit_and_family_unchanged(self):
        c, i, o, m = fixture()
        i["oa_x"], o["oa_y"] = "T", "T"
        i["eaf_x"], o["eaf_y"] = .2, .8
        pairs = build_pairs(c, i, o, m)
        default, _ = expanded_scan(pairs)
        inferred, audit = expanded_scan(pairs, palindrome="infer")
        self.assertEqual((default.status == "palindromic_dropped").sum(), 2)
        self.assertTrue(default.palindrome_policy.eq("drop").all())
        self.assertEqual((inferred.status == "estimated").sum(), 2)
        self.assertTrue(inferred.palindrome_policy.eq("infer").all())
        self.assertTrue(audit.palindrome_policy.eq("infer").all())
        self.assertTrue(inferred.loc[inferred.status == "estimated", "allele_flipped"].eq(True).all())
        self.assertAlmostEqual(inferred.loc[(inferred.assay == "OID1") & (inferred.status == "estimated"), "beta"].iloc[0], .1)
        self.assertEqual(default.family_n_global.tolist(), inferred.family_n_global.tolist())
        self.assertEqual(default.cis_family_n.tolist(), inferred.cis_family_n.tolist())

    def test_palindrome_inference_preserves_missing_and_ambiguous_exclusions(self):
        for fx, fy, status in [(np.nan, .2, "palindromic_missing_eaf"),
                              (.2, np.nan, "palindromic_missing_eaf"),
                              (.49, .49, "palindromic_ambiguous_eaf")]:
            c, i, o, m = fixture()
            i["oa_x"], o["oa_y"] = "T", "T"
            i["eaf_x"], o["eaf_y"] = fx, fy
            full, _ = expanded_scan(build_pairs(c, i, o, m), palindrome="infer")
            rows = full[(full.phenotype == "left_T2") & full.has_cis_instrument]
            self.assertTrue(rows.status.eq(status).all())
            self.assertTrue(rows.p.isna().all())
            self.assertTrue(full.cis_family_n.eq(4).all())

    def test_missing_beta_and_empty_outcomes_are_not_null_effects(self):
        c, i, o, m = fixture()
        o["beta_y"] = np.nan
        full, _ = expanded_scan(build_pairs(c, i, o, m))
        self.assertEqual((full.status == "missing_statistics").sum(), 2)
        self.assertTrue(full.p.isna().all())
        empty_full, _ = expanded_scan(build_pairs(c, i, pd.DataFrame(), m))
        self.assertEqual(len(empty_full), 6)
        self.assertEqual((empty_full.status == "missing_outcome").sum(), 4)
        self.assertTrue(empty_full.p.isna().all())

    def test_missing_selected_lead_differs_from_no_published_cis(self):
        c, i, o, m = fixture()
        full, _ = expanded_scan(build_pairs(c, i.iloc[[0]], o, m))
        self.assertTrue(full.loc[full.assay == "OID2", "status"].eq("missing_selected_instrument").all())
        self.assertTrue(full.loc[full.assay == "OID3", "status"].eq("no_cis_instrument").all())

    def test_duplicate_instruments_or_incomplete_plan_fail(self):
        c, i, o, m = fixture()
        with self.assertRaisesRegex(ValueError, "duplicate key"):
            build_pairs(c, pd.concat([i, i.iloc[[0]]], ignore_index=True), o, m)
        pairs = build_pairs(c, i, o, m)
        with self.assertRaisesRegex(ValueError, "Cartesian"):
            expanded_scan(pairs.iloc[:-1])
        pairs.loc[0, "pos_y"] = 999
        full, _ = expanded_scan(pairs)
        self.assertEqual(full.iloc[0].status, "coordinate_mismatch")


if __name__ == "__main__":
    unittest.main(verbosity=2)
