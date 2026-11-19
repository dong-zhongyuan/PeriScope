"""One-lead MR over the complete assay-by-phenotype Cartesian product.

The statistical unit is exposure_id (a measured assay), never gene. Metadata is
retained using assay_meta__, instrument__, phenotype_meta__, and outcome__
prefixes, alongside canonical analysis fields. Missing instruments/statistics
remain rows with NA estimates and occupy slots in the declared testing family.
No-cis means no supplied selected instrument, not evidence of a null effect.
"""
from collections import defaultdict
import numpy as np
import pandas as pd

from .mr import harmonize, wald
from .multiplicity import adjust_table, adjust_pvalues


def _present(value):
    return value is not None and not pd.isna(value) and str(value).strip() != ""


def _first(row, names):
    for name in names:
        if name in row and _present(row[name]):
            return row[name]
    return None


def _build(value):
    if not _present(value):
        return None
    return {"grch37": "GRCh37", "hg19": "GRCh37", "37": "GRCh37",
            "grch38": "GRCh38", "hg38": "GRCh38", "38": "GRCh38"}.get(str(value).lower())


def _chrom(value):
    if not _present(value):
        return None
    value = str(value).strip().upper().removeprefix("CHR").removesuffix(".0")
    # UKB-PPP ST9 uses numeric 23 for X; Oxford encodes X as 0X and
    # may zero-pad autosomes. Preserve raw source fields in the metadata.
    value = {"0X": "X", "23": "X", "0Y": "Y", "24": "Y"}.get(value, value)
    if value.isdecimal() and 1 <= int(value) <= 22:
        value = str(int(value))
    if value == "M":
        value = "MT"
    return value if value in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"} else None


def _position(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if np.isfinite(number) and number > 0 and number.is_integer() else None


def _compatible(ea, oa, ya, yo):
    if not all(_present(v) for v in (ea, oa, ya, yo)):
        return False
    ea, oa, ya, yo = [str(v).strip().upper() for v in (ea, oa, ya, yo)]
    if any(v not in {"A", "C", "G", "T"} for v in (ea, oa, ya, yo)) or ea == oa or ya == yo:
        return False
    complement = str.maketrans("ACGT", "TGCA")
    return ((ea, oa) in {(ya, yo), (yo, ya),
                         (ya.translate(complement), yo.translate(complement)),
                         (yo.translate(complement), ya.translate(complement))})


def _exposure_allele_status(inst):
    alleles = [inst.get("ea_x"), inst.get("oa_x")]
    if not all(_present(a) for a in alleles):
        return "missing_exposure_alleles"
    if any(str(a).strip().upper() not in {"A", "C", "G", "T"} for a in alleles):
        return "unsupported_allele"
    return None


def _require_table(df, columns, name, unique=None):
    if df.columns.duplicated().any():
        raise ValueError(f"{name}: duplicate columns")
    missing = set(columns) - set(df.columns)
    if missing:
        raise ValueError(f"{name}: missing columns {sorted(missing)}")
    for col in columns:
        if col in {"exposure_id", "phenotype", "assay", "family"} and (df[col].isna() | df[col].astype(str).str.strip().eq("")).any():
            raise ValueError(f"{name}: missing {col}")
    if unique and df.duplicated(unique).any():
        raise ValueError(f"{name}: duplicate key {unique}; resolve before scanning")


def _prefixed(row, prefix):
    return {prefix + str(key): value for key, value in row.items()}


def _boolean(value, name):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if str(value).lower() in {"true", "1", "1.0"}:
        return True
    if str(value).lower() in {"false", "0", "0.0"}:
        return False
    raise ValueError(f"{name} must be an explicit boolean")


def build_pairs(assay_coverage, instruments, outcomes, phenotypes):
    """Build a full plan without dropping assays lacking a supplied cis lead.

    coverage requires exposure_id, assay, gene; phenotype catalog requires
    phenotype, family. Instruments have at most one row per exposure_id and
    require gene, assay, variant, chrom, pos37 and/or pos38, ea_x, oa_x, beta_x,
    se_x. Source metadata and additional columns are retained without inference.

    Nonempty outcomes require phenotype, build_y, ea_y, oa_y, beta_y, se_y and
    coordinates (chrom_y/chrom/chr/chromosome plus pos_y/pos/position/pos37/
    pos38). The explicit build must agree with a build-specific position column.
    variant/variant_y/rsid/variant_id is optional and never sufficient to match.
    Lookup requires same build, chromosome, position, and compatible SNP alleles.
    Multiple compatible rows at a locus are an error, not an invitation to pick
    the smallest P. Original source row numbers are zero-based positional indexes.
    Selected indels remain unsupported_allele rows with candidate counts; they
    never trigger fallback to another exposure lead or ambiguous-outcome errors.
    """
    _require_table(assay_coverage, ["exposure_id", "assay", "gene"], "assay_coverage", ["exposure_id"])
    _require_table(phenotypes, ["phenotype", "family"], "phenotypes", ["phenotype"])
    if assay_coverage.empty or phenotypes.empty:
        raise ValueError("assay and phenotype catalogs must be nonempty")
    if not instruments.empty:
        _require_table(instruments, ["exposure_id", "assay", "gene", "variant", "chrom", "ea_x", "oa_x", "beta_x", "se_x"],
                       "instruments", ["exposure_id"])
        if not {"pos37", "pos38", "pos_x"}.intersection(instruments.columns):
            raise ValueError("instruments: missing positions")
        if not set(instruments.exposure_id).issubset(set(assay_coverage.exposure_id)):
            raise ValueError("instrument exposure_id absent from assay catalog")
    inst_by_id = {row["exposure_id"]: row for row in instruments.to_dict("records")}
    if not outcomes.empty:
        _require_table(outcomes, ["phenotype", "build_y", "ea_y", "oa_y", "beta_y", "se_y"], "outcomes")
        if not {"chrom_y", "chrom", "chr", "chromosome"}.intersection(outcomes.columns):
            raise ValueError("outcomes: missing chromosome field; rsID-only matching is forbidden")
        if not {"pos_y", "pos", "position", "pos37", "pos38"}.intersection(outcomes.columns):
            raise ValueError("outcomes: missing position field; rsID-only matching is forbidden")
        if not set(outcomes.phenotype).issubset(set(phenotypes.phenotype)):
            raise ValueError("outcome phenotype absent from requested catalog")
    indexed, by_rsid = defaultdict(list), defaultdict(list)
    canonical = []
    for i, raw in enumerate(outcomes.to_dict("records")):
        build = _build(raw.get("build_y"))
        generic_pos = _first(raw, ["pos_y", "pos", "position"])
        assembly_pos = _first(raw, ["pos37"] if build == "GRCh37" else ["pos38"] if build == "GRCh38" else [])
        if generic_pos is not None and assembly_pos is not None and _position(generic_pos) != _position(assembly_pos):
            raise ValueError(f"outcome row {i}: conflicting generic/build-specific positions")
        pos = _position(generic_pos if generic_pos is not None else assembly_pos)
        chrom = _chrom(_first(raw, ["chrom_y", "chrom", "chr", "chromosome"]))
        variant = _first(raw, ["variant_y", "variant", "rsid", "variant_id"])
        row = {**raw, "build_y": build, "chrom_y": chrom, "pos_y": pos, "variant_y": variant,
               "outcome_row_index": i, "_raw": raw}
        canonical.append(row)
        if build and chrom and pos:
            indexed[(raw["phenotype"], build, chrom, pos)].append(row)
        if _present(variant):
            by_rsid[(raw["phenotype"], str(variant))].append(row)

    rows = []
    for cov in assay_coverage.to_dict("records"):
        inst = inst_by_id.get(cov["exposure_id"])
        published_cis = _boolean(cov["has_cis_primary"], "has_cis_primary") if "has_cis_primary" in cov else None
        if inst and published_cis is False:
            raise ValueError(f"instrument conflicts with no-cis coverage: {cov['exposure_id']}")
        if inst:
            for key in ("assay", "gene"):
                if _present(cov.get(key)) and _present(inst.get(key)) and str(cov[key]) != str(inst[key]):
                    raise ValueError(f"coverage/instrument {key} mismatch for {cov['exposure_id']}")
        coords = {}
        if inst:
            chrom = _chrom(inst.get("chrom"))
            coords = {b: (chrom, _position(inst.get(p))) for b, p in (("GRCh37", "pos37"), ("GRCh38", "pos38"))}
            if "pos_x" in inst and _build(inst.get("build_x")):
                coords[_build(inst["build_x"])] = (chrom, _position(inst["pos_x"]))
            coords = {b: value for b, value in coords.items() if all(v is not None for v in value)}
        for meta in phenotypes.to_dict("records"):
            rec = {**_prefixed(cov, "assay_meta__"), **_prefixed(meta, "phenotype_meta__"),
                   "exposure_id": cov["exposure_id"], "gene": cov["gene"], "assay": cov["assay"],
                   "phenotype": meta["phenotype"], "family": meta["family"],
                   "trait": meta.get("trait", meta["phenotype"]),
                   "pair_status": "missing_selected_instrument" if published_cis is True and inst is None else "no_cis_instrument",
                   "has_cis_instrument": inst is not None,
                   "n_planned_instruments": int(inst is not None),
                   "outcome_candidate_count": 0, "outcome_candidate_rows": "",
                   "variant_match": None}
            for key in ("source_y", "units_y", "n_y", "modality"):
                if key in meta:
                    rec[key] = meta[key]
            if not inst:
                rows.append(rec)
                continue
            rec.update(_prefixed(inst, "instrument__"))
            for key, value in inst.items():
                if key not in {"exposure_id", "gene", "assay", "family", "phenotype", "trait"}:
                    rec[key] = value
            allele_status = _exposure_allele_status(inst)
            if not coords:
                rec["pair_status"] = allele_status or "invalid_exposure_coordinates"
                rows.append(rec)
                continue
            default_build = "GRCh37" if "GRCh37" in coords else next(iter(coords))
            rec.update(build_x=default_build, chrom_x=coords[default_build][0], pos_x=coords[default_build][1])
            candidates = []
            for build, (chrom, pos) in coords.items():
                candidates.extend(indexed.get((meta["phenotype"], build, chrom, pos), []))
            rec["outcome_candidate_count"] = len(candidates)
            rec["outcome_candidate_rows"] = ";".join(str(r["outcome_row_index"]) for r in candidates)
            if allele_status:
                # Count every same-position candidate, but do not select between
                # incompatible rows or try a different lead for unsupported alleles.
                rec["pair_status"] = allele_status
                if len(candidates) == 1:
                    chosen = candidates[0]
                    rec.update(_prefixed(chosen["_raw"], "outcome__"))
                    rec.update({key: value for key, value in chosen.items()
                                if key.endswith("_y") or key in {"outcome_row_index", "n_y_status"}})
                    rec["n_y"] = chosen.get("n_y", np.nan)
                    rec.update(build_x=chosen["build_y"], chrom_x=coords[chosen["build_y"]][0],
                               pos_x=coords[chosen["build_y"]][1])
                rows.append(rec)
                continue
            compatible = [r for r in candidates if _compatible(inst["ea_x"], inst["oa_x"], r["ea_y"], r["oa_y"])]
            if len(compatible) > 1:
                raise ValueError(f"multiple compatible outcome rows: {cov['exposure_id']} / {meta['phenotype']}; rows {rec['outcome_candidate_rows']}")
            if compatible:
                chosen = compatible[0]
                rec["pair_status"] = "ready"
                rec["variant_match"] = "same_build_position_and_alleles"
            elif len(candidates) == 1:
                chosen = candidates[0]
                rec["pair_status"] = "allele_mismatch"
            elif len(candidates) > 1:
                # Reject ambiguity rather than losing the candidate source records.
                raise ValueError(f"multiple incompatible rows at requested locus: {cov['exposure_id']} / {meta['phenotype']}")
            else:
                rsid_rows = by_rsid.get((meta["phenotype"], str(inst["variant"])), [])
                if len(rsid_rows) > 1:
                    raise ValueError(f"multiple rsID records without coordinate match: {cov['exposure_id']} / {meta['phenotype']}")
                if not rsid_rows:
                    rec["pair_status"] = "missing_outcome"
                    rows.append(rec)
                    continue
                chosen = rsid_rows[0]
                rec["pair_status"] = "coordinate_mismatch" if chosen["build_y"] in coords else "build_mismatch_or_unknown"
                rec["variant_match"] = "rsid_only_rejected"
                rec["outcome_candidate_count"] = 1
                rec["outcome_candidate_rows"] = str(chosen["outcome_row_index"])
            rec.update(_prefixed(chosen["_raw"], "outcome__"))
            for key, value in chosen.items():
                if key.endswith("_y") or key in {"outcome_row_index", "n_y_status"}:
                    rec[key] = value
            # A matched row's NA n is meaningful (e.g. unknown chrX sample size),
            # so never fall back to the phenotype catalog's autosomal n.
            rec["n_y"] = chosen.get("n_y", np.nan)
            if chosen["build_y"] in coords:
                rec.update(build_x=chosen["build_y"], chrom_x=coords[chosen["build_y"]][0], pos_x=coords[chosen["build_y"]][1])
            rec["rsid_agreement"] = (str(inst["variant"]) == str(chosen["variant_y"])) if _present(chosen["variant_y"]) else None
            rows.append(rec)
    return pd.DataFrame(rows)


def expanded_scan(pairs, *, palindrome="drop"):
    """Return (full_result_table, full_row_audit) for an assay-complete plan.

    One selected lead per assay is assumed. No instrument selection, LD pruning,
    outcome-based filtering, or multi-instrument pooling occurs here. All columns
    from build_pairs survive in both outputs. Global/family correction includes
    all planned assay-by-phenotype slots, including NA outcomes and no-cis assays.
    q_bh_cis_family/q_by_cis_family also use the full selected-cis assay-by-IDP
    family (all modalities together); its membership depends only on exposure
    instrument availability, never outcome values or availability. Missing
    outcomes remain slots. No-cis rows have NA for these additional corrections.
    palindrome='drop' retains the conservative default. Explicit 'infer' is a
    separate sensitivity using harmonize's existing EAF tolerance 0.1 and MAF
    ceiling 0.42; absent or ambiguous frequencies remain unusable. Every output
    row records the policy, and the full planned family is identical either way.
    """
    if palindrome not in {"drop", "infer"}:
        raise ValueError("palindrome must be 'drop' or 'infer'")
    _require_table(pairs, ["exposure_id", "assay", "gene", "phenotype", "family", "pair_status", "has_cis_instrument"],
                   "pairs", ["exposure_id", "phenotype"])
    if pairs.empty or len(pairs) != pairs.exposure_id.nunique() * pairs.phenotype.nunique():
        raise ValueError("incomplete assay-by-phenotype Cartesian product")
    if (pairs.groupby("exposure_id").assay.nunique() != 1).any() or (pairs.groupby("phenotype").family.nunique() != 1).any():
        raise ValueError("inconsistent assay identity or phenotype family")
    cis_flags = pairs.has_cis_instrument.map(lambda value: _boolean(value, "has_cis_instrument"))
    if (pd.DataFrame({"exposure_id": pairs.exposure_id, "cis": cis_flags}).groupby("exposure_id").cis.nunique() != 1).any():
        raise ValueError("cis instrument availability must be constant across phenotypes")
    if ((pairs.pair_status == "ready") & ~cis_flags).any():
        raise ValueError("ready pair has no selected cis instrument")
    rows = []
    for raw in pairs.to_dict("records"):
        rec = {**raw, "status": raw["pair_status"], "usable": False, "method": None,
               "beta": np.nan, "se": np.nan, "p": np.nan, "beta_y_aligned": np.nan,
               "n_instruments": 0, "harmonization_status": None, "palindrome_policy": palindrome}
        if raw["pair_status"] != "ready":
            rows.append(rec)
            continue
        required = ["ea_y", "oa_y", "beta_y", "se_y", "ea_x", "oa_x", "beta_x", "se_x"]
        if any(not _present(raw.get(key)) for key in required):
            rec["status"] = "missing_statistics"
            rows.append(rec)
            continue
        # Recheck coordinates so externally constructed/tampered plans cannot
        # bypass the build_pairs identity safeguards by setting pair_status=ready.
        bx, by = _build(raw.get("build_x")), _build(raw.get("build_y"))
        cx, cy = _chrom(raw.get("chrom_x")), _chrom(raw.get("chrom_y"))
        px, py = _position(raw.get("pos_x")), _position(raw.get("pos_y"))
        if bx is None or by is None or bx != by:
            rec["status"] = "build_mismatch_or_unknown"
        elif None in (cx, cy, px, py):
            rec["status"] = "invalid_coordinates"
        elif (cx, px) != (cy, py):
            rec["status"] = "coordinate_mismatch"
        else:
            h = harmonize(raw["ea_x"], raw["oa_x"], raw["ea_y"], raw["oa_y"], raw["beta_y"],
                          eaf_x=raw.get("eaf_x"), eaf_y=raw.get("eaf_y"), exposure_build=bx,
                          outcome_build=by, palindrome=palindrome)
            rec.update(harmonization_status=h["status"], beta_y_aligned=h["beta_y"],
                       eaf_y_aligned=h["eaf_y"], allele_flipped=h["flipped"], eaf_difference=h["eaf_difference"])
            if not h["usable"]:
                rec["status"] = h["status"]
            else:
                try:
                    rec.update(wald(raw["beta_x"], raw["se_x"], h["beta_y"], raw["se_y"]))
                    rec.update(status="estimated", usable=True)
                except (TypeError, ValueError, OverflowError) as exc:
                    rec.update(status="invalid_statistics", error=str(exc))
        rows.append(rec)
    audit = pd.DataFrame(rows)
    full = adjust_table(audit, family_cols=("family",))
    cis_flags = cis_flags.to_numpy(dtype=bool)
    full["cis_family_n"] = int(cis_flags.sum())
    full["cis_family_n_observed"] = int(np.isfinite(full.loc[cis_flags, "p"]).sum())
    full["cis_family_n_assays"] = int(full.loc[cis_flags, "exposure_id"].nunique())
    full["cis_family_n_phenotypes"] = int(full.phenotype.nunique())
    full["q_bh_cis_family"] = np.nan
    full["q_by_cis_family"] = np.nan
    for name, method in (("q_bh_cis_family", "bh"), ("q_by_cis_family", "by")):
        full.loc[cis_flags, name] = adjust_pvalues(full.loc[cis_flags, "p"].to_numpy(dtype=float), method)
    return full, audit


# Keep the misspelling from the initial task description usable for callers.
exanded_scan = expanded_scan
