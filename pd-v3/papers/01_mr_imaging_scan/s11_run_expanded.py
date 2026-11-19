"""Run complete assay-by-IDP scans with parallel palindrome policies.

The all-assay and selected-cis family denominators depend only on exposure
availability and the declared phenotype catalog, never outcome P values.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from common import ROOT, PAPER, provenance
from pdv3.expanded import build_pairs, expanded_scan
from pdv3.provenance import sha256


def normalize_chrom(value):
    value = str(value).upper().removeprefix('CHR').removesuffix('.0').lstrip('0')
    return 'X' if value == '23' else value


def attach_frequencies(outcomes, frequencies):
    """Exact forward allele pairing to Oxford's A2 frequencies; no LD proxies.

    Frequency rows are identified by GRCh37 coordinates and exact allele pair.
    Reversed pairs transform AF to 1-AF. Strand-complement guesses are not used
    here, particularly for palindromes; they would make inference circular.
    """
    lookup = {}
    for r in frequencies.to_dict('records'):
        key = (normalize_chrom(r['chr']), int(r['pos37']))
        lookup.setdefault(key, []).append(r)
    rows = []
    for raw in outcomes.to_dict('records'):
        rec = dict(raw, eaf_y=np.nan, eaf_attachment_status='frequency_not_available',
                   eaf_source=None, eaf_reference_cohort='BIG40 release2 genotype frequency resource')
        if raw['build_y'] != 'GRCh37':
            rec['eaf_attachment_status'] = 'build_not_GRCh37'
            rows.append(rec)
            continue
        candidates = lookup.get((normalize_chrom(raw['chrom_y']), int(raw['pos_y'])), [])
        matches = []
        for f in candidates:
            if (raw['ea_y'], raw['oa_y']) == (f['ea_y'], f['oa_y']):
                matches.append((f['eaf_y'], f['source_y'], 'exact_alleles'))
            elif (raw['ea_y'], raw['oa_y']) == (f['oa_y'], f['ea_y']):
                matches.append((1 - f['eaf_y'], f['source_y'], 'reversed_alleles'))
        if len(matches) == 1 and pd.notna(matches[0][0]):
            rec.update(eaf_y=matches[0][0], eaf_source=matches[0][1], eaf_attachment_status=matches[0][2])
        elif len(matches) > 1:
            rec['eaf_attachment_status'] = 'ambiguous_frequency_rows'
        elif candidates:
            rec['eaf_attachment_status'] = 'frequency_allele_mismatch_or_missing'
        rows.append(rec)
    return pd.DataFrame(rows)


def write(frame, path, inputs, extra):
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {'method': 'gzip', 'mtime': 0} if path.name.endswith('.gz') else None
    frame.to_csv(path, index=False, compression=compression)
    provenance(path, inputs, extra)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phenotypes', type=Path, default=ROOT / 'data/big40/phase1_17_phenotypes.csv')
    parser.add_argument('--outcomes', type=Path, default=ROOT / 'data/big40/extracted/phase1_17_outcomes.csv')
    parser.add_argument('--acquisition-manifest', type=Path, default=ROOT / 'data/big40/PHASE1_17_MANIFEST.json')
    parser.add_argument('--output', type=Path, default=PAPER / 'results/expanded_20261116/big40_17')
    args = parser.parse_args()
    instrument_path = ROOT / 'data/ukb_ppp/extracted/cis_primary_selected.csv'
    coverage_path = ROOT / 'data/ukb_ppp/extracted/assay_coverage.csv'
    af_path = ROOT / 'data/big40_variants/extracted/cis_instrument_frequencies.csv'
    source_inputs = [coverage_path, instrument_path, args.outcomes, args.phenotypes, af_path,
                     args.acquisition_manifest, ROOT / 'data/ukb_ppp/MANIFEST.json',
                     ROOT / 'data/big40_variants/MANIFEST.json']
    if not all(p.exists() for p in source_inputs):
        raise SystemExit(f'Missing inputs: {[str(p) for p in source_inputs if not p.exists()]}')
    acquisition = json.loads(args.acquisition_manifest.read_text())
    if acquisition.get('complete') is not True:
        raise SystemExit('Acquisition manifest is incomplete; do not label query failures as null results')
    if acquisition.get('output_sha256') != sha256(args.outcomes):
        raise SystemExit('Outcome extraction differs from acquisition manifest')
    cov, inst = pd.read_csv(coverage_path), pd.read_csv(instrument_path)
    meta, y = pd.read_csv(args.phenotypes), pd.read_csv(args.outcomes)
    af = pd.read_csv(af_path)
    y = attach_frequencies(y, af)
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    extra = {'design': 'complete assay-by-declared-IDP plan', 'outcome_selection_by_results': False,
             'coverage_assays': len(cov), 'selected_cis_assays': len(inst), 'phenotypes': len(meta)}
    write(cov, output / 'assay_coverage.csv', source_inputs, extra)
    write(inst, output / 'selected_cis_instruments.csv', source_inputs, extra)
    write(meta, output / 'phenotypes.csv', source_inputs, extra)
    write(y, output / 'outcome_locus_statistics_with_frequencies.csv.gz', source_inputs, extra)
    print('Building complete plan...', flush=True)
    pairs = build_pairs(cov, inst, y, meta)
    write(pairs, output / 'planned_pairs.csv.gz', source_inputs, extra)
    summaries = {}
    for policy in ['drop', 'infer']:
        print(f'Running palindrome={policy}: {len(pairs):,} planned pairs', flush=True)
        full, audit = expanded_scan(pairs, palindrome=policy)
        full['effect_rank_within_assay_modality'] = full.groupby(['exposure_id', 'family']).beta.transform(
            lambda s: s.abs().rank(method='min', ascending=False))
        full['p_rank_entire_plan'] = full.p.rank(method='min')
        name = f'full_family_palindrome_{policy}.csv.gz'
        write(full, output / name, source_inputs, dict(extra, palindrome_policy=policy))
        # Full table contains every audit field plus adjustments; no redundant second 100MB table.
        coverage = (full.groupby(['phenotype', 'trait', 'family', 'status'], dropna=False)
                    .size().rename('n_assays').reset_index())
        write(coverage, output / f'coverage_by_phenotype_{policy}.csv', source_inputs, extra)
        historic = full[full.gene.isin(['CD22', 'CD40', 'TFRC'])]
        write(historic, output / f'three_historical_proteins_{policy}.csv', source_inputs,
              dict(extra, full_family_file=name, selection='three proteins named before this scan'))
        compact_cols = ['exposure_id', 'gene', 'assay', 'phenotype', 'trait', 'family', 'status',
                        'beta', 'se', 'p', 'log_p', 'q_bh_global', 'q_by_global', 'q_bh_family',
                        'q_bh_cis_family', 'q_by_cis_family', 'F', 'n_y', 'palindrome_policy']
        compact = full[[c for c in compact_cols if c in full]].sort_values('p', na_position='last', kind='stable')
        write(compact, output / f'all_associations_ordered_by_p_{policy}.csv.gz', source_inputs, extra)
        summaries[policy] = dict(planned_pairs=len(full), estimated=int(full.status.eq('estimated').sum()),
                                 status_counts=full.status.value_counts().to_dict(),
                                 full_family_size=int(full.family_n_global.iloc[0]),
                                 selected_cis_family_size=int(full.cis_family_n.iloc[0]),
                                 exposure_assays=len(cov), selected_cis_assays=len(inst), phenotypes=len(meta),
                                 finite_result_assays=full.loc[full.status.eq('estimated'), 'exposure_id'].nunique())
        print(json.dumps(summaries[policy]), flush=True)
    result = dict(**extra, policies=summaries, interpretation='Exploratory MR estimates; no biological support decisions',
                  sample_overlap='UKB exposure and imaging may overlap; zero-covariance SE approximation',
                  n_x='ST9 per-assay sample size not reported; left missing',
                  n_y='Autosomal exact IDP n where verified; X and unverifiable QSM n remain missing',
                  direction='MR sign is imaging change per protein analysis unit, not clinical benefit',
                  allele_frequency='BIG40 reference frequencies aid orientation; QSM-specific frequencies unavailable',
                  global_family='all measured assays x declared IDPs; unavailable estimates occupy missing slots',
                  cis_family='selected cis assays x declared IDPs; selection uses exposure data only',
                  command=sys.argv)
    summary = output / 'scan_summary.json'
    summary.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    provenance(summary, source_inputs, extra)


if __name__ == '__main__':
    main()
