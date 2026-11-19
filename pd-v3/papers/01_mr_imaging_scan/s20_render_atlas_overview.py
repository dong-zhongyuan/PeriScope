"""Render a descriptive overview from completed, verified MR presentation tables."""
import argparse
import os
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from s19_prepare_atlas_presentation import completed_source

FAMILY_LABELS = {
    'volume': 'Volume', 'QSM': 'QSM', 'T2star': 'T2*',
    'regional and tissue volume': 'Reg. volume', 'cortical area': 'Cort. area',
    'cortical thickness': 'Thickness', 'WM tract diffusivity': 'WM diff.',
    'WM tract FA': 'WM FA', 'WM tract MO': 'WM MO', 'WM tract ICVF': 'WM ICVF',
    'WM tract OD': 'WM OD', 'WM tract ISOVF': 'WM ISOVF',
    'cortical grey-white contrast': 'GW contrast',
    'regional and tissue intensity': 'Reg. intensity',
    'tfMRI activation': 'Task fMRI', 'rfMRI connectivity': 'Rest connect.',
    'rfMRI node amplitude': 'Rest ampl.', 'regional T2*': 'Reg. T2*',
    'white matter hyperintensity volume': 'WMH volume', 'QC': 'QC',
}


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify_hashes(hashes):
    for relative, expected in hashes.items():
        if sha(ROOT / relative) != expected:
            raise ValueError(f'Input identity changed: {relative}')


def load_verified(mode):
    folder = HERE / 'results/full_atlas_20261116' / mode
    run = completed_source(folder, mode)  # No output before this guard.
    source = folder / 'presentation'
    manifest_path = source / 'MANIFEST.json'
    manifest = json.loads(manifest_path.read_text())
    if (manifest.get('complete') is not True or manifest.get('mode') != mode
            or manifest.get('formal_full_atlas_result') != (mode == 'full')
            or manifest.get('historical_identity_used') is not False
            or manifest.get('full_original_family_denominators_preserved') is not True):
        raise ValueError('Presentation is not an eligible completed source')
    inputs = dict(manifest['inputs'])
    inputs[str(manifest_path.relative_to(ROOT))] = sha(manifest_path)
    for name, digest in manifest['outputs'].items():
        inputs[str((source / name).relative_to(ROOT))] = digest
    if mode == 'full':
        verification_path = folder / 'COMPLETE_RESULT_VERIFICATION.json'
        verification = json.loads(verification_path.read_text())
        if (verification.get('complete') is not True
                or verification.get('formal_full_atlas_verified') is not True
                or verification.get('mode') != mode):
            raise ValueError('Full result has not passed independent complete-table verification')
        for policy in ('drop', 'infer'):
            if verification['policies'][policy]['rows'] != run['planned_pairs']:
                raise ValueError('Complete-result verification has the wrong scope')
        for relative, digest in manifest['inputs'].items():
            if verification['inputs_sha256'].get(relative) != digest:
                raise ValueError(f'Presentation and independent verification differ: {relative}')
        inputs.update(verification['inputs_sha256'])
        inputs[str(verification_path.relative_to(ROOT))] = sha(verification_path)
        verifier = HERE / 'verify_complete_atlas.py'
        inputs[str(verifier.relative_to(ROOT))] = verification['verifier_sha256']
    inputs.update(manifest['code_sha256'])
    verify_hashes(inputs)
    tables = {name: pd.read_csv(source / name, float_precision='round_trip')
              for name in manifest['outputs']}
    return folder, run, manifest, inputs, tables


def derive_tables(tables, run, manifest):
    """Only deterministic display transforms; source values are retained in full."""
    ranks = {p: tables[f'all_assay_ranking_{p}.csv'].copy() for p in ('drop', 'infer')}
    display = tables['overview_display_assays.csv'].copy()
    expected_display = ranks['drop'].loc[ranks['drop'].n_estimated.gt(0)].head(24)
    if display.exposure_id.tolist() != expected_display.exposure_id.tolist():
        raise ValueError('Display rows differ from the declared primary top 24')
    if display.display_row.tolist() != list(range(1, len(display) + 1)):
        raise ValueError('Display row order is invalid')
    source_meta = pd.read_csv(HERE / f'results/full_atlas_20261116/{run["mode"]}/phenotypes.csv')
    family_order = source_meta.family.drop_duplicates().tolist()
    if set(family_order) - FAMILY_LABELS.keys():
        raise ValueError('Unreviewed category labels; add an explicit label before rendering')
    family_counts = source_meta.groupby('family').size().to_dict()
    categories = pd.DataFrame([dict(family=f, label=FAMILY_LABELS[f],
        n_phenotypes=family_counts[f], is_qc=f == 'QC', family_order=i + 1)
        for i, f in enumerate(family_order)])
    # A single genomic order is shared by both policies; instrument-bearing but
    # unestimated assays keep their place and their unestimated source records.
    order = ranks['infer'].loc[ranks['infer'].chrom.notna(),
        ['exposure_id', 'chrom', 'pos37']].sort_values(
        ['chrom', 'pos37', 'exposure_id'], kind='stable').reset_index(drop=True)
    order['genome_order'] = np.arange(1, len(order) + 1)
    chromosome_rows = order.groupby('chrom').genome_order.agg(['min', 'max', 'mean']).reset_index()
    chromosome_rows['label'] = chromosome_rows.chrom.map(lambda x: 'X' if int(x) == 23 else str(int(x)))
    all_rank = []
    matrix = []
    coverage = []
    for policy in ('drop', 'infer'):
        rank = ranks[policy]
        if len(rank) != run['assays'] or rank.exposure_id.duplicated().any():
            raise ValueError('Ranking lost assay identities')
        if not rank.n_planned.eq(manifest['audits'][policy]['biological_phenotypes']).all():
            raise ValueError('Ranking phenotype coverage differs from source audit')
        if not (rank.n_estimated + rank.n_not_estimated).eq(rank.n_planned).all():
            raise ValueError('Ranking coverage is incomplete')
        rank = rank.merge(order[['exposure_id', 'genome_order']], on='exposure_id', how='left', validate='one_to_one')
        estimated = rank.n_estimated.gt(0)
        if rank.loc[estimated, ['genome_order', 'best_log_p', 'best_beta', 'best_se', 'F_x']].isna().any().any():
            raise ValueError('Estimated display record lacks required source values')
        rank['minus_log10_p'] = -rank.best_log_p / np.log(10)
        rank['best_z'] = rank.best_beta / rank.best_se
        rank['policy'] = policy
        all_rank.append(rank)
        profile = tables[f'overview_category_profiles_{policy}.csv'].copy()
        bio_families = [f for f in family_order if f != 'QC']
        expected_keys = {(e, f) for e in display.exposure_id for f in bio_families}
        actual_keys = list(zip(profile.exposure_id, profile.imaging_family))
        if (len(actual_keys) != len(expected_keys) or set(actual_keys) != expected_keys
                or profile.is_qc.any()):
            raise ValueError('Matrix is not the complete non-QC display Cartesian product')
        source_profiles = tables[f'all_assay_category_profiles_{policy}.csv.gz']
        expected = source_profiles.loc[~source_profiles.is_qc & source_profiles.exposure_id.isin(display.exposure_id)]
        fields = source_profiles.columns.tolist()
        keys = ['exposure_id', 'imaging_family']
        pd.testing.assert_frame_equal(profile[fields].sort_values(keys).reset_index(drop=True),
                                      expected[fields].sort_values(keys).reset_index(drop=True))
        profile['best_z'] = profile.best_beta / profile.best_se
        valid = profile.n_estimated.gt(0)
        if ((profile.loc[valid, 'best_z'] < profile.loc[valid, 'z_min'] - 1e-12)
                | (profile.loc[valid, 'best_z'] > profile.loc[valid, 'z_max'] + 1e-12)).any():
            raise ValueError('Best pair lies outside its category z range')
        profile['policy'] = policy
        matrix.append(profile)
        cov = tables[f'coverage_by_phenotype_{policy}.csv'].copy()
        if cov.duplicated(['phenotype', 'status']).any():
            raise ValueError('Duplicate phenotype coverage status')
        if (set(cov.phenotype) != set(source_meta.phenotype)
                or not cov.groupby('phenotype').n_assays.sum().eq(run['assays']).all()):
            raise ValueError('Coverage does not retain every planned assay-phenotype slot')
        totals = cov.groupby('status').n_assays.sum().to_dict()
        if totals != run['policies'][policy]['status_counts']:
            raise ValueError('Coverage status totals differ from completed MR')
        counts = cov.loc[cov.status.eq('estimated')].set_index('phenotype').n_assays
        for record in source_meta[['phenotype', 'family']].to_dict('records'):
            coverage.append(dict(**record, policy=policy,
                n_estimated=int(counts.get(record['phenotype'], 0)),
                n_planned=run['assays'], is_qc=record['family'] == 'QC'))
    all_rank = pd.concat(all_rank, ignore_index=True)
    matrix = pd.concat(matrix, ignore_index=True)
    coverage = pd.DataFrame(coverage)
    coverage['estimated_fraction'] = coverage.n_estimated / coverage.n_planned
    # One row per assay with an estimate in either policy. No pseudo-zero for
    # absent estimates, and no duplicate copy of primary estimates under infer.
    primary = all_rank.loc[all_rank.policy.eq('drop') & all_rank.n_estimated.gt(0)].copy()
    primary['availability'] = 'Primary'
    added = all_rank.loc[all_rank.policy.eq('infer') & all_rank.n_estimated.gt(0)
                        & ~all_rank.exposure_id.isin(primary.exposure_id)].copy()
    added['availability'] = 'EAF-only'
    strength = pd.concat([primary, added], ignore_index=True)
    if strength.exposure_id.duplicated().any() or (strength.F_x <= 0).any():
        raise ValueError('Invalid strength display records')
    strength['absolute_best_z'] = strength.best_z.abs()
    return {'A_all_assays.csv': all_rank, 'A_chromosomes.csv': chromosome_rows,
            'B_profiles.csv': matrix, 'B_display_rows.csv': display,
            'C_instrument_strength.csv': strength, 'D_phenotype_coverage.csv': coverage,
            'categories.csv': categories}


def run_render(mode):
    folder, run, source_manifest, inputs, tables = load_verified(mode)
    output = folder / 'overview_figure'
    if output.exists():
        raise FileExistsError(f'Output exists; inspect instead of overwriting: {output}')
    derived = derive_tables(tables, run, source_manifest)
    code_paths = [Path(__file__), HERE / 'figures_r/common.R',
                  *sorted((HERE / 'figures_atlas_r').glob('*.R')),
                  HERE / 'figures_atlas_r/layout.csv', HERE / 'figures_atlas_r/DESIGN.md']
    code = {str(p.relative_to(ROOT)): sha(p) for p in code_paths}
    verify_hashes(inputs)
    (output / 'plot_data').mkdir(parents=True)
    for name, table in derived.items():
        table.to_csv(output / 'plot_data' / name, index=False, lineterminator='\n')
    config = dict(mode=mode, formal_full_atlas_result=mode == 'full',
        assays=run['assays'], phenotypes=run['phenotypes'], planned_pairs=run['planned_pairs'],
        policies={p: dict(estimated_assays=source_manifest['audits'][p]['ranked_assays'],
            unavailable_assays=run['assays'] - source_manifest['audits'][p]['ranked_assays'])
            for p in ('drop', 'infer')},
        display_assays=len(derived['B_display_rows.csv']),
        matrix_z_limit=max(1, float(np.ceil(derived['B_profiles.csv'][['z_min', 'z_max']].abs().max().max()))),
        no_new_mr_estimation=True, no_new_multiplicity_correction=True,
        historical_additions=False, category_qc_in_matrix=False, category_qc_in_coverage=True)
    (output / 'plot_data/CONFIG.json').write_text(json.dumps(config, indent=2) + '\n')
    manifest = dict(complete=False, mode=mode, created_at_utc=datetime.now(timezone.utc).isoformat(),
        formal_full_atlas_result=mode == 'full', inputs_sha256=inputs, code_sha256=code,
        plot_rows={name: len(frame) for name, frame in derived.items()},
        plot_data_sha256={p.name: sha(p) for p in (output / 'plot_data').iterdir()},
        interpretation='Assay/category minima are descriptive best pairs, not combined tests. Matrix ranges are z ranges, not confidence intervals. No target nomination.',
        visual_review_complete=False)
    manifest_path = output / 'MANIFEST.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    with (output / 'render.log').open('w') as log:
        completed = subprocess.run([os.environ.get('PDV3_RSCRIPT', 'Rscript'),
            str(HERE / 'figures_atlas_r/render.R'), str(output)], stdout=log, stderr=subprocess.STDOUT)
    manifest['render_exit_code'] = completed.returncode
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    if completed.returncode:
        raise RuntimeError(f'R renderer failed; retained log: {output / "render.log"}')
    verify_hashes(inputs)
    verify_hashes(code)
    manifest.update(complete=True, completed_at_utc=datetime.now(timezone.utc).isoformat(),
                    outputs_sha256={str(p.relative_to(output)): sha(p)
                        for p in sorted(output.rglob('*')) if p.is_file() and p != manifest_path})
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(dict(output=str(output), complete=True, mode=mode,
                         formal_full_atlas_result=mode == 'full', plot_rows=manifest['plot_rows']), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['full', 'validation45'], default='full')
    run_render(parser.parse_args().mode)
