"""Prepare source-backed ranking and figure tables after full MR completion."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.atlas_presentation import summarize_table


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def completed_source(folder, mode):
    path = folder / 'RUN.json'
    if not path.exists():
        raise RuntimeError('Completed full MR run is unavailable; presentation disabled')
    run = json.loads(path.read_text())
    required = (2941, 3963, 11655183) if mode == 'full' else (2941, 45, 132345)
    actual = (run.get('assays'), run.get('phenotypes'), run.get('planned_pairs'))
    if run.get('complete') is not True or run.get('status') != 'complete' or run.get('mode') != mode or actual != required or run.get('completed_traits') != required[1]:
        raise RuntimeError('Full declared MR is not complete; presentation disabled')
    return run


def prepare(mode, output=None):
    folder = ROOT / f'papers/01_mr_imaging_scan/results/full_atlas_20261116/{mode}'
    run = completed_source(folder, mode)  # Must precede any output creation.
    output = output or folder / 'presentation'
    if output.exists():
        raise FileExistsError(f'Presentation already exists; inspect before replacing: {output}')
    inputs = {folder / 'RUN.json': sha(folder / 'RUN.json')}
    for relative, digest in run['input_sha256'].items():
        inputs[ROOT / relative] = digest
    for relative, digest in run['source_code_sha256'].items():
        if sha(ROOT / relative) != digest:
            raise ValueError(f'Original statistical code identity changed: {relative}')
    inputs[folder / 'phenotypes.csv'] = run['phenotypes_sha256']
    for policy in ('drop', 'infer'):
        inputs[folder / f'full_family_palindrome_{policy}.csv.gz'] = run['policies'][policy]['output_sha256']
    for path, digest in inputs.items():
        if sha(path) != digest:
            raise ValueError(f'Input SHA mismatch: {path}')
    cov = pd.read_csv(ROOT / 'data/ukb_ppp/extracted/assay_coverage.csv')
    meta = pd.read_csv(folder / 'phenotypes.csv')
    instruments = pd.read_csv(ROOT / 'data/ukb_ppp/extracted/cis_primary_selected.csv')
    if len(cov) != run['assays'] or len(meta) != run['phenotypes']:
        raise ValueError('Plan dimensions disagree with completed run')
    annotations = instruments[['exposure_id', 'variant', 'chrom', 'pos37', 'pos38',
        'pav_consequence', 'pav_coding_consequence', 'pav_annotation_scope', 'F_x']]
    tables = {}
    audits = {}
    for policy in ('drop', 'infer'):
        rank, profiles, coverage, audit = summarize_table(
            folder / f'full_family_palindrome_{policy}.csv.gz', cov, meta, mode,
            run['policies'][policy])
        rank = rank.merge(annotations, on='exposure_id', how='left', validate='one_to_one', sort=False)
        tables[f'all_assay_ranking_{policy}.csv'] = rank
        tables[f'all_assay_category_profiles_{policy}.csv.gz'] = profiles
        tables[f'coverage_by_phenotype_{policy}.csv'] = coverage
        audits[policy] = audit
    primary = tables['all_assay_ranking_drop.csv']
    # This only determines an overview's display rows. It nominates no targets.
    display = primary.loc[primary.n_estimated.gt(0)].head(24).copy()
    display.insert(0, 'display_row', range(1, len(display) + 1))
    tables['overview_display_assays.csv'] = display
    for policy in ('drop', 'infer'):
        profiles = tables[f'all_assay_category_profiles_{policy}.csv.gz']
        selected = profiles.loc[~profiles.is_qc & profiles.exposure_id.isin(display.exposure_id)]
        tables[f'overview_category_profiles_{policy}.csv'] = selected.merge(
            display[['exposure_id', 'display_row']], on='exposure_id', validate='many_to_one')
    for path, digest in inputs.items():
        if sha(path) != digest:
            raise ValueError(f'Input changed during presentation preparation: {path}')
    output.mkdir(parents=True)
    for name, table in tables.items():
        table.to_csv(output / name, index=False, lineterminator='\n',
            compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)
    manifest = dict(complete=True, mode=mode, created_at_utc=datetime.now(timezone.utc).isoformat(),
        formal_full_atlas_result=mode == 'full', audits=audits,
        inputs={str(p.relative_to(ROOT)): h for p, h in inputs.items()},
        code_sha256={str(p.relative_to(ROOT)): sha(p) for p in
            (Path(__file__), ROOT / 'pdv3/atlas_presentation.py')},
        outputs={name: sha(output / name) for name in tables},
        ranking='Descriptive minimum log P over non-QC phenotypes; ties by phenotype ID, then assay ID. Every assay retained.',
        display_selection='First 24 estimable assays in primary descriptive ranking only; no historical additions; display is not nomination.',
        interpretation='Best P/q refers to one source pair, not an assay/category combined test; category sizes and phenotype dependence differ.',
        no_new_mr_estimation=True, no_new_multiplicity_correction=True,
        historical_identity_used=False, full_original_family_denominators_preserved=True,
        figures_rendered=False)
    (output / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return output, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['full', 'validation45'], default='full')
    args = parser.parse_args()
    output, manifest = prepare(args.mode)
    print(json.dumps(dict(output=str(output), mode=args.mode, audits=manifest['audits']), indent=2))


if __name__ == '__main__':
    main()
