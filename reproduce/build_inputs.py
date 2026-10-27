"""Rebuild six model inputs from public source files and bundled frozen annotations.

python reproduce/build_inputs.py --output /large-disk/periscope-inputs --rscript Rscript
Use --raw-cache /path/to/raw to reuse SHA-256-verified public downloads. No existing
matrices or old build reports are needed. Use a new output directory, or --resume
to reuse completed stages after checking every recorded output hash.
Python dependencies: reproduce/requirements-inputs.txt; R: SeuratObject and Matrix.
Cell-state assignments and donor partitions are frozen inputs included in the
repository. This command does not rerun the historical state-clustering procedure.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import urllib.request

REPO = Path(__file__).resolve().parents[1]
RESOURCES = REPO / 'reproduce/input_resources'
STAGES = REPO / '_work/rescreen_20261002'


def sha256(path, compressed=False):
    h = hashlib.sha256()
    with (gzip.open(path, 'rb') if compressed else Path(path).open('rb')) as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def acquire(spec, raw, cache, offline):
    dest = raw / spec['path']
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file():
        if dest.stat().st_size != spec['bytes'] or sha256(dest) != spec['sha256']:
            raise ValueError(f'Existing source checksum mismatch: {dest}')
        print(f'Verified existing source: {spec["path"]}', flush=True)
        return 'verified-existing'
    cached = cache / spec['path'] if cache else None
    if cached and cached.is_file():
        if cached.stat().st_size != spec['bytes'] or sha256(cached) != spec['sha256']:
            raise ValueError(f'Cached source checksum mismatch: {cached}')
        dest.symlink_to(cached)
        print(f'Verified cached source: {spec["path"]}', flush=True)
        return 'verified-cache'
    if offline:
        raise FileNotFoundError(f'Offline cache missing {spec["path"]}')
    part = dest.with_name(dest.name + '.part')
    req = urllib.request.Request(spec['url'], headers={'User-Agent': 'PeriScope-reproduction/1.0'})
    print(f'Downloading {spec["url"]}', flush=True)
    with urllib.request.urlopen(req, timeout=120) as response, part.open('wb') as f:
        shutil.copyfileobj(response, f, length=8 * 1024 * 1024)
    if part.stat().st_size != spec['bytes'] or sha256(part) != spec['sha256']:
        raise ValueError(f'Download checksum mismatch: {part}')
    part.rename(dest)
    return 'downloaded-and-verified'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--raw-cache', type=Path, help='Read-only source cache with the manifest-relative paths')
    ap.add_argument('--offline', action='store_true', help='Require every source file in --raw-cache')
    ap.add_argument('--rscript', default='Rscript')
    ap.add_argument('--resume', action='store_true', help='Verify and reuse completed stages in this output')
    args = ap.parse_args()
    out = args.output.expanduser().resolve()
    if out.exists() and not args.resume:
        raise FileExistsError(f'Choose a new output directory: {out}')
    if args.resume and not (out / 'completed_stages.json').is_file():
        raise FileNotFoundError('Resume requires an output containing completed_stages.json')
    if REPO == out or REPO in out.parents:
        raise ValueError('Use an output directory outside the checkout (large generated data)')
    manifest = json.loads((REPO / 'reproduce/input_sources.json').read_text())
    # Fail before acquiring the large files if the frozen metadata is incomplete.
    for name, spec in manifest['resources'].items():
        path = RESOURCES / name
        compressed = name == 'hgnc_complete_set.txt'
        if compressed:
            path = path.with_suffix('.txt.gz')
        if sha256(path, compressed) != spec['sha256']:
            raise ValueError(f'Frozen resource checksum mismatch: {path}')
    if not shutil.which(args.rscript):
        raise FileNotFoundError(f'Rscript not found: {args.rscript}')
    subprocess.run([args.rscript, '-e', 'library(SeuratObject); library(Matrix)'], check=True)
    out.mkdir(parents=True, exist_ok=args.resume)
    raw = out / 'raw'
    cache = args.raw_cache.expanduser().resolve() if args.raw_cache else None
    source_status = {spec['path']: acquire(spec, raw, cache, args.offline) for spec in manifest['raw']}
    (out / 'source_verification.json').write_text(json.dumps(source_status, indent=2) + '\n')
    env = os.environ.copy()
    env['PYTHONPATH'] = str(REPO / 'src') + os.pathsep + env.get('PYTHONPATH', '')
    env.setdefault('OMP_NUM_THREADS', '4')
    env.setdefault('OPENBLAS_NUM_THREADS', '4')
    env.setdefault('NUMBA_NUM_THREADS', '4')
    env.setdefault('MKL_NUM_THREADS', '4')
    logs = out / 'logs'
    logs.mkdir(exist_ok=args.resume)
    receipts = out / 'stage_receipts'
    receipts.mkdir(exist_ok=args.resume)
    completed = json.loads((out / 'completed_stages.json').read_text()) if args.resume else []

    def run(label, command, outputs=()):
        receipt = receipts / (label + '.json')
        run_manifest = out / 'runs' / label / 'run_manifest.json'
        if label in completed and label != 'verify-six-inputs':
            recorded = json.loads((receipt if receipt.exists() else run_manifest).read_text())
            if recorded.get('status') != 'completed' or not recorded.get('outputs'):
                raise ValueError(f'Missing completed-stage output inventory: {label}')
            for item in recorded['outputs']:
                path = Path(item['path'])
                if out not in path.resolve().parents or path.stat().st_size != item['bytes'] or sha256(path) != item['sha256']:
                    raise ValueError(f'Completed output changed: {path}')
            print(f'REUSE verified {label}', flush=True)
            return
        if args.resume and run_manifest.exists():
            if json.loads(run_manifest.read_text()).get('status') == 'completed':
                raise ValueError(f'Completed manifest lacks stage receipt: {run_manifest}')
            run_manifest.unlink()  # This run's failed stage; allow its original entrypoint to retry.
        print(f'START {label}', flush=True)
        with (logs / (label + '.log')).open('w') as log:
            result = subprocess.run([str(x) for x in command], cwd=REPO, env=env,
                                    stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f'{label} failed; see {logs / (label + ".log")}')
        if run_manifest.exists():
            recorded = json.loads(run_manifest.read_text())
        else:
            recorded = {'status': 'completed', 'outputs': [
                {'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha256(path)} for path in outputs]}
        receipt.write_text(json.dumps(recorded, indent=2) + '\n')
        if label not in completed:
            completed.append(label)
        (out / 'completed_stages.json').write_text(json.dumps(completed, indent=2) + '\n')
        print(f'DONE {label}', flush=True)

    def python_stage(label, script, *options, outputs=()):
        run(label, [sys.executable, '-u', REPO / script, *options], outputs=outputs)

    processed = {}
    for cohort in ('gse178265', 'gse223138', 'mb2025'):
        interim = out / f'{cohort}_unfiltered.h5ad'
        annotated = out / f'{cohort}_annotated.h5ad'
        if cohort == 'gse178265':
            options = ['--raw-dir', raw / 'GSE178265/geo_original', '--meta-json', RESOURCES / 'GSE178265_human_meta.json']
        elif cohort == 'gse223138':
            options = ['--raw-tar', raw / 'GSE223138/geo_original/GSE223138_RAW.tar',
                       '--extract-dir', out / 'gse223138_extracted', '--meta-json', RESOURCES / 'GSE223138_human_meta.json']
        else:
            extract = out / 'mb2025_extracted'
            if args.resume and extract.exists() and 'export-mb2025' not in completed:
                expected = {'matrix.mtx', 'features.tsv', 'barcodes.tsv', 'metadata.csv'}
                if extract.is_symlink() or any(not p.is_file() or p.is_symlink() or p.name not in expected for p in extract.iterdir()):
                    raise ValueError(f'Unexpected files in failed export; inspect {extract}')
                for path in extract.iterdir():
                    path.unlink()  # Only this run's incomplete export, never the source RDS.
                extract.rmdir()
            run('export-mb2025', [args.rscript, REPO / 'reproduce/export_mb2025.R',
                                 raw / 'moquin_beaudry_2025/scPD_AllCells_Annotated.rds', extract],
                outputs=[extract / name for name in ('matrix.mtx', 'features.tsv', 'barcodes.tsv', 'metadata.csv')])
            options = ['--extract-dir', extract]
        python_stage('ingest-' + cohort, f'scripts/preprocess_{cohort}.py', *options,
                     '--out-h5ad', interim, '--report-json', out / f'{cohort}_ingest.json',
                     '--run-dir', out / 'runs' / ('ingest-' + cohort))
        python_stage('annotate-' + cohort, f'scripts/annotate_{cohort}.py',
                     '--in-h5ad', interim, '--out-h5ad', annotated,
                     '--report-json', out / f'{cohort}_annotation.json',
                     '--run-dir', out / 'runs' / ('annotate-' + cohort))
        processed[cohort] = str(annotated)
    tensors = out / 'inputs'
    config = {
        'inputs': {key: {'path': processed[cohort]} for key, cohort in
                   [('brain', 'gse178265'), ('blood223', 'gse223138'), ('bloodmb', 'mb2025')]},
        'purify_dir': str(RESOURCES / 'purification'), 'out_dir': str(tensors),
        'cite_raw_dir': str(raw / 'citeseq_hao'),
    }
    config['inputs']['split'] = {'path': str(RESOURCES / 'GSE178265_sn_donor_split_v1.json')}
    config_path = out / 'input_config.json'
    config_path.write_text(json.dumps(config, indent=2) + '\n')
    env['PERISCOPE_INPUT_CONFIG'] = str(config_path)
    stage_outputs = {
        'tissue-features': [tensors / 'tissue_features.json'],
        'cite-bridge': [tensors / name for name in ('bridge_data.npz', 'pairing_audit.json', 'protected_feature_sources.json')],
        'model-inputs': [tensors / name for name in (*[n for n in manifest['expected_arrays'] if n != 'bridge_data.npz'], 'build_report.json')],
    }
    for label, script in [('tissue-features', 'input_utils.py'), ('cite-bridge', 'preprocess_citeseq_hao.py'), ('model-inputs', 'run_build.py')]:
        run(label, [sys.executable, '-u', STAGES / script], outputs=stage_outputs[label])
    python_stage('verify-six-inputs', 'reproduce/verify_inputs.py', '--inputs', tensors,
                 '--report', out / 'input_verification.json', outputs=[out / 'input_verification.json'])
    print(f'PASS: all six inputs rebuilt and every array matched. Output: {tensors}', flush=True)


if __name__ == '__main__':
    main()
