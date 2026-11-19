"""Acquire a declared 17-IDP discovery panel and retain all selected cis loci.

Input filters use exposure coordinates only. Complete gzip files are retained;
all rows are parsed and validated, including EOF/CRC. No outcome-P selection.
"""
import argparse
import concurrent.futures
import gzip
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

from common import ROOT
from s07_acquire_big40 import BASE, acquire, parse_index, sha256, atomic_json, utc

IDS = [f'{i:04d}' for i in range(11, 19)] + ['0342'] + [f'{i:04d}' for i in range(1438, 1446)]
DEFAULT_EXPOSURES = ROOT / 'data/ukb_ppp/extracted/cis_primary_selected.csv'
DEFAULT_DATA = ROOT / 'data/big40'


def chrom(series):
    return (series.astype(str).str.upper().str.removeprefix('CHR').str.removesuffix('.0')
            .str.lstrip('0').replace({'23': 'X'}))


def extract_loci(path, instruments_path, output, meta):
    """Parse full source once; write exact coordinates without allele/P filtering."""
    instruments = pd.read_csv(instruments_path)
    targets = set(zip(chrom(instruments.chrom), instruments.pos37.astype(int)))
    positions = {p for c, p in targets}
    pieces, rows, errors = [], 0, {}
    chr_counts, first_row, last_row = {}, None, None
    source = meta['source_y']
    with gzip.open(path, 'rt', encoding='ascii') as handle:
        expected = ['chr', 'rsid', 'pos', 'a1', 'a2', 'beta', 'se', 'pval(-log10)']
        actual = handle.readline().split()
        if actual != expected:
            raise ValueError(f'Unexpected header for {meta["phenotype"]}: {actual}')
        for chunk in pd.read_csv(handle, sep=r'\s+', names=expected, chunksize=300000,
                                 dtype={'chr': str, 'rsid': str, 'a1': str, 'a2': str}):
            ch = chrom(chunk['chr'])
            for name, mask in {
                'bad_coordinate': ~ch.isin([str(i) for i in range(1, 23)] + ['X']) | (chunk.pos <= 0),
                'bad_beta': ~np.isfinite(chunk.beta),
                'bad_se': ~np.isfinite(chunk.se) | (chunk.se <= 0),
                'bad_logp': ~np.isfinite(chunk['pval(-log10)']) | (chunk['pval(-log10)'] < 0),
                'bad_alleles': chunk.a1.isna() | chunk.a2.isna() | chunk.a1.eq(chunk.a2),
            }.items():
                errors[name] = errors.get(name, 0) + int(mask.sum())
            for c, count in ch.value_counts().items():
                chr_counts[c] = chr_counts.get(c, 0) + int(count)
            if first_row is None:
                first_row = chunk.iloc[0].to_dict()
            last_row = chunk.iloc[-1].to_dict()
            subset = chunk.loc[chunk.pos.isin(positions)].copy()
            subset['chrom_y'] = ch.loc[subset.index]
            subset = subset.loc[[p in targets for p in zip(subset.chrom_y, subset.pos)]]
            subset['source_data_row'] = subset.index.to_numpy() + 1
            pieces.append(subset)
            rows += len(chunk)
    if rows == 0 or any(errors.values()):
        raise ValueError(f'Invalid source rows in {source}: {errors}; rows={rows}')
    matched = pd.concat(pieces, ignore_index=True)
    out = pd.DataFrame({
        'phenotype': meta['phenotype'], 'trait': meta['trait'], 'family': meta['family'],
        'variant_y': matched.rsid, 'chrom_y': matched.chrom_y, 'pos_y': matched.pos,
        'ea_y': matched.a2, 'oa_y': matched.a1, 'beta_y': matched.beta, 'se_y': matched.se,
        'p_y': np.power(10., -matched['pval(-log10)']),
        'neg_log10_p_y': matched['pval(-log10)'], 'eaf_y': np.nan,
        'n_y': np.where(matched.chrom_y.eq('X'), np.nan, meta['n_autosomal']),
        'n_y_status': np.where(matched.chrom_y.eq('X'), 'X_PAR_nonPAR_not_assigned', 'catalog_autosomal_n'),
        'build_y': 'GRCh37', 'units_y': 'SD of deconfounded imaging phenotype',
        'source_y': source, 'source_data_row': matched.source_data_row,
    })
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix('.csv.tmp')
    out.to_csv(temporary, index=False)
    temporary.replace(output)
    validation = {
        'gzip_crc_and_eof_verified': True, 'source_rows': rows, 'source_header': actual,
        'source_chromosome_counts': chr_counts, 'source_error_counts': errors,
        'selected_target_loci': len(targets), 'matched_source_rows': len(out),
        'matched_target_loci': len(set(zip(out.chrom_y, out.pos_y))),
        'outcome_p_filter': None, 'allele_filter_at_extraction': None,
        'first_source_row': first_row, 'last_source_row': last_row,
        'source_sha256': sha256(path), 'instrument_sha256': sha256(instruments_path),
        'extracted_path': str(output.relative_to(ROOT)), 'extracted_sha256': sha256(output),
    }
    atomic_json(output.with_suffix('.validation.json'), validation)
    return validation


def work(meta, instruments_path, data_dir):
    output = data_dir / 'extracted/phase1_17' / (meta['idp'] + '.cis_loci.csv')
    target = data_dir / 'download/stats' / (meta['idp'] + '.txt.gz')
    record = acquire(meta['source_y'], target,
                     lambda p: extract_loci(p, instruments_path, output, meta))
    validation_path = output.with_suffix('.validation.json')
    validation = json.loads(validation_path.read_text()) if validation_path.exists() else None
    if not (validation and output.exists() and
            validation['source_sha256'] == record['sha256'] and
            validation['instrument_sha256'] == sha256(instruments_path) and
            validation['extracted_sha256'] == sha256(output)):
        validation = extract_loci(target, instruments_path, output, meta)
    print(f'EXTRACTED {meta["idp"]}: {validation["matched_source_rows"]} rows at '
          f'{validation["matched_target_loci"]} of {validation["selected_target_loci"]} loci', flush=True)
    return {'phenotype': meta, 'download': record, 'extraction': validation}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workers', type=int, default=3)
    p.add_argument('--exposures', type=Path, default=DEFAULT_EXPOSURES)
    p.add_argument('--data-dir', type=Path, default=DEFAULT_DATA)
    args = p.parse_args()
    if not 1 <= args.workers <= 4:
        p.error('workers must be 1..4')
    index = args.data_dir / 'download/IDPs_v4.html'
    catalog = {row['idp']: row for row in parse_index(index)}
    metas = []
    for idp in IDS:
        row = catalog[idp]
        metas.append(dict(phenotype=f'BIG40:release2:discovery:{idp}', idp=idp,
                          trait=row['short_name'], family='T2star' if 'T2star' in row['short_name'] else 'volume',
                          n_autosomal=int(row['n_discovery']), build_y='GRCh37',
                          source_y=BASE + f'release2/stats/{idp}.txt.gz',
                          units_y='SD of deconfounded imaging phenotype',
                          source_catalog_sha256=sha256(index)))
    plan = args.data_dir / 'phase1_17_phenotypes.csv'
    pd.DataFrame(metas).to_csv(plan, index=False)
    started = utc()
    outputs = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs = [pool.submit(work, m, args.exposures, args.data_dir) for m in metas]
        for job in concurrent.futures.as_completed(jobs):
            outputs.append(job.result())
    outputs.sort(key=lambda row: row['phenotype']['idp'])
    combined = pd.concat([pd.read_csv(ROOT / r['extraction']['extracted_path']) for r in outputs], ignore_index=True)
    output = args.data_dir / 'extracted/phase1_17_outcomes.csv'
    combined.to_csv(output, index=False)
    record = {'started_at_utc': started, 'finished_at_utc': utc(), 'complete': len(outputs) == 17,
              'source': 'BIG40 release2 discovery', 'selection': 'original 16 IDPs plus midbrain volume 0342',
              'outcome_driven_selection': False, 'n_phenotypes': len(outputs), 'extracted_rows': len(combined),
              'plan_sha256': sha256(plan), 'exposure_sha256': sha256(args.exposures),
              'output': str(output.relative_to(ROOT)), 'output_sha256': sha256(output),
              'script_sha256': sha256(__file__), 'command': sys.argv, 'datasets': outputs}
    atomic_json(args.data_dir / 'PHASE1_17_MANIFEST.json', record)
    print(json.dumps({k: v for k, v in record.items() if k != 'datasets'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
