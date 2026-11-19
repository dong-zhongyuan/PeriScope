"""Unfiltered public CD40 / bilateral pallidum regions, with source checks.

This is a distinct SCALLOP exposure, not a substitute labelled UKB-PPP. The
region is fixed around the historical ST9 CD40 lead, before regional results.
"""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.tabix_http import save_json

DATA = ROOT / 'data/scallop_cd40'
OUT = ROOT / 'papers/01_mr_imaging_scan/results/public_cd40_20261116'
CENTER = 44747947
START, END = CENTER - 1000000, CENTER + 1000000


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def write(df, name):
    path = OUT / name
    df.to_csv(path, index=False, compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)
    return {'path': str(path.relative_to(ROOT)), 'rows': len(df), 'sha256': sha(path)}


def key(chrom, pos, a, b):
    if a not in 'ACGT' or b not in 'ACGT' or len(a) != 1 or len(b) != 1 or a == b:
        return None
    return f'{chrom}:{pos}:' + '_'.join(sorted([a, b]))


def extract_protein():
    path = DATA / 'download/CD40.txt.gz'
    receipt = json.loads((DATA / 'DOWNLOAD_MANIFEST.json').read_text())
    assert path.stat().st_size == receipt['size_bytes']
    assert sha(path) == receipt['sha256']
    rows, counts = [], Counter()
    with gzip.open(path, 'rt') as f:
        header = next(f).split()
        assert header == 'MarkerName Allele1 Allele2 Freq1 FreqSE Effect StdErr P-value Direction TotalSampleSize'.split()
        for line in f:
            v = line.split()
            if len(v) != len(header): raise ValueError('Malformed SCALLOP row')
            counts['source_rows'] += 1
            chrom, pos, alleles = v[0].split(':')
            if chrom != '20' or not START <= int(pos) <= END: continue
            r = dict(zip(header, v))
            a, b = alleles.upper().split('_')
            r.update(chrom=chrom, pos37=int(pos), marker_a=a, marker_b=b,
                     ea_x=v[1].upper(), oa_x=v[2].upper(), beta_x=float(v[5]),
                     se_x=float(v[6]), p_x=float(v[7]), eaf_x=float(v[3]), n_x=float(v[9]),
                     variant_key=key(chrom, int(pos), a, b), build='GRCh37')
            if not np.isfinite([r['beta_x'], r['se_x'], r['eaf_x'], r['n_x']]).all() or r['se_x'] <= 0:
                raise ValueError('Invalid regional protein statistics')
            if not 0 <= r['p_x'] <= 1 or not 0 <= r['eaf_x'] <= 1 or r['n_x'] <= 0:
                raise ValueError('Invalid regional protein frequency / p / n')
            # Marker alleles define genomic identity; METAL's output alleles
            # need only be compatible. Coloc ABFs use squared Z, not direction.
            comp = str.maketrans('ACGT', 'TGCA')
            if r['variant_key'] and set([r['ea_x'], r['oa_x']]) not in [set([a,b]), set([a.translate(comp), b.translate(comp)])]:
                raise ValueError('METAL alleles disagree with marker identity')
            rows.append(r)
    counts.update(regional_rows=len(rows), gzip_crc_and_eof_verified=1)
    receipt.update(gzip_crc='verified by reading the full source to EOF',source_rows=counts['source_rows'])
    save_json(DATA / 'DOWNLOAD_MANIFEST.json',receipt)
    return pd.DataFrame(rows), dict(counts), receipt


def extract_imaging(idp, manifest):
    source = next(r for r in manifest['datasets'] if r['phenotype']['idp'] == idp)
    path = ROOT / f'data/big40/download/stats/{idp}.txt.gz'
    expected = source['download']['sha256']
    assert sha(path) == expected
    rows, total = [], 0
    with gzip.open(path, 'rt') as f:
        assert next(f).split() == 'chr rsid pos a1 a2 beta se pval(-log10)'.split()
        for line in f:
            total += 1
            v = line.split()
            if len(v) != 8: raise ValueError('Malformed Oxford row')
            if v[0] != '20' or not START <= int(v[2]) <= END: continue
            r = dict(chrom='20', pos37=int(v[2]), variant_y=v[1], oa_y=v[3], ea_y=v[4],
                     beta_y=float(v[5]), se_y=float(v[6]), neg_log10_p_y=float(v[7]),
                     n_y=source['phenotype']['n_autosomal'], build='GRCh37',
                     variant_key=key('20', int(v[2]), v[3], v[4]))
            if not np.isfinite([r['beta_y'],r['se_y'],r['neg_log10_p_y']]).all() or r['se_y'] <= 0 or r['neg_log10_p_y'] < 0:
                raise ValueError('Invalid regional imaging statistics')
            rows.append(r)
    assert total == 17103079
    return pd.DataFrame(rows), dict(source_rows=total,regional_rows=len(rows),
                                   gzip_crc_and_eof_verified=True,source_sha256=expected,
                                   source_url=source['phenotype']['source_y'])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inst = pd.read_csv(ROOT / 'data/ukb_ppp/extracted/cis_primary_selected.csv')
    lead = inst[inst.gene.eq('CD40')]
    assert len(lead) == 1 and lead.iloc[0].variant == 'rs4810485' and int(lead.iloc[0].pos37) == CENTER
    manifest = json.loads((ROOT / 'data/big40/PHASE1_17_MANIFEST.json').read_text())
    x, x_audit, receipt = extract_protein()
    x_audit['output'] = write(x, 'protein_region_all_rows.csv.gz')
    x_snps = x[x.variant_key.notna()].copy()
    if x_snps.variant_key.duplicated().any(): raise ValueError('Duplicate protein SNP identity')
    audits = {}
    for idp in ['0017', '0018']:
        y, audit = extract_imaging(idp, manifest)
        audit['raw_output'] = write(y, f'{idp}_imaging_region_all_rows.csv.gz')
        y_snps = y[y.variant_key.notna()].copy()
        if y_snps.variant_key.duplicated().any(): raise ValueError('Duplicate imaging SNP identity')
        union = x_snps.merge(y_snps, on='variant_key', suffixes=('_protein','_imaging'), how='outer', indicator=True, validate='one_to_one')
        audit['match_counts'] = union['_merge'].value_counts().to_dict()
        audit['unsupported_protein_rows'] = int(x.variant_key.isna().sum())
        audit['unsupported_imaging_rows'] = int(y.variant_key.isna().sum())
        audit['union_output'] = write(union, f'{idp}_variant_identity_audit.csv.gz')
        shared = union[union['_merge'].eq('both')].copy()
        assert shared.pos37_protein.eq(shared.pos37_imaging).all()
        shared['pos37'] = shared.pos37_protein.astype(int)
        shared['idp'] = idp
        shared['palindromic'] = shared.apply(lambda r: set([r.marker_a,r.marker_b]) in [{'A','T'},{'C','G'}], axis=1)
        shared['coloc_input'] = True
        # No palindromic sign assumption is needed for ABFs. Direction plots
        # must omit these SNPs unless separately frequency-resolved.
        shared['direction_status'] = np.where(shared.palindromic,'not_used_for_signed_direction','requires_allele_alignment')
        assert len(shared) > 100
        audit['coloc_input'] = write(shared.sort_values(['pos37','variant_key']), f'{idp}_shared_snps.csv.gz')
        audit['n_x_min'] = float(shared.n_x.min()); audit['n_x_max'] = float(shared.n_x.max())
        audit['historical_lead'] = shared[shared.variant_y.eq('rs4810485')][['variant_key','ea_x','oa_x','beta_x','se_x','n_x','ea_y','oa_y','beta_y','se_y']].to_dict('records')
        if len(audit['historical_lead']) != 1: raise ValueError('Historical lead not uniquely matched')
        audits[idp] = audit
        print(idp, audit['match_counts'], flush=True)
    save_json(OUT / 'INPUT_AUDIT.json', {
        'complete': True, 'region': {'chromosome':20,'build':'GRCh37','center':CENTER,'start':START,'end':END,'center_reason':'Historical CD40 rs4810485 selected before regional analysis'},
        'exposure':'SCALLOP CVD1 CD40; distinct from UKB-PPP',
        'exposure_build_evidence': ['https://doi.org/10.1038/s43856-024-00530-x (Methods, SCALLOP explicitly hg19)', 'Exact genomic allele / coordinate matching to original GRCh37 Oxford source; historical rs4810485 identity checked'],
        'effect_scale_evidence':'https://pmc.ncbi.nlm.nih.gov/articles/PMC7611474/ (genetic analyses: inverse-normal transformed and/or unit variance); Oxford BIG40 deconfounded SD',
        'effect_allele_evidence':'https://github.com/statgen/METAL/blob/master/metal/Main.cpp (Effect is for allele1)',
        'source_download':receipt,'protein':x_audit,'outcomes':audits,
        'matching':'Same GRCh37 chromosome, position, exact unordered genomic marker alleles; biallelic SNVs only; no rsID-only or P-based matching',
        'p_value_filter':None,'palindromes':'Included only in sign-invariant ABF analysis after genomic identity match; not used for signed direction',
        'source_code_sha256':sha(Path(__file__))})


if __name__ == '__main__': main()
