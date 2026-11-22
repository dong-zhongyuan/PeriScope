#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Orientation of TFRC's instrument on the clinical soluble transferrin receptor scale (Figure 4d,i; Methods).

rs3817672 encodes TfR1 Gly142Ser. Affinity-proteomic assays read its Ser142 (T) allele as lowering plasma TFRC
(UKB-PPP beta = -0.524 SD); the clinical soluble transferrin receptor assay reads it as raising the receptor. This
writes, for the TFRC variants that carry the nigral iron signal and the independent clinical-assay leads, the
clinical-assay effect next to the nigral QSM and T2* effects, all on the same effect allele. Figure 4i scales the
affinity-proteomic Wald ratios of TFRC by beta_olink / beta_clinical at rs3817672 from this table.

Inputs
  soluble transferrin receptor GWAS, GCST90451684 (Allara et al. 2024; CROATIA-VIS, FinDonor, InCHIANTI, INTERVAL,
    KORA; n = 45,330; GRCh37), cut to TFRC +/- 500 kb:
      curl -s https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90451001-GCST90452000/GCST90451684/\
GCST90451684.tsv.gz | zcat | awk -F'\\t' 'NR==1 || ($1==3 && $2>=195250000 && $2<=196310000)'
  nigral QSM (0015 left, 0016 right) and T2* (0033 left, 0034 right) from the Oxford QSM/T2* discovery GWAS
    (Wang et al. 2022), cut to the same window by extract_qsm_tfrc.py; T2* signs are reversed so that positive
    means more iron.
  1000 Genomes phase 3 EUR reference (plink bed) for r2 < 0.01 clumping of the clinical-assay leads.
Usage: python tfrc_orientation.py <stfr_window.tsv> <qsm_region_dir> <EUR_bfile_prefix> <out_dir>
"""
import csv
import gzip
import json
import os
import sys

import numpy as np

STFR, QDIR, REF, OUT = sys.argv[1:5]
TFRC = (195754053, 195809060)            # GRCh37 gene body
WIN = (TFRC[0] - 500000, TFRC[1] + 500000)
IMG = [('QSM_L_SN', '0015', 1), ('QSM_R_SN', '0016', 1), ('T2s_L_SN', '0033', -1), ('T2s_R_SN', '0034', -1)]
SHOW = ['rs3817672', 'rs406271', 'rs111625920', 'rs9867473']
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

X = {}
for r in csv.DictReader(open(STFR), delimiter='\t'):
    pos = int(r['base_pair_location'])
    if WIN[0] <= pos <= WIN[1] and r['rs_id'].startswith('rs'):
        X[r['rs_id']] = {'pos': pos, 'ea': r['effect_allele'].upper(), 'oa': r['other_allele'].upper(),
                         'b': float(r['beta']), 's': float(r['standard_error']),
                         'f': float(r['effect_allele_frequency']), 'p': float(r['p_value'])}
Y = {}
for name, code, sign in IMG:
    d = {}
    with gzip.open('%s/%s_chr3_tfrc.tsv.gz' % (QDIR, code), 'rt') as fh:
        fh.readline()
        for line in fh:
            c, rs, pos, a1, a2, b, s, nlp = line.split()
            d[rs] = (a2.upper(), a1.upper(), sign * float(b), 10 ** -float(nlp))      # a2 is the effect allele
    Y[name] = d


def aligned(rs, name):
    o = Y[name].get(rs)
    if o is None:
        return None
    oe, oo, b, p = o
    ea, oa = X[rs]['ea'], X[rs]['oa']
    if {oe, oo} == {ea, oa} and COMP.get(ea) != oa:
        return (b if oe == ea else -b), p
    return None


idx = {}
with open(REF + '.bim') as fh:
    for i, line in enumerate(fh):
        f = line.split()
        if f[0] == '3' and f[1] in X:
            idx[f[1]] = i
nsamp = sum(1 for _ in open(REF + '.fam'))
NB = (nsamp + 3) // 4
bed = open(REF + '.bed', 'rb')
assert bed.read(3) == b'\x6c\x1b\x01'
_c = {}


def dos(rs):
    if rs not in _c:
        bed.seek(3 + idx[rs] * NB)
        b = np.frombuffer(bed.read(NB), dtype=np.uint8)
        codes = ((b[:, None] >> np.array([0, 2, 4, 6], dtype=np.uint8)) & 3).ravel()[:nsamp]
        _c[rs] = np.choose(codes, [2.0, np.nan, 1.0, 0.0])
    return _c[rs]


def r2(a, b):
    x, y = dos(a), dos(b)
    m = ~(np.isnan(x) | np.isnan(y))
    return float(np.corrcoef(x[m], y[m])[0, 1] ** 2)


leads = []
for rs in sorted((r for r in X if X[r]['p'] < 5e-8 and r in idx), key=lambda r: X[r]['p']):
    if all(r2(rs, k) < 0.01 for k in leads):
        leads.append(rs)

rows = []
for rs in SHOW + [r for r in leads if r not in SHOW]:
    x = X[rs]
    row = {'variant': rs, 'pos_grch37': x['pos'], 'ea': x['ea'], 'eaf': x['f'],
           'sTfR_clinical_b': x['b'], 'sTfR_p': x['p']}
    for name, _, _ in IMG:
        h = aligned(rs, name)
        row[name + '_b'] = '' if h is None else h[0]
        row[name + '_p'] = '' if h is None else h[1]
    rows.append(row)
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, 'tfrc_variant_panel.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
json.dump({'window_grch37': WIN, 'n_variants': len(X), 'clinical_assay_leads_r2_0.01': leads},
          open(os.path.join(OUT, 'tfrc_orientation.run.json'), 'w'), indent=1)
for r in rows:
    print('%-11s %s  sTfR %+.4f (p %.1e)  QSM L %s  R %s  T2* L %s  R %s' % (
        r['variant'], r['ea'], r['sTfR_clinical_b'], r['sTfR_p'],
        *['%+.4f' % r[n + '_b'] if r[n + '_b'] != '' else 'NA' for n in ('QSM_L_SN', 'QSM_R_SN', 'T2s_L_SN', 'T2s_R_SN')]))
