"""Independently verify overview display records and physical vector exports."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[3]


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return pd.read_csv(path, float_precision='round_trip')


def exact(left, right, keys):
    pd.testing.assert_frame_equal(left.sort_values(keys).reset_index(drop=True),
        right.sort_values(keys).reset_index(drop=True), check_dtype=False,
        check_exact=True)


def check(output):
    manifest = json.loads((output / 'MANIFEST.json').read_text())
    assert manifest['complete'] and manifest['render_exit_code'] == 0
    assert manifest['formal_full_atlas_result'] == (manifest['mode'] == 'full')
    for relative, expected in {**manifest['inputs_sha256'], **manifest['code_sha256']}.items():
        assert sha(ROOT / relative) == expected, relative
    for relative, expected in manifest['outputs_sha256'].items():
        assert sha(output / relative) == expected, relative
    source = output.parent / 'presentation'
    data = output / 'plot_data'
    config = json.loads((data / 'CONFIG.json').read_text())
    all_assays = read(data / 'A_all_assays.csv')
    profiles = read(data / 'B_profiles.csv')
    display = read(data / 'B_display_rows.csv')
    strength = read(data / 'C_instrument_strength.csv')
    coverage = read(data / 'D_phenotype_coverage.csv')
    categories = read(data / 'categories.csv')
    maxima = {}
    original_ranks = {}
    for policy in ('drop', 'infer'):
        rank = read(source / f'all_assay_ranking_{policy}.csv')
        original_ranks[policy] = rank
        actual = all_assays.loc[all_assays.policy.eq(policy)]
        exact(actual[rank.columns], rank, ['exposure_id'])
        valid = actual.n_estimated.gt(0)
        logp = np.array([-float(v) / math.log(10) for v in actual.loc[valid, 'best_log_p']])
        z = np.array([float(b) / float(s) for b, s in
                      zip(actual.loc[valid, 'best_beta'], actual.loc[valid, 'best_se'])])
        maxima[policy] = dict(minus_log10_p=float(np.max(np.abs(logp - actual.loc[valid, 'minus_log10_p']))),
                              best_z=float(np.max(np.abs(z - actual.loc[valid, 'best_z']))))
        assert max(maxima[policy].values()) < 1e-12
        assert actual.loc[~valid, ['minus_log10_p', 'best_z']].isna().all().all()
        profile = read(source / f'overview_category_profiles_{policy}.csv')
        exact(profiles.loc[profiles.policy.eq(policy), profile.columns], profile,
              ['exposure_id', 'imaging_family'])
        cov = read(source / f'coverage_by_phenotype_{policy}.csv')
        expected_cov = []
        for phenotype, group in cov.groupby('phenotype'):
            assert group.family.nunique() == 1 and group.is_qc.nunique() == 1
            expected_cov.append(dict(phenotype=phenotype, family=group.family.iloc[0],
                policy=policy, n_estimated=int(group.loc[group.status.eq('estimated'), 'n_assays'].sum()),
                n_planned=int(group.n_assays.sum()), is_qc=bool(group.is_qc.iloc[0])))
        expected_cov = pd.DataFrame(expected_cov)
        expected_cov['estimated_fraction'] = expected_cov.n_estimated / expected_cov.n_planned
        exact(coverage.loc[coverage.policy.eq(policy)], expected_cov, ['phenotype'])
    top = original_ranks['drop'].loc[original_ranks['drop'].n_estimated.gt(0)].head(24).copy()
    top.insert(0, 'display_row', range(1, len(top) + 1))
    exact(display, top, ['display_row'])
    order = original_ranks['infer'].dropna(subset=['chrom'])
    order = sorted(order.to_dict('records'), key=lambda d:(d['chrom'], d['pos37'], d['exposure_id']))
    expected_order = {row['exposure_id']: i + 1 for i, row in enumerate(order)}
    for row in all_assays.to_dict('records'):
        if row['exposure_id'] in expected_order:
            assert row['genome_order'] == expected_order[row['exposure_id']]
        else:
            assert pd.isna(row['genome_order'])
    chosen = {}
    for policy in ('drop', 'infer'):
        for row in original_ranks[policy].to_dict('records'):
            if row['n_estimated'] > 0 and row['exposure_id'] not in chosen:
                chosen[row['exposure_id']] = (row, policy)
    assert len(strength) == len(chosen) and not strength.exposure_id.duplicated().any()
    for row in strength.to_dict('records'):
        expected, policy = chosen[row['exposure_id']]
        assert row['policy'] == policy
        assert row['availability'] == ('Primary' if policy == 'drop' else 'EAF-only')
        for field in ['F_x', 'best_log_p', 'best_beta', 'best_se', 'best_phenotype']:
            assert row[field] == expected[field]
        assert abs(row['absolute_best_z'] - abs(expected['best_beta'] / expected['best_se'])) < 1e-12
    assert not profiles.is_qc.any()
    assert coverage.is_qc.eq(coverage.family.eq('QC')).all()
    for row in categories.to_dict('records'):
        counts = coverage.loc[coverage.family.eq(row['family'])].groupby('policy').size()
        assert counts.eq(row['n_phenotypes']).all()
    layout = read(output / 'qa/layout_mm.csv')
    texts = pd.concat([read(p) for p in (output / 'qa').glob('text_*.csv')])
    assert texts.pt.min() >= 6 - 1e-6
    for i, a in layout.iterrows():
        assert a.x_mm >= 1 and a.y_mm >= 1 and a.x_mm+a.w_mm <= 173 and a.y_mm+a.h_mm <= 218
        for _, b in layout.iloc[i + 1:].iterrows():
            assert not (a.x_mm < b.x_mm+b.w_mm and b.x_mm < a.x_mm+a.w_mm
                        and a.y_mm < b.y_mm+b.h_mm and b.y_mm < a.y_mm+a.h_mm)
    stem = 'Full_atlas_overview' if manifest['mode'] == 'full' else 'Validation45_overview'
    exports = {}
    pages = [(stem, 174, 225)] + [(f'panels/{r.panel}', r.w_mm, r.h_mm) for r in layout.itertuples()]
    for name, width, height in pages:
        pdf = output / f'{name}.pdf'
        raw = pdf.read_bytes()
        reader = PdfReader(pdf)
        assert len(reader.pages) == 1
        page = reader.pages[0]
        actual_mm = [float(page.mediabox.width)*25.4/72, float(page.mediabox.height)*25.4/72]
        assert max(abs(a-b) for a,b in zip(actual_mm, [width,height])) < .01
        fonts = sorted({v.decode() for v in re.findall(rb'/BaseFont\s*/([^\s/>]+)', raw)})
        assert fonts and all('Arial' in f for f in fonts)
        assert raw.count(b'/FontFile2') >= 2
        assert not re.findall(rb'/Subtype\s*/Image', raw)
        png = Image.open(output / f'{name}.png')
        assert all(abs(a-b/25.4*600) <= 1 for a,b in zip(png.size,[width,height]))
        assert all(abs(x-600) < 1 for x in png.info['dpi'])
        if name == stem and manifest['mode'] != 'full':
            assert 'VALIDATION ONLY' in page.extract_text()
        exports[name] = dict(pdf_sha256=sha(pdf), page_mm=actual_mm, png_dimensions=png.size,
                             png_dpi=png.info['dpi'], fonts=fonts, vector_only=True)
    report = dict(complete=True, verified_at_utc=datetime.now(timezone.utc).isoformat(),
        mode=manifest['mode'], formal_full_atlas_verified=manifest['mode']=='full',
        plotted_rows=manifest['plot_rows'], independent_display_transform_max_error=maxima,
        all_assay_source_columns_exact=True, category_profile_source_columns_exact=True,
        display_selection_exact=True, whole_phenotype_coverage_exact=True,
        no_duplicate_assays_in_strength_panel=True, no_estimate_imputed_for_missing=True,
        minimum_text_pt=float(texts.pt.min()), exports=exports,
        checker_sha256=sha(Path(__file__)), visual_review_complete=False)
    path = output / 'qa/EXPORT_VERIFICATION.json'
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(complete=True, mode=report['mode'], plotted_rows=report['plotted_rows'],
                         exports_verified=len(exports), visual_review_complete=False),indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    check(parser.parse_args().output.resolve())
