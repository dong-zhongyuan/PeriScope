"""Independent dense audit and adversarial guards for atlas presentation."""
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.atlas_presentation import summarize_table, qc_flags
from s19_prepare_atlas_presentation import completed_source

FOLDER = ROOT / 'papers/01_mr_imaging_scan/results/full_atlas_20261116/validation45'
OUTPUT = FOLDER / 'presentation'


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return pd.read_csv(path, float_precision='round_trip')


def dense_profile(raw, cov):
    """Use whole-table pandas grouping, independent of the streaming accumulator."""
    rows = raw.loc[raw.status.eq('estimated')].copy()
    rows['z'] = rows.beta / rows.se
    expected = cov[['exposure_id', 'gene', 'assay']].copy()
    count = raw.groupby('exposure_id').size()
    observed = rows.groupby('exposure_id').size()
    expected['n_planned'] = expected.exposure_id.map(count)
    expected['n_estimated'] = expected.exposure_id.map(observed).fillna(0).astype(int)
    expected['n_not_estimated'] = expected.n_planned - expected.n_estimated
    for key, values in [('z_min', rows.groupby('exposure_id').z.min()),
                        ('z_max', rows.groupby('exposure_id').z.max())]:
        expected[key] = expected.exposure_id.map(values)
    best = rows.sort_values(['exposure_id', 'log_p', 'phenotype']).drop_duplicates('exposure_id')
    fields = ['phenotype', 'trait', 'family', 'beta', 'se', 'ci_lower', 'ci_upper',
              'p', 'log_p', 'q_bh_global', 'q_by_global', 'q_bh_family', 'F']
    best = best[['exposure_id', *fields]].rename(columns={x: 'best_' + x for x in fields})
    return expected.merge(best, on='exposure_id', how='left', validate='one_to_one')


def assert_frame(actual, expected, keys):
    a = actual.sort_values(keys).reset_index(drop=True)[expected.columns]
    b = expected.sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b, check_dtype=False, check_exact=True)


def main():
    manifest = json.loads((OUTPUT / 'MANIFEST.json').read_text())
    if manifest['mode'] != 'validation45' or manifest['formal_full_atlas_result']:
        raise AssertionError('Validation incorrectly described as formal full atlas')
    for relative, digest in {**manifest['inputs'], **manifest['code_sha256']}.items():
        assert sha(ROOT / relative) == digest
    for name, digest in manifest['outputs'].items():
        assert sha(OUTPUT / name) == digest
    cov = read(ROOT / 'data/ukb_ppp/extracted/assay_coverage.csv')
    instruments = read(ROOT / 'data/ukb_ppp/extracted/cis_primary_selected.csv')
    reports = {}
    for policy in ('drop', 'infer'):
        raw = read(FOLDER / f'full_family_palindrome_{policy}.csv.gz')
        expected = dense_profile(raw, cov).sort_values(
            ['best_log_p', 'exposure_id'], na_position='last').reset_index(drop=True)
        ranks = pd.Series(pd.NA, index=expected.index, dtype='Int64')
        valid = expected.n_estimated.gt(0)
        ranks.loc[valid] = np.arange(1, valid.sum() + 1)
        expected.insert(0, 'descriptive_rank', ranks)
        actual = read(OUTPUT / f'all_assay_ranking_{policy}.csv')
        assert actual.exposure_id.tolist() == expected.exposure_id.tolist()
        assert_frame(actual, expected, ['exposure_id'])
        columns = ['exposure_id', 'variant', 'chrom', 'pos37', 'pos38', 'pav_consequence',
                   'pav_coding_consequence', 'pav_annotation_scope', 'F_x']
        annotation = cov[['exposure_id']].merge(instruments[columns], on='exposure_id', how='left')
        assert_frame(actual, annotation, ['exposure_id'])
        profiles = []
        for family, rows in raw.groupby('family', sort=False):
            group = dense_profile(rows, cov)
            group.insert(3, 'imaging_family', family)
            group.insert(4, 'is_qc', False)
            profiles.append(group)
        expected_profiles = pd.concat(profiles, ignore_index=True)
        actual_profiles = read(OUTPUT / f'all_assay_category_profiles_{policy}.csv.gz')
        assert_frame(actual_profiles, expected_profiles, ['imaging_family', 'exposure_id'])
        expected_coverage = raw.groupby(['phenotype', 'family', 'status']).size().reset_index(name='n_assays')
        expected_coverage.insert(2, 'is_qc', False)
        assert_frame(read(OUTPUT / f'coverage_by_phenotype_{policy}.csv'), expected_coverage,
                     ['phenotype', 'family', 'status'])
        if policy == 'drop':
            display = expected.loc[valid].head(24)
            assert read(OUTPUT / 'overview_display_assays.csv').exposure_id.tolist() == display.exposure_id.tolist()
        selected = expected_profiles.loc[expected_profiles.exposure_id.isin(display.exposure_id)].copy()
        selected['display_row'] = selected.exposure_id.map(dict(zip(display.exposure_id, range(1, 25))))
        assert_frame(read(OUTPUT / f'overview_category_profiles_{policy}.csv'), selected,
                     ['imaging_family', 'exposure_id'])
        reports[policy] = dict(planned_rows=len(raw), all_assays=len(actual),
            all_category_profiles=len(actual_profiles), coverage_rows=len(expected_coverage),
            independent_values_exact=True, all_annotations_exact=True)
    guards = {}
    with tempfile.TemporaryDirectory(prefix='atlas_presentation_check_') as tmp:
        tmp = Path(tmp)
        fixture_cov = pd.DataFrame(dict(exposure_id=['a', 'b', 'c', 'd'],
            gene=['CD22', 'NEW', 'NONE', 'TIE'], assay=['A', 'B', 'C', 'D']))
        fixture_meta = pd.DataFrame(dict(phenotype=['bio_z', 'qc', 'bio_a'],
            family=['bio', 'QC', 'bio'], is_qc=[False, True, False]))
        rows = []
        for trait in fixture_meta.to_dict('records'):
            for i, record in enumerate(fixture_cov.to_dict('records')):
                lp = {'bio_z': [-1., -900., 0., -5.], 'qc': [-10000., -2., 0., -3.],
                      'bio_a': [-2., -1000., 0., -5.]}[trait['phenotype']][i]
                rows.append(dict(**record, phenotype=trait['phenotype'], trait=trait['phenotype'],
                    family=trait['family'], status='missing_outcome' if i == 2 else 'estimated',
                    beta=1., se=1., ci_lower=-.96, ci_upper=2.96, p=np.exp(lp), log_p=lp,
                    q_bh_global=.4, q_by_global=.6, q_bh_family=.3, F=100.,
                    family_n_global=12, family_n=4 if trait['is_qc'] else 8, cis_family_n=12))
        fixture = pd.DataFrame(rows)
        path = tmp / 'pairs.csv'
        fixture.to_csv(path, index=False)
        expected_counts = dict(family_n=12, cis_family_n=12,
            status_counts=dict(Counter(fixture.status)))
        ranked, profiled, _, _ = summarize_table(path, fixture_cov, fixture_meta, 'full', expected_counts)
        assert ranked.exposure_id.tolist() == ['b', 'd', 'a', 'c']
        assert ranked.set_index('exposure_id').loc['b', 'best_phenotype'] == 'bio_a'
        assert ranked.set_index('exposure_id').loc['d', 'best_phenotype'] == 'bio_a'
        assert ranked.set_index('exposure_id').loc['a', 'best_log_p'] == -2.
        assert pd.isna(ranked.set_index('exposure_id').loc['c', 'descriptive_rank'])
        assert len(profiled) == 8 and profiled.loc[profiled.is_qc, 'n_planned'].eq(1).all()
        guards.update(qc_excluded_from_biological_ranking=True, qc_profile_preserved=True,
            underflow_sorted_by_log_p=True, tie_sorted_by_phenotype_id=True,
            unestimated_assay_retained=True, historical_name_receives_no_priority=True)
        for label, bad in [
            ('assay_order', fixture.iloc[[1, 0, *range(2, len(fixture))]]),
            ('missing_rows', fixture.iloc[:-1]),
            ('extra_rows', pd.concat([fixture, fixture.iloc[:1]])),
            ('denominator', fixture.assign(family_n_global=11)),
        ]:
            bad.to_csv(path, index=False)
            try:
                summarize_table(path, fixture_cov, fixture_meta, 'full', expected_counts)
            except ValueError:
                guards[label + '_rejected'] = True
            else:
                raise AssertionError(f'Malformed input accepted: {label}')
        bad_meta = fixture_meta.copy()
        bad_meta['is_qc'] = np.nan
        try:
            qc_flags(bad_meta, 'full')
        except ValueError:
            guards['missing_qc_designation_rejected'] = True
        else:
            raise AssertionError('Unknown QC flags silently accepted')
        for label, state in [('incomplete_run', dict(complete=False)),
            ('wrong_dimensions', dict(complete=True, status='complete', mode='full',
                assays=2941, phenotypes=45, completed_traits=45, planned_pairs=132345))]:
            (tmp / 'RUN.json').write_text(json.dumps(state))
            try:
                completed_source(tmp, 'full')
            except RuntimeError:
                guards[label + '_rejected'] = True
            else:
                raise AssertionError(f'Completion guard failed: {label}')
    full = FOLDER.parent / 'full'
    before = (full / 'presentation').exists()
    try:
        completed_source(full, 'full')
    except RuntimeError:
        guards['actual_incomplete_full_atlas_rejected'] = True
    else:
        raise AssertionError('Reassess: the full run is now available')
    assert (full / 'presentation').exists() == before
    result = dict(complete=True, created_at_utc=datetime.now(timezone.utc).isoformat(),
        scope='45-trait implementation validation plus adversarial guards; not completed full-atlas results',
        source_manifest_sha256=sha(OUTPUT / 'MANIFEST.json'),
        verifier_sha256=sha(Path(__file__)), independent_dense_checks=reports, guards=guards,
        full_atlas_presentation_created=False, figures_rendered=False)
    (OUTPUT / 'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
