"""Read a completed atlas without refitting MR or changing its corrections."""
from collections import Counter

import numpy as np
import pandas as pd


BEST_FIELDS = [
    'phenotype', 'trait', 'family', 'beta', 'se', 'ci_lower', 'ci_upper',
    'p', 'log_p', 'q_bh_global', 'q_by_global', 'q_bh_family', 'F',
]
READ_FIELDS = list(dict.fromkeys([
    'exposure_id', 'gene', 'assay', 'status', 'family_n_global',
    'family_n', 'cis_family_n', *BEST_FIELDS,
]))


def qc_flags(meta, mode):
    """Only the frozen catalog determines whether a trait is QC."""
    flags = []
    for row in meta.to_dict('records'):
        value = str(row.get('is_qc', '')).lower()
        if value in ('true', 'false'):
            flag = value == 'true'
        elif mode == 'validation45' or (
            str(row['phenotype']).startswith('OxfordQSM2022:')
            and row['family'] in ('QSM', 'T2star')
        ):
            flag = False
        else:
            raise ValueError(f"Missing frozen QC designation: {row['phenotype']}")
        if flag != (row['family'] == 'QC'):
            raise ValueError('QC flag and frozen family disagree')
        flags.append(flag)
    return np.asarray(flags, dtype=bool)


class Profile:
    def __init__(self, size):
        self.planned = 0
        self.estimated = np.zeros(size, dtype=np.int64)
        self.best_log_p = np.full(size, np.inf)
        self.best_id = np.full(size, '', dtype=object)
        self.best_rows = [None] * size
        self.z_min = np.full(size, np.inf)
        self.z_max = np.full(size, -np.inf)

    def add(self, block):
        self.planned += 1
        valid = block.status.eq('estimated').to_numpy()
        self.estimated += valid
        log_p = block.log_p.to_numpy(dtype=float)
        phenotype = str(block.phenotype.iloc[0])
        replace = valid & ((log_p < self.best_log_p) | (
            (log_p == self.best_log_p) & (phenotype < self.best_id)))
        for i in np.flatnonzero(replace):
            self.best_rows[i] = block.iloc[i][BEST_FIELDS].to_dict()
        self.best_log_p[replace] = log_p[replace]
        self.best_id[replace] = phenotype
        z = np.full(len(block), np.nan)
        z[valid] = block.loc[valid, 'beta'] / block.loc[valid, 'se']
        self.z_min[valid] = np.minimum(self.z_min[valid], z[valid])
        self.z_max[valid] = np.maximum(self.z_max[valid], z[valid])

    def frame(self, coverage):
        result = coverage[['exposure_id', 'gene', 'assay']].reset_index(drop=True).copy()
        result['n_planned'] = self.planned
        result['n_estimated'] = self.estimated
        result['n_not_estimated'] = self.planned - self.estimated
        result['z_min'] = np.where(self.estimated > 0, self.z_min, np.nan)
        result['z_max'] = np.where(self.estimated > 0, self.z_max, np.nan)
        best = pd.DataFrame([
            row if row is not None else {key: np.nan for key in BEST_FIELDS}
            for row in self.best_rows
        ]).add_prefix('best_')
        return pd.concat([result, best], axis=1)


def summarize_table(source, coverage, meta, mode, expected):
    """Keep every assay; rank biological profiles descriptively by minimum log P.

    QC rows stay in source corrections and in coverage/category tables. They
    cannot determine the biological ranking. A best-row P/q is copied from an
    individual pair, never interpreted as a combined assay/category test.
    """
    size = len(coverage)
    if coverage.exposure_id.duplicated().any() or meta.phenotype.duplicated().any():
        raise ValueError('Duplicate assay or phenotype in the declared plan')
    flags = qc_flags(meta, mode)
    biological = Profile(size)
    categories = {family: Profile(size) for family in meta.family.drop_duplicates()}
    family_sizes = meta.groupby('family').size().to_dict()
    observed_counts = Counter()
    coverage_rows = []
    reader = pd.read_csv(source, usecols=READ_FIELDS, chunksize=size,
                         float_precision='round_trip')
    total = 0
    try:
        for i, trait in enumerate(meta.to_dict('records')):
            try:
                block = next(reader)
            except StopIteration as error:
                raise ValueError('Incomplete full-family table') from error
            block = block.reset_index(drop=True)
            if len(block) != size or block.exposure_id.tolist() != coverage.exposure_id.tolist():
                raise ValueError('Assay order/coverage differs from declared Cartesian product')
            for key in ('gene', 'assay'):
                if block[key].fillna('').tolist() != coverage[key].fillna('').tolist():
                    raise ValueError(f'Assay annotation mismatch: {key}')
            for key in ('phenotype', 'family'):
                if not block[key].eq(trait[key]).all():
                    raise ValueError(f'Phenotype block mismatch: {key}')
            denominators = {
                'family_n_global': size * len(meta),
                'family_n': size * family_sizes[trait['family']],
                'cis_family_n': expected['cis_family_n'],
            }
            for key, value in denominators.items():
                if not block[key].eq(value).all():
                    raise ValueError(f'Correction denominator mismatch: {key}')
            valid = block.status.eq('estimated')
            numeric = block.loc[valid, ['beta', 'se', 'p', 'log_p', 'q_bh_global', 'q_by_global', 'q_bh_family']]
            if not np.isfinite(numeric.to_numpy()).all() or (numeric.se <= 0).any():
                raise ValueError('Invalid numerical value in an estimated pair')
            if (numeric.log_p > 0).any() or ((numeric[['p', 'q_bh_global', 'q_by_global', 'q_bh_family']] < 0) | (numeric[['p', 'q_bh_global', 'q_by_global', 'q_bh_family']] > 1)).any().any():
                raise ValueError('P/q outside valid range')
            if block.status.isna().any():
                raise ValueError('A planned pair has no explicit estimation status')
            counts = block.status.value_counts().to_dict()
            observed_counts.update(counts)
            for status, count in sorted(counts.items()):
                coverage_rows.append(dict(phenotype=trait['phenotype'], family=trait['family'],
                    is_qc=bool(flags[i]), status=status, n_assays=count))
            categories[trait['family']].add(block)
            if not flags[i]:
                biological.add(block)
            total += len(block)
        try:
            next(reader)
        except StopIteration:
            pass
        else:
            raise ValueError('Unexpected rows beyond the declared Cartesian product')
    finally:
        reader.close()
    if total != expected['family_n'] or dict(observed_counts) != expected['status_counts']:
        raise ValueError('Complete-table count/status audit differs from the completed run')
    ranking = biological.frame(coverage).sort_values(
        ['best_log_p', 'exposure_id'], na_position='last', kind='stable').reset_index(drop=True)
    ranks = pd.Series(pd.NA, index=ranking.index, dtype='Int64')
    available = ranking.n_estimated.gt(0)
    ranks.loc[available] = np.arange(1, available.sum() + 1)
    ranking.insert(0, 'descriptive_rank', ranks)
    profiles = []
    for family, state in categories.items():
        frame = state.frame(coverage)
        frame.insert(3, 'imaging_family', family)
        frame.insert(4, 'is_qc', family == 'QC')
        profiles.append(frame)
    return ranking, pd.concat(profiles, ignore_index=True), pd.DataFrame(coverage_rows), {
        'planned_rows': total, 'status_counts': dict(observed_counts),
        'biological_phenotypes': int((~flags).sum()), 'qc_phenotypes': int(flags.sum()),
        'all_assays_retained': len(ranking), 'ranked_assays': int(available.sum()),
        'qc_in_ranking': False,
    }
