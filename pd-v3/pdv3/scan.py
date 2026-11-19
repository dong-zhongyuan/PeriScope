"""Complete planned protein/outcome pairs, including failed and missing estimates.

Input is one row per planned instrument/outcome pair. Source adapters must build
that complete plan before joining acquired outcome statistics onto it.
"""
import numpy as np
import pandas as pd

from .mr import harmonize, wald, ivw
from .multiplicity import adjust_table

REQUIRED = ['gene', 'assay', 'variant', 'phenotype', 'family', 'ea_x', 'oa_x',
            'ea_y', 'oa_y', 'beta_x', 'se_x', 'beta_y', 'se_y', 'build_x',
            'build_y', 'chrom_x', 'chrom_y', 'pos_x', 'pos_y']


def scan_pairs(pairs, *, independent=False):
    """Return one result per planned gene/phenotype and a full instrument audit.

    Multiple instruments are pooled only when independence is explicitly attested.
    No strength/significance cutoff is applied here. Duplicate/metadata errors fail
    rather than allowing a many-to-many join or mixing exposure cohorts silently.
    """
    if not isinstance(independent, bool):
        raise ValueError('independent must be an explicit boolean')
    missing = set(REQUIRED) - set(pairs.columns)
    if missing:
        raise ValueError(f'Missing columns: {sorted(missing)}')
    if pairs.empty:
        raise ValueError('Empty analysis plan')
    if pairs[['gene', 'assay', 'variant', 'phenotype', 'family']].isna().any().any():
        raise ValueError('Plan identifiers must not be missing')
    if pairs.duplicated(['gene', 'phenotype', 'variant']).any():
        raise ValueError('Duplicate variant within gene/outcome; separate assays/cohorts into runs')
    observed = pairs[['gene', 'phenotype']].drop_duplicates()
    if len(observed) != pairs.gene.nunique() * pairs.phenotype.nunique():
        raise ValueError('Incomplete gene/outcome plan; retain missing outcome rows explicitly')
    audit_rows = []
    for row in pairs.to_dict('records'):
        a = dict(row, usable=False, status='missing_outcome', beta_y_aligned=np.nan)
        if any(pd.isna(row[k]) for k in ['beta_y', 'se_y', 'ea_y', 'oa_y']):
            audit_rows.append(a)
            continue
        if row['build_x'] not in ('GRCh37', 'GRCh38') or row['build_y'] != row['build_x']:
            a['status'] = 'build_mismatch_or_unknown'
        elif any(pd.isna(row[k]) for k in ['chrom_x', 'chrom_y', 'pos_x', 'pos_y']):
            a['status'] = 'missing_coordinates'
        elif any(not np.isfinite(float(row[k])) or float(row[k]) <= 0 or not float(row[k]).is_integer()
                 for k in ['pos_x', 'pos_y']):
            a['status'] = 'invalid_coordinates'
        elif str(row['chrom_x']).removesuffix('.0') != str(row['chrom_y']).removesuffix('.0') or row['pos_x'] != row['pos_y']:
            a['status'] = 'coordinate_mismatch'
        else:
            h = harmonize(row['ea_x'], row['oa_x'], row['ea_y'], row['oa_y'], row['beta_y'],
                          eaf_x=row.get('eaf_x'), eaf_y=row.get('eaf_y'),
                          exposure_build=row['build_x'], outcome_build=row['build_y'])
            a.update(status=h['status'], usable=h['usable'], beta_y_aligned=h['beta_y'])
            if a['usable']:
                try:
                    a.update(wald(row['beta_x'], row['se_x'], h['beta_y'], row['se_y']))
                except ValueError as exc:
                    a.update(usable=False, status='invalid_statistics', error=str(exc))
        audit_rows.append(a)
    audit = pd.DataFrame(audit_rows)
    results = []
    for (gene, phenotype), group in audit.groupby(['gene', 'phenotype'], sort=True):
        if group.family.nunique() != 1:
            raise ValueError(f'Inconsistent family: {gene}/{phenotype}')
        if group.assay.nunique() != 1:
            raise ValueError(f'Mixed assays: {gene}/{phenotype}; use separate runs')
        for metadata in ['source_x', 'units_x', 'units_y', 'source_y']:
            if metadata in group and group[metadata].dropna().nunique() > 1:
                raise ValueError(f'Mixed {metadata}: {gene}/{phenotype}')
        good = group[group.usable]
        rec = dict(gene=gene, phenotype=phenotype, family=group.family.iloc[0],
                   n_planned_instruments=len(group), n_instruments=len(good),
                   status='no_usable_instrument', method=None, beta=np.nan, se=np.nan, p=np.nan)
        for metadata in ['source_x', 'source_y', 'units_x', 'units_y', 'n_x', 'n_y']:
            if metadata in group:
                vals = group[metadata].dropna().unique()
                rec[metadata] = vals[0] if len(vals) == 1 else np.nan
        if len(good) == 1:
            r = good.iloc[0]
            rec.update(wald(r.beta_x, r.se_x, r.beta_y_aligned, r.se_y))
            rec.update(status='estimated', method='Wald', assay=r.assay, variant=r.variant)
        elif len(good) > 1 and independent:
            rec.update(ivw(good.beta_x, good.se_x, good.beta_y_aligned, good.se_y, independent=True))
            rec.update(status='estimated', method='IVW_independent')
        elif len(good) > 1:
            rec['status'] = 'independence_not_attested'
        results.append(rec)
    return adjust_table(pd.DataFrame(results), family_cols=('family',)), audit
