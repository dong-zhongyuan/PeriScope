"""Adapter for the untouched PeriScope 2026-10-03 MR input snapshot.

Historical cohort selection/F>=10 are reproduced explicitly, not proposed as a
new discovery rule. No saved MR estimates are read while constructing inputs.
"""
from pathlib import Path
import numpy as np
import pandas as pd

from .scan import scan_pairs


def prepare_legacy(snapshot):
    src = Path(snapshot) / 'results'
    exposure = pd.read_csv(src / 'combined_exposure_instruments.csv', float_precision='round_trip')
    outcome = pd.read_csv(src / 'brain_sentinel_statistics.csv', float_precision='round_trip')
    if exposure.assay.duplicated().any() or outcome.duplicated(['exposure_assay', 'outcome_accession']).any():
        raise ValueError('Duplicate input keys')
    meta_cols = ['outcome_accession', 'outcome_build', 'outcome_n', 'outcome_trait']
    meta = outcome[meta_cols].drop_duplicates()
    if meta.outcome_accession.duplicated().any():
        raise ValueError('Inconsistent outcome metadata')
    grid = exposure.merge(meta, how='cross')
    x = grid.merge(outcome.drop(columns=meta_cols[1:]), left_on=['assay', 'outcome_accession'],
                   right_on=['exposure_assay', 'outcome_accession'], how='left',
                   suffixes=('_x', '_y'), validate='one_to_one', indicator=True)
    acquired = x['_merge'].eq('both')
    if not x.loc[acquired, 'rsid'].eq(x.loc[acquired, 'rsid_outcome']).all():
        raise ValueError('Legacy extraction has nonmatching variant identifiers')
    pairs = pd.DataFrame({
        'gene': x.gene, 'assay': x.assay, 'variant': x.rsid,
        'phenotype': x.outcome_accession,
        'family': np.where(x.outcome_trait.str.contains('T2star'), 'T2star', 'volume'),
        'trait': x.outcome_trait,
        'ea_x': x.effect_allele_x, 'oa_x': x.other_allele_x,
        'ea_y': x.effect_allele_y, 'oa_y': x.other_allele_y,
        'beta_x': x.beta_exposure, 'se_x': x.se_exposure,
        'beta_y': x.beta_outcome, 'se_y': x.se_outcome,
        'eaf_x': x.eaf_exposure, 'eaf_y': x.eaf_outcome,
        'chrom_x': x.chrom, 'chrom_y': x.chromosome,
        'pos_x': np.where(x.outcome_build.eq('GRCh37'), x.pos37, x.pos),
        'pos_y': x.position, 'build_x': x.outcome_build, 'build_y': x.outcome_build,
        'source_x': x.source, 'source_y': x.source_url_y,
        'n_x': np.nan, 'n_y': x.outcome_n,
        'units_x': 'reported protein analysis unit', 'units_y': 'standardized imaging phenotype SD',
        'F_source': x.F_exposure, 'PAV_linked': x.PAV_linked,
        'acquisition_status': np.where(x['_merge'].eq('both'), 'acquired', 'not_in_extraction'),
        'coordinate_source': 'legacy pos(GRCh38)/pos37 map; no new liftover',
    })
    # Audit each assay separately because distinct assays for a gene are not IVs.
    temporary = pairs.copy()
    temporary['gene'] = temporary['assay']
    _, audit = scan_pairs(temporary)
    audit['gene'] = pairs.gene.to_numpy()
    # Legacy code did not test coordinates; verify them separately before proceeding.
    if audit.status.isin(['coordinate_mismatch', 'missing_coordinates', 'build_mismatch_or_unknown']).any():
        raise ValueError('Legacy extraction fails coordinate verification')
    audit['legacy_eligible'] = audit.usable & (audit.F_source >= 10)
    ukb = audit.source_x.eq('Sun2023_ST9_discovery')
    ukb_genes = set(audit.loc[ukb & audit.legacy_eligible, 'gene'])
    ages = exposure[exposure.source.eq('Gudjonsson2022_AGES_serum') & ~exposure.gene.isin(ukb_genes)]
    ages = ages.sort_values(['F_exposure', 'assay'], ascending=[False, True]).drop_duplicates('gene')
    audit['legacy_selected'] = audit.legacy_eligible & (ukb | audit.assay.isin(ages.assay))
    audit['selection_reason'] = np.select(
        [audit.legacy_selected, ~audit.legacy_eligible, audit.gene.isin(ukb_genes)],
        ['historical_selected', 'not_legacy_eligible', 'AGES_UKB_covered'], default='AGES_alternative_assay')
    selected = pairs.loc[audit.legacy_selected].copy()
    if selected.duplicated(['gene', 'phenotype']).any():
        raise ValueError('Legacy selected family is not single-instrument')
    return selected, audit
