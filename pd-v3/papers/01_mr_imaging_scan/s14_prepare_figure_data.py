"""Prepare exact, inspectable plotting units for the R small-multiple figure."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from common import ROOT, PAPER
from pdv3.provenance import sha256

REGIONS = ['Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Hippocampus', 'Amygdala', 'Accumbens', 'SN']
SHORT = dict(zip(REGIONS, ['THA', 'CAU', 'PUT', 'PAL', 'HIP', 'AMY', 'ACC', 'SN']))


def region_meta(trait):
    low = trait.lower()
    region = next((r for r in REGIONS if r.lower() in low), None)
    if 'midbrain' in low:
        region = 'Midbrain'
    if 'wmh' in low:
        region = 'WMH adjusted' if 'wmhvout' in low else 'WMH'
    hemisphere = 'L' if 'left' in low else 'R' if 'right' in low else ''
    return region, hemisphere


def prepare(folder):
    folder = Path(folder).resolve()
    out = folder / 'figures/plot_data'
    out.mkdir(parents=True, exist_ok=True)
    source = folder / 'full_family_palindrome_drop.csv.gz'
    cols = ['exposure_id', 'gene', 'assay', 'phenotype', 'trait', 'family', 'status',
            'beta', 'se', 'se_delta', 'p', 'q_bh_global', 'F', 'chrom', 'pos37',
            'eaf_x', 'eaf_y', 'pav_coding_consequence', 'ci_lower', 'ci_upper']
    all_rows = pd.read_csv(source, usecols=cols, low_memory=False)
    assert len(all_rows) == 2941 * 45
    all_rows['region'], all_rows['hemisphere'] = zip(*all_rows.trait.map(region_meta))
    all_rows['modality'] = all_rows.family.map({'volume': 'Volume', 'QSM': 'QSM', 'T2star': 'T2*'})
    all_rows['region_short'] = all_rows.region.map(SHORT).fillna(all_rows.region)
    all_rows['idp_label'] = (all_rows.hemisphere + ' ' + all_rows.region_short).str.strip()
    estimated = all_rows[all_rows.status.eq('estimated')].copy()
    assert len(estimated) == 68175
    estimated['z'] = estimated.beta / estimated.se
    estimated['se_increase_percent'] = 100 * (estimated.se_delta / estimated.se - 1)
    meta = pd.read_csv(folder / 'phenotypes.csv')
    def write(frame, name):
        frame.to_csv(out / f'{name}.csv', index=False)
    write(meta, 'phenotypes')
    # A: every estimable assay-by-IDP unit, ordered by its cis instrument coordinate.
    positions = estimated[['exposure_id', 'chrom', 'pos37']].drop_duplicates()
    assert not positions.exposure_id.duplicated().any()
    positions['chr_n'] = positions.chrom.replace({'X': 23}).astype(float).astype(int)
    positions = positions.sort_values(['chr_n', 'pos37', 'exposure_id']).reset_index(drop=True)
    positions['assay_order'] = np.arange(1, len(positions) + 1)
    landscape = estimated.merge(positions[['exposure_id', 'assay_order', 'chr_n']], on='exposure_id', validate='many_to_one')
    write(landscape[['exposure_id', 'gene', 'assay', 'phenotype', 'modality', 'z', 'assay_order', 'chr_n']], 'A_landscape')
    write(positions.groupby('chr_n').assay_order.agg(['min', 'max', 'median']).reset_index(), 'A_chromosomes')
    # B: display top 24 assays by their best raw P, plus all three historical assays.
    ranking = estimated.groupby('exposure_id', sort=False).p.min().sort_values(kind='stable')
    selected = list(ranking.index[:24])
    for exposure in estimated.loc[estimated.gene.isin(['CD22', 'CD40', 'TFRC']), 'exposure_id'].unique():
        if exposure not in selected:
            selected.append(exposure)
    matrix = estimated[estimated.exposure_id.isin(selected)].copy()
    matrix['assay_row'] = matrix.exposure_id.map({e: i + 1 for i, e in enumerate(selected)})
    # Order by modality, then anatomy and hemisphere; full source IDs remain in CSV.
    order = pd.DataFrame({'phenotype': meta.phenotype})
    order = order.merge(estimated[['phenotype', 'modality', 'region', 'hemisphere', 'idp_label']].drop_duplicates(), on='phenotype')
    order['m_order'] = order.modality.map({'Volume': 0, 'QSM': 1, 'T2*': 2})
    order['r_order'] = order.region.map({r: i for i, r in enumerate(REGIONS + ['Midbrain', 'WMH', 'WMH adjusted'])})
    order = order.sort_values(['m_order', 'r_order', 'hemisphere']).reset_index(drop=True)
    order['idp_order'] = np.arange(1, 46)
    matrix = matrix.merge(order[['phenotype', 'idp_order']], on='phenotype', validate='many_to_one')
    assert len(matrix) == len(selected) * 45
    write(matrix, 'B_matrix')
    write(order, 'B_columns')
    historical = estimated[estimated.gene.isin(['CD22', 'CD40', 'TFRC'])].copy()
    assert len(historical) == 135
    write(historical.merge(order[['phenotype', 'idp_order']], on='phenotype'), 'S1_historical')
    focus = historical[historical.region.eq('SN')].copy()
    focus['outcome_label'] = focus.modality + ' ' + focus.hemisphere
    assert len(focus) == 12
    write(focus, 'C_black_substantia_nigra')
    regional = estimated[estimated.gene.eq('TFRC') & estimated.region.isin(REGIONS) & estimated.modality.isin(['QSM', 'T2*'])]
    assert len(regional) == 32
    write(regional, 'D_tfrc_regions')
    # E: hemisphere-specific paired QSM/T2* estimates (not independent replication).
    base = estimated[estimated.region.isin(REGIONS) & estimated.modality.isin(['QSM', 'T2*'])]
    cross = base.pivot(index=['exposure_id', 'gene', 'region', 'region_short', 'hemisphere'], columns='modality', values='z').reset_index()
    assert len(cross) == 1515 * 8 * 2 and not cross[['QSM', 'T2*']].isna().any().any()
    cross = cross.rename(columns={'QSM': 'z_qsm', 'T2*': 'z_t2star'})
    write(cross, 'E_cross_modality')
    # F: within-assay left-minus-right z is descriptive, not an independent-sample test.
    bilateral = estimated[estimated.hemisphere.isin(['L', 'R'])].pivot(
        index=['exposure_id', 'gene', 'region', 'region_short', 'modality'], columns='hemisphere', values='z').reset_index()
    assert len(bilateral) == 1515 * 20
    bilateral['left_minus_right_z'] = bilateral.L - bilateral.R
    write(bilateral, 'F_bilateral')
    write(estimated[['exposure_id', 'phenotype', 'modality', 'se_increase_percent']], 'G_delta_se')
    # H: all 269 palindromic tools once, using one source phenotype; no outcome selection.
    infer_source = folder / 'full_family_palindrome_infer.csv.gz'
    inf = pd.read_csv(infer_source, usecols=['exposure_id', 'gene', 'phenotype', 'status', 'eaf_x', 'eaf_y'], low_memory=False)
    pal_ids = set(all_rows.loc[all_rows.status.eq('palindromic_dropped'), 'exposure_id'])
    pal = inf[inf.exposure_id.isin(pal_ids) & inf.phenotype.eq('BIG40:release2:discovery:0011')].copy()
    assert len(pal) == 269 and pal.status.eq('estimated').sum() == 227
    pal['maf_protein'] = np.minimum(pal.eaf_x, 1 - pal.eaf_x)
    pal['maf_imaging'] = np.minimum(pal.eaf_y, 1 - pal.eaf_y)
    pal = pal.sort_values(['maf_protein', 'exposure_id']).reset_index(drop=True)
    pal['tool_order'] = np.arange(1, 270)
    pal['retained'] = pal.status.eq('estimated')
    write(pal, 'H_palindromes')
    record = {'inputs': {str(p.relative_to(ROOT)): sha256(p) for p in [source, infer_source, folder / 'phenotypes.csv']},
              'outputs': {p.name: sha256(p) for p in sorted(out.glob('*.csv'))},
              'matrix_assays': selected, 'matrix_selection': 'top 24 by best raw P across declared 45 IDPs, plus three historical assays; display only',
              'complete_assay_idp_units_in_A': len(estimated), 'all_planned_slots': len(all_rows),
              'descriptive_only': ['cross-modality pairing', 'left-minus-right Wald z'],
              'preparation_script_sha256': sha256(__file__)}
    (out / 'MANIFEST.json').write_text(json.dumps(record, indent=2) + '\n')
    return out


if __name__ == '__main__':
    print(prepare(PAPER / 'results/expanded_20261116/combined_45'))
