"""Observed held-out CITE-seq identity context; never a fitted nomination gate.

Calculate within donor and annotated fine subtype, rather than correlating
across mixed blood lineages. All biological antigens receive the same audit.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata

A = Path('/public/home/mengxl/dzy/pd_product_assets')
I = A / 'interim/rescreen_20261002'
O = A / 'results/optimization_20261003'


def correlation(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a -= a.mean()
    b -= b.mean()
    denominator = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / denominator) if denominator > 1e-12 else np.nan


def main():
    registry = json.loads((O / 'candidate_registry.json').read_text())
    z = np.load(I / 'citeseq.npz')
    X, Y, counts = z['X'], z['Y'], z['adt_counts']
    gene_mu = z['gene_mu']
    genes = z['genes'].astype(str).tolist()
    gi = {g: i for i, g in enumerate(genes)}
    controls = np.median(np.column_stack([
        Y[:, r['channels']].mean(1) for r in registry if r['is_control']
    ]), axis=1)
    state, donor, subtype = z['state'], z['donor'], z['subtype']
    state_names = z['state_vocab'].astype(str)
    test = z['is_test'].astype(bool)
    rows = []
    for st in ['cDC', 'classical_mono', 'nonclassical_mono', 'pDC']:
        for d in np.unique(donor[(state == list(state_names).index(st)) & test]):
            ix_state = (state == list(state_names).index(st)) & (donor == d) & test
            for sub in np.unique(subtype[ix_state]):
                ix = np.where(ix_state & (subtype == sub))[0]
                if len(ix) < 25:
                    continue
                control_rank = rankdata(controls[ix])
                control_design = np.column_stack([np.ones(len(ix)), control_rank])
                for target in registry:
                    if target['is_control']:
                        continue
                    mapped = target['genes']
                    available = [g for g in mapped if g in gi]
                    y = Y[np.ix_(ix, target['channels'])].mean(1)
                    y_rank = rankdata(y)
                    y_resid = y_rank - control_design @ np.linalg.lstsq(control_design, y_rank, rcond=None)[0]
                    row = dict(target=target['target'], genes=';'.join(mapped),
                        mapping=target['mapping'], blood_state=st,
                        donor=str(z['donor_names'][d]), citeseq_subtype=str(sub),
                        n_cells=len(ix), split='held_out_test',
                        n_RNA_genes_available=len(available), RNA_genes_available=';'.join(available),
                        median_ADT_CLR=float(np.median(y)), ADT_CLR_IQR=float(np.quantile(y, .75)-np.quantile(y, .25)),
                        fraction_ADT_counts_positive=float((counts[np.ix_(ix, target['channels'])].sum(1) > 0).mean()),
                        ADT_isotype_rho=correlation(y_rank.copy(), control_rank.copy()))
                    if available:
                        idx = [gi[g] for g in available]
                        x = X[np.ix_(ix, idx)] + gene_mu[idx]
                        aggregate = x.mean(1)
                        x_rank = rankdata(aggregate)
                        x_resid = x_rank - control_design @ np.linalg.lstsq(control_design, x_rank, rcond=None)[0]
                        row.update(RNA_fraction_any_detected=float((x > 1e-6).any(1).mean()),
                            RNA_fraction_all_detected=float((x > 1e-6).all(1).mean()),
                            RNA_mean_log1p_CP10k=float(aggregate.mean()),
                            ADT_mapped_RNA_rho=correlation(y_rank.copy(), x_rank.copy()),
                            ADT_mapped_RNA_partial_rho_isotype=correlation(y_resid, x_resid))
                    rows.append(row)
        print('ASSAY STATE COMPLETE', st, len(rows), flush=True)
    df = pd.DataFrame(rows)
    assert set(df.donor) == {'P6', 'P7', 'P8'}
    df.to_csv(O / 'assay_identity_by_heldout_subtype.csv', index=False)
    values = ['RNA_fraction_any_detected', 'RNA_fraction_all_detected',
        'RNA_mean_log1p_CP10k', 'ADT_mapped_RNA_rho',
        'ADT_mapped_RNA_partial_rho_isotype', 'ADT_isotype_rho',
        'median_ADT_CLR', 'ADT_CLR_IQR', 'fraction_ADT_counts_positive']
    # Equal subtype weight inside each donor, then equal donor weight.
    donor_means = df.groupby(['target', 'blood_state', 'donor'], as_index=False)[values].mean()
    donor_means.to_csv(O / 'assay_identity_by_heldout_donor.csv', index=False)
    means = donor_means.groupby(['target', 'blood_state'])[values].mean().add_prefix('heldout_observed_')
    means['heldout_identity_n_donors'] = donor_means.groupby(['target', 'blood_state']).size()
    means.reset_index().to_csv(O / 'assay_identity_context.csv', index=False)
    (O / 'assay_identity_context_definition.json').write_text(json.dumps(dict(
        status='complete', n_antigens=int(df.target.nunique()), n_rows=len(df),
        source=str(I / 'citeseq.npz'),
        registry_sha256=hashlib.sha256((O / 'candidate_registry.json').read_bytes()).hexdigest(),
        unit='Observed held-out donors P6-P8, within fine subtype; >=25 cells; equal subtype weight per donor, equal donor weight.',
        complex='For defined complexes, aggregate available mapped-gene log1p RNA; also retain any/all subunit detection.',
        partial='Rank-residual correlation after adjusting both mapped-RNA and antigen ADT ranks for the median isotype rank.',
        use='Post-screen interpretation context for every antigen, without modifying frozen nomination criteria.',
        interpretation='Low RNA detection can reflect capture sensitivity; observed RNA-ADT concordance is distinct from model protein-prediction accuracy.'
    ), indent=2))
    print('ASSAY IDENTITY CONTEXT', len(df), flush=True)


if __name__ == '__main__':
    main()
