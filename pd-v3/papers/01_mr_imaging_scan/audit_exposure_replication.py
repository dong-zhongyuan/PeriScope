"""Audit all selected exposure instruments across Sun ST9 cohorts.

No imaging outcomes are read and no instruments are reselected. This is a
diagnostic within the replacement full MR analysis, not an additional MR scan.
"""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.tabix_http import save_json


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    folder = ROOT / 'data/ukb_ppp/extracted'
    source = folder / 'cis_primary_selected.csv'
    raw_source = folder / 'ST9_raw.csv'
    x = pd.read_csv(source)
    raw = pd.read_csv(raw_source).set_index('original_row')
    assert len(x) == 1954 and x.exposure_id.is_unique
    columns = {
        'beta_x': 'BETA (discovery, wrt. A1)',
        'se_x': 'SE (discovery)',
        'eaf_x': 'A1FREQ (discovery)',
        'beta_x_replication': 'BETA (replication)',
        'se_x_replication': 'SE (replication)',
        'eaf_x_replication': 'A1FREQ (replication)',
    }
    errors = {}
    for canonical, original in columns.items():
        a = x[canonical].to_numpy()
        b = raw.loc[x.original_row, original].to_numpy(dtype=float)
        assert np.isfinite(a).all() and np.isfinite(b).all()
        errors[canonical] = float(np.max(np.abs(a-b)))
        assert errors[canonical] == 0
    assert (x.se_x > 0).all() and (x.se_x_replication > 0).all()
    assert x.beta_x.ne(0).all() and x.beta_x_replication.ne(0).all()
    x['F_discovery'] = (x.beta_x/x.se_x)**2
    x['F_replication'] = (x.beta_x_replication/x.se_x_replication)**2
    x['same_direction'] = x.beta_x * x.beta_x_replication > 0
    x['beta_replication_minus_discovery'] = x.beta_x_replication-x.beta_x
    x['difference_se_disjoint_cohorts'] = np.hypot(x.se_x, x.se_x_replication)
    x['difference_z'] = x.beta_replication_minus_discovery/x.difference_se_disjoint_cohorts
    x['difference_p_descriptive'] = 2*norm.sf(np.abs(x.difference_z))
    x['replication_to_discovery_beta_ratio'] = x.beta_x_replication/x.beta_x
    x['mr_effect_rescaling_factor_if_using_replication'] = x.beta_x/x.beta_x_replication
    x['replication_ci_includes_zero'] = np.abs(x.beta_x_replication) <= norm.ppf(.975)*x.se_x_replication
    x['eaf_replication_minus_discovery'] = x.eaf_x_replication-x.eaf_x
    out = ROOT/'papers/01_mr_imaging_scan/results/exposure_replication_20261116'
    out.mkdir(exist_ok=True)
    target = out/'all_selected_instruments.csv.gz'
    x.to_csv(target, index=False, compression={'method': 'gzip', 'mtime': 0})
    qs = [0, .01, .05, .25, .5, .75, .95, .99, 1]
    quantities = ['F_discovery', 'F_replication', 'replication_to_discovery_beta_ratio',
                  'mr_effect_rescaling_factor_if_using_replication',
                  'eaf_replication_minus_discovery']
    summary = dict(
        complete=True, selected_assays=len(x), outcome_data_used=False,
        instrument_reselection=False, direction_concordant=int(x.same_direction.sum()),
        replication_F_below_10_diagnostic_only=int(x.F_replication.lt(10).sum()),
        replication_F_below_1_diagnostic_only=int(x.F_replication.lt(1).sum()),
        replication_nominal_ci_includes_zero=int(x.replication_ci_includes_zero.sum()),
        pearson_beta=float(pearsonr(x.beta_x, x.beta_x_replication).statistic),
        spearman_beta=float(spearmanr(x.beta_x, x.beta_x_replication).statistic),
        quantiles={c: {str(q): float(v) for q,v in x[c].quantile(qs).items()} for c in quantities},
        raw_ST9_numeric_max_errors=errors,
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [source,raw_source]},
        output_sha256=sha(target), source_code_sha256=sha(Path(__file__)),
        source_publication='https://pmc.ncbi.nlm.nih.gov/articles/PMC10567551/',
        cohort_evidence='Discovery: randomly selected European ancestry baseline cohort; replication: remaining samples, multiple ancestries and enriched selection/batches. Both are within UKB.',
        sample_sizes='Headline discovery n=34557 and replication n=17806 are not assigned as per-SNP N.',
        interpretation='Exposure consistency diagnostic; not an independent imaging replication or an identified winner-curse correction. Differences combine selection, ancestry, ascertainment, batch and sampling effects.',
        single_instrument_identity='For fixed aligned outcome and nonzero beta_x, first-order absolute Wald z equals absolute outcome beta/se. Replacing beta_x changes effect and SE but is not another independent P-value test.',
    )
    save_json(out/'AUDIT.json', summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['quantiles','inputs','cohort_evidence','interpretation']},indent=2))


if __name__ == '__main__':
    main()
