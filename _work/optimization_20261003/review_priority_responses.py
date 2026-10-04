"""Inspect every B/A program's per-seed disease projection without reranking."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

A = Path('/public/home/mengxl/dzy/pd_product_assets')
O = A / 'results/optimization_20261003'
J = O / 'joint_analysis'
e = pd.read_csv(J / 'all_program_evidence.csv')
e = e[e.tier.isin(['A', 'B'])]
programs = json.loads((J / 'programs.json').read_text())
brain = np.load(A / 'interim/rescreen_20261002/brain_locked_test.npz')
brain_X = brain['X']
genes = brain['genes'].astype(str).tolist()
states = brain['state_vocab'].astype(str).tolist()
gi = {g: i for i, g in enumerate(genes)}
rows = []
for rec in e.to_dict('records'):
    state = rec['brain_state']
    means, conditions = [], []
    for donor in np.unique(brain['donor'][brain['state'] == states.index(state)]):
        ix = np.where((brain['donor'] == donor) & (brain['state'] == states.index(state)))[0]
        if len(ix) < 25:
            continue
        means.append(brain_X[ix].mean(0))
        conditions.append(brain['cond2'][ix[0]])
    means, conditions = np.stack(means), np.array(conditions)
    contrast = means[conditions == 1].mean(0) - means[conditions == 0].mean(0)
    ii = np.array([gi[g] for g in programs[rec['combo']]])
    observed = contrast[ii]
    for seed in range(47, 52):
        z = np.load(O / 'curves' / f"{rec['axis']}__seed{seed}.npz")
        assert list(z['genes'].astype(str)) == genes
        pi = list(z['proteins'].astype(str)).index(rec['protein'])
        delta = z['C'][pi, -1, ii] - z['C'][pi, 0, ii]
        projection = float(delta @ observed / max(float(observed @ observed), 1e-16))
        cosine = float(delta @ observed / max(float(np.linalg.norm(delta) * np.linalg.norm(observed)), 1e-16))
        rows.append(dict(combo=rec['combo'], protein=rec['protein'], blood_state=rec['blood_state'],
                         brain_state=state, seed=seed, disease_projection=projection, cosine=cosine,
                         response_to_disease_norm_ratio=float(np.linalg.norm(delta) / np.linalg.norm(observed)),
                         high_minus_low_program_mean=float(delta.mean())))
d = pd.DataFrame(rows)
d.to_csv(J / 'priority_response_seed_review.csv', index=False)
review = []
for rec in e.to_dict('records'):
    v = d[d.combo.eq(rec['combo'])]
    assert np.isclose(v.disease_projection.mean(), rec['confirmation_disease_projection_mean'], atol=1e-8)
    majority = 1 if (v.disease_projection > 0).mean() >= .8 else -1 if (v.disease_projection < 0).mean() >= .8 else 0
    mean_sign = int(np.sign(v.disease_projection.mean()))
    review.append(dict(combo=rec['combo'], protein=rec['protein'], tier=rec['tier'],
        above_isotype_controls=rec['above_isotype_controls'],
        positive_projection_seeds=int((v.disease_projection > 0).sum()),
        negative_projection_seeds=int((v.disease_projection < 0).sum()),
        projection_mean=float(v.disease_projection.mean()),
        projection_median=float(v.disease_projection.median()),
        projection_min=float(v.disease_projection.min()), projection_max=float(v.disease_projection.max()),
        majority_and_mean_projection_agree=majority != 0 and majority == mean_sign,
        response_to_disease_norm_ratio=rec['confirmation_response_to_disease_norm_ratio'],
        independent_p=rec['independent_p'], independent_q=rec['independent_q'],
        protein_panel_q=rec['q_protein_decoy']))
r = pd.DataFrame(review)
r.to_csv(J / 'priority_quality_review.csv', index=False)
(J / 'priority_quality_review_definition.json').write_text(json.dumps(dict(
    scope='Diagnostic review of all nominated programs; nomination rules and tier assignments unchanged.',
    response='Measured 95th minus 5th percentile dose, confirmation seeds 47-51.',
    disease='Identical held-out donor contrasts and RNA coordinates to program_effect_sizes.py.',
    checks='Reproduce stored mean projection, compare its sign with the seed-majority sign, inspect isotype controls and reported uncertainty.'), indent=2))
print(r.to_string(index=False))
