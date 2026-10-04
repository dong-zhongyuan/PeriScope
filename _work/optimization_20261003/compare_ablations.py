"""Compare the new complete model and its two matched training controls, using identical queries."""
from pathlib import Path
import json,time,hashlib
import numpy as np
import pandas as pd
O=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
while not ((O/'ablations/completed.json').exists() and (O/'programs.json').exists() and (O/'models/seed_42/screened.json').exists()):
    time.sleep(15)
programs=json.loads((O/'programs.json').read_text());full=pd.read_csv(O/'evaluation/brain_direction_seed42.csv');rows=[];effects=[]
for kind in ['random_pairing','coupling_off']:
    root=O/'ablations'/kind
    control=pd.read_csv(root/'evaluation/brain_direction_seed42.csv')
    joined=full.merge(control,on=['blood_state','brain_state','mode'],suffixes=('_full','_control'))
    for r in joined.to_dict('records'):
        r['control']=kind;r['same_direction']=bool(r['rho_full']*r['rho_control']>0)
        r['absolute_direction_retention']=abs(r['rho_control'])/max(abs(r['rho_full']),1e-12)
        rows.append(r)
    for axis in sorted({c.split('__',1)[1].rsplit('__',1)[0] for c in programs}):
        z=np.load(O/'curves'/f'{axis}__seed42.npz');zc=np.load(root/'curves'/f'{axis}__seed42.npz')
        assert np.array_equal(z['u'],zc['u']) and np.array_equal(z['genes'],zc['genes']) and np.array_equal(z['proteins'],zc['proteins'])
        gi={g:i for i,g in enumerate(z['genes'])};ps=z['proteins'].astype(str).tolist()
        C=z['C'];C_control=zc['C']  # Decompress each unchanged response array once per axis.
        for name,members in programs.items():
            protein,rest=name.split('__',1)
            if rest.rsplit('__',1)[0]!=axis:continue
            pi=ps.index(protein);ii=[gi[g] for g in members]
            a=(C[pi,-1]-C[pi,0])[ii];b=(C_control[pi,-1]-C_control[pi,0])[ii]
            na=float(np.linalg.norm(a));nb=float(np.linalg.norm(b));dot=float(a@b)
            effects.append(dict(combo=name,control=kind,full_response_norm=na,control_response_norm=nb,
                response_norm_retention=nb/max(na,1e-16),response_cosine=dot/max(na*nb,1e-16),
                full_program_mean_high_minus_low=float(a.mean()),control_program_mean_high_minus_low=float(b.mean())))
pd.DataFrame(rows).to_csv(O/'ablation_brain_direction_comparison.csv',index=False)
pd.DataFrame(effects,columns=['combo','control','full_response_norm','control_response_norm','response_norm_retention','response_cosine','full_program_mean_high_minus_low','control_program_mean_high_minus_low']).to_csv(O/'ablation_program_response_comparison.csv',index=False)
protein=[]
for kind,root in [('full',O),('random_pairing',O/'ablations/random_pairing'),('coupling_off',O/'ablations/coupling_off')]:
    ac=pd.read_csv(root/'evaluation/protein_accuracy_seed42.csv');ac['model']=kind;protein.append(ac)
pd.concat(protein,ignore_index=True).to_csv(O/'ablation_protein_accuracy.csv',index=False)
(O/'ablation_comparison_definition.json').write_text(json.dumps(dict(status='complete',seed=42,programs_sha256=hashlib.sha256((O/'programs.json').read_bytes()).hexdigest(),
    model_selection='identical internal donor validation and training budget for full and controls',
    pairing='blood training donor/state sampling independent of brain stratum; no held-out blood cells used',
    coupling='direct coupling potential removed, blood initialization and inferred protein inputs retained',
    queries='identical measured subtype doses and paired blood/latent draws',n_brain_direction_comparisons=len(rows),n_program_comparisons=len(effects)),indent=2))
print('ABLATION COMPARISON COMPLETE',len(rows),len(effects),flush=True)
