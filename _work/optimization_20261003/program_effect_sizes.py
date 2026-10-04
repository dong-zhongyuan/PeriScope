"""Compare confirmed dose responses with actual held-out brain PD-control contrasts in the same RNA units."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
O=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003');I=O.parent.parent/'interim/rescreen_20261002'
def main():
    programs=json.loads((O/'programs.json').read_text());brain=np.load(I/'brain_locked_test.npz');X=brain['X'];genes=brain['genes'].astype(str).tolist();gi={g:i for i,g in enumerate(genes)};states=brain['state_vocab'].astype(str).tolist()
    contrasts={};donor_support={}
    for state in ['astro','microglia_homeostatic','microglia_mhc2']:
        means=[];conditions=[];donors=[]
        for d in np.unique(brain['donor'][brain['state']==states.index(state)]):
            ix=np.where((brain['donor']==d)&(brain['state']==states.index(state)))[0]
            if len(ix)<25:continue
            means.append(X[ix].mean(0));conditions.append(brain['cond2'][ix[0]]);donors.append(brain['donor_names'][d])
        conditions=np.array(conditions)
        if len(means) and min(sum(conditions==0),sum(conditions==1))>=2:
            means=np.stack(means);contrasts[state]=means[conditions==1].mean(0)-means[conditions==0].mean(0)
            donor_support[state]=dict(n_PD=int(sum(conditions==1)),n_control=int(sum(conditions==0)),donors=donors)
    rows=[]
    axes=sorted({c.split('__',1)[1].rsplit('__',1)[0] for c in programs})
    for axis in axes:
        state=axis.split('__x__')[1]
        if state not in contrasts:continue
        C=np.stack([np.load(O/'curves'/f'{axis}__seed{s}.npz')['C'] for s in range(47,52)])
        ps=np.load(O/'curves'/f'{axis}__seed47.npz')['proteins'].astype(str).tolist()
        for name,members in programs.items():
            protein,rest=name.split('__',1)
            if rest.rsplit('__',1)[0]!=axis:continue
            ii=np.array([gi[g] for g in members]);delta=C[:,ps.index(protein),-1][:,ii]-C[:,ps.index(protein),0][:,ii];observed=contrasts[state][ii]
            denominator=float(observed@observed);norm=np.linalg.norm(delta,axis=1);observed_norm=np.linalg.norm(observed)
            projection=delta@observed/max(denominator,1e-16);cosine=(delta@observed)/np.maximum(norm*observed_norm,1e-16)
            positive=float(np.mean(projection>0));negative=float(np.mean(projection<0));agreement=max(positive,negative)
            action='decrease' if positive>=.8 else 'increase' if negative>=.8 else 'unresolved'
            rows.append(dict(combo=name,confirmation_disease_projection_mean=float(projection.mean()),confirmation_disease_projection_sd=float(projection.std(ddof=1)),
                confirmation_response_disease_cosine=float(cosine.mean()),confirmation_response_to_disease_norm_ratio=float(np.mean(norm/max(observed_norm,1e-16))),
                confirmation_projection_direction_fraction=agreement,predicted_restorative_protein_change=action,
                brain_observed_program_mean_PD_minus_control=float(observed.mean()),confirmation_high_minus_low_program_mean=float(delta.mean()),
                effect_size_n_PD=donor_support[state]['n_PD'],effect_size_n_control=donor_support[state]['n_control']))
    cols=['combo','confirmation_disease_projection_mean','confirmation_disease_projection_sd','confirmation_response_disease_cosine','confirmation_response_to_disease_norm_ratio','confirmation_projection_direction_fraction','predicted_restorative_protein_change','brain_observed_program_mean_PD_minus_control','confirmation_high_minus_low_program_mean','effect_size_n_PD','effect_size_n_control']
    pd.DataFrame(rows,columns=cols).to_csv(O/'program_effect_sizes.csv',index=False)
    (O/'program_effect_size_definition.json').write_text(json.dumps(dict(dose='measured 95th minus 5th percentile, paired source cells and latent noise',
        response='confirmation models 47-51 only',observed='donor-mean held-out brain PD minus control in the same log1p(CP10k) feature coordinates',
        projection='dot(model high-low, observed PD-control) / squared norm(observed PD-control), within frozen program genes',
        action='positive projection in >=4/5 models: decrease protein; negative projection in >=4/5: increase; otherwise unresolved',
        interpretation='model-based directional hypothesis; positive projection moves along the observed disease contrast',donor_support=donor_support),indent=2))
    print('EFFECT SIZES',len(rows),flush=True)
if __name__=='__main__':main()
