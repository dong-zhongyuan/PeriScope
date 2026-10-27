"""Recompute protein accuracy, donor-level blood readouts and brain directions."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from ccwm import CCWM,CCWMConfig
from training_utils import predict,mmd

A=Path('/public/home/mengxl/dzy/pd_product_assets');I=A/'interim/rescreen_20261002';O=A/'results/optimization_20261003'

def main(seed,output_root=None):
    global O
    if output_root is not None:O=Path(output_root)
    torch.set_num_threads(2);dev='cuda'
    from model_registry import checkpoint, sha256, atomic_json
    model_path=checkpoint(O,seed)
    ck=torch.load(model_path,map_location='cpu',weights_only=False)
    m=CCWM(CCWMConfig(**ck['cfg']));m.load_state_dict(ck['state_dict']);m.eval().to(dev)
    cg=np.load(I/'citeseq.npz');bl=np.load(I/'blood.npz');br=np.load(I/'brain_locked_test.npz')
    states=bl['state_vocab'].astype(str).tolist();genes=bl['genes'].astype(str);registry=json.loads((O/'candidate_registry.json').read_text())
    out=O/'figure_inputs_20261004';out.mkdir(exist_ok=True)
    citeX=cg['X'];ci=np.where(cg['is_test'])[0];yp=[]
    with torch.no_grad():
        for start in range(0,len(ci),512):
            mu,_=m.enc_blood(torch.from_numpy(citeX[ci[start:start+512]]).to(dev));yp.append(m.protein_head(mu).cpu().numpy())
    ph=np.concatenate(yp);Y=cg['Y'][ci];rows=[];fine_rows=[]
    pairs=[]; donors=cg['donor']; donor_names=cg['donor_names']; cell_states=cg['state']; subtypes=cg['subtype']
    for r in registry:
        if r['target'] not in ['CD22','CD40','CD71']:continue
        for j,idx in enumerate(ci):
            pairs.append(dict(target=r['target'],donor=str(donor_names[donors[idx]]),blood_state=states[cell_states[idx]],subtype=str(subtypes[idx]),cell_row=int(idx),observed=float(Y[j,r['channels']].mean()),predicted=float(ph[j,r['channels']].mean())))
    pd.DataFrame(pairs).to_csv(out/'protein_cell_pairs.csv',index=False)
    (out/'protein_pairs_provenance.json').write_text(json.dumps(dict(seed=seed,checkpoint=str(model_path),sha256=sha256(model_path),input=str(I/'citeseq.npz'),input_sha256=sha256(I/'citeseq.npz'),unit='same held-out CITE-seq NPZ row; no pseudopairs',rows=len(pairs)),indent=2))
    for st in ['cDC','classical_mono','nonclassical_mono','pDC']:
        for d in np.unique(cg['donor'][ci]):
            ix=(cg['state'][ci]==states.index(st))&(cg['donor'][ci]==d)
            if ix.sum()<25:continue
            for r in registry:
                y=Y[ix][:,r['channels']].mean(1);p=ph[ix][:,r['channels']].mean(1)
                rho=float(spearmanr(y,p).statistic) if np.std(y)>1e-8 and np.std(p)>1e-8 else np.nan
                rows.append(dict(seed=seed,target=r['target'],blood_state=st,donor=cg['donor_names'][d],n=int(ix.sum()),rho=rho,mse=float(np.mean((p-y)**2)),is_control=r['is_control']))
            for subtype in np.unique(cg['subtype'][ci][ix]):
                fine=ix&(cg['subtype'][ci]==subtype)
                if fine.sum()<25:continue
                for r in registry:
                    y=Y[fine][:,r['channels']].mean(1);p=ph[fine][:,r['channels']].mean(1)
                    rho=float(spearmanr(y,p).statistic) if np.std(y)>1e-8 and np.std(p)>1e-8 else np.nan
                    fine_rows.append(dict(seed=seed,target=r['target'],blood_state=st,citeseq_subtype=str(subtype),donor=cg['donor_names'][d],n=int(fine.sum()),rho=rho,mse=float(np.mean((p-y)**2)),is_control=r['is_control']))
    pd.DataFrame(rows).to_csv(out/f'protein_accuracy_model_state_seed{seed}.csv',index=False)
    pd.DataFrame(fine_rows).to_csv(out/f'protein_accuracy_seed{seed}.csv',index=False)
    print('EVALUATED',seed,flush=True)


def export_programs():
    from model_registry import sha256,valid_curve,curve_provenance
    out=O/'figure_inputs_20261004';out.mkdir(exist_ok=True)
    e=pd.read_csv(O/'mr_parallel_gate_20261003/priority_all_programs.csv')
    programs=json.loads((O/'joint_analysis/programs.json').read_text())
    scores=pd.read_csv(O/'joint_analysis/program_seed_scores.csv')
    discovery=scores[scores.seed.le(46)].groupby('combo').signed_program_slope.mean()
    e['discovery_signed_slope']=e.combo.map(discovery)
    # One representative for each available target/brain state; ordered before confirmation rendering.
    chosen=e.sort_values(['protein','brain_state','discovery_signed_slope','combo'],ascending=[True,True,False,True]).drop_duplicates(['protein','brain_state']).head(5)
    chosen.to_csv(out/'five_examples.csv',index=False)
    examples=set(chosen.combo); curves=[];gene_rows=[];null_rows=[];obs=[];provenance={}
    for axis,es in e.groupby('axis'):
        f=O/'curves'/f'{axis}__slopes.npz';z=np.load(f);D=z['D'];gn=z['genes'].astype(str);prots=z['proteins'].astype(str).tolist();seeds=z['seeds'].astype(int);gidx={g:i for i,g in enumerate(gn)}
        provenance[str(f)]=sha256(f)
        for row in es.itertuples():
            pi=prots.index(row.protein);gi=np.array([gidx[g] for g in programs[row.combo]]);sgn=1 if row.response_direction=='up' else -1
            disc=D[seeds<=46,pi].mean(0);conf=D[seeds>=47,pi].mean(0)
            if row.combo in examples:
                for j,g in enumerate(gn):gene_rows.append(dict(combo=row.combo,gene=g,discovery=disc[j],confirmation=conf[j],positive_fraction=(D[:,pi,j]>0).mean(),member=g in programs[row.combo]))
                for seed in seeds:
                    f=O/'curves'/f'{axis}__seed{seed}.npz';assert valid_curve(f,curve_provenance(O,int(seed)));c=np.load(f);cv=c['C'][pi][:,gi].mean(1);doses=c['u'][pi]
                    for k,q in enumerate([.05,.25,.5,.75,.95]):curves.append(dict(combo=row.combo,seed=int(seed),quantile=q,dose=doses[k],response=cv[k]-cv[2],phase='Discovery' if seed<=46 else 'Confirmation'))
                signed=sgn*conf;denom=np.abs(signed).mean()+1e-12;v=signed[gi].mean()/denom
                rng=np.random.RandomState(7);ix=np.stack([rng.choice(len(gn),len(gi),replace=False) for _ in range(5000)])
                vals=signed[ix].mean(1)/denom;p=(1+(vals>=v).sum())/5001
                assert np.isclose(p,row.p_gene_decoy,rtol=1e-10,atol=1e-12),(row.combo,p,row.p_gene_decoy)
                null_rows.extend(dict(combo=row.combo,draw=i,value=float(x)) for i,x in enumerate(vals))
                obs.append(dict(combo=row.combo,observed=v,p=p,q=row.q_gene_decoy))
    for name,rows in [('program_dose_curves',curves),('gene_responses',gene_rows),('gene_null_draws',null_rows),('gene_null_observed',obs)]:pd.DataFrame(rows).to_csv(out/(name+'.csv'),index=False)
    (out/'program_export_provenance.json').write_text(json.dumps(dict(inputs=provenance,programs_sha256=sha256(O/'joint_analysis/programs.json'),models_sha256=sha256(O/'model_registry.json'),selection='highest discovery mean signed slope in each priority target/brain state; five first target/state groups in lexical order; confirmation not used for within-group selection',null_seed=7,null_draws=5000),indent=2))
    print('EXPORTED PROGRAMS',len(examples),flush=True)

if __name__=='__main__':
    main(42)
    export_programs()
