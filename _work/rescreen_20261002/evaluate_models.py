"""Recompute protein accuracy, donor-level blood readouts and brain directions."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from ccwm import CCWM,CCWMConfig
from training_utils import predict,mmd

A=Path('/public/home/mengxl/dzy/pd_product_assets');I=A/'interim/rescreen_20261002';O=A/'results/rescreen_20261002'

def main(seed,output_root=None):
    global O
    if output_root is not None:O=Path(output_root)
    torch.set_num_threads(2);dev='cuda'
    ck=torch.load(O/f'models/seed_{seed}/model.pt',map_location='cpu',weights_only=False)
    m=CCWM(CCWMConfig(**ck['cfg']));m.load_state_dict(ck['state_dict']);m.eval().to(dev)
    cg=np.load(I/'citeseq.npz');bl=np.load(I/'blood.npz');br=np.load(I/'brain_locked_test.npz')
    states=bl['state_vocab'].astype(str).tolist();genes=bl['genes'].astype(str);registry=json.loads((O/'candidate_registry.json').read_text())
    out=O/'evaluation';out.mkdir(exist_ok=True)
    citeX=cg['X'];ci=np.where(cg['is_test'])[0];yp=[]
    with torch.no_grad():
        for start in range(0,len(ci),512):
            mu,_=m.enc_blood(torch.from_numpy(citeX[ci[start:start+512]]).to(dev));yp.append(m.protein_head(mu).cpu().numpy())
    ph=np.concatenate(yp);Y=cg['Y'][ci];rows=[];fine_rows=[]
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
    X=bl['X'];blood_ids=sorted({s[1] for s in json.loads((O/f'models/seed_{seed}/training_inputs.json').read_text())['strata']})
    brain_ids=sorted({s[2] for s in json.loads((O/f'models/seed_{seed}/training_inputs.json').read_text())['strata']})
    donor_rows=[];direction=[];readouts=[];readout_labels=[];brainX=br['X']
    # Original gradient solve, with inferred proteins in every inference path.
    grad=torch.autograd.grad
    def first_order(*a,**kw):kw['create_graph']=False;return grad(*a,**kw)
    torch.autograd.grad=first_order
    for bs in blood_ids:
        bydonor={}
        for d in np.unique(bl['donor'][bl['state']==bs]):
            ix=np.where((bl['state']==bs)&(bl['donor']==d))[0]
            if len(ix)<25:continue
            with torch.no_grad():
                buf=[]
                for start in range(0,len(ix),512):
                    mu,_=m.enc_blood(torch.from_numpy(X[ix[start:start+512]]).to(dev));buf.append(m.protein_head(mu).cpu().numpy())
                mean=np.concatenate(buf).mean(0)
            for r in registry:donor_rows.append(dict(seed=seed,target=r['target'],blood_state=states[bs],donor=bl['donor_names'][d],condition=int(bl['cond2'][ix[0]]),cohort=int(bl['cohort'][ix[0]]),n=len(ix),predicted_ADT=float(mean[r['channels']].mean()),training_donor=bool(not(bl['is_val'][ix[0]] or bl['is_test'][ix[0]]))))
            if not(bl['is_val'][ix[0]] or bl['is_test'][ix[0]]):continue
            sel=np.random.default_rng(6600+int(d)).choice(ix,min(96,len(ix)),replace=False)
            bydonor[d]=(sel,int(bl['cond2'][sel[0]]))
        for brain in brain_ids:
            observed=[];conditions=[]
            for d in np.unique(br['donor'][br['state']==brain]):
                ix=np.where((br['donor']==d)&(br['state']==brain))[0]
                if len(ix)<25:continue
                observed.append(brainX[ix].mean(0));conditions.append(int(br['cond2'][ix[0]]))
            conditions=np.array(conditions)
            if not observed or not all(np.sum(conditions==c)>=2 for c in (0,1)):continue
            observed=np.stack(observed);diff=observed[conditions==1].mean(0)-observed[conditions==0].mean(0)
            for mode in ['fixed_control','fixed_PD','condition_aware']:
                means=[];labels=[]
                for d,(sel,condition) in bydonor.items():
                    cond=0 if mode=='fixed_control' else 1 if mode=='fixed_PD' else condition
                    c=torch.zeros(len(sel),m.cfg.cond_dim,device=dev);c[:,cond]=1;c[:,2+brain]=1;c[:,-1]=1
                    p=predict(m,torch.from_numpy(X[sel]).to(dev),c).cpu().numpy().mean(0)
                    means.append(p);labels.append(condition);readouts.append(p);readout_labels.append(f'{states[bs]}|{states[brain]}|{mode}|{bl["donor_names"][d]}')
                if not means:continue
                means=np.stack(means);labels=np.array(labels)
                if not all(np.sum(labels==c)>=2 for c in (0,1)):continue
                delta=means[labels==1].mean(0)-means[labels==0].mean(0)
                direction.append(dict(seed=seed,blood_state=states[bs],brain_state=states[brain],mode=mode,
                    rho=float(spearmanr(delta,diff).statistic),n_blood_donors=len(labels),n_brain_donors=len(conditions),
                    n_brain_PD=int(sum(conditions==1)),n_brain_control=int(sum(conditions==0)),
                    predicted_direction_norm=float(np.linalg.norm(delta)),observed_direction_norm=float(np.linalg.norm(diff))))
    pd.DataFrame(donor_rows).to_csv(out/f'blood_donor_predictions_seed{seed}.csv',index=False)
    pd.DataFrame(direction).to_csv(out/f'brain_direction_seed{seed}.csv',index=False)
    np.savez_compressed(out/f'brain_donor_readouts_seed{seed}.npz',values=np.stack(readouts),labels=readout_labels,genes=genes)
    torch.autograd.grad=grad
    print('EVALUATED',seed,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True);ap.add_argument('--output-root');args=ap.parse_args();main(args.seed,args.output_root)
