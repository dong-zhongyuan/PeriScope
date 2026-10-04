"""Update only repaired PeriScope seeds in the existing external sensitivity."""
from pathlib import Path
import sys,json,os
import numpy as np
import pandas as pd
import torch
from model_registry import ROOT, sha256, atomic_json, checkpoint

P=ROOT.parent/'pig_external_validation_20261003'
W=Path('/public/home/mengxl/dzy/pd_product/_work/pig_external_validation_20261003')
sys.path.insert(0,str(W))
import predict_external as pred
from common import load_model,condition
from evaluate_external import response,metrics

def main():
    torch.set_num_threads(2)
    device=os.environ.get('PIG_DEVICE','cuda:1');pred.DEVICE=device
    registry=json.loads((ROOT/'model_registry.json').read_text())
    fixed={n:sha256(P/n) for n in ['adapter_locked.json','frozen_primary_predictions.npz','primary_summary.csv','response_metrics.csv']}
    z=dict(np.load(P/'frozen_ten_seed_predictions.npz'))
    old=z['linear_profiles'].copy()
    inp=dict(np.load(P/'adapted_pig_blood_reference_inputs.npz'))
    ref=dict(np.load(P/'human_reference_cells.npz'));states=ref['state_vocab'].tolist()
    assert np.array_equal(z['sample_ids'],inp['sample_ids']) and np.array_equal(z['genes'],inp['genes'])
    bslist=['cDC','classical_mono','nonclassical_mono']
    def one(model,axis,c,j):
        bs,rs=axis.split('__')
        eps=torch.as_tensor(np.random.default_rng(42000+bslist.index(bs)).normal(size=(128,model.cfg.z_dim)),device=device,dtype=torch.float32)
        cond=condition(128,c,states.index(rs),len(states),device)
        out,_,_=pred.predict_periscope(model,inp[bs][j],cond,eps)
        return pred.to_linear_mean(out,model.brain_center.cpu().numpy())
    # Confirm the unchanged inference route against an existing seed before replacing anything.
    model=load_model(43,device)
    check=one(model,str(z['axes'][0]),0,0)
    si=list(z['seeds']).index(43)
    assert np.allclose(check,old[si,0,0,0],rtol=1e-5,atol=2e-5),float(np.max(np.abs(check-old[si,0,0,0])))
    del model
    for seed in [47,48,49]:
        checkpoint(ROOT,seed);model=load_model(seed,device);si=list(z['seeds']).index(seed)
        for ai,axis in enumerate(z['axes']):
            for c in [0,1]:
                for j in range(len(z['sample_ids'])):z['linear_profiles'][si,c,ai,j]=one(model,str(axis),c,j)
        del model
        print('PIG_SEED_UPDATED',seed,flush=True)
    unchanged=[i for i,s in enumerate(z['seeds']) if s not in [47,48,49]]
    assert np.array_equal(old[unchanged],z['linear_profiles'][unchanged])
    assert np.isfinite(z['linear_profiles']).all()
    tmp=P/'frozen_ten_seed_predictions.pending.npz';np.savez_compressed(tmp,**z);tmp.replace(P/'frozen_ten_seed_predictions.npz')
    lock=json.loads((P/'predictions_locked.json').read_text())
    lock.update(ten_seed_predictions_sha256=sha256(P/'frozen_ten_seed_predictions.npz'),
        checkpoint_sha256={s:v['sha256'] for s,v in registry['models'].items()},
        ten_seed_readout='Original PeriScope readout with authoritative repaired ten-seed weights; separate from calibrated primary comparison',
        ten_seed_refresh_script_sha256=sha256(__file__),model_registry_sha256=sha256(ROOT/'model_registry.json'))
    atomic_json(P/'predictions_locked.json',lock)
    names=z['sample_ids'].astype(str);mask=ref['mask']
    meta=pd.read_csv(P/'inputs/matched_animals.csv').set_index('sample_id').loc[names]
    group=meta.group.eq('LPS').to_numpy();ctrl=~group
    brain=pd.read_csv(P/'inputs/brain_one2one_sum_provided_transcript_TPM.csv',index_col=0).loc[z['genes'],names].to_numpy().T
    obs=response(brain,mask,ctrl);rows=[]
    for si,seed in enumerate(z['seeds']):
        for c in [0,1]:
            pp=np.stack([response(z['linear_profiles'][si,c,ai],mask,ctrl) for ai in range(len(z['axes']))]).mean(0)
            rows.append(dict(seed=int(seed),condition=c,**metrics(pp,obs,group)))
    pd.DataFrame(rows).to_csv(P/'ten_seed_sensitivity.csv',index=False)
    done=json.loads((P/'evaluation_completed.json').read_text());done['predictions_lock_sha256']=sha256(P/'predictions_locked.json')
    done['ten_seed_model_registry_sha256']=sha256(ROOT/'model_registry.json');atomic_json(P/'evaluation_completed.json',done)
    assert all(sha256(P/n)==h for n,h in fixed.items())
    atomic_json(P/'ten_seed_refresh_verification.json',{'status':'passed','replaced_seeds':[47,48,49],
        'unchanged_seeds':[int(z['seeds'][i]) for i in unchanged],'primary_comparison_unchanged':True,
        'model_registry_sha256':sha256(ROOT/'model_registry.json')})

if __name__=='__main__':main()
