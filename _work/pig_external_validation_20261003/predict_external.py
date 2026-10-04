"""Frozen external predictions from measured pig blood; pig brain is not opened."""
import os,time
import numpy as np,pandas as pd
from scipy.special import logsumexp
from common import *
import sys
sys.path.insert(0,"/public/home/mengxl/dzy/pd_product/_work/benchmark_shared_20261003")
import current_benchmark as shared_b
from shared_pipeline import Predictors, METHODS
shared_b.DEVICE=os.environ.get("PIG_DEVICE","cuda:1")

DEVICE=os.environ.get('PIG_DEVICE','cuda:1')

def to_linear_mean(pred,center):
    raw=np.maximum(pred+center,0)
    if not np.isfinite(raw).all():raise ValueError('Nonfinite external prediction')
    # Exact stable counterpart of mean(expm1(raw)), followed by a common
    # library rescaling. Prevent overflow without clipping model magnitudes.
    rr=raw.astype(np.float64)
    with np.errstate(divide='ignore',invalid='ignore'):
        lab=rr+np.log(-np.expm1(-rr))
    lm=logsumexp(lab,axis=0)-np.log(len(raw))
    return (np.exp(lm-logsumexp(lm))*10000).astype(np.float32)

def predict_periscope(model,x,cond,eps):
    xx=torch.as_tensor(x,device=DEVICE,dtype=torch.float32)
    with torch.no_grad():
        mu,lv=model.enc_blood(xx);u=model.protein_head(mu)
        z=mu+torch.exp(.5*lv)*eps if model.cfg.sample_latent_a else mu
    with torch.enable_grad():_,q=model.solve_equilibrium(z,z.clone(),u,cond)
    with torch.no_grad():pred=model.dec_brain(q.detach()).cpu().numpy()
    return pred,u.cpu().numpy(),mu.cpu().numpy()

def main():
    torch.set_num_threads(2)
    lock=json.loads((O/'adapter_locked.json').read_text());assert lock['protocol_sha256']==sha(O/'protocol.txt')
    assert lock['checkpoint_sha256']==sha(A/'results/optimization_20261003/models/seed_42/model.pt')
    ref=dict(np.load(O/'human_reference_cells.npz'));genes=ref['genes'];mask=ref['mask'];states=list(ref['state_vocab'])
    meta=pd.read_csv(O/'inputs/matched_animals.csv');names=meta.sample_id.tolist();control=meta.group.eq('Control').to_numpy()
    pig=pd.read_csv(O/'inputs/blood_one2one_sum_provided_transcript_TPM.csv',index_col=0).loc[genes,names].to_numpy().T
    assert np.isfinite(pig[:,mask]).all()
    pig[:,~mask]=0 # Outside mask never supplies a measured value to the adapter.
    pr=lock['selected'];refpig=pig[control].mean(0)
    fcs=np.stack([logfactor(row,refpig,mask,pr['alpha'],pr['cap_log']) for row in pig])
    qc=[]
    for j,name in enumerate(names):
        qc.append(dict(sample_id=name,group=meta.group.iloc[j],fraction_observed_genes_at_cap=float(np.mean(np.isclose(np.abs(fcs[j,mask]),pr['cap_log']))),median_absolute_log_factor=float(np.median(abs(fcs[j,mask])))))
    pd.DataFrame(qc).to_csv(O/'blood_input_change_diagnostics.csv',index=False)
    adapted={}
    for bs in BS:
        linear=np.expm1(np.maximum(ref[bs]+ref['blood_center'],0))
        adapted[bs]=np.stack([np.log1p(adapt(linear,fc,pr['gain']))-ref['blood_center'] for fc in fcs]).astype(np.float32)
    np.savez_compressed(O/'adapted_pig_blood_reference_inputs.npz',**adapted,genes=genes,sample_ids=np.array(names))
    methods=METHODS
    selected=os.environ.get('PIG_METHODS',','.join(methods)).split(',')
    preds=np.load(O/'frozen_primary_predictions.npz')['linear_profiles'].copy() if len(selected)<len(methods) else np.empty((2,len(methods),len(AXES),len(names),len(genes)),np.float32)
    diagnostics=pd.read_csv(O/'prediction_diagnostics.csv').query('method not in @selected').to_dict('records') if len(selected)<len(methods) else [];model=load_model(42,DEVICE);center=model.brain_center.cpu().numpy()
    baseline_sources={}
    for ai,axis in enumerate(AXES):
        bs,rs=axis.split('__');k=dict(np.load(B/'cache'/f'{axis}.npz'))
        assert np.array_equal(k['genes'],genes)
        baseline_sources[str(B/'cache'/f'{axis}.npz')]=sha(B/'cache'/f'{axis}.npz')
        eps=torch.as_tensor(np.random.default_rng(42000+BS.index(bs)).normal(size=(128,model.cfg.z_dim)),device=DEVICE,dtype=torch.float32)
        for c in [0,1]:
            cond=condition(128,c,states.index(rs),len(states),DEVICE)
            factory=Predictors(model,k,axis,42)
            for j,name in enumerate(names):
                x=adapted[bs][j]
                with torch.no_grad():
                    mu,lv=model.enc_blood(torch.as_tensor(x,device=DEVICE));u=model.protein_head(mu)
                    z=mu+torch.exp(.5*lv)*eps if model.cfg.sample_latent_a else mu
                for mi,method in enumerate(methods):
                    if method not in selected:continue
                    pred=factory.predict(method,z.cpu().numpy(),u.cpu().numpy(),c)
                    parfile=B/'shared_readout/parameters/seed_42'/f'{method}__{axis}.npz'
                    cal=dict(np.load(parfile));baseline_sources[str(parfile)]=sha(parfile)
                    pred=pred.astype(float)+float(cal['alpha'])*cal['offset']+float(cal['residual_scale'])*cal['noise']
                    preds[c,mi,ai,j]=to_linear_mean(pred,center)
                    diagnostics.append(dict(condition=c,axis=axis,sample_id=name,method=method,
                        fraction_negative_implied_log_expression=float(np.mean(pred+center<0)),
                        max_implied_log_expression=float(np.max(pred+center)),
                        fraction_above_log10001=float(np.mean(pred+center>np.log(10001))),
                        predicted_mean_linear=float(preds[c,mi,ai,j].mean())))
            print('PREDICTED',axis,c,flush=True)
    np.savez_compressed(O/'frozen_primary_predictions.npz',linear_profiles=preds,methods=np.array(methods),axes=np.array(AXES),sample_ids=np.array(names),genes=genes)
    pd.DataFrame(diagnostics).to_csv(O/'prediction_diagnostics.csv',index=False)
    checkpoints=json.loads((O/'predictions_locked.json').read_text())['checkpoint_sha256']
    locked={'adapter_sha256':sha(O/'adapter_locked.json'),'protocol_sha256':sha(O/'protocol.txt'),
        'primary_predictions_sha256':sha(O/'frozen_primary_predictions.npz'),
        'ten_seed_predictions_sha256':sha(O/'frozen_ten_seed_predictions.npz'),
        'checkpoint_sha256':checkpoints,'baseline_sources_sha256':baseline_sources,
        'prediction_script_sha256':sha(Path(__file__)),'pig_brain_expression_accessed':False,
        'primary_seed':42,'primary_condition':0,'shared_representation':True,'readout':'human-validation-selected shared calibration','ten_seed_readout':'original PeriScope; unchanged sensitivity, not calibrated main results','n_reference_cells_per_source':128,
        'completed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    dump(O/'predictions_locked.json',locked);print('ALL_EXTERNAL_PREDICTIONS_LOCKED',flush=True)

if __name__=='__main__':main()
