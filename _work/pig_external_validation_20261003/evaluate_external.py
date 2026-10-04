"""Open independent measured pig brain only after all predictions are locked."""
import itertools,time
import numpy as np,pandas as pd
from scipy.stats import spearmanr,rankdata
from common import *

def corr(x,y,rank=True):
    if np.ptp(x)<1e-12 or np.ptp(y)<1e-12:return None
    if rank:x=rankdata(x);y=rankdata(y)
    return float(np.corrcoef(x,y)[0,1])

def response(linear,mask,ctrl):
    values=normalize(linear,mask)
    return np.log2((values+.1)/(values[ctrl].mean(0)+.1))

def effect(x,group):return x[group].mean(0)-x[~group].mean(0)

def metrics(p,y,group):
    a,b=effect(p,group),effect(y,group)
    v=np.mean(b*b);use=(abs(a)>1e-12)|(abs(b)>1e-12)
    return dict(spearman=corr(a,b),pearson=corr(a,b,False),
        response_RMSE=float(np.sqrt(np.mean((a-b)**2))),
        zero_response_RMSE=float(np.sqrt(v)),
        improvement_over_zero=float(1-np.mean((a-b)**2)/v) if v>0 else None,
        predicted_RMS=float(np.sqrt(np.mean(a*a))),observed_RMS=float(np.sqrt(v)),
        magnitude_ratio=float(np.sqrt(np.mean(a*a)/v)) if v>0 else None,
        sign_agreement=float(np.mean(np.sign(a[use])==np.sign(b[use]))),n_genes=len(a))

def exact_group_p(p,y,group):
    a=effect(p,group);r=corr(a,effect(y,group))
    if r is None:return dict(exact_group_p_positive=None,exact_group_p_two_sided=None)
    null=[]
    for subset in itertools.combinations(range(len(group)),int(group.sum())):
        g=np.zeros(len(group),bool);g[list(subset)]=True
        null.append(corr(a,effect(y,g)))
    v=np.asarray(null,float)
    return dict(exact_group_p_positive=float(np.mean(v>=r-1e-12)),
                exact_group_p_two_sided=float(np.mean(abs(v)>=abs(r)-1e-12)))

def paired_test(p,y,group):
    pp=p.copy();yy=y.copy()
    for g in [group,~group]:pp[g]-=pp[g].mean(0);yy[g]-=yy[g].mean(0)
    pn=np.linalg.norm(pp,axis=1);yn=np.linalg.norm(yy,axis=1)
    if np.any(pn<1e-12) or np.any(yn<1e-12):return {'within_group_pairing_cosine':None,'within_group_pairing_p':None},np.array([])
    cost=(pp/pn[:,None])@(yy/yn[:,None]).T
    cc=np.flatnonzero(~group);ll=np.flatnonzero(group)
    csum=np.array([cost[cc,list(v)].sum() for v in itertools.permutations(cc)])
    lsum=np.array([cost[ll,list(v)].sum() for v in itertools.permutations(ll)])
    null=((csum[:,None]+lsum[None,:])/len(group)).ravel()
    obs=float(np.trace(cost)/len(group))
    return dict(within_group_pairing_cosine=obs,within_group_pairing_p=float(np.mean(null>=obs-1e-12))),null

def bh(p):
    p=np.asarray(p,float);ix=np.argsort(p);v=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty_like(v);out[ix]=np.minimum(v,1);return out

def main():
    lock=json.loads((O/'predictions_locked.json').read_text())
    assert lock['primary_predictions_sha256']==sha(O/'frozen_primary_predictions.npz')
    assert lock['ten_seed_predictions_sha256']==sha(O/'frozen_ten_seed_predictions.npz')
    assert lock['adapter_sha256']==sha(O/'adapter_locked.json')
    assert lock['protocol_sha256']==sha(O/'protocol.txt')
    brain_path=O/'inputs/brain_one2one_sum_provided_transcript_TPM.csv'
    z=dict(np.load(O/'frozen_primary_predictions.npz'));genes=z['genes'].astype(str);names=z['sample_ids'].astype(str);methods=z['methods'].astype(str);axes=z['axes'].astype(str)
    ref=dict(np.load(O/'human_reference_cells.npz'));mask=ref['mask'];mapped_genes=genes[mask]
    meta=pd.read_csv(O/'inputs/matched_animals.csv').set_index('sample_id').loc[names].reset_index()
    group=meta.group.eq('LPS').to_numpy();ctrl=~group
    brain=pd.read_csv(brain_path,index_col=0).loc[genes,names].to_numpy().T
    assert np.isfinite(brain[:,mask]).all();assert group.sum()==4 and ctrl.sum()==6
    obs=response(brain,mask,ctrl)
    pd.DataFrame(obs.T,index=mapped_genes,columns=names).to_csv(O/'observed_brain_log2_response.csv',index_label='gene')
    audit=pd.read_csv(O/'inputs/existing_TPM_feature_audit.csv')
    hc=set(audit[(audit.tissue=='Brain')&audit.high_confidence].model_gene)
    markers=pd.read_csv(O/'human_training_marker_sets.csv')
    subsets={'all_one2one':np.ones(mask.sum(),bool),'high_confidence':np.isin(mapped_genes,list(hc))}
    for rs in RS:subsets[rs+'_markers']=np.isin(mapped_genes,markers[markers.receiver.eq(rs)].gene)
    primary={};allrows=[];axisrows=[];nulls={};group_genes=pd.DataFrame({'gene':mapped_genes,'observed_brain_effect':effect(obs,group)})
    for c in [0,1]:
        for mi,method in enumerate(methods):
            pp=np.stack([response(z['linear_profiles'][c,mi,ai],mask,ctrl) for ai in range(len(axes))])
            aggregate=pp.mean(0)
            if c==0:
                primary[method]=aggregate
                pd.DataFrame(aggregate.T,index=mapped_genes,columns=names).to_csv(O/f'predicted_response_{method}.csv',index_label='gene')
                group_genes[method]=effect(aggregate,group)
            for subset,keep in subsets.items():
                row=dict(method=method,condition=c,subset=subset,**metrics(aggregate[:,keep],obs[:,keep],group))
                if c==0 and subset=='all_one2one':
                    row.update(exact_group_p(aggregate,obs,group));pr,null=paired_test(aggregate,obs,group);row.update(pr);nulls[method]=null
                allrows.append(row)
            for ai,axis in enumerate(axes):
                axisrows.append(dict(method=method,condition=c,axis=axis,**metrics(pp[ai],obs,group)))
    blood=pd.read_csv(O/'inputs/blood_one2one_sum_provided_transcript_TPM.csv',index_col=0).loc[genes,names].to_numpy().T
    primary['no_blood_zero']=np.zeros_like(obs)
    for method in ['no_blood_zero']:
        pp=primary[method];row=dict(method=method,condition=0,subset='all_one2one',**metrics(pp,obs,group));row.update(exact_group_p(pp,obs,group));pr,null=paired_test(pp,obs,group);row.update(pr);nulls[method]=null;allrows.append(row);group_genes[method]=effect(pp,group)
    result=pd.DataFrame(allrows)
    primary_mask=result.condition.eq(0)&result.subset.eq('all_one2one')&result.method.isin(methods)
    valid=primary_mask&result.exact_group_p_positive.notna()
    result.loc[valid,'BH_q_shared_models']=bh(result.loc[valid,'exact_group_p_positive'])
    result.to_csv(O/'response_metrics.csv',index=False)
    pd.DataFrame(axisrows).to_csv(O/'all_six_axis_metrics.csv',index=False)
    group_genes.to_csv(O/'gene_level_group_effects.csv',index=False)
    np.savez_compressed(O/'within_group_pairing_nulls.npz',**nulls)
    # Paired-animal bootstrap: same resampled animal indices for blood predictions and brain truth.
    rng=np.random.default_rng(20261003);cc=np.flatnonzero(ctrl);ll=np.flatnonzero(group);bg=np.r_[np.zeros(6,bool),np.ones(4,bool)]
    boot=[]
    for rep in range(1000):
        ix=np.r_[rng.choice(cc,6,replace=True),rng.choice(ll,4,replace=True)]
        for method,pp in primary.items():
            m=metrics(pp[ix],obs[ix],bg);boot.append(dict(replicate=rep,method=method,spearman=m['spearman'],improvement_over_zero=m['improvement_over_zero']))
    boot=pd.DataFrame(boot);boot.to_csv(O/'paired_animal_bootstrap.csv',index=False)
    cis=[]
    for method,q in boot.groupby('method'):
        rr=q.spearman.dropna();im=q.improvement_over_zero.dropna()
        cis.append(dict(method=method,spearman_ci_low=float(rr.quantile(.025)) if len(rr) else None,
            spearman_ci_high=float(rr.quantile(.975)) if len(rr) else None,
            improvement_ci_low=float(im.quantile(.025)),improvement_ci_high=float(im.quantile(.975))))
    pd.DataFrame(cis).to_csv(O/'primary_bootstrap_intervals.csv',index=False)
    wide=boot.pivot(index='replicate',columns='method',values='spearman');diff=[]
    for method in wide.columns:
        if method in ['periscope','no_blood_zero']:continue
        v=(wide.periscope-wide[method]).dropna()
        diff.append(dict(comparator=method,delta_spearman_ci_low=float(v.quantile(.025)),delta_spearman_ci_high=float(v.quantile(.975)),bootstrap_fraction_periscope_higher=float((v>0).mean())))
    pd.DataFrame(diff).to_csv(O/'periscope_comparator_bootstrap.csv',index=False)
    seeds=dict(np.load(O/'frozen_ten_seed_predictions.npz'));seedrows=[]
    for si,seed in enumerate(seeds['seeds']):
        for c in [0,1]:
            pp=np.stack([response(seeds['linear_profiles'][si,c,ai],mask,ctrl) for ai in range(len(axes))]).mean(0)
            seedrows.append(dict(seed=int(seed),condition=c,**metrics(pp,obs,group)))
    pd.DataFrame(seedrows).to_csv(O/'ten_seed_sensitivity.csv',index=False)
    loo=[]
    for j,name in enumerate(names):
        keep=np.arange(len(names))!=j
        for method,pp in primary.items():loo.append(dict(removed_animal=name,removed_group=meta.group.iloc[j],removed_duration_min=int(meta.duration_min.iloc[j]),method=method,**metrics(pp[keep],obs[keep],group[keep])))
    pd.DataFrame(loo).to_csv(O/'leave_one_animal_out.csv',index=False)
    # Prespecified zero-gain sensitivity only if needed.
    extra=O/'unshrunk_sensitivity_predictions.npz'
    if extra.exists():
        zz=np.load(extra);rows=[]
        for c in [0,1]:
            pp=np.stack([response(zz['linear_profiles'][c,ai],mask,ctrl) for ai in range(len(axes))]).mean(0)
            rows.append(dict(condition=c,interpretation='unshrunk sensitivity; not human-selected adapter',**metrics(pp,obs,group)))
        pd.DataFrame(rows).to_csv(O/'unshrunk_sensitivity_metrics.csv',index=False)
    summary=result[result.condition.eq(0)&result.subset.eq('all_one2one')].merge(pd.DataFrame(cis),on='method',how='left')
    summary.to_csv(O/'primary_summary.csv',index=False)
    dump(O/'evaluation_completed.json',{'status':'completed_external_response_prediction_evaluation',
        'n_animals':10,'n_LPS':4,'n_control':6,'n_genes':int(mask.sum()),
        'brain_input_sha256':sha(brain_path),'predictions_lock_sha256':sha(O/'predictions_locked.json'),
        'script_sha256':sha(Path(__file__)),'model_retrained':False,'candidate_gates_changed':False,
        'evaluation_completed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        'primary_condition':'fixed_human_control','primary_seed':42,
        'interpretation':'Independent systemic-intervention bulk brain response test using reference-anchored blood input. No claim of measured pig subtype expression or target-specific intervention.'})
    print(summary.to_string(index=False),flush=True)

if __name__=='__main__':main()
