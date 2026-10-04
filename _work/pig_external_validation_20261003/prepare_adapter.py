"""Train-data-only bulk response adapter; never opens pig brain expression."""
import itertools,time
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from common import *

def main():
    torch.set_num_threads(2);O.mkdir(parents=True,exist_ok=True)
    protocol_sha=sha(O/'protocol.txt')
    b=dict(np.load(I/'blood.npz'));genes=b['genes'].astype(str);states=list(b['state_vocab'].astype(str))
    pig=pd.read_csv(O/'inputs/blood_one2one_sum_provided_transcript_TPM.csv',index_col=0)
    assert pig.index.tolist()==genes.tolist()
    mask=pig.notna().all(axis=1).to_numpy();assert mask.sum()==3730
    train=~(b['is_val']|b['is_test']);val=b['is_val'];mu=b['gene_mu']
    assert np.all(b['gene_sd']==1)
    assert set(b['donor'][train]).isdisjoint(set(b['donor'][val]))
    donors=np.unique(b['donor'][train|val]);whole={};sub={};meta={};counts={}
    # Training and validation only, no existing human test donors.
    for d in donors:
        ix=np.flatnonzero((b['donor']==d)&(train|val))
        raw=np.expm1(np.maximum(b['X'][ix]+mu,0)).astype(np.float32)
        whole[int(d)]=raw.mean(0,dtype=np.float64)
        st=b['state'][ix]
        for bs in BS:
            ss=st==states.index(bs);counts[(int(d),bs)]=int(ss.sum())
            if ss.sum()>=25:sub[(int(d),bs)]=raw[ss].mean(0,dtype=np.float64)
        meta[int(d)]={'donor':str(b['donor_names'][d]),'condition':int(b['cond2'][ix[0]]),
                      'cohort':int(b['cohort'][ix[0]]),'split':'train' if train[ix[0]] else 'validation'}
    del raw
    cases=[];case_rows=[]
    for d in donors:
        d=int(d);mm=meta[d]
        refs=[k for k in donors if int(k)!=d and meta[int(k)]['split']=='train' and meta[int(k)]['condition']==0 and meta[int(k)]['cohort']==mm['cohort']]
        if not refs:continue
        bref=np.mean([whole[int(k)] for k in refs],0)
        for bs in BS:
            sr=[int(k) for k in refs if (int(k),bs) in sub]
            if (d,bs) not in sub or not sr:continue
            ref=np.mean([sub[(k,bs)] for k in sr],0)
            cases.append((d,bs,whole[d],bref,ref,sub[(d,bs)],mm['split']))
            case_rows.append(dict(**mm,blood_state=bs,n_reference_donors=len(sr),target_cells=counts[d,bs]))
    pd.DataFrame(case_rows).to_csv(O/'human_adapter_cases.csv',index=False)
    scores=[]
    for alpha,gain,cap in itertools.product([.1,1.,10.],[0.,.25,.5,.75,1.],np.log([2.,4.,8.])):
        by_donor={}
        for d,bs,target,bref,ref,truth,split in cases:
            if split!='train':continue
            f=logfactor(target,bref,mask,alpha,cap)
            pred=adapt(ref,f,gain)
            mse=float(np.mean((np.log1p(pred[mask])-np.log1p(truth[mask]))**2))
            by_donor.setdefault(d,[]).append(mse)
        scores.append(dict(alpha=alpha,gain=gain,cap_log=float(cap),human_CV_mse=float(np.mean([np.mean(v) for v in by_donor.values()])),n_donors=len(by_donor)))
    score=pd.DataFrame(scores).sort_values(['human_CV_mse','gain','alpha','cap_log'])
    score.to_csv(O/'human_adapter_training_grid.csv',index=False)
    selected=score.iloc[0].to_dict();alpha=selected['alpha'];gain=selected['gain'];cap=selected['cap_log']
    checks=[]
    for d,bs,target,bref,ref,truth,split in cases:
        f=logfactor(target,bref,mask,alpha,cap);pred=adapt(ref,f,gain)
        truth_delta=np.log1p(truth[mask])-np.log1p(ref[mask]);pred_delta=np.log1p(pred[mask])-np.log1p(ref[mask])
        checks.append(dict(donor=meta[d]['donor'],split=split,blood_state=bs,
            adapted_MSE=float(np.mean((pred_delta-truth_delta)**2)),
            no_change_MSE=float(np.mean(truth_delta**2)),
            direction_spearman=float(spearmanr(pred_delta,truth_delta).statistic) if np.ptp(pred_delta)>1e-10 else None))
    checks=pd.DataFrame(checks);checks.to_csv(O/'human_adapter_validation.csv',index=False)
    # Reference cells: equal quota for each available human control training donor.
    rng=np.random.default_rng(20261003);refout={};refrows=[]
    for bs in BS:
        st=states.index(bs);ii=np.flatnonzero(train&(b['cond2']==0)&(b['state']==st))
        ds=[int(d) for d in np.unique(b['donor'][ii]) if np.sum(b['donor'][ii]==d)>=25]
        selected_donors=np.resize(rng.permutation(ds),128)
        ix=np.array([rng.choice(ii[b['donor'][ii]==d]) for d in selected_donors])
        refout[bs]=b['X'][ix];refout[bs+'_indices']=ix
        refrows.extend([dict(blood_state=bs,reference_index=int(i),human_donor=str(b['donor_names'][b['donor'][i]])) for i in ix])
    refout.update(genes=genes,mask=mask,blood_center=mu,state_vocab=b['state_vocab'])
    np.savez_compressed(O/'human_reference_cells.npz',**refout)
    pd.DataFrame(refrows).to_csv(O/'human_reference_cells.csv',index=False)
    del b
    # Brain training-derived markers only; never pig brain expression.
    br=dict(np.load(I/'brain_train.npz'));tr=~br['is_val'];profiles={}
    for state in RS+['neural_other']:
        rows=[]
        for d in np.unique(br['donor'][tr]):
            ii=np.flatnonzero(tr&(br['donor']==d)&(br['state']==states.index(state)))
            if len(ii)<25:continue
            sums=np.zeros(len(genes));n=0
            for lo in range(0,len(ii),512):
                raw=np.expm1(np.maximum(br['X'][ii[lo:lo+512]]+br['gene_mu'],0))
                sums+=raw.sum(0,dtype=np.float64);n+=len(raw)
            rows.append(sums/n)
        profiles[state]=np.mean(rows,0)
    marker=[]
    for state in RS:
        enrich=np.log2((profiles[state]+.1)/(profiles['neural_other']+.1))
        eligible=mask&(profiles[state]>.1)&(enrich>0)
        ids=np.flatnonzero(eligible);ids=ids[np.argsort(-enrich[ids])[:100]]
        marker.extend([dict(receiver=state,gene=genes[i],index=int(i),human_training_enrichment_log2=float(enrich[i])) for i in ids])
    pd.DataFrame(marker).to_csv(O/'human_training_marker_sets.csv',index=False)
    selection={'protocol_sha256':protocol_sha,'selected':selected,'features':int(mask.sum()),
        'train_donor_macro_MSE':float(checks[checks.split.eq('train')].groupby('donor').adapted_MSE.mean().mean()),
        'validation_donor_macro_MSE':float(checks[checks.split.eq('validation')].groupby('donor').adapted_MSE.mean().mean()),
        'validation_no_change_MSE':float(checks[checks.split.eq('validation')].groupby('donor').no_change_MSE.mean().mean()),
        'validation_direction_median':float(checks[checks.split.eq('validation')].direction_spearman.median()) if gain>0 else None,
        'reference_sha256':sha(O/'human_reference_cells.npz'),'marker_sha256':sha(O/'human_training_marker_sets.csv'),
        'pig_brain_expression_accessed':False,'checkpoint_sha256':sha(A/'results/optimization_20261003/models/seed_42/model.pt'),
        'script_sha256':sha(Path(__file__)),'selected_gain_zero':bool(gain==0),
        'interpretation':'Control-referenced bulk-to-human-cell response transfer, not measured pig subtype deconvolution.'}
    dump(O/'adapter_locked.json',selection);print(json.dumps(selection,indent=2),flush=True)

if __name__=='__main__':main()
