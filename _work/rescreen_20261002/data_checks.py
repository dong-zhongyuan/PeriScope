import json,hashlib
from pathlib import Path
import numpy as np
from input_utils import A,ROOT,REPORT,log_count_blocks,metadata

O=A/'results/rescreen_20261002';info={};datasets={}
for name in ['brain_train','brain_locked_test','brain_calibration','blood','citeseq']:
    z=np.load(ROOT/(name+'.npz'));datasets[name]=z
    X=z['X'];genes=z['genes'].astype(str)
    assert len(set(genes))==len(genes) and X.shape[1]==len(genes) and 'CR1' in genes
    assert len(set(z['barcodes']))==len(X)
    assert np.isfinite(X).all() and (X+z['gene_mu']).min()>-1e-5
    assert np.allclose(z['gene_sd'],1)
    info[name]={'shape':X.shape,'genes_sha256':hashlib.sha256('\n'.join(genes).encode()).hexdigest(),
        'donors':z['donor_names'][np.unique(z['donor'])].tolist(),'state_counts':dict(zip(z['state_vocab'][np.unique(z['state'])].tolist(),np.unique(z['state'],return_counts=True)[1].tolist())),
        'raw_log_min':float((X+z['gene_mu']).min())}
assert len({r['genes_sha256'] for r in info.values()})==1
b=datasets['brain_train'];t=datasets['brain_locked_test'];cal=datasets['brain_calibration']
arms=[set(b['donor_names'][np.unique(b['donor'][m])]) for m in [~b['is_val'],b['is_val']]]
arms += [set(z['donor_names'][np.unique(z['donor'])]) for z in [t,cal]]
assert all(not arms[i]&arms[j] for i in range(4) for j in range(i))
bl=datasets['blood'];c=datasets['citeseq']
assert np.array_equal(bl['gene_mu'],c['gene_mu'])
for z in [bl,c]:
    assert not np.any(z['is_test']&z['is_val'])
    ds=[set(z['donor'][m]) for m in [~(z['is_test']|z['is_val']),z['is_val'],z['is_test']]]
    assert all(not ds[i]&ds[j] for i in range(3) for j in range(i))
raw=np.log1p(c['adt_counts'].astype(np.float64));raw-=raw.mean(1,keepdims=True)
assert np.max(np.abs(raw-c['Y']))<1e-4
info['ADT_reconstruction_max_error']=float(np.max(np.abs(raw-c['Y'])))
info['ADT_CLR_max_abs_sum']=float(np.max(np.abs(c['Y'].sum(1))))
for source,dataset in [('brain',b),('blood223',bl),('bloodmb',bl)]:
    path=REPORT['inputs'][source]['path'];obs,_=metadata(path)
    names=dataset['barcodes'].astype(str)
    prefix='' if source=='brain' else ('b223:' if source=='blood223' else 'bmb:')
    eligible=np.where(np.char.startswith(names,prefix))[0] if prefix else np.arange(len(names))
    chosen=eligible[np.linspace(0,len(eligible)-1,20,dtype=int)]
    pos={v:i for i,v in enumerate(obs.index.astype(str))}
    rr=np.array([pos[names[i][len(prefix):]] for i in chosen])
    ix,x,_=next(log_count_blocks(path,rows=rr,selected=genes))
    error=float(np.max(np.abs(x.toarray()-(dataset['X'][chosen]+dataset['gene_mu']))))
    assert error<1e-5,(source,error)
    info[source+'_raw_count_crosscheck']=error
info['status']='passed'
(O/'data_checks.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2))
