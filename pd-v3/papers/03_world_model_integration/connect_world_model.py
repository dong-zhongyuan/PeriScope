"""Load the actual PeriScope implementation, verify all registered weights and real-input parity."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,importlib.util,sys,os,csv
import numpy as np,torch
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from pdv3.vendor.periscope.ccwm import CCWM,CCWMConfig
ASSETS=Path('/public/home/mengxl/dzy/pd_product_assets')
OLD=ASSETS/'results/optimization_20261003'
OUT=ROOT/'papers/03_world_model_integration/results/model_connection_20261117'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(n,j):(OUT/n).write_text(json.dumps(j,indent=2)+'\n')
def now():return datetime.now(timezone.utc).isoformat()
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'RUN.json').exists():raise RuntimeError('Inspect existing run before restarting')
 torch.set_num_threads(2)
 run=dict(status='running',pid=os.getpid(),started_at_utc=now(),complete=False);dump('RUN.json',run)
 registry_path=OLD/'model_registry.json';registry=json.loads(registry_path.read_text())
 source=Path('/public/home/mengxl/dzy/pd_product/_work/optimization_20261003/ccwm.py')
 assert sha(source)==sha(ROOT/'pdv3/vendor/periscope/ccwm.py')
 spec=importlib.util.spec_from_file_location('periscope_original_checked',source)
 original=importlib.util.module_from_spec(spec);sys.modules[spec.name]=original;spec.loader.exec_module(original)
 data_path=ASSETS/'interim/rescreen_20261002/blood.npz';data_hash=sha(data_path)
 with np.load(data_path,allow_pickle=False) as blood:
  X=blood['X'];states=blood['state_vocab'].astype(str).tolist();genes=blood['genes'].astype(str)
  train=~(blood['is_test']|blood['is_val']);pdmask=blood['cond2']==1
  # Same original PD training block; this is an implementation comparison, not new training or accuracy selection.
  ids=np.flatnonzero(train & pdmask)[:32];xb=torch.from_numpy(X[ids]).float()
  donor=np.asarray(blood['donor'][ids]).astype(str)
  np.savez_compressed(OUT/'real_input_subset.npz',X=X[ids],source_row=ids,donor=donor,genes=genes,state_vocab=np.array(states))
  assert len(ids)==32 and np.isfinite(X[ids]).all()
 rows=[];weight_refs=[]
 for seed in registry['seeds']:
  rec=registry['models'][str(seed)];path=Path(rec['path']);assert sha(path)==rec['sha256']
  ck=torch.load(path,map_location='cpu',weights_only=False)
  a=original.CCWM(original.CCWMConfig(**ck['cfg'])).eval();b=CCWM(CCWMConfig(**ck['cfg'])).eval()
  a.load_state_dict(ck['state_dict'],strict=True);b.load_state_dict(ck['state_dict'],strict=True)
  train_path=path.parent/'training_inputs.json';training=json.loads(train_path.read_text())
  assert ck['cfg']['n_genes']==len(genes)
  pairs=sorted({(int(bs),int(br)) for condition,bs,br,*_ in training['strata'] if condition==1})
  for bs,br in pairs:
   c=torch.zeros(len(xb),ck['cfg']['cond_dim']);c[:,1]=1;c[:,2+br]=1;c[:,-1]=1
   results=[]
   for model in [a,b]:
    with torch.no_grad():
     mu,lv=model.enc_blood(xb);u=model.protein_head(mu)
     if model.cfg.sample_latent_a:
      gen=torch.Generator(device='cpu').manual_seed(7701);mu=mu+torch.exp(.5*lv)*torch.randn(mu.shape,generator=gen)
    _,z,tr=model.solve_equilibrium(mu,mu.clone(),u,c,True)
    with torch.no_grad():y=model.dec_brain(z.detach())
    assert torch.isfinite(y).all();results.append((y,u,tr))
   diff=float((results[0][0]-results[1][0]).abs().max());udiff=float((results[0][1]-results[1][1]).abs().max())
   assert diff==0 and udiff==0
   rows.append(dict(seed=seed,blood_state=states[bs],brain_state=states[br],n_real_inputs=len(ids),maximum_RNA_difference=diff,maximum_ADT_difference=udiff))
  weight_refs.append(dict(seed=seed,path=str(path),sha256=rec['sha256'],training_inputs_sha256=sha(train_path),cfg=ck['cfg']))
  run.update(current_seed=seed,verified_seeds=len(weight_refs));dump('RUN.json',run);print('VERIFIED_SEED',seed,len(pairs),flush=True)
 with (OUT/'implementation_comparison.csv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
 assert sha(data_path)==data_hash and sha(source)==sha(ROOT/'pdv3/vendor/periscope/ccwm.py')
 dump('SOURCE_IDENTITIES.json',dict(model_registry_sha256=sha(registry_path),original_code_sha256=sha(source),blood_data_sha256=data_hash,weights=weight_refs))
 run.update(status='complete',complete=True,finished_at_utc=now(),seeds=len(weight_refs),conditions_checked=len(rows),maximum_RNA_difference=max(x['maximum_RNA_difference'] for x in rows),scientific_claim='Implementation equivalence only; not new validation of biological or causal predictions')
 dump('RUN.json',run)
 dump('MANIFEST.json',{'artifacts_sha256':{str(p.relative_to(ROOT)):sha(p) for p in OUT.iterdir() if p.is_file() and p.name not in ['MANIFEST.json','controller.log']}})
if __name__=='__main__':main()
