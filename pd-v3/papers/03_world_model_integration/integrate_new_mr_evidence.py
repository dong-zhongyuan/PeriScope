"""Attach the replacement MR atlas to all original model targets without a new MR or selection score."""
from pathlib import Path
from datetime import datetime,timezone
import gzip,hashlib,json,os,sys
import pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from pdv3.native_imaging import HashReader
ATLAS=ROOT/'papers/01_mr_imaging_scan/results/full_atlas_20261116/full';CROSS=ROOT/'papers/03_world_model_integration/results/cross_domain_inputs_20261117';OUT=ROOT/'papers/03_world_model_integration/results/new_mr_integration_20261117'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def dump(n,j):
 p=OUT/n;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(p)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'RUN.json').exists():raise RuntimeError('Existing integration; inspect before restarting')
 run=dict(status='running',pid=os.getpid(),complete=False,started_at_utc=datetime.now(timezone.utc).isoformat(),policies={});dump('RUN.json',run)
 verification=json.loads((ATLAS/'COMPLETE_RESULT_VERIFICATION.json').read_text());pv=json.loads((ATLAS/'presentation/COMPLETE_PRESENTATION_VERIFICATION.json').read_text());assert verification['formal_full_atlas_verified'] and pv['formal_full_presentation_verified']
 assert json.loads((CROSS/'RUN.json').read_text())['complete']
 mapping=pd.read_csv(CROSS/'all_model_target_assay_links.csv',keep_default_na=False);assert mapping.model_target_index.nunique()==217
 selected=set(mapping.loc[mapping.exposure_id!='','exposure_id']);assert len(selected)==146 and all(',' not in s and '"' not in s for s in selected)
 pheno=pd.read_csv(ATLAS/'phenotypes.csv');assert len(pheno)==3963 and pheno.phenotype.is_unique
 inputs={str(x.relative_to(ROOT)):sha(x) for x in [Path(__file__),ROOT/'pdv3/native_imaging.py',CROSS/'RUN.json',CROSS/'all_model_target_assay_links.csv',CROSS/'all_targets_dose_response_across_seeds.csv',CROSS/'original_four_state_crosswalk.csv',ATLAS/'phenotypes.csv',ATLAS/'COMPLETE_RESULT_VERIFICATION.json',ATLAS/'presentation/COMPLETE_PRESENTATION_VERIFICATION.json']}
 summary=[]
 for policy in ['drop','infer']:
  path=ATLAS/f'full_family_palindrome_{policy}.csv.gz';target=OUT/f'all_linked_assay_imaging_rows_{policy}.csv.gz';temp=target.with_suffix(target.suffix+'.tmp');count=0;kept=0;payload=hashlib.sha256()
  with path.open('rb') as raw, temp.open('wb') as dest:
   reader=HashReader(raw)
   with gzip.GzipFile(fileobj=reader,mode='rb') as src,gzip.GzipFile(fileobj=dest,mode='wb',mtime=0,filename='') as out:
    header=src.readline();assert header.startswith(b'exposure_id,');out.write(header);payload.update(header)
    for line in src:
     count+=1;identity=line.split(b',',1)[0].decode('ascii')
     if identity in selected:out.write(line);payload.update(line);kept+=1
   assert reader.n==path.stat().st_size
  expected=verification['inputs_sha256'][str(path.relative_to(ROOT))];assert reader.h.hexdigest()==expected;inputs[str(path.relative_to(ROOT))]=expected
  assert count==2941*3963 and kept==len(selected)*3963
  temp.replace(target)
  linked=pd.read_csv(target,usecols=['exposure_id','phenotype','status','palindrome_policy'],keep_default_na=False)
  assert not linked.duplicated(['exposure_id','phenotype']).any() and set(linked.exposure_id)==selected and linked.palindrome_policy.eq(policy).all()
  assert linked.groupby('exposure_id').size().eq(3963).all() and set(linked.phenotype)==set(pheno.phenotype)
  coverage=linked.groupby(['exposure_id','status'],sort=False).size().unstack(fill_value=0);coverage.to_csv(OUT/f'all_linked_assay_status_counts_{policy}.csv')
  ranking=ATLAS/'presentation'/f'all_assay_ranking_{policy}.csv';assert sha(ranking)==pv['inputs_sha256'][str(ranking.relative_to(ROOT))];inputs[str(ranking.relative_to(ROOT))]=sha(ranking)
  rank=pd.read_csv(ranking,float_precision='round_trip',keep_default_na=False);assert len(rank)==2941 and rank.exposure_id.is_unique
  joined=mapping.merge(rank,on='exposure_id',how='left',validate='many_to_one',indicator=True)
  assert joined.loc[joined.exposure_id!='','_merge'].eq('both').all() and joined.loc[joined.exposure_id=='','_merge'].eq('left_only').all()
  joined['palindrome_policy']=policy;joined['MR_model_connection']='author_gene_label_only_no_unit_conversion';joined=joined.drop(columns='_merge');summary.append(joined)
  run['policies'][policy]=dict(complete=True,full_source_rows=count,linked_rows=kept,linked_assays=len(selected),all_phenotypes=3963,output_sha256=sha(target),selected_uncompressed_text_sha256=payload.hexdigest(),numeric_values='Source CSV lines preserved byte-for-byte; no MR estimates or q values recalculated',status_counts=linked.status.value_counts().to_dict());dump('RUN.json',run);print('MR_LINKED',policy,kept,flush=True)
 frame=pd.concat(summary,ignore_index=True);frame.to_csv(OUT/'all_model_target_MR_summary.csv',index=False)
 response=pd.read_csv(CROSS/'all_targets_dose_response_across_seeds.csv',float_precision='round_trip')
 assert len(response)==2604 and not response.duplicated(['axis','model_target_index']).any()
 together=frame.merge(response,on=['model_target_index','model_target','is_control'],how='left',validate='many_to_many')
 assert len(together)==len(frame)*12 and together.n_seeds.eq(10).all();together.to_csv(OUT/'all_model_targets_MR_and_expression_responses.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 pheno.to_csv(OUT/'all_3963_imaging_phenotypes.csv',index=False)
 for f,h in inputs.items():
  if 'full_family_palindrome_' not in f:assert sha(ROOT/f)==h,f
 dump('SOURCE_IDENTITIES.json',dict(input_sha256=inputs,original_model_targets=217,original_model_axes=12,original_seeds=10))
 run.update(status='complete',complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),model_targets_retained=217,matched_model_targets=int(mapping.loc[mapping.exposure_id!='','model_target_index'].nunique()),model_target_assay_rows=len(mapping),target_policy_summary_rows=len(frame),target_policy_axis_summary_rows=len(together),new_MR_branch=False,candidate_acceptance_criterion=None,measurement_scope='Genetic plasma-protein-to-MRI evidence and cellular-ADT-to-within-state-RNA model responses kept as distinct observed quantities. No product of incompatible effects, state-proportion or MRI model predictions. Both astro subtypes remain unresolved by the original astro model.')
 dump('RUN.json',run);dump('MANIFEST.json',{'outputs_sha256':{str(f.relative_to(ROOT)):sha(f) for f in OUT.iterdir() if f.is_file() and f.name not in ['MANIFEST.json','controller.log','JOB.json']}})
if __name__=='__main__':main()
