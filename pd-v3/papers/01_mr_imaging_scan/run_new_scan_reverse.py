"""Original reverse-direction follow-up: source imaging peaks, reference LD, AGES protein outcomes."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,sys,os,itertools
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from pdv3.mr import harmonize,wald,ivw
DATA=ROOT/'data/new_scan_regional_20261117';LD=DATA/'reverse_LD';OUT=ROOT/'papers/01_mr_imaging_scan/results/new_scan_reverse_20261117'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(n,j):(OUT/n).write_text(json.dumps(j,indent=2)+'\n')
def key(c,p,a,b):return f'{int(c)}:{int(p)}:'+'_'.join(sorted([a,b]))
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'RUN.json').exists():raise RuntimeError('Inspect existing reverse follow-up before rerunning')
 run=dict(status='running',pid=os.getpid(),complete=False,started_at_utc=datetime.now(timezone.utc).isoformat());dump('RUN.json',run)
 assert json.loads((LD/'RUN.json').read_text())['complete'] and json.loads((DATA/'RUN.json').read_text())['complete']
 records=json.loads((DATA/'FIRST_DOWNLOAD_BATCH.json').read_text())['records'];inputs={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'pdv3/mr.py',DATA/'FIRST_DOWNLOAD_BATCH.json',LD/'RUN.json']}
 genotypes={};ldrows={};reference_samples=None
 for d in sorted(LD.glob('chr*')):
  m=json.loads((d/'MANIFEST.json').read_text());assert m['complete']
  for n,h in m['outputs_sha256'].items():assert sha(d/n)==h;inputs[str((d/n).relative_to(ROOT))]=h
  inputs[str((d/'MANIFEST.json').relative_to(ROOT))]=sha(d/'MANIFEST.json')
  coverage=pd.read_csv(d/'coverage.csv')
  with np.load(d/'effect_allele_dosages.npz',allow_pickle=False) as z:
   samples=z['samples'];assert len(samples)==503
   if reference_samples is None:reference_samples=samples
   else:assert np.array_equal(samples,reference_samples)
   for k,g in zip(z['variant_keys'],z['dosages']):genotypes[k]=g.copy()
  for row in coverage.to_dict('records'):ldrows[row['variant_key']]=row
 all_audits=[];pairs=[];estimates=[];summary=[]
 for rec in records:
  acc=rec['accession'];folder=DATA/acc
  for n in ['all_source_imaging_peaks.csv','imaging_peak_outcomes.csv.gz','ACQUISITION.json']:inputs[str((folder/n).relative_to(ROOT))]=sha(folder/n)
  x=pd.read_csv(folder/'all_source_imaging_peaks.csv',float_precision='round_trip');y=pd.read_csv(folder/'imaging_peak_outcomes.csv.gz',float_precision='round_trip')
  x['variant_key']=[key(v.chromosome,v.position_grch37,v.effect_allele,v.other_allele) for v in x.itertuples()];y['variant_key']=[key(v.chromosome,v.base_pair_location,v.effect_allele,v.other_allele) for v in y.itertuples()]
  assert x.variant_key.is_unique and y.variant_key.is_unique
  correlations={}
  for a,b in itertools.combinations(x.to_dict('records'),2):
   ka,kb=a['variant_key'],b['variant_key'];within=a['chromosome']==b['chromosome'] and abs(a['position_grch37']-b['position_grch37'])<=10000000
   correlation=np.nan;n=0
   if ka in genotypes and kb in genotypes:
    ga,gb=genotypes[ka],genotypes[kb];valid=(ga>=0)&(gb>=0);n=int(valid.sum())
    if n>=3 and np.var(ga[valid])>0 and np.var(gb[valid])>0:correlation=float(np.corrcoef(ga[valid],gb[valid])[0,1])
   correlations[frozenset([ka,kb])]=correlation
   pairs.append(dict(accession=acc,variant_a=ka,variant_b=kb,n_complete_reference_samples=n,r=correlation,r2=correlation**2,within_10Mb_same_chromosome=within))
  # Select using imaging discovery association and LD alone, before looking at protein availability/effects.
  kept=[];clump=[]
  for a in x.sort_values(['discovery_neglog10_p','variant_key'],ascending=[False,True]).to_dict('records'):
   k=a['variant_key'];status='retained_by_exposure_only_LD_pruning';blocker='';r2=np.nan
   if not (len(a['effect_allele'])==len(a['other_allele'])==1 and a['effect_allele'] in 'ACGT' and a['other_allele'] in 'ACGT'):status='unsupported_alleles'
   elif a['discovery_neglog10_p'] < -np.log10(5e-8):status='outside_discovery_instrument_P_threshold'
   elif k not in genotypes:status='missing_reference_variant'
   elif np.all(genotypes[k]<0) or np.var(genotypes[k][genotypes[k]>=0])==0:status='monomorphic_or_no_reference_calls'
   else:
    for b in kept:
     if a['chromosome']==b['chromosome'] and abs(a['position_grch37']-b['position_grch37'])<=10000000:
      corr=correlations[frozenset([k,b['variant_key']])]
      if not np.isfinite(corr):status='LD_not_estimable';blocker=b['variant_key'];break
      if corr**2>0.001:status='LD_pruned';blocker=b['variant_key'];r2=corr**2;break
   if status=='retained_by_exposure_only_LD_pruning':kept.append(a)
   clump.append(dict(**a,accession=acc,gene=rec['gene'],ld_selection_status=status,ld_blocker=blocker,ld_blocker_r2=r2))
  kept_keys={a['variant_key'] for a in kept};cy=y.set_index('variant_key');audit=pd.DataFrame(clump);audit['protein_outcome_available']=audit.variant_key.isin(set(y.variant_key));all_audits.append(audit)
  for policy in ['drop','infer']:
   for effect_set in ['discovery','replication']:
    valid=[]
    for a in clump:
     row={k:a[k] for k in ['accession','gene','phenotype','variant_key','source_variant_id','chromosome','position_grch37','ld_selection_status']};row.update(palindrome_policy=policy,imaging_effect_set=effect_set,selected_after_LD=a['variant_key'] in kept_keys,protein_outcome_available=a['variant_key'] in cy.index)
     bx=float(a[effect_set+'_beta']);sx=float(a[effect_set+'_se']);row.update(beta_imaging=bx,se_imaging=sx)
     if a['variant_key'] not in cy.index:row['status']='protein_outcome_unavailable';estimates.append(row);continue
     b=cy.loc[a['variant_key']];h=harmonize(a['effect_allele'],a['other_allele'],b.effect_allele,b.other_allele,b.beta,a['effect_allele_frequency'],b.effect_allele_frequency,exposure_build='GRCh37',outcome_build='GRCh37',palindrome=policy)
     row.update(protein_beta_source=float(b.beta),protein_se=float(b.standard_error),protein_effect_allele=b.effect_allele,imaging_effect_allele=a['effect_allele'],harmonization_status=h['status'],protein_beta_aligned=h['beta_y'])
     if not h['usable']:row['status']=h['status'];estimates.append(row);continue
     if not np.isfinite([bx,sx]).all() or bx==0 or sx<=0:row['status']='invalid_imaging_effect_for_ratio';estimates.append(row);continue
     row.update(wald(bx,sx,h['beta_y'],b.standard_error));row['status']='estimated';estimates.append(row)
     if row['selected_after_LD']:valid.append(row)
    common=dict(accession=acc,gene=rec['gene'],phenotype=rec['best_phenotype'],palindrome_policy=policy,imaging_effect_set=effect_set,source_peaks=len(x),LD_selected_before_outcome_join=len(kept),usable_instruments=len(valid),direction='imaging_to_AGES_protein',p_scope='nominal_followup_P_not_a_support_decision')
    if len(valid)==0:common.update(status='not_estimated_no_usable_instrument',method='none',n_instruments=0)
    elif len(valid)==1:common.update(status='estimated',**wald(valid[0]['beta_imaging'],valid[0]['se_imaging'],valid[0]['protein_beta_aligned'],valid[0]['protein_se']))
    else:common.update(status='estimated',**ivw([v['beta_imaging'] for v in valid],[v['se_imaging'] for v in valid],[v['protein_beta_aligned'] for v in valid],[v['protein_se'] for v in valid],independent=True))
    summary.append(common)
 pd.concat(all_audits,ignore_index=True).to_csv(OUT/'all_peaks_exposure_only_LD_audit.csv',index=False);pd.DataFrame(pairs).to_csv(OUT/'all_peak_pair_reference_LD.csv',index=False);pd.DataFrame(estimates).to_csv(OUT/'all_peak_reverse_ratios_and_coverage.csv',index=False);pd.DataFrame(summary).to_csv(OUT/'reverse_direction_estimates.csv',index=False)
 for f,h in inputs.items():assert sha(ROOT/f)==h,f
 dump('SOURCE_IDENTITIES.json',dict(input_sha256=inputs))
 run.update(status='complete',complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),planned_gene_policy_effectset_rows=12,source_peak_rows=48,reference='1000 Genomes Phase3 EUR503 GRCh37; external, not cohort matched',LD_method='Pairwise complete observed genotype dosage Pearson r squared; greedy pruning by discovery P then variant key, r2 > .001 within 10 Mb. Not claimed to be a PLINK run.',LD_parameter_reference='https://mrcieu.github.io/TwoSampleMR/reference/clump_data.html',LD_selection_before_outcome_join=True,source_peak_selection='Reported Oxford peak table, not a fresh whole-genome peak discovery. Complete original source peak rows retained.',imaging_replication_scope='Replication effect sensitivity for the same discovery-selected instruments; not a separate population and no replication-P filter',protein_units='Reported AGES protein beta units, not assumed equivalent to UKB Olink or ADT CLR',MR_arithmetic='Unchanged pdv3.mr Wald and independence-asserted IVW after reference LD pruning; fixed and multiplicative random SE, heterogeneity, delta SE for individual ratios',limitations=['External small LD reference cannot prove cohort independence','UKB discovery-selected instrument effects may have winner bias; replication effects shown separately','Horizontal pleiotropy, tissue mediation and therapeutic direction not established','Nominal follow-up P values with no candidate acceptance threshold','No proxies substituted for missing outcomes or missing reference variants'])
 dump('RUN.json',run);dump('MANIFEST.json',dict(outputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in OUT.iterdir() if p.is_file() and p.name not in ['MANIFEST.json','controller.log','JOB.json']}));print(pd.DataFrame(summary)[['gene','palindrome_policy','imaging_effect_set','usable_instruments','status']].to_string(index=False),flush=True)
if __name__=='__main__':main()
