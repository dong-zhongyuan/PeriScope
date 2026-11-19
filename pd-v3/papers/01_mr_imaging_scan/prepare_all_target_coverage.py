"""All registered target coverage and explicit regional/reverse follow-up worklist."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,json,hashlib
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'papers/01_mr_imaging_scan/results/all_target_coverage_20261118'
DATA=ROOT/'data/all_target_followup_20261118'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
INPUTS={}
def read(rel):
 p=ROOT/rel;INPUTS[rel]=sha(p)
 with p.open() as f:return list(csv.DictReader(f))
def write(name,rows):
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
 OUT.mkdir(exist_ok=True)
 assays=read('data/ukb_ppp/extracted/assay_coverage.csv')
 selected={r['exposure_id']:r for r in read('data/ukb_ppp/extracted/cis_primary_selected.csv')}
 ranks={p:{r['exposure_id']:r for r in read(f'papers/01_mr_imaging_scan/results/full_atlas_20261116/full/presentation/all_assay_ranking_{p}.csv')} for p in ['drop','infer']}
 ages=read('data/ages_catalog/ukb_all_assays_to_ages.csv');amap=defaultdict(list)
 for r in ages:
  if r['accession']:amap[r['exposure_id']].append(r)
 scallop={r['exposure_id']:r for r in read('data/scallop_catalog/ukb_all_2941_assays.csv')}
 links=read('papers/03_world_model_integration/results/cross_domain_inputs_20261117/all_model_target_assay_links.csv')
 lmap=defaultdict(list)
 for r in links:
  if r['exposure_id']:lmap[r['exposure_id']].append(r)
 perturb=read('papers/03_world_model_integration/results/current_measurement_links_20261117/all_217_targets_by_three_perturbation_sources.csv')
 pmap=defaultdict(list)
 for r in perturb:
  if r['source_has_matched_model_antigen']=='True':pmap[r['model_target']].append(r['study'])
 batch=ROOT/'data/new_scan_regional_20261117/FIRST_DOWNLOAD_BATCH.json'
 INPUTS[str(batch.relative_to(ROOT))]=sha(batch)
 done={(r['exposure_id'],r['accession'],r['best_phenotype']) for r in json.loads(batch.read_text())['records']}
 coverage=[];work=[]
 for a in assays:
  eid=a['exposure_id'];d=ranks['drop'][eid];i=ranks['infer'][eid];s=selected.get(eid,{})
  pubs=sorted({g for l in lmap[eid] for g in pmap[l['model_target']]})
  common=dict(exposure_id=eid,gene=a['gene'],assay=a['assay'],uniprot=a['uniprot'],
   has_selected_cis_instrument=eid in selected,selected_variant=s.get('variant',''),
   full_MR_planned_phenotypes=3963,main_MR_computation='complete_with_missing_slots_retained',
   non_QC_planned_phenotypes=3947,non_QC_estimated_drop=d['n_estimated'],non_QC_estimated_infer=i['n_estimated'],
   best_phenotype_drop=d['best_phenotype'],best_global_q_drop=d['best_q_bh_global'],
   best_phenotype_infer=i['best_phenotype'],best_global_q_infer=i['best_q_bh_global'],
   ages_exact_gene_assays=len(amap[eid]),ages_accessions=';'.join(r['accession'] for r in amap[eid]),
   scallop_status=scallop[eid]['status'],scallop_same_gene_uniprot_files=scallop[eid]['n_same_gene_uniprot_public_files'],
   model_targets=';'.join(r['model_target'] for r in lmap[eid]),
   world_model_response='completed_all_ten_seeds' if lmap[eid] else 'not_in_frozen_model_antigen_universe',
   public_perturbation_source_overlaps=';'.join(pubs),
   matched_intervention_validation='not_completed',
   MRI_to_original_four_cellstates_calibration='not_completed',
   wetlab_validation='not_completed',
   first_batch_completed_pairs=sum(x[0]==eid for x in done),
   full_suite_complete=False)
  coverage.append(common)
  for policy,r in [('drop',d),('infer',i)]:
   for hit in amap[eid] or [None]:
    phenotype=r['best_phenotype'];has_pair=bool(phenotype) and eid in selected and hit is not None
    key=(eid,hit['accession'] if hit else '',phenotype)
    status='already_completed_first_batch' if key in done else ('pending_region_and_reverse' if has_pair else 'no_estimated_imaging_pair' if not phenotype else 'no_exact_AGES_source')
    work.append(dict(exposure_id=eid,gene=a['gene'],policy=policy,phenotype=phenotype,
     phenotype_description=r['best_trait'],accession=hit['accession'] if hit else '',
     ages_assay=hit['ages_assay'] if hit else '',chromosome=s.get('chrom',''),center_grch37=s.get('pos37',''),
     center_grch38=s.get('pos38',''),status=status,source_identity_scope='exact_gene_label_only',
     ordering_scope='one_best_non_QC_imaging_pair_per_policy_for_each_assay;operational_continuation_of_first_batch',
     no_result_based_assay_exclusion=True))
 write('all_2941_assay_coverage.csv',coverage);write('all_assays_first_pair_followup_worklist.csv',work)
 # All 217 model entries, including controls and unmatched channels, remain separately addressable.
 model=[]
 for target in sorted({r['model_target'] for r in links}):
  found=[r for r in links if r['model_target']==target];first=found[0]
  model.append(dict(model_target=target,model_target_index=first['model_target_index'],
   model_genes=first['model_genes'],is_control=first['is_control'],
   matched_UKB_assays=';'.join(sorted({r['exposure_id'] for r in found if r['exposure_id']})),
   dose_response='complete_all_seeds_and_axes',perturbation_source_overlap=';'.join(sorted(pmap[target])),
   matched_intervention_validation='not_completed',full_suite_complete=False))
 assert len(coverage)==len({r['exposure_id'] for r in coverage})==2941
 assert len(model)==217 and len(work)==2*len(ages)
 write('all_217_model_target_coverage.csv',model)
 result=dict(complete=True,scope='Coverage/worklist generation only, not completion of the scientific goal',
 all_assays=2941,all_model_targets=217,full_MR_pairs_per_policy=2941*3963,
 no_cis_instrument=sum(not r['has_selected_cis_instrument'] for r in coverage),
 ages_exact_gene_assays=sum(r['ages_exact_gene_assays']>0 for r in coverage),
 ages_unique_accessions=len({r['accession'] for r in ages if r['accession']}),
 worklist_rows=len(work),worklist_status_counts=dict(Counter(r['status'] for r in work)),
 unique_pending_source_phenotype_pairs=len({(r['accession'],r['phenotype']) for r in work if r['status']=='pending_region_and_reverse'}),
 model_matched_assays=sum(bool(r['model_targets']) for r in coverage),
 model_targets_with_listed_perturbation_overlap=sum(bool(r['perturbation_source_overlap']) for r in model),
 full_suite_complete=False,inputs=INPUTS,code_sha256=sha(Path(__file__)),
 outputs={p.name:sha(p) for p in OUT.glob('*.csv')})
 (OUT/'MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['inputs','outputs']},indent=2))
if __name__=='__main__':main()
