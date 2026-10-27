from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
ROOT=Path(__file__).resolve().parents[4];W=Path(__file__).resolve().parents[1];D=ROOT/'delivery';O=D/'optimization_20261003';OUT=W/'data';OUT.mkdir(exist_ok=True);sources={}
def read(p,**kw):
 p=Path(p);sources[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p,**kw)
def save(d,n):d.to_csv(OUT/(n+'.csv'),index=False)
e=read(O/'mr_parallel_gate_20261003/all_programs_parallel_gate.csv');priority=e[e.priority_pass].sort_values(['protein','brain_state','blood_state','program_definition','response_direction']).copy()
# Preserve the recipient-matched Ma 2025 effects used in the revised Figure 6B.
# This is the frozen diagnostic export, joined by complete program identity.
ma=read(D/'tables/Figure6_Ma_2025_effects.csv',index_col=0)
assert ma.index.is_unique and set(ma.index)==set(priority.combo)
priority=priority.join(ma[['cell_state_sensitivity_ma_effect']],on='combo',validate='one_to_one')
assert priority.cell_state_sensitivity_ma_effect.notna().all()
save(priority,'priority')
assert len(e)==6849 and len(priority)==16
s=read(O/'joint_analysis/program_seed_scores.csv');save(s[s.combo.isin(priority.combo)],'priority_seeds')
acc=pd.concat([read(p) for p in sorted((O/'evaluation').glob('protein_accuracy_seed*.csv'))]);save(acc.groupby(['target','blood_state','citeseq_subtype','donor','n','is_control'],as_index=False)[['rho','mse']].mean(),'protein_accuracy')
coverage=acc[['donor','citeseq_subtype','blood_state','n']].drop_duplicates();assert not coverage.duplicated(['donor','citeseq_subtype']).any();save(coverage,'protein_coverage')
tr=[]
for seed in range(42,52):
 t=read(O/f'models/seed_{seed}/train_log.csv');t['seed']=seed;tr.append(t)
save(pd.concat(tr),'training')
b=read(D/'benchmark_20261002/baseline_metrics.csv');b=b[~b.method.eq('cellot')];raw=b[b.readout.eq('raw_shared')].copy();cal=b[b.readout.ne('raw_shared')].copy();c=read(D/'benchmark_20261002/cellot_restored/cellot_metrics.csv');c['readout']='original_cellot';raw=pd.concat([raw,c],ignore_index=True);c2=read(D/'benchmark_20261002/cellot_restored/residual_metrics.csv');c2['readout']='original_cellot_observation_residual';cal=pd.concat([cal,c2],ignore_index=True)
# Restored summary exports may already contain CellOT; keep one exact method/seed/axis/readout entry.
for t,n in [(raw,'benchmark_raw'),(cal,'benchmark_calibrated')]:
 assert not t.duplicated(['method','seed','axis']).any(),n
 assert t.method.nunique()==7 and len(t)==168,(n,t.method.value_counts())
 save(t,n)
bridge=read(D/'benchmark_20261002/bridge_benchmark_summary.csv');save(bridge,'protein_bridge')
pig=read(D/'pig_external_validation_20261003/primary_summary.csv');pig=pig[~pig.method.eq('cellot')];pc=read(D/'pig_external_validation_20261003/cellot_restored/primary_summary.csv');pig=pd.concat([pig,pc[pc.method.eq('cellot')]],ignore_index=True);save(pig,'pig_summary')
save(read(D/'pig_external_validation_20261003/gene_level_group_effects.csv'),'pig_genes')
# External robustness is evaluated over all ten seeds and all animal omissions.
save(read(D/'pig_external_validation_20261003/ten_seed_sensitivity.csv'),'pig_seed_sensitivity')
lo=read(D/'pig_external_validation_20261003/leave_one_animal_out.csv');lo=lo[~lo.method.eq('cellot')]
cp=D/'pig_external_validation_20261003/cellot_restored/leave_one_animal_out.csv'
if cp.exists():
 cc=read(cp);lo=pd.concat([lo,cc[cc.method.eq('cellot')]],ignore_index=True)
else:
 cc=read(D/'pig_external_validation_20261003/leave_one_animal_out.csv');lo=pd.concat([lo,cc[cc.method.eq('cellot')]],ignore_index=True)
assert not lo.duplicated(['removed_animal','method']).any()
save(lo,'pig_animal_sensitivity')
# Complete comparison tables retain raw/calibrated errors and external response amplitude.
(D/'tables').mkdir(exist_ok=True)
raw.to_csv(D/'tables/Figure3_benchmark_raw.csv',index=False)
cal.to_csv(D/'tables/Figure3_benchmark_calibrated.csv',index=False)
pig.to_csv(D/'tables/Figure3_pig_all_metrics.csv',index=False)

# No pooling across raw and calibrated predictions or across different biological axes.
abl=[]
for model,root in [('Full',O),('Random pairing',O/'ablations/random_pairing'),('Coupling off',O/'ablations/coupling_off')]:
 t=read(root/'evaluation/brain_direction_seed42.csv');t['model']=model;abl.append(t)
save(pd.concat(abl),'ablation_directions')
mr=read(D/'blood_brain_MR_20261003/all_brain_MR_estimates.csv');mr.to_csv(D/'tables/Figure5_all_MR_estimates.csv',index=False);save(mr[mr.gene.isin(['CD22','CD40','TFRC'])],'mr_priority_all_outcomes')
# Programme-level flow is restricted to the already defined, broader 6-criterion candidate set.
f=e[e.six_criteria_pass].copy();f['source']=f.blood_state;f['target']=f.MR_gate;f['third']=np.where(f.priority_pass,'MR + perturbation',np.where(f.MR_gate.eq('pending'),'MR pending','MR not supported'));f['count']=1;f['group']='Same program identities';save(f[['combo','source','target','third','count','group']],'evidence_flow')
save(e[['protein','MR_gate','MR_gene','MR_q']].drop_duplicates('protein'),'mr_coverage')
save(e[e.C7.eq('pass') & e.six_criteria_pass],'spatial_candidates')
leaf=read(O/'spatial_subtype_author_20261003/all_leaf_program_spatial_tests.csv');ss=read(O/'spatial_subtype_author_20261003/spatial_cell_state_donor_program_scores.csv');c7=e[e.C7.eq('pass') & e.six_criteria_pass];leaf=leaf[leaf.combo.isin(c7.combo)&leaf.leaf_C7.eq('pass')];assert len(leaf)==6
heat=ss.merge(leaf[['combo','spatial_subtype']],left_on=['combo','state'],right_on=['combo','spatial_subtype'],validate='many_to_one');save(heat,'spatial_heat_scores');save(leaf,'spatial_heat_annotation')
save(read(O/'spatial_subtype_author_20261003/state_identifiability_QC.csv'),'spatial_quality')
save(priority[['combo','protein','blood_state','brain_state','response_direction','program_definition','n_genes','q_gene_decoy','matched_brain_q','cell_state_sensitivity_kamath_effect','cell_state_sensitivity_gse157783_effect','matched_program_abundance_hypothesis']],'priority_matrix')
drug=read(O/'six_criteria_20261003/verified_drug_target_mechanisms.csv');drug=drug[drug.protein.isin(priority.protein)&drug.direct_interaction.eq(1)];drug=drug.drop_duplicates(['protein','molecule_chembl_id','action']);save(drug,'priority_drugs')
coord=read(ROOT/'pd_product_figure_design/template_semantic_review/data/spatial_coordinates.csv');fr=[]
for p in sorted((O/'spatial_subtype_author_20261003/samples').glob('GSM*/cell_state_RNA_fractions.csv')):
 t=read(p);t=t.rename(columns={t.columns[0]:'spot_id'});t['sample_gsm']=p.parent.name;fr.append(t)
fr=pd.concat(fr);assert not fr.spot_id.duplicated().any();sp=coord.merge(fr,on=['spot_id','sample_gsm'],how='left',validate='one_to_one');assert len(sp)==11859 and sp.astro_fibrous.notna().sum()==len(fr);save(sp,'spatial_fractions')
programs=json.loads((O/'joint_analysis/programs.json').read_text());sources[str((O/'joint_analysis/programs.json').relative_to(ROOT))]=hashlib.sha256((O/'joint_analysis/programs.json').read_bytes()).hexdigest()
from collections import Counter
reuse=Counter(g for gs in programs.values() for g in set(gs));save(pd.DataFrame([dict(gene=g,n_programs=n) for g,n in reuse.items()]),'gene_reuse')
files=O/'figure_inputs_20261004'
if files.exists():
 for p in files.glob('*.csv'):save(read(p),p.stem)
# The discovery annotation panel shows the actual top-200 query enrichment, not a newly independent validation of the GO-defined program.
go=read(O/'enrichment_terms.csv');go=go[(go.combo=='CD22__cDC__x__astro__up') & go.fold.eq('A')];save(go,'discovery_go')
(OUT/'source_manifest.json').write_text(json.dumps(sources,indent=2))

go=read(O/"enrichment_terms.csv");go=go[go.combo.eq("CD40__cDC__x__astro__up")&go.fold.eq("A")&go.library.eq("GO_Biological_Process_2023")].sort_values("q").head(6);save(go,"cd40_functional_go")

# Stable cohort panels reuse unchanged source data and stored ggplot objects.
extra=[ROOT/'pd_product_figure_design/panel_grid/objects'/f'{x}.rds' for x in ['1A','1B','1C','1D','1E','1G']]
extra += [ROOT/'pd_product_figure_design'/p for p in [
 'template_semantic_review/data/input_1A.csv',
 'template_semantic_review/data/real_cell_embedding.csv',
 'template_semantic_review/inputs/brain_diagnosis_lookup.json',
 'full_template_reuse/inputs/all_cohort_composition.csv']]
for p in extra:sources[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(OUT/'source_manifest.json').write_text(json.dumps(sources,indent=2))
print('PREPARED',len(sources),'source files; priority programs',len(priority),'spatial scores',len(heat))
