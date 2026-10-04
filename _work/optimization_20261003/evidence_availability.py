"""Separate untested cell states from evaluated, nonpassing disease support."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd

ap=argparse.ArgumentParser();ap.add_argument('--result-dir',required=True);ap.add_argument('--out-dir',required=True)
args=ap.parse_args();R=Path(args.result_dir);O=Path(args.out_dir);O.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(R/'all_program_evidence.csv');brain=json.loads((R/'brain_spatial_competitive.json').read_text())
ext=pd.read_csv(R/'independent_validation_stats.csv').set_index('combo');rows=[]
for r in df.to_dict('records'):
    sid='gsea::'+r['combo'];ma=brain['ma_sets'].get(sid,{});ka=brain['kamath_sets'].get(sid,{})
    mstat=ma.get('status','not_tested');kstat=ka.get('status','not_tested')
    evalbrain=mstat=='ok' and kstat=='ok';evalext=r['combo'] in ext.index
    scientific_fail=[];unassessed=[]
    for label,flag in [('protein_prediction','protein_prediction_pass'),('dose_range','measured_dose_pass'),
                       ('confirmation_response','reproducible_pass'),('specificity','specificity_pass')]:
        if not r[flag]:scientific_fail.append(label)
    if not r['mapping_pass']:unassessed.append('unresolved_antigen_identity')
    if not evalbrain:
        unassessed.extend(label+':'+status for label,status in [('Ma',mstat),('Kamath',kstat)] if status!='ok')
    elif not r['brain_pass']:scientific_fail.append('brain_disease_support')
    if not evalext:unassessed.append('GSE157783:insufficient_state_or_gene_coverage')
    elif evalbrain and not r['independent_direction_pass']:scientific_fail.append('independent_disease_direction')
    status='nominated' if r['tier'] in ['A','B'] else 'incomplete_coverage' if unassessed else 'tested_below_requirements'
    rows.append(dict(r,availability_status=status,ma_test_status=mstat,kamath_test_status=kstat,
        independent_test_available=evalext,unassessed_evidence=';'.join(unassessed),
        evaluated_nonpassing_evidence=';'.join(scientific_fail),
        all_nonbrain_gates_pass=bool(r['technical_pass'] and r['reproducible_pass'] and r['specificity_pass']),
        coverage_limited_nonbrain_pass=bool(unassessed and not scientific_fail)))
t=pd.DataFrame(rows);t.to_csv(O/'evidence_with_availability.csv',index=False)
follow=t[t.coverage_limited_nonbrain_pass].sort_values(['q_gene_decoy','protein_percentile'],ascending=[True,False])
follow.to_csv(O/'coverage_limited_followup_programs.csv',index=False)
summary=dict(source=str(R),status_counts=t.availability_status.value_counts().to_dict(),
    coverage_limited_nonbrain_pass_programs=len(follow),
    note='Missing validation coverage is not negative biology. Follow-up rows pass nonbrain gates; available brain evidence may still be weak or inconsistent. Original tiers retained, with no promotion to A or B.')
(O/'availability_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
