"""Restore an existing spatial assessment after a byte-identical six-gate rerun."""
from pathlib import Path
import argparse,json,hashlib
p=argparse.ArgumentParser();p.add_argument('--result-dir',required=True);a=p.parse_args();P=Path(a.result_dir)
C=P/'seven_criteria_20261003/completed.json'
if not C.exists():raise SystemExit(0)
c=json.loads(C.read_text());source=P/'six_criteria_20261003/all_6849_programs_six_criteria.csv'
if hashlib.sha256(source.read_bytes()).hexdigest()!=c['source_six_criteria_sha256']:
 print('Spatial assessment not reapplied: six-criterion input changed; rerun spatial program evaluation.')
 raise SystemExit(0)
S=P/'spatial_subtype_author_20261003'
if not (S/'verification.json').exists() or json.loads((S/'verification.json').read_text())['status']!='passed':raise SystemExit(0)
if hashlib.sha256((S/'protocol.json').read_bytes()).hexdigest()!=c['protocol_sha256']:raise SystemExit('Spatial protocol changed; do not reuse assessment.')
pointer=json.loads((P/'active_nomination.json').read_text());pointer.update(spatial_evidence_directory=S.name,spatial_assessment_directory='seven_criteria_20261003')
if c['spatial_layer_adoptable']:
 pointer.update(active_directory='seven_criteria_20261003',policy_id='20261003_seven_criteria_spatial_subtype',candidate_table='A7_antigens.csv',approved_subset='A7_approved_drug_subset.csv',pending_spatial_table='A6_pending_spatial_resolution.csv',previous_six_criteria_directory='six_criteria_20261003')
else:pointer.update(active_directory='six_criteria_20261003',policy_id='20261003_six_biological_and_drug_criteria',candidate_table='A6_antigens.csv',approved_subset='A6_approved_drug_subset.csv',spatial_supported_candidate_table='A7_antigens.csv',spatial_pending_program_table='A6_pending_spatial_resolution.csv',spatial_layer_status='Partial validation; global six-criterion screen retained')
if not c['spatial_layer_adoptable']:
 for obsolete in ['pending_spatial_table','previous_six_criteria_directory']:pointer.pop(obsolete,None)
(P/'active_nomination.json').write_text(json.dumps(pointer,indent=2))
print('Completed spatial assessment restored; original six-criterion input verified.')
