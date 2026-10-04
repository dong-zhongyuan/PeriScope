from pathlib import Path
import json,csv,hashlib,datetime
W=Path(__file__).resolve().parent;O=W.parents[2]/'delivery/blood_brain_MR_20261003';R=W/'raw/indexed_checked'
plan=json.loads((W/'indexed_plan.json').read_text());ex=list(csv.DictReader((O/'combined_exposure_instruments.csv').open()));ed={r['assay']:r for r in ex};rows=[];audit=[];files={};counts=[]
for x in plan['outcomes']:
 p=R/(x['accession']+'.sentinels.json');d=json.loads(p.read_text());assert len(d['audit'])==57 and not d.get('errors'),(p,len(d['audit']),d.get('errors'))
 if x['id']!='0011':assert d['engine']=='checked_http_range'
 for r in d['rows']:
  a=r['exposure_assay'];e=ed[a];build='GRCh37' if x['id']=='0011' else 'GRCh38';epos=int(e['pos37'] if build=='GRCh37' else e['pos']);assert int(r['base_pair_location'])==epos
  assert r['variant_id']==e['rsid']
  rows.append(dict(exposure_assay=a,outcome_accession=x['accession'],outcome_id=x['id'],outcome_trait=x['trait'],outcome_n=x['sample_size'],outcome_build=build,chromosome=r['chromosome'],position=r['base_pair_location'],rsid_outcome=r['variant_id'],effect_allele=r['effect_allele'],other_allele=r['other_allele'],beta_outcome=r['beta'],se_outcome=r['standard_error'],p_outcome=r['p_value'],eaf_outcome=r['effect_allele_frequency'],source_url=d['source_url']))
 for a in d['audit']:audit.append(dict(outcome_accession=x['accession'],**a))
 files[str(p.relative_to(W))]=hashlib.sha256(p.read_bytes()).hexdigest();counts.append(dict(outcome=x['trait'],accession=x['accession'],requested=57,matched=len(d['rows'])))
assert len({(r['exposure_assay'],r['outcome_accession']) for r in rows})==len(rows)
for name,data in [('brain_sentinel_statistics.csv',rows),('brain_extraction_audit.csv',audit)]:
 with (O/name).open('w') as f:w=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in data))));w.writeheader();w.writerows(data)
summary=dict(complete=True,n_outcomes=16,n_exposure_assays=57,n_requested_pairs=912,n_matched_pairs=len(rows),n_AGES_pairs=sum(r['exposure_assay'].startswith('GCST') for r in rows),n_UKB_pairs=sum(not r['exposure_assay'].startswith('GCST') for r in rows),matched=counts,unfiltered_sources=True,query_failures=0,sources_sha256=files,finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'brain_data_acquisition.json').write_text(json.dumps(summary,indent=2));print(summary['n_matched_pairs'],summary['n_AGES_pairs'],summary['n_UKB_pairs'])
