from pathlib import Path
import pandas as pd,json,hashlib
W=Path(__file__).resolve().parent;ROOT=W.parents[2];O=ROOT/'delivery/blood_brain_MR_20261003';O.mkdir(exist_ok=True)
c=pd.read_csv(ROOT/'delivery/mr_feasibility_20261003/59_target_UKB_published_cis_coverage.csv').fillna('')
gene_ant={}
for _,r in c.iterrows():
 for g in r.genes.split(';'):
  gene_ant.setdefault(g,[]).append(r.protein)
d=pd.read_excel(W/'AGES_SD3.xlsx',header=3);a=pd.read_excel(W/'AGES_SD16.xlsx',header=2)
d=d[d['Protein (Entrez symbol)'].isin(gene_ant)&d['cis/trans'].eq('cis')].copy()
d['F']=(d.beta/d.se)**2;d=d[d.F>=10]
d['abs_z']=abs(d.beta/d.se)
best=d.sort_values('abs_z',ascending=False).drop_duplicates('SOMAmer').merge(a,left_on='SOMAmer',right_on='Study tag',validate='one_to_one')
best['candidate_antigens']=best['Protein (Entrez symbol)'].map(lambda g:';'.join(gene_ant[g]))
best.to_csv(O/'AGES_published_cis_instruments.csv',index=False)
coverage=[]
for _,r in c.iterrows():
 gs=r.genes.split(';');avail=sorted(set(best.loc[best['Protein (Entrez symbol)'].isin(gs),'Protein (Entrez symbol)']))
 coverage.append(dict(antigen=r.protein,genes=r.genes,AGES_cis_genes=';'.join(avail),n_AGES_assays=int(best['Protein (Entrez symbol)'].isin(gs).sum()),any_component_has_cis=bool(avail),all_components_have_cis=set(gs)<=set(avail)))
pd.DataFrame(coverage).to_csv(O/'59_target_independent_exposure_coverage.csv',index=False)
protocol={'created_before_outcome_extraction':True,'exposure_cohort':'AGES-Reykjavik, Iceland; Gudjonsson2022','specimen':'serum','exposure_selection':'Published cis pQTLs; strongest marginal absolute beta/SE per SOMAmer, F>=10. All measured candidate components retained; multiple aptamers stay separate sensitivity assays. No selection by brain associations.','units':'Verify original paper normalization before MR interpretation','harmonization':'Obtain exposure other allele and exact sample size from original full GWAS. Supplement has EA only. GRCh37 supplemental positions; harmonized GWAS usually GRCh38, never cross-match coordinates without build conversion.','source_hashes':{n:hashlib.sha256((W/n).read_bytes()).hexdigest() for n in ['AGES_SD3.xlsx','AGES_SD16.xlsx','brain_panel_lock.json']},'n_cis_assays':len(best),'n_candidate_antigens':sum(x['any_component_has_cis'] for x in coverage),'brain_panel':json.loads((W/'brain_panel_lock.json').read_text())}
(O/'data_protocol.json').write_text(json.dumps(protocol,indent=2))
print(best[['Protein (Entrez symbol)','Study Accession','rsID','EA','F','candidate_antigens']].to_string(index=False));print('N',len(best),sum(x['any_component_has_cis'] for x in coverage))
