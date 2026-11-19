"""Extract the author's deposited guide panel and full reported perturbation table without significance filtering."""
from pathlib import Path
import csv,gzip,json,zipfile,hashlib,xml.etree.ElementTree as E
from collections import Counter
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'data/perturbation_sources_20261117/GSE178317'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=D/'41593_2022_1131_MOESM4_ESM.xlsx';before=sha(p);z=zipfile.ZipFile(p);assert z.testzip() is None
 uri='http://schemas.openxmlformats.org/spreadsheetml/2006/main';q='{'+uri+'}'
 strings=[''.join(e.itertext()) for e in E.fromstring(z.read('xl/sharedStrings.xml'))]
 rels={e.attrib['Id']:e.attrib['Target'] for e in E.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
 sheets={}
 for s in E.fromstring(z.read('xl/workbook.xml')).find(q+'sheets'):
  target=rels[s.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
  sheets[s.attrib['name']]=target.lstrip('/') if target.startswith('/') else 'xl/'+target
 def rows(name):
  with z.open(sheets[name]) as f:
   for _,e in E.iterparse(f,events=['end']):
    if e.tag==q+'row':
     row={}
     for c in e:
      v=c.find(q+'v');value=v.text if v is not None else ''
      if c.attrib.get('t')=='s':value=strings[int(value)]
      row[''.join(x for x in c.attrib['r'] if x.isalpha())]=value
     yield int(e.attrib['r']),row
     e.clear()
 table5=list(rows('Supplementary Table 5'));assert table5[0][1]=={'A':'target.gene','B':'sgRNA.name','C':'sgRNA.sequence'}
 panel=[dict(source_sheet_row=n,target_gene=r.get('A',''),guide_name=r.get('B',''),guide_sequence=r.get('C','')) for n,r in table5[1:]]
 with (D/'author_guide_panel.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(panel[0]));w.writeheader();w.writerows(panel)
 assert len(panel)==81
 targets=set(x['target_gene'] for x in panel);registry=json.loads((ROOT/'papers/03_world_model_integration/results/cross_domain_inputs_20261117/original_all_217_model_targets.json').read_text())
 links=[]
 for x in registry:
  hits=sorted(set(x['genes'])&targets)
  links.append(dict(model_target=x['target'],model_genes=';'.join(x['genes']),is_control=x['is_control'],exact_author_guide_gene_overlap=';'.join(hits),intervention_cell_matches=False,recipient_subtype_resolved=False))
 with (D/'all_217_target_coverage.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(links[0]));w.writeheader();w.writerows(links)
 pertarget=Counter();p_ranges={};reported_genes=set();total=0
 with gzip.open(D/'all_author_reported_RNA_effects.csv.gz','wt',newline='') as out:
  w=csv.writer(out);w.writerow(['source_sheet_row','TargetGene','Gene','log2CPM','log2FC','Pvalue','FDR'])
  for n,r in rows('Supplementary Table 9'):
   if n==1:assert [r[k] for k in 'ABCDEF']==['TargetGene','Gene','log2CPM','log2FC','Pvalue','FDR'];continue
   vals=[r.get(k,'') for k in 'ABCDEF'];w.writerow([n]+vals);total+=1;pertarget[vals[0]]+=1;reported_genes.add(vals[1])
   pv=float(vals[4]);lo,hi=p_ranges.get(vals[0],(1.,0.));p_ranges[vals[0]]=(min(lo,pv),max(hi,pv))
 # Preserve author normalized cluster values as reported; they are not raw proportions or the project's named states.
 t8=list(rows('Supplementary Table 8'));cols=sorted(t8[0][1],key=lambda c:(len(c),c))
 with (D/'author_normalized_cluster_values.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['source_sheet_row','author_cluster','target_gene','reported_normalized_value'])
  for n,r in t8[1:]:
   for col in cols:w.writerow([n,r.get('A',''),t8[0][1][col],r.get(col,'')])
 assert sha(p)==before and total==343707
 report=dict(complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),workbook_sha256=before,workbook_zip_crc_all_members=True,guide_rows=len(panel),author_target_labels=sorted(targets),non_NTC_targets=len(targets-{'NTC'}),model_overlap={x['model_target']:x['exact_author_guide_gene_overlap'] for x in links if x['exact_author_guide_gene_overlap']},reported_RNA_rows=total,reported_genes=len(reported_genes),rows_per_target=dict(pertarget),source_Pvalue_range=p_ranges,extraction_added_no_P_or_effect_filter=True,source_reported_gene_universe_is_not_claimed_complete_transcriptome=True,cell_level_assignment_available=False,biological_model_validation_complete=False,script_sha256=sha(Path(__file__)))
 (D/'AUTHOR_TABLE_EXTRACTION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
