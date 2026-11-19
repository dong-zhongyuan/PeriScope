"""Independently verify author worksheet values and all 160 CSF1R model comparisons."""
from pathlib import Path
from datetime import datetime,timezone
import json,zipfile,xml.etree.ElementTree as E,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];D=R/'data/perturbation_sources_20261117/GSE178317';O=R/'papers/03_world_model_integration/results/perturbation_author_readouts_20261117'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 sources=json.loads((O/'SOURCE_IDENTITIES.json').read_text());assert sha(D/'41593_2022_1131_MOESM4_ESM.xlsx')==sources['workbook_sha256']
 z=zipfile.ZipFile(D/'41593_2022_1131_MOESM4_ESM.xlsx');q='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}';strings=[''.join(e.itertext()) for e in E.fromstring(z.read('xl/sharedStrings.xml'))]
 expected={}
 with z.open('xl/worksheets/sheet9.xml') as f:
  for _,row in E.iterparse(f,events=['end']):
   if row.tag!=q+'row':continue
   vals={}
   for c in row:
    v=c.find(q+'v')
    if v is not None:vals[''.join(x for x in c.attrib['r'] if x.isalpha())]=strings[int(v.text)] if c.attrib.get('t')=='s' else v.text
   if vals.get('A')=='CSF1R':
    gene=vals['B'];assert gene not in expected;expected[gene]=[int(row.attrib['r'])]+[float(vals[c]) for c in 'CDEF']
   row.clear()
 actual=pd.read_csv(O/'all_8610_author_CSF1R_RNA_effects.csv.gz',float_precision='round_trip').set_index('Gene')
 assert len(expected)==8610 and set(actual.index)==set(expected)
 for gene,values in expected.items():assert np.array_equal(actual.loc[gene,['source_sheet_row','log2CPM','log2FC','Pvalue','FDR']].to_numpy(dtype=float),values)
 with np.load(O/'frozen_model_decrease_responses.npz',allow_pickle=False) as z:
  genes=z['gene'].astype(str);deltas=z['low_minus_high'];axes=z['axis'].astype(str);seeds=z['seed'];common=np.array([g in expected for g in genes]);assert common.sum()>1
  rows=pd.read_csv(O/'all_direction_comparisons.csv',float_precision='round_trip');assert len(rows)==160
  maximum=0.
  for i in range(80):
   for scope in ['all_shared_genes','exclude_perturbed_gene']:
    use=common if scope=='all_shared_genes' else common&(genes!='CSF1R')
    observed=np.array([expected[g][2] for g in genes[use]]);predicted=deltas[i,use];assert np.ptp(observed)>0 and np.ptp(predicted)>0 and np.isfinite(observed).all() and np.isfinite(predicted).all()
    row=rows[(rows.axis==axes[i])&(rows.seed==seeds[i])&(rows.gene_scope==scope)];assert len(row)==1;row=row.iloc[0];assert row.n_shared_genes==use.sum()
    pr=np.corrcoef(observed,predicted)[0,1];sr=pd.Series(observed).rank().corr(pd.Series(predicted).rank())
    nz=(observed!=0)&(predicted!=0);agree=np.mean(np.sign(observed[nz])==np.sign(predicted[nz]))
    maximum=max(maximum,abs(pr-row.pearson),abs(sr-row.spearman),abs(agree-row.sign_agreement))
  assert maximum<1e-12
 s=pd.read_csv(O/'direction_comparison_summary.csv',float_precision='round_trip').set_index(['axis','gene_scope'])
 for k,part in rows.groupby(['axis','gene_scope']):
  v=s.loc[k];assert part.seed.nunique()==v.seeds==10
  for field,actual_value in [('spearman_median',part.spearman.median()),('spearman_min',part.spearman.min()),('spearman_max',part.spearman.max()),('pearson_median',part.pearson.median()),('sign_agreement_median',part.sign_agreement.median())]:assert abs(v[field]-actual_value)<1e-12
 report=dict(complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),all_8610_source_worksheet_rows_exact=True,all_160_metrics_verified=True,shared_genes=int(common.sum()),maximum_metric_difference=maximum,source_statistics_are_author_reported_not_new_DE_fit=True,model_validation_complete=False,verifier_sha256=sha(Path(__file__)))
 (O/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
