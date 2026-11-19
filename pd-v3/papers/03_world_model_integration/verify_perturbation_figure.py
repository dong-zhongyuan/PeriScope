"""Verify the actual exported perturbation figure and every plotted summary against source tables."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,math
import numpy as np,pandas as pd
from pypdf import PdfReader
from PIL import Image
ROOT=Path(__file__).resolve().parents[2];B=ROOT/'papers/03_world_model_integration/results/perturbation_readout_20261117';O=B/'figures'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 pdf=O/'Public_perturbation_direction_comparison.pdf';reader=PdfReader(pdf);assert len(reader.pages)==1;p=reader.pages[0]
 mm=[float(p.mediabox.width)*25.4/72,float(p.mediabox.height)*25.4/72];assert abs(mm[0]-174)<.5 and abs(mm[1]-225)<.5 and len(p.images)==0
 fonts=[]
 for obj in p['/Resources']['/Font'].values():
  f=obj.get_object();fonts.append(str(f['/BaseFont']));desc=f.get('/FontDescriptor')
  if not desc and '/DescendantFonts' in f:desc=f['/DescendantFonts'][0].get_object().get('/FontDescriptor')
  assert desc and any(k in desc.get_object() for k in ['/FontFile','/FontFile2','/FontFile3'])
 assert all('Arial' in f for f in fonts);sizes=[]
 def visit(text,cm,tm,font,size):
  if text.strip():sizes.append(size*math.hypot(tm[2]*cm[0]+tm[3]*cm[2],tm[2]*cm[1]+tm[3]*cm[3]))
 text=p.extract_text(visitor_text=visit);assert min(sizes)>=6-1e-6
 for token in ['CD86','IFNGR1','MERTK','unmatched','one iPSC clone']:assert token in text,token
 with Image.open(O/'Public_perturbation_direction_comparison.png') as im:
  pixels=im.size;dpi=im.info['dpi'];assert all(abs(x-600)<.1 for x in dpi) and pixels==(4110,5314)
 effects=pd.read_csv(B/'target_gene_readouts.csv',float_precision='round_trip')
 actual=pd.read_csv(O/'plot_data/target_effects.csv',float_precision='round_trip');assert len(effects)==len(actual)==168
 for c in ['target_gene_delta','n_cells','n_control_cells','target_raw_count','control_raw_count']:assert np.allclose(actual[c],effects[c],rtol=1e-12,atol=1e-14,equal_nan=True)
 data=pd.read_csv(B/'all_direction_comparisons.csv.gz',float_precision='round_trip')
 data=data[(data.kind=='gene')&(data.gene_scope=='exclude_perturbed_gene')]
 k=['target_gene','axis','sample'];expected=data[data.assignment=='author_inclusive'].groupby(k).spearman.median()
 m=pd.read_csv(O/'plot_data/well_correlations.csv',float_precision='round_trip').set_index(k);assert len(m)==168 and set(m.index)==set(expected.index)
 assert np.allclose(m.spearman,expected.reindex(m.index),rtol=1e-12,atol=1e-14,equal_nan=True)
 k+=['seed'];a=data[data.assignment=='author_inclusive'].set_index(k).spearman;b=data[data.assignment=='single_positive_guide'].set_index(k).spearman
 pair=pd.read_csv(O/'plot_data/assignment_sensitivity.csv',float_precision='round_trip').set_index(k);assert len(pair)==1680 and set(pair.index)==set(a.index)==set(b.index)
 assert np.allclose(pair.spearman_author,a.reindex(pair.index),rtol=1e-12,atol=1e-14,equal_nan=True) and np.allclose(pair.spearman_single,b.reindex(pair.index),rtol=1e-12,atol=1e-14,equal_nan=True)
 report=dict(complete=True,verified_at_utc=datetime.now(timezone.utc).isoformat(),actual_pdf_mm=mm,effective_font_range_pt=[min(sizes),max(sizes)],embedded_fonts=fonts,pdf_raster_images=len(p.images),png_dimensions=pixels,png_dpi=dpi,target_source_rows=168,target_points=int(np.isfinite(effects.target_gene_delta).sum()),undefined_target_groups=int((~np.isfinite(effects.target_gene_delta)).sum()),correlation_tiles=168,sensitivity_rows=1680,sensitivity_points=int(pair[["spearman_author","spearman_single"]].notna().all(axis=1).sum()),all_plotted_summaries_verified=True,undefined_sensitivity_pairs=int(pair.isna().any(axis=1).sum()),pdf_sha256=sha(pdf),verifier_sha256=sha(Path(__file__)),visual_inspection_still_required=True)
 (O/'qa/EXPORT_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
