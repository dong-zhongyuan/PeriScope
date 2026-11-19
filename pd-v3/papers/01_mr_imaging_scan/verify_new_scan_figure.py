"""Check the generated regional figure, including source values and actual PDF exports."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,math
import pypdf
from pypdf import PdfReader
from PIL import Image
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'papers/01_mr_imaging_scan/results/new_scan_regions_20261117';OUT=BASE/'figures'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 # This generated vector PDF contains tens of thousands of points. Bound the
 # parser at10million tokens for this process; leave other limits unchanged.
 pypdf.overwrite_configuration(stream_parsing_work_maximum_count=10000000)
 pdf=OUT/'New_scan_public_regions.pdf';reader=PdfReader(pdf);assert len(reader.pages)==1;page=reader.pages[0];mm=[float(page.mediabox.width)*25.4/72,float(page.mediabox.height)*25.4/72];assert abs(mm[0]-174)<.5 and abs(mm[1]-225)<.5 and len(page.images)==0
 fonts=[]
 for x in page['/Resources']['/Font'].values():
  f=x.get_object();fonts.append(str(f['/BaseFont']));d=f.get('/FontDescriptor')
  if not d and '/DescendantFonts' in f:d=f['/DescendantFonts'][0].get_object().get('/FontDescriptor')
  assert d and any(k in d.get_object() for k in ['/FontFile','/FontFile2','/FontFile3'])
 assert all('Arial' in f for f in fonts);sizes=[]
 def visit(t,c,m,f,s):
  if t.strip():sizes.append(s*math.hypot(m[2]*c[0]+m[3]*c[2],m[2]*c[1]+m[3]*c[3]))
 text=page.extract_text(visitor_text=visit);assert min(sizes)>=6-1e-6
 for idp in ['1763','1441','2046']:assert idp in text
 with Image.open(OUT/'New_scan_public_regions.png') as im:pixel=im.size;dpi=im.info['dpi'];assert all(abs(d-600)<.1 for d in dpi)
 plot={n:pd.read_csv(OUT/'plot_data'/f'{n}.csv',float_precision='round_trip') for n in ['regional','shared','prior','reverse']};audit=json.loads((BASE/'INPUT_AUDIT.json').read_text());assert len(plot['shared'])==sum(x['shared_snps'] for x in audit['datasets']) and len(plot['prior'])==384 and len(plot['reverse'])==12
 assert len(plot['regional'])==sum(x['protein_source_rows']+x['imaging_source_rows'] for x in audit['datasets'])
 original=pd.read_csv(ROOT/'papers/01_mr_imaging_scan/results/new_scan_reverse_20261117/reverse_direction_estimates.csv',float_precision='round_trip');z=plot['reverse'];assert z.gene.tolist()==original.gene.tolist()
 for c in ['beta','ci_lower','ci_upper','n_instruments']:assert np.allclose(z[c],original[c],rtol=1e-12,atol=1e-14)
 for dataset in audit['datasets']:
  acc,gene=dataset['accession'],dataset['gene'];d=ROOT/'data/new_scan_regional_20261117'/acc;center=dataset['region']['center']
  for source,name in [('AGES protein','regional_statistics.csv.gz'),('Oxford imaging','imaging_region_GRCh37.csv.gz')]:
   frame=pd.read_csv(d/name,float_precision='round_trip');actual=plot['regional'][(plot['regional'].gene==gene)&(plot['regional'].source==source)]
   assert len(actual)==len(frame);assert np.allclose(actual.z,abs(frame.beta/frame.standard_error),rtol=1e-12,atol=1e-12);assert np.allclose(actual.relative_kb,(frame.base_pair_location-center)/1000,rtol=0,atol=1e-12)
 q=dict(complete=True,verified_at_utc=datetime.now(timezone.utc).isoformat(),actual_pdf_mm=mm,effective_font_range_pt=[min(sizes),max(sizes)],embedded_fonts=fonts,pdf_raster_images=len(page.images),png_dimensions=pixel,png_dpi=dpi,plot_data_rows={k:len(v) for k,v in plot.items()},all_regional_plotted_values_checked=True,phenotype_identifiers_in_actual_pdf=True,parser_budget=10000000,parser_budget_reason='Default2million token budget exhausted by dense generated vector scatter; bounded larger process-local budget, PDF unchanged',pdf_sha256=sha(pdf),verifier_sha256=sha(Path(__file__)))
 (OUT/'qa/EXPORT_VERIFICATION.json').write_text(json.dumps(q,indent=2)+'\n');print(json.dumps(q),flush=True)
if __name__=='__main__':main()
