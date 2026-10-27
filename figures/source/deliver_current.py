# cp deliver_grid.py, then edit destinations, current manifest and 7-figure assembly.
from pathlib import Path
import os
import fitz,json,csv,shutil,gc,tempfile
from PIL import Image,ImageChops
from lxml import etree
W=Path(__file__).resolve().parents[1];ROOT=W.parents[2];D=ROOT/'delivery';B=W/'build';OUT=D/'panels';OUT.mkdir(exist_ok=True)
rows=json.loads((W/'panel_manifest.json').read_text());layout=json.loads((W/'layout.json').read_text());stage=Path(tempfile.mkdtemp(prefix='current-panels-',dir=W));checks=[]
only=int(os.environ.get('FIGURE','0'))
selected=set(filter(None,os.environ.get('PD_EXPORT_ONLY','').split(',')))
if only:checks=[x for x in json.loads((W/'export_validation.json').read_text())if not x['id'].startswith(str(only))]
if selected:checks=[x for x in json.loads((W/'export_validation.json').read_text())if x['id'] not in selected]
for id,r in sorted(rows.items()):
 if only and not id.startswith(str(only)):continue
 if selected and id not in selected:continue
 source=B/f'{id}.pdf';shutil.copy2(source,stage/f'{id}.pdf')
 with fitz.open(source) as doc:
  p=doc[0];root=etree.fromstring(p.get_svg_image(text_as_path=True).encode());root.set('width',f'{r["width_in"]}in');root.set('height',f'{r["height_in"]}in');(stage/f'{id}.svg').write_bytes(etree.tostring(root,encoding='UTF-8',xml_declaration=True))
  pix=p.get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False);assert (pix.width,pix.height)==(r['width_px'],r['height_px']),(id,pix.width,pix.height)
  im=Image.frombytes('RGB',(pix.width,pix.height),pix.samples);im.save(stage/f'{id}.png',dpi=(600,600));im.save(stage/f'{id}.tiff',dpi=(600,600),compression='tiff_lzw')
  with Image.open(stage/f'{id}.tiff') as t:assert ImageChops.difference(t,im).getbbox() is None
  for el in root.iter():
   for k,v in el.attrib.items():
    if k.endswith('href'):assert v.startswith(('#','data:'))
  checks.append(dict(id=id,width_px=pix.width,height_px=pix.height,dpi=600,png_tiff_identical=True,native_size=True))
  del pix,im;gc.collect()
 print('EXPORTED',id,flush=True)
assert len(checks)==49
for p in stage.iterdir():p.replace(OUT/p.name)
stage.rmdir()
# Removed numbered panels cannot remain as misleading previous versions.
for p in OUT.iterdir():
 if p.suffix.lower() in ['.pdf','.png','.svg','.tiff'] and p.stem[:1].isdigit() and p.stem not in rows:p.unlink()
with (OUT/'目录.csv').open('w',newline='',encoding='utf-8-sig') as f:
 wr=csv.DictWriter(f,fieldnames=list(next(iter(rows.values())).keys()));wr.writeheader();wr.writerows(rows.values())
(OUT/'说明.txt').write_text('49张独立子图。基础模块1500像素，600 dpi，2.5英寸；部分横向宽度按版面微调，4H/G各2500×1500像素。均为原生画布重排。普通文字12 pt，轴标题13 pt，字母16 pt。PDF/SVG为矢量，PNG/TIFF为600 dpi。主图拼版没有缩放子图。\n')
for n in ['preview','pdf']:(D/n).mkdir(exist_ok=True)
for fig in layout:
 if only and fig['figure']!=only:continue
 if selected and str(fig['figure']) not in {x[0] for x in selected}:continue
 doc=fitz.open();p=doc.new_page(width=fig['width'],height=fig['height'])
 for q in fig['all_panels']:
  with fitz.open(OUT/f'{q["id"]}.pdf') as src:
   rect=fitz.Rect(q['x'],q['y'],q['x']+q['width'],q['y']+q['height']);assert abs(rect.width-src[0].rect.width)<.001 and abs(rect.height-src[0].rect.height)<.001
   p.show_pdf_page(rect,src,0,keep_proportion=True)
 path=D/'pdf'/f'Figure{fig["figure"]}.pdf';tmp=path.with_suffix('.pending.pdf');doc.save(tmp,garbage=4,deflate=True);tmp.replace(path)
 pix=p.get_pixmap(matrix=fitz.Matrix(100/72,100/72));pix.save(str(D/'preview'/f'Figure{fig["figure"]}.png'));doc.close()
 print('ASSEMBLED',fig['figure'],flush=True)
(W/'export_validation.json').write_text(json.dumps(checks,indent=2))
(D/'panels'/'native_dimensions.json').write_text(json.dumps(rows,indent=2))
