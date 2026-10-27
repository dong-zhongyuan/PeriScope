from pathlib import Path
import subprocess,fitz,json,os,shutil
W=Path(__file__).resolve().parents[1];ROOT=W.parents[2];OUT=W/'ppt/rendered';OUT.mkdir(exist_ok=True)
profile=W/'ppt/.build/lo_layout_profile';profile.mkdir(exist_ok=True);layouts=json.loads((W/'layout.json').read_text());checks=[]
only=int(os.environ.get('FIGURE','0'))
if only:checks=[x for x in json.loads((W/'audit_layout_20261004/lo_render_checks.json').read_text())if x['figure']!=only]
for row in layouts:
 if only and row['figure']!=only:continue
 n=row['figure'];src=ROOT/'delivery/pptx'/f'Figure{n}.pptx'
 subprocess.run([os.environ.get('PD_SOFFICE') or shutil.which('soffice') or '/Applications/LibreOffice.app/Contents/MacOS/soffice','-env:UserInstallation='+profile.resolve().as_uri(),'--headless','--convert-to','pdf','--outdir',str(OUT),str(src)],check=True,timeout=90)
 with fitz.open(OUT/f'Figure{n}.pdf') as d:
  assert len(d)==1;page=d[0];assert abs(page.rect.width-row['width'])<.1 and abs(page.rect.height-row['height'])<.1
  page.get_pixmap(matrix=fitz.Matrix(1.15,1.15)).save(OUT/f'Figure{n}.png')
  page.get_pixmap(matrix=fitz.Matrix(1.5,1.5),clip=fitz.Rect(0,0,720,360)).save(OUT/f'workflow{n}.png')
  checks.append(dict(figure=n,width=page.rect.width,height=page.rect.height,text_count=len(page.get_text().split()),source_mtime=src.stat().st_mtime))
 print('PPT_RENDERED',n,flush=True)
(W/'audit_layout_20261004/lo_render_checks.json').write_text(json.dumps(checks,indent=2))
shutil.rmtree(profile)
