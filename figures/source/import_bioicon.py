"""Convert the licensed vector paths to native editable scene primitives."""
from pathlib import Path
import fitz,json,hashlib
W=Path(__file__).resolve().parents[1];p=W/'assets/bioicons/antibody7.svg'
d=fitz.open(p);pdf=fitz.open('pdf',d.convert_to_pdf());page=pdf[0];out=[]
def color(c):return '#'+''.join(f'{round(v*255):02X}' for v in c[:3]) if c else 'none'
for dr in page.get_drawings():
    pts=[]
    for item in dr['items']:
        if item[0]=='l':
            if not pts:pts.append(tuple(item[1]))
            pts.append(tuple(item[2]))
        elif item[0]=='c':
            p0,p1,p2,p3=item[1:];
            if not pts:pts.append(tuple(p0))
            for j in range(1,17):
                t=j/16;pts.append(tuple((1-t)**3*p0[k]+3*(1-t)**2*t*p1[k]+3*(1-t)*t*t*p2[k]+t**3*p3[k] for k in range(2)))
        elif item[0]=='re':
            r=item[1];pts.extend([(r.x0,r.y0),(r.x1,r.y0),(r.x1,r.y1),(r.x0,r.y1)])
    if pts:out.append(dict(pts=pts,fill=color(dr['fill']),stroke=color(dr['color']),lw=dr['width'],closed=dr['closePath']))
box=[min(p[0]for e in out for p in e['pts']),min(p[1]for e in out for p in e['pts']),max(p[0]for e in out for p in e['pts']),max(p[1]for e in out for p in e['pts'])]
(W/'assets/bioicons/antibody_paths.json').write_text(json.dumps(dict(box=box,paths=out)))
record=dict(asset='antibody-7.svg',author='Servier',source='https://github.com/duerrsimon/bioicons/blob/main/static/icons/cc-by-3.0/Blood_Immunology/Servier/antibody-7.svg',license='CC BY 3.0 Unported',license_url='https://creativecommons.org/licenses/by/3.0/',modification='Converted curves to dense polylines; recoloured heavy/light chains; arranged with membrane and target schematics. Individual paths remain editable in PPTX and Blender.',sha256=hashlib.sha256(p.read_bytes()).hexdigest())
(W/'assets/bioicons/used_asset_license.json').write_text(json.dumps(record,indent=2))
print('Imported',len(out),'native paths',box)
