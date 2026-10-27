from pathlib import Path
import fitz,json
from PIL import Image,ImageDraw
w=Path(__file__).resolve().parents[1];fs=[w/'build'/f'{k}.pdf' for k in sorted(json.loads((w/'panel_manifest.json').read_text()))];good=[]
for f in fs:
 try:
  doc=fitz.open(f);pix=doc[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5));im=Image.frombytes('RGB',(pix.width,pix.height),pix.samples);good.append((f.stem,im))
 except Exception as e:print(f.name,str(e))
for i in range(0,len(good),12):
 out=Image.new('RGB',(1600,1500),'#E7EBED');d=ImageDraw.Draw(out)
 for j,(name,im) in enumerate(good[i:i+12]):
  im.thumbnail((390,455));x=j%4*400;y=j//4*500;out.paste(im,(x+(400-im.width)//2,y+25));d.text((x+5,y+5),name,fill='black')
 out.save(w/'build'/f'contact_{i//12+1}.png')
print('rendered',len(good))
