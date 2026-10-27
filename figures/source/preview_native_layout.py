from pathlib import Path
import json,fitz
W=Path(__file__).resolve().parents[1];O=W/'audit_layout_20261004';B=W/'build'
for fig in json.loads((W/'layout.json').read_text()):
 doc=fitz.open();p=doc.new_page(width=fig['width'],height=fig['height'])
 for q in fig['all_panels']:
  with fitz.open(B/f'{q["id"]}.pdf') as src:
   rect=fitz.Rect(q['x'],q['y'],q['x']+q['width'],q['y']+q['height']);assert abs(rect.width-src[0].rect.width)<.01 and abs(rect.height-src[0].rect.height)<.01
   p.show_pdf_page(rect,src,0,keep_proportion=True)
 p.get_pixmap(matrix=fitz.Matrix(1.25,1.25)).save(O/f'Figure{fig["figure"]}.png');doc.close()
print('Native preview rendered')
