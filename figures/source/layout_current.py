from pathlib import Path
import json,fitz
W=Path(__file__).resolve().parents[1];B=W/'build';U=180
# Positions are in modules. Every panel remains at its native PDF dimensions.
positions={
1:{'a':(0,0),'b':(3.5,0),'c':(0,2),'d':(2,2),'e':(4,2),'f':(0,4),'g':(0,6),'h':(2,6)},
2:{'a':(0,0),'b':(4,0),'c':(4,1),'d':(4,3),'e':(0,2),'f':(4,5),'g':(0,6)},
3:{'a':(0,0),'d':(4,0),'b':(0,2),'e':(4,2),'c':(0,4),'f':(3.5,4)},
4:{'a':(0,0),'e':(4,0),'b':(0,2),'f':(4,2),'c':(0,4),'d':(4,4),'h':(8/3,5),'g':(13/3,5)},
5:{'a':(0,0),'c':(4,0),'b':(0,2),'e':(4,3),'d':(0,5),'f':(3,5)},
6:{'a':(0,0),'b':(4,0),'c':(0,2),'e':(4,3),'d':(0,4),'f':(4,4),'g':(0,6),'h':(3,6)},
7:{'a':(0,0),'e':(4,0),'b':(0,2),'d':(3,2),'c':(0,6),'f':(3,6)}}

scenes=json.loads((W/'workflow_scenes.json').read_text());manifest=json.loads((W/'panel_manifest.json').read_text());layout=[]
for n,xy in positions.items():
 items=[]
 for l,(x,y) in xy.items():
  id=f'{n}{l}'
  with fitz.open(B/f'{id}.pdf') as d:width,height=d[0].rect.width,d[0].rect.height
  assert abs(width/72*600-round(width/72*600))<1e-4 and abs(height/72*600-round(height/72*600))<1e-4,(id,width,height)
  q=dict(panel=l,id=id,x=x*U,y=y*U,width=width,height=height,svg=f'../../../../../../delivery/panels/{id}.svg',fallback=f'../../../../../../delivery/panels/{id}.png',description=f'Current real-data panel {id}; native dimensions, no image scaling')
  # Resolve these project-relative paths at use time so relocation needs no rewrite.
  q['svg']=f'delivery/panels/{id}.svg';q['fallback']=f'delivery/panels/{id}.png'
  items.append(q)
  manifest[id]=dict(id=id,units_x=width/U,units_y=height/U,width_in=width/72,height_in=height/72,width_px=round(width/72*600),height_px=round(height/72*600),base_unit_px=1500)
 width=max(q['x']+q['width']for q in items);height=max(q['y']+q['height']for q in items)
 # 4H uses the empty sixth facet of 4C and continues into the area below 4D; 4G completes that row.
 overlap_area=0
 for i,q in enumerate(items):
  for r in items[i+1:]:
   a=fitz.Rect(q['x'],q['y'],q['x']+q['width'],q['y']+q['height']);b=fitz.Rect(r['x'],r['y'],r['x']+r['width'],r['y']+r['height'])
   overlap=a&b
   if not overlap.is_empty:
    assert n==4 and {q['panel'],r['panel']}=={'c','h'},(n,q['panel'],r['panel'])
    overlap_area+=overlap.width*overlap.height
 area=sum(q['width']*q['height'] for q in items)-overlap_area;assert abs(area-width*height)<.1,(n,'unfilled panel rectangle',area,width*height)
 layout.append(dict(figure=n,width=width,height=height,workflow=next(q for q in items if q['panel']=='a'),panels=[q for q in items if q['panel']!='a'],all_panels=items))
(W/'panel_manifest.json').write_text(json.dumps(manifest,indent=2));(W/'layout.json').write_text(json.dumps(layout,indent=2));print('Native layouts',[(x['figure'],x['width']/72,x['height']/72)for x in layout])
