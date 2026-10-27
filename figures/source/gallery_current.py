from pathlib import Path
import os
import json,csv,html,shutil,fitz
from PIL import Image,ImageDraw,ImageFont
W=Path(__file__).resolve().parents[1];R=W.parents[2];D=R/'delivery';V=D/'preview';(V/'comparisons').mkdir(exist_ok=True);(V/'panel_previews').mkdir(exist_ok=True)
m=json.loads((W/'panel_manifest.json').read_text());used=json.loads((W/'template_usage.json').read_text());font_path=Path(os.environ.get('PD_FONT_REGULAR','/System/Library/Fonts/Supplemental/Arial.ttf'));font=ImageFont.truetype(str(font_path),22) if font_path.exists() else ImageFont.load_default()
reg={x['id']:x for x in csv.DictReader((R/'scMOBA_template_library/registry.csv').open())};audit={x['template']:x for x in []};summaries=[]
old=W/'copy_review/reference_panels'
for id,r in sorted(m.items()):
 with fitz.open(D/'panels'/f'{id}.pdf') as d:
  p=d[0];pix=p.get_pixmap(matrix=fitz.Matrix(min(1.5,1000/p.rect.width),min(1.5,1000/p.rect.width)));pix.save(str(V/'panel_previews'/f'{id}.png'))
 if id.endswith('a'):
  source=old/f'{min(int(id[0]),6)}a.png';label='Earlier editable workflow';template='workflow'
 elif id=='6c':source=old/'5c.png';label='Earlier 5C: annotated heatmap';template='old5C'
 elif id=='7d':source=W/'references/user_wheel_structure.png';label='User wheel structure (illustrative data)';template='user_evidence_wheel'
 elif id=='7f':source=W/'references/UpSet_baseline_same_real_data.png';label='R UpSet baseline (same real data)';template='new_upset'
 else:
  template=used[id]['template'];source=R/'scMOBA_template_library'/reg[template]['preview'];label=f'Template {template}'
 assert source.exists(),source
 left=Image.open(source).convert('RGB');right=Image.open(V/'panel_previews'/f'{id}.png').convert('RGB');out=Image.new('RGB',(1800,1100),'white');draw=ImageDraw.Draw(out);draw.line((900,0,900,1100),fill='#c7cdd0',width=2)
 draw.text((24,14),label,font=font,fill='#17212b');draw.text((924,14),f'{id.upper()}  |  {r["width_in"]:g} x {r["height_in"]:g} in  |  native 12 pt',font=font,fill='#17212b')
 for i,im in enumerate([left,right]):
  im.thumbnail((850,1010));out.paste(im,(i*900+(900-im.width)//2,65+(1020-im.height)//2))
 out.save(V/'comparisons'/f'{id}.png');summaries.append(dict(id=id,template=template,source=str(source),width_in=r['width_in'],height_in=r['height_in']))
sections=[];titles={1:'数据与框架',2:'细胞内蛋白预测',3:'方法比较与消融',4:'剂量响应与程序',5:'模型扰动、MR 与稳健性',6:'独立脑队列与亚型空间证据',7:'候选靶点与药理信息'}
for n in range(1,8):
 cards=[]
 for id,r in sorted(m.items()):
  if not id.startswith(str(n)):continue
  cards.append(f'<article><h3>{id.upper()} · {r["width_in"]:g} × {r["height_in"]:g} 英寸</h3><a href="../panels/{id}.pdf"><img loading="lazy" src="panel_previews/{id}.png"></a><p><a href="../panels/{id}.svg">SVG</a> · <a href="../panels/{id}.tiff">TIFF</a> · <a href="comparisons/{id}.png">模板左右对比</a></p></article>')
 sections.append(f'<section id="fig{n}"><h2>Figure {n} · {titles[n]}</h2><p><a href="../pptx/Figure{n}.pptx">可编辑 PPT</a> · <a href="../pdf/Figure{n}.pdf">整图 PDF</a> · <a href="Figure{n}.png">整图预览</a></p><details><summary>展开主图拼版</summary><img class="main" loading="lazy" src="Figure{n}.png"></details><div class="grid">'+''.join(cards)+'</div></section>')
(V/'index.html').write_text('''<!doctype html><html lang="zh"><meta charset="utf-8"><title>PeriScope 当前图件</title><style>body{font:16px/1.6 -apple-system,sans-serif;background:#f3f5f6;color:#20292d;margin:28px auto;max-width:1500px;padding:0 24px}nav{position:sticky;top:0;background:#fff;padding:12px;z-index:2}a{color:#126478}h1{font-size:28px}h2{padding-top:20px}h3{font-size:17px}section{margin:30px 0}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}article{background:#fff;padding:18px;border:1px solid #d9e0e3}article img{width:100%;height:500px;object-fit:contain}.main{width:100%}summary{cursor:pointer}@media(max-width:850px){.grid{grid-template-columns:1fr}} </style><h1>PeriScope · 七张主图 / 四十九张独立子图</h1><p>基础模块 1500 px（600 dpi、2.5 英寸）。普通文字 12 pt、轴标题 13 pt、字母 16 pt。点击子图打开原生 PDF；页面缩略图仅用于浏览。</p><nav>'''+ ' · '.join(f'<a href="#fig{n}">Figure {n}</a>'for n in range(1,8))+'</nav>'+''.join(sections)+'</html>')
(W/'comparison_sources.json').write_text(json.dumps(summaries,ensure_ascii=False,indent=2))
counts={}
for u in used.values():counts[u['template']]=counts.get(u['template'],0)+1
(W/'template_counts.json').write_text(json.dumps(counts,indent=2));print('Gallery and 49 explicit source comparisons ready')
