# Reuse the copied drawing primitives and per-element renderer; current scenes redesigned.
from pathlib import Path
import json,copy,math,os
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.pdfbase.pdfmetrics import stringWidth
from workflow_primitives import Scene,BLUE as B,PURPLE as P,GREEN as G,GOLD as Y,INK
from PIL import Image
import pandas as pd
W=Path(__file__).resolve().parents[1];ROOT=W.parents[2];PANEL_BUILD=W/'build';ASSETS=ROOT/'pd_product_figure_design/workflow_assets';out=[]
E=pd.read_csv(ROOT/'delivery/optimization_20261003/mr_parallel_gate_20261003/all_programs_parallel_gate.csv');N=len(E);NP=int(E.priority_pass.sum());NT=E.loc[E.priority_pass,'protein'].nunique();NC=int(E.six_criteria_pass.sum());NA=E.loc[E.six_criteria_pass,'protein'].nunique()
def setup(n):
 s=Scene(n,720,360,['Reference_1C','Reference_1D','Reference_2A','Reference_3A','Reference_4A','Reference_6A','Reference_7A','Reference_S1D']);s.e[0].update(size=16,x=4,y=3,w=20,h=20);return s
def tx(s,x,y,w,h,t,b=False,align='center'):s.text(x,y,w,h,t,13 if b else 12,b,align)
def arrow(s,*p):s.arrow(*p,lw=.85,head=5)
def box(s,x,y,w,h,t,c=B):s.rect(x,y,w,h,c,INK,.65);tx(s,x+5,y+3,w-10,h-6,t,True)
def region(s,x,y,w,h,c='#F5F8FA'):s.rect(x,y,w,h,c,'#D8DFE3',.55)
def mat(s,x,y,c=B,cols=5,rows=4,cell=8):
 start=len(s.e);s.matrix(x,y,cols,rows,cell,[c,'#F5F7F8','#A9BFC8'])
 for e in s.e[start:]:
  if 'lw' in e:e['lw']=.65

def asset(s,key,cx,cy,h=65):
 with Image.open(ASSETS/f'{key}.png') as im:
  bb=im.convert('RGBA').getchannel('A').point(lambda a:255 if a>8 else 0).getbbox();sc=h/(bb[3]-bb[1]);s.asset(key,cx-(bb[0]+bb[2])/2*sc,cy-(bb[1]+bb[3])/2*sc,im.width*sc,im.height*sc)
def lock(s,x,y):s.lock(x,y,.35)
def seeds(s,x,y,vals,fill=G,cross=None):
 for i,v in enumerate(vals):
  xx=x+i*28;s.rect(xx,y,24,20,fill,'#95A8B2',.5);tx(s,xx,y,24,20,str(v))
  if v==cross:s.line((xx,y),(xx+24,y+20),stroke='#B85865',lw=1.1)
def finish(s):
 for e in s.e:
  if e['type']=='text':e['size']=16 if e['group']=='frame' else (13 if e['bold']else 12)
 out.append(dict(figure=s.n,width=s.w,height=s.h,reference_panels=s.refs,elements=s.e))
# Compact scientific process glyphs. They carry no empirical plotted values.
def network(s,x,y,w=100,h=54):
 for a in range(2):
  for i in range(3):
   for j in range(3):s.line((x+a*w/2,y+i*h/2),(x+(a+1)*w/2,y+j*h/2),stroke='#C2CCD3',lw=.6)
 for a in range(3):
  for i in range(3):s.ellipse(x+a*w/2-4,y+i*h/2-4,8,8,[B,P,G][a],INK,.5)
def vector(s,x,y,c=B,n=8):
 for i in range(n):s.rect(x+i*9,y,7,10,c if i%3 else'#DCE6E9','none',0)
def people(s,x,y,n=3,c=B,step=17):
 for i in range(n):
  xx=x+i*step;s.ellipse(xx+4,y,7,7,c,'#78949E',.5);s.rect(xx+2,y+9,11,16,c,'#78949E',.5)
def genes(s,x,y,w=49,h=62,c=G,rows=5,offset=0):
 s.rect(x,y,w,h,'#FFFFFF','#90A3AC',.55)
 for i in range(rows):
  s.rect(x+5,y+7+i*(h-12)/rows,w-10,3,c if(i+offset)%3 else'#DEE6E8','none',0)
def dna(s,x,y,w=55,h=52):
 a=[];b=[]
 for i in range(21):
  t=i/20;xx=x+t*w;yy=y+h/2+math.sin(t*math.tau)*h*.36;a.append((xx,yy));b.append((xx,2*y+h-yy))
 s.path(a,stroke='#609BB9',lw=1.4);s.path(b,stroke='#A995BF',lw=1.4)
 for i in range(2,20,3):s.line(a[i],b[i],stroke='#C9A66C',lw=.85)
def matchcloud(s,x,y):
 # Two cohorts represented by point sets and matching links.
 aa=[(0,0),(12,15),(-7,21),(16,-8),(-13,-9)];bb=[(59,2),(71,17),(54,23),(77,-5),(50,-10)]
 for i,(a,b) in enumerate(zip(aa,bb)):
  s.line((x+a[0],y+a[1]),(x+b[0],y+b[1]),stroke='#C4CCD3',lw=.5)
  for pt,c in [(a,B),(b,P)]:s.ellipse(x+pt[0]-3,y+pt[1]-3,6,6,c,'#66858F',.4)
def nullsets(s,x,y):
 genes(s,x,y,37,49,G,5)
 arrow(s,(x+44,y+24),(x+62,y+24))
 for j in range(3):genes(s,x+69+j*28,y+(j%2)*9,25,40,B,5,j)
def metric_glyph(s,key,x,y):
 if key=='MMD':matchcloud(s,x+17,y+16)
 elif key=='MSE':
  vector(s,x+2,y+1,B,5);vector(s,x+2,y+22,G,5);s.line((x+57,y+4),(x+64,y+4),(x+64,y+29),(x+57,y+29),stroke=INK,lw=.7)
 elif key=='Direction':
  for i,v in enumerate([1,-1,1,-1]):
   xx=x+6+i*18;s.line((xx,y+15),(xx,y+15-v*12),stroke=['#B66B79','#5B99B0'][v<0],lw=3)
  s.line((x,y+15),(x+72,y+15),stroke='#A0ACB3',lw=.5)
 elif key=='Blood':
  vector(s,x+4,y+1,B,5);lock(s,x+58,y-3);arrow(s,(x+23,y+16),(x+23,y+36));vector(s,x+4,y+38,P,5)
def barcode(s,x,y,w=70,h=18):
 for i in range(25):s.rect(x+i*w/25,y,(1+(i%3==0))*w/70,h,INK,'none',0)
def spatialgrid(s,x,y):
 for r in range(4):
  for c in range(6):s.ellipse(x+c*10+(r%2)*4,y+r*9,4,4,[B,G,P,'#DBBEC6'][(r+c)%4],'none',0)

# 1. Compact native geometry; font sizes are unchanged.
s=setup(1);s.w=630;s.group='observations'
region(s,20,28,594,91)
for cx,key,label in [(78,'pbmc','Blood RNA'),(289,'cite_cell','Paired RNA + ADT'),(518,'brain','Brain RNA reference')]:
 asset(s,key,cx,64,48);tx(s,cx-88,91,176,23,label,True)
mat(s,135,45,B,4,4,7);mat(s,347,45,G,4,4,7);mat(s,574,45,P,4,4,7)
s.group='blood_encoder';arrow(s,(88,121),(88,144));network(s,45,151,85,29);tx(s,20,185,159,22,'Blood RNA encoder',True);vector(s,44,211,B,10);lock(s,144,207);tx(s,20,231,159,20,'64-D blood state')
arrow(s,(289,121),(289,144));network(s,245,151,88,29);tx(s,206,185,166,22,'Protein head',True);vector(s,246,211,G,10);tx(s,207,231,166,20,'228 ADT channels');arrow(s,(140,165),(231,165))
s.group='equilibrium';region(s,197,267,202,77,'#FBF8F0');tx(s,205,267,184,20,'Coupled potential',True);asset(s,'energy',233,317,39)
s.ellipse(292,301,15,15,B,INK,.5);s.ellipse(370,301,15,15,P,INK,.5);s.line((308,308),(368,308),stroke='#85AE90',lw=2);tx(s,272,321,55,19,'Blood');tx(s,350,321,53,19,'Brain');tx(s,304,288,68,16,'Coupling')
arrow(s,(88,252),(88,308),(194,308));arrow(s,(289,252),(289,264))
s.group='brain_readout';arrow(s,(402,308),(436,308));tx(s,422,260,181,20,'Brain RNA decoder',True);network(s,469,287,93,25);arrow(s,(512,313),(512,320));mat(s,450,326,P,5,3,5);tx(s,495,322,110,24,'Predicted RNA')
s.group='reference_matching';matchcloud(s,470,163);arrow(s,(518,120),(518,145));arrow(s,(602,335),(615,335),(615,176),(562,176));tx(s,421,191,187,23,'Distribution matching',True);tx(s,425,217,179,35,'Disease × cell state\n4,741 shared genes')
finish(s)
# 2. Same-cell measurements, frozen prediction and donor-resolved evaluation.
s=setup(2);s.group='heldout_inputs';tx(s,24,30,199,25,'Held-out CITE-seq donors',True)
for i,t in enumerate(['P6','P7','P8']):
 asset(s,'cite_cell',55+i*70,92,52);tx(s,31+i*70,126,48,22,t,True)
region(s,25,168,194,156);tx(s,32,175,180,23,'Same-cell measurements',True);mat(s,53,219,B,6,5,9);mat(s,142,219,Y,6,5,9);tx(s,42,276,78,28,'RNA');tx(s,132,276,72,28,'Observed ADT')
s.group='protein_prediction';arrow(s,(221,238),(242,238),(242,161),(286,161));region(s,258,78,193,246,'#F4F8F4');lock(s,425,86);tx(s,267,90,168,22,'Frozen model',True)
network(s,298,140,103,42);tx(s,278,194,152,20,'RNA encoder');arrow(s,(353,216),(353,233));box(s,277,236,155,29,'Protein head',G);vector(s,308,281,G,10);tx(s,272,298,164,24,'228 ADT channels')
s.group='prediction_assessment';arrow(s,(434,250),(468,250),(468,219),(484,219));tx(s,491,30,205,26,'Observed / predicted ADT',True)
mat(s,518,80,Y,6,5,10);mat(s,619,80,G,6,5,10);tx(s,505,136,91,24,'Observed');tx(s,607,136,91,24,'Predicted')
box(s,489,185,207,54,'Matched cells\nwithin donor and subtype',G);arrow(s,(592,241),(592,266));metric_glyph(s,'MSE',515,280);matchcloud(s,623,296);tx(s,494,327,110,22,'MSE');tx(s,602,327,100,22,'Correlation')
finish(s)
# 3. Biological routes and comparator methods are distinct; pig is an evaluation domain.
s=setup(3);s.group='routes';tx(s,16,28,168,26,'Eight biological routes',True)
for j,t in enumerate(['cDC','cMono','nMono','pDC']):
 y=82+j*40;box(s,22,y,63,25,t,B)
 for yy in [108,210]:s.line((86,y+12),(119,yy),stroke='#A8B9C2',lw=.65)
box(s,123,94,62,28,'Astro',P);box(s,123,196,62,28,'MG-AP',P);tx(s,17,265,176,46,'Shared data splits\nSeeds 42, 43, 44')
s.group='methods';arrow(s,(189,164),(207,164));region(s,210,28,283,302,'#F8F5FA');tx(s,216,33,271,24,'Seven methods',True)
# Each comparator has its own structural glyph, retained as separate editable objects.
region(s,219,65,265,56,'#E3EFE3');asset(s,'energy',250,88,30)
s.ellipse(281,81,13,13,B,INK,.5);s.ellipse(318,81,13,13,P,INK,.5);s.line((295,88),(317,88),stroke='#719D85',lw=2)
tx(s,336,73,139,32,'PeriScope',True)
# OT-Ridge: transport links followed by a linear readout.
s.group='ot_ridge'
for i in range(3):
 s.ellipse(239,135+i*9,5,5,B,INK,.4);s.ellipse(273,136+i*9,5,5,P,INK,.4);s.line((244,137+i*9),(273,138+i*9),stroke='#A8B9C2',lw=.65)
s.line((293,160),(317,135),stroke='#719D85',lw=1.5);tx(s,220,169,124,23,'OT-Ridge')
s.group='conditional_mlp';network(s,382,136,62,25);tx(s,351,169,137,23,'Conditional MLP')
s.group='conditional_icnn';asset(s,'energy',281,220,35);tx(s,219,243,133,23,'Conditional ICNN')
s.group='label_mean'
for j in range(3):
 for k in range(4):s.rect(378+k*10,206+j*8,7,5,[B,P,G][j],'none',0)
arrow(s,(422,217),(435,217));vector(s,441,212,G,4);tx(s,351,243,137,23,'Label mean')
s.group='label_marginal'
for j,hh in enumerate([9,18,29,23,11]):s.rect(237+j*8,303-hh,6,hh,B,'none',0)
for xx,yy in [(298,282),(310,296),(319,278),(305,288)]:s.ellipse(xx,yy,5,5,P,INK,.4)
tx(s,219,304,133,23,'Label marginal')
s.group='cellot'
for i in range(4):
 y=276+(i%3)*10;s.ellipse(370+(i//3)*10,y,5,5,B,INK,.4);s.ellipse(443+(i//3)*10,y+3,5,5,P,INK,.4)
 s.line((377+(i//3)*10,y+2),(441+(i//3)*10,y+5),stroke='#B0A4C1',lw=.7)
tx(s,351,304,137,23,'CellOT')
s.group='evaluation';arrow(s,(495,164),(510,164));tx(s,512,29,195,26,'Evaluation',True)
for y,key,label in [(68,'MMD','MMD'),(118,'MSE','Mean MSE'),(167,'Direction','Disease rho'),(217,'Blood','Blood-only rho')]:
 metric_glyph(s,key,512,y);tx(s,597,y+1,109,36,label,True)
asset(s,'pig',553,307,38);tx(s,590,286,117,40,'Pig brain response',True)
s.group='protocol';tx(s,28,336,677,20,'Raw readout · Validation-calibrated readout')
finish(s)
# 4. Dose inputs, equilibrium, response tensor, discovery and confirmation.
s=setup(4);s.group='dose_input';tx(s,24,30,194,24,'Protein input perturbation',True);asset(s,'cite_cell',72,98,60)
for j,q in enumerate([5,25,50,75,95]):
 for k in range(4):s.rect(116+j*18,67+k*14,12,10,[B,G,Y,P,'#DEAEAC'][j]if k==1 else'#E1E8EC','none',0)
 tx(s,111+j*18,128,22,20,str(q))
s.rect(113,79,92,15,'none','#B64E65',.85);tx(s,28,170,188,47,'Measured ADT quantiles\nCLR re-centering')
s.group='solver';arrow(s,(221,104),(249,104));region(s,254,51,195,162,'#FAF7EF');tx(s,263,61,177,24,'Update brain state',True);asset(s,'energy',305,134,61);tx(s,354,111,86,65,'Blood fixed\nBrain free');box(s,273,174,158,28,'Decode brain RNA',P)
s.group='response';arrow(s,(450,131),(483,131));tx(s,489,31,210,25,'Gene-response tensor',True)
for off in [14,7,0]:
 s.rect(549+off,81-off,86,81,'#FFFFFF','#81959F',.65)
mat(s,559,90,P,8,7,9);tx(s,490,178,207,37,'4,741 genes × 5 doses\n10 current model seeds')
s.group='discovery_confirmation';s.line((27,240),(696,240),stroke='#D6DFE2',lw=.65);tx(s,24,250,194,23,'Discovery',True);seeds(s,52,284,range(42,47),P);arrow(s,(218,295),(264,295));genes(s,312,267,44,55,G);lock(s,363,271);tx(s,262,327,173,25,f'{N:,} frozen programs',True);arrow(s,(443,295),(480,295));tx(s,478,250,220,23,'Confirmation',True);seeds(s,518,284,range(47,52),G)
finish(s)
# 5. Parallel model response and cis-protein MR evidence, with a seed-omission diagram.
s=setup(5);s.group='model_perturbation';region(s,23,29,323,194,'#F4F8F9');tx(s,32,36,304,23,'Model perturbation',True)
nullsets(s,65,87);tx(s,37,144,93,36,'Frozen program genes');tx(s,153,144,178,36,'Size-matched sets\n5,000 draws')
arrow(s,(191,185),(191,198));tx(s,36,198,298,22,'BH across 6,849 programs')
s.group='protein_mr';region(s,372,29,324,194,'#F8F5FA');tx(s,380,36,308,23,'Blood protein → brain MRI',True)
dna(s,395,86,59,50);arrow(s,(464,112),(488,112));asset(s,'protein',531,110,44);arrow(s,(572,112),(594,112));asset(s,'brain',646,110,46)
tx(s,383,147,87,27,'cis-pQTL');tx(s,487,147,88,27,'Blood protein');tx(s,597,147,94,27,'Brain MRI');tx(s,380,194,308,23,'UKB · AGES')
s.group='joint_evidence';arrow(s,(184,225),(184,251),(287,251));arrow(s,(534,225),(534,251),(438,251));box(s,289,233,147,34,'Joint support',G)
s.group='robustness';arrow(s,(362,269),(362,288),(239,288));tx(s,25,298,89,34,'Seed\nomission')
for j,v in enumerate(range(47,52)):tx(s,116+j*23,272,22,19,str(v))
for r in range(5):
 for c in range(5):
  x=117+c*23;y=297+r*9;s.rect(x,y,19,6,G,'none',0)
  if r==c:s.line((x,y-1),(x+19,y+7),stroke='#B85865',lw=.8)
arrow(s,(239,320),(292,320));genes(s,305,299,31,39,G,5);tx(s,349,301,152,35,'Recompute gene-set q');arrow(s,(503,320),(533,320));vector(s,545,313,G,13);tx(s,540,332,147,21,'Program retention')
finish(s)
# 6. Matched fine states, frozen gene scores, and spatial subtype evidence.
s=setup(6);s.group='frozen_program';tx(s,22,30,160,23,'Frozen programs',True)
for off in [10,5,0]:s.rect(54+off,93-off,75,101,G,'#718991',.65)
for j in range(5):s.line((66,109+j*16),(115,109+j*16),stroke='#83AAA0',lw=3)
lock(s,141,86);tx(s,25,211,151,42,'Fixed gene lists\nMatched brain states')
s.group='brain_donors';region(s,217,29,479,135,'#F8F5FA');asset(s,'brain',251,82,43);tx(s,283,37,408,24,'Kamath · GSE157783 · Ma 2025',True)
people(s,287,85,3,B);people(s,344,85,3,'#E4B0BC');arrow(s,(399,102),(424,102));mat(s,437,77,P,6,5,8);arrow(s,(495,101),(520,101));vector(s,534,83,G,8);vector(s,534,103,P,8)
tx(s,282,127,128,24,'Donor × subtype');tx(s,422,127,92,24,'Pseudobulk');tx(s,523,127,166,24,'Program scores')
s.group='spatial_subtypes';region(s,217,193,479,140,'#F3F8FA');asset(s,'spatial_section',249,245,40);tx(s,282,202,410,23,'Spatial RNA · 10 donors')
barcode(s,283,247,61,21);arrow(s,(349,258),(365,258));spatialgrid(s,375,242);arrow(s,(437,257),(451,257));mat(s,462,240,G,6,4,9);arrow(s,(520,257),(540,257));genes(s,553,235,32,43,P,5);genes(s,592,235,32,43,G,5);genes(s,631,235,32,43,B,5)
tx(s,276,286,81,30,'Barcode');tx(s,356,286,97,30,'Coordinates');tx(s,446,284,92,35,'RNA fractions');tx(s,547,284,134,35,'Subtype programs')
s.group='branches';arrow(s,(166,143),(194,143),(194,97),(215,97));arrow(s,(166,143),(194,143),(194,263),(215,263))
finish(s)
# 7. Native molecular diagrams replace generic token icons.
from molecular_glyphs import antibody,membrane,receptor,protein,binder,molecule
s=setup(7);s.group='candidate_filter';region(s,24,31,196,140,'#F3F8F5')
for off in [10,5,0]:genes(s,37+off,63-off,40,62,G,6)
mat(s,98,54,G,6,5,8);tx(s,85,104,122,30,f'{NC} programs\n{NA} antigens');tx(s,30,140,185,22,'Six evidence criteria',True)
s.group='protein_mr';region(s,242,31,198,140,'#F7F3FA');dna(s,256,63,76,47);arrow(s,(338,88),(356,88));asset(s,'brain',395,83,47);tx(s,249,119,85,19,'cis-pQTL');tx(s,350,119,83,19,'Brain MRI');tx(s,250,142,182,22,'Blood protein MR',True)
s.group='integrated_targets';arrow(s,(220,101),(237,101));arrow(s,(440,101),(456,101));
for cx,label in [(493,'CD22'),(576,'CD40'),(659,'CD71')]:
 membrane(s,cx-30,116,64);receptor(s,cx,116,label,.76);tx(s,cx-34,31,68,23,label,True)
tx(s,466,151,224,22,f'{NP} programs · {NT} targets')
s.group='drug_mechanisms';s.line((24,186),(696,186),stroke='#D4DDE2',lw=.65)
# All three mechanism illustrations share the same 100 pt vertical frame.
for x,title in [(24,'Antibody'),(256,'Binder'),(488,'Inhibitor')]:
 region(s,x,203,208,131,'#FAFBFC');tx(s,x+9,307,190,24,title,True)
s.group='antibody_binding';membrane(s,52,280,143);protein(s,154,264,45,38);antibody(s,101,240,62);s.path([(127,229),(135,238),(141,246)],stroke='#C79465',lw=1.4)
s.group='binder_interface';membrane(s,283,280,142);protein(s,388,261,45,42);binder(s,324,246,49,45);s.path([(350,244),(356,245),(364,250)],stroke='#C79465',lw=1.4)
s.group='small_molecule';membrane(s,515,280,144);protein(s,580,258,71,42,cleft=True);molecule(s,581,229,10);s.line((581,239),(581,248),stroke='#9A5568',lw=1.5)
finish(s)

# Original per-object renderer retained; text reflow records real overflow instead of shrinking.
def lines_for(e):
 font='Helvetica-Bold'if e.get('bold')else'Helvetica';size=e['size'];w=e['w'];lines=[]
 for paragraph in e['text'].split('\n'):
  cur=''
  for word in paragraph.split(' '):
   trial=(cur+' '+word).strip()
   if stringWidth(trial,font,size)>w and cur:lines.append(cur);cur=word
   else:cur=trial
  lines.append(cur)
 return lines
issues=[]
for sc in out:
 if os.environ.get('FIGURE') and sc['figure']!=int(os.environ['FIGURE']):continue
 c=canvas.Canvas(str(PANEL_BUILD/f'{sc["figure"]}a.pdf'),pagesize=(sc['width'],sc['height']));c.setTitle(f'{sc["figure"]}a')
 for e in sc['elements']:
  typ=e['type'];fill=e.get('fill','none');stroke=e.get('stroke','none')
  c.setLineWidth(e.get('lw',.7));c.setFillColor(HexColor(fill)if fill!='none'else HexColor('#FFFFFF'));c.setStrokeColor(HexColor(stroke)if stroke!='none'else HexColor('#FFFFFF'))
  if typ in ('rect','ellipse'):
   fn=c.rect if typ=='rect'else c.ellipse
   args=(e['x'],sc['height']-e['y']-e['h'],e['w'],e['h'])if typ=='rect'else(e['x'],sc['height']-e['y']-e['h'],e['x']+e['w'],sc['height']-e['y'])
   fn(*args,fill=int(fill!='none'),stroke=int(stroke!='none'))
  elif typ=='path':
   p=c.beginPath();p.moveTo(e['pts'][0][0],sc['height']-e['pts'][0][1]);[p.lineTo(x,sc['height']-y)for x,y in e['pts'][1:]]
   if e.get('closed'):p.close()
   c.drawPath(p,fill=int(fill!='none'),stroke=int(stroke!='none'))
  elif typ=='image':
   paths=list(ASSETS.glob('*'+e['asset']+'*.png'))
   if not paths:raise RuntimeError(e['asset'])
   c.drawImage(str(paths[0]),e['x'],sc['height']-e['y']-e['h'],e['w'],e['h'],preserveAspectRatio=True,anchor='c',mask='auto')
  elif typ=='text':
   font='Helvetica-Bold'if e.get('bold')else'Helvetica';size=e['size'];c.setFont(font,size);c.setFillColor(HexColor(e['color']));lines=lines_for(e);lh=size*1.17;total=lh*len(lines)
   if total>e['h']+4 and e['text'].strip():issues.append((sc['figure'],e['text'],round(total-e['h'],1)))
   base=sc['height']-e['y']-(e['h']-total)/2-size*.85
   for j,t in enumerate(lines):
    x=e['x']if e.get('align')=='left'else e['x']+e['w']/2
    (c.drawString if e.get('align')=='left'else c.drawCentredString)(x,base-j*lh,t)
 c.showPage();c.save()
for s in out:
 for e in s['elements']:
  if e['type']=='text':e['text']='\n'.join(lines_for(e))
(W/'workflow_scenes.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
(W/'workflow_text_fit.json').write_text(json.dumps(issues,ensure_ascii=False,indent=2));print('Text fit flags',issues)
