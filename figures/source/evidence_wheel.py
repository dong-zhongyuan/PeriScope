# User-supplied reference, preserved verbatim below the provenance line.
# fig5_e_wheel.py — Fig5b hero: 证据轮（cp AD figure3/src/panel_n_circular.py 改）
# 扇区 = 7 靶点; 9 环 = 证据层（淡=未过 / 实=过）; 外环柱 = 过层数; 圆心 = 层图例
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from pycirclize import Circos

# cp of the user-supplied wheel; data, sectors, labels, fonts and native size edited.
from pathlib import Path
import pandas as pd
import json
W = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'savefig.bbox':None})
m = pd.read_csv(W/'data/evidence_wheel.csv')
TARGETS = list(m.protein.unique())
layers = ['C1','C2','C3','C4','C5','C6','MR_gate']
names = ['Blood antigen','Measured dose','Perturbation','Brain disease','Accessibility','Clinical drug','Protein MR']
colors = ['#599CAF','#D0AE67','#C66C7E','#879C66','#9D88AE','#6D91AF','#4A9D8B']
all_cmap = [mcolors.ListedColormap([mcolors.to_rgba(c,.12),mcolors.to_rgba(c,1)]) for c in colors]
fig = plt.figure(figsize=(7.5,10))
ax = fig.add_axes([.01,.19,.98,.735],polar=True)
sector_group = {t:int(m.protein.eq(t).sum()) for t in TARGETS}
circo_group = Circos(sector_group, start=8, end=352, space=7, endspace=False)
for sector in circo_group.sectors:
    z = m[m.protein.eq(sector.name)]
    start_radius = 72
    for i,(l,c) in enumerate(zip(layers,all_cmap)):
        radius = start_radius-i*4.6
        track_group=sector.add_track((radius-4.2,radius))
        track_group.heatmap(z[l].eq('pass').astype(int).to_numpy().reshape(1,-1),vmin=0,vmax=1,cmap=c,rect_kws=dict(ec='white',lw=.27))
        for j,status in enumerate(z[l]):
            if status=='unresolved':track_group.rect(j,j+1,fc='#F5F5F5',ec='#9EA8AE',lw=.2,hatch='//')
    x=np.arange(len(z))+.5
    Sum=z[layers].eq('pass').sum(axis=1)
    track_sum=sector.add_track((77,93))
    track_sum.bar(x,Sum.to_numpy(),vmin=0,vmax=7,width=.76,color=['#B34F69' if v==7 else '#597B85' for v in Sum],ec='none')
    track_sum.line(x,Sum.to_numpy(),vmin=0,vmax=7,color='#526770',lw=.45,ls='--')
    sector.text(sector.name,r=102,size=13,weight='bold')
circo_group.plotfig(ax=ax)
# Ring key in the actual centre; labels fit without reducing the 12 pt font.
for i,(name,c) in enumerate(zip(names,colors)):
    yy=.621-i*.022
    fig.patches.append(mpatches.Rectangle((.379,yy-.005),.013,.01,transform=fig.transFigure,color=c,clip_on=False))
    fig.text(.402,yy,f'{i+1}  {name}',va='center',ha='left',fontsize=12)
# Explicit visual encodings, not an additional screening threshold.
handles=[mpatches.Patch(facecolor='#6B8F85',label='Met'),mpatches.Patch(facecolor='#EDF1EF',edgecolor='#BAC5C1',label='Not supported'),mpatches.Patch(facecolor='#F5F5F5',edgecolor='#9EA8AE',hatch='//',label='Unresolved'),mpatches.Patch(facecolor='#B34F69',label='All 7 met')]
fig.legend(handles=handles,loc='lower center',ncol=2,fontsize=12,frameon=False,bbox_to_anchor=(.5,.105),handlelength=1,handletextpad=.5,columnspacing=1.3)
fig.text(4/540,1-5/720,'D',fontsize=16,weight='bold',va='top')
fig.savefig(W/'build/7d.pdf')
# Register a native 3 x 4 module panel without resizing other panels.
p=W/'panel_manifest.json';manifest=json.loads(p.read_text());manifest['7d']=dict(id='7d',units_x=3,units_y=4,width_in=7.5,height_in=10,width_px=4500,height_px=6000,base_unit_px=1500);p.write_text(json.dumps(manifest,indent=2))
p=W/'template_usage.json';u=json.loads(p.read_text());u['7d']=dict(template='user_evidence_wheel',source='cp references/user_fig5_e_wheel.py -> source/evidence_wheel.py',n=len(m),question='Evidence alignment over every frozen program for the three nominated targets');p.write_text(json.dumps(u,indent=2))
plt.close(fig)
print('DRAW 7d',len(m),'unaltered program identities')
