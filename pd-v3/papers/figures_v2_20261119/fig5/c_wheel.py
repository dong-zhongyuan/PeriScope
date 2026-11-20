# Figure 5c - evidence wheel. Reuses the original PeriScope wheel code (figure_redesign_20261004/source/evidence_wheel.py,
# itself a copy of the user's fig5_e_wheel.py): sectors = the two final targets, angular columns = all of their frozen
# programs, rings = the agreed gates only (outside in: C1-C5 on the same program, protein MR, coloc on the same pair,
# direction-matched drug); outer bars = gates met, all met highlighted. Edits vs the original: data, ring set, 7 pt Arial,
# native 2x2 block (85.5 mm), no panel letter (added at assembly), no index numbers in the ring key.
# Data: six_criteria_20261003/all_6849_programs_six_criteria.csv (C1-C5 per program, as in build_target_funnel_final.py);
# target_funnel_final_20261119/antigen_gate_matrix.csv (gate6 = MR, gate7 = coloc, gate8_expanded = drug).
import os, sys, json, subprocess
from pathlib import Path
ROOT = Path(os.environ['FIGV2_ROOT']); sys.path.insert(0, str(ROOT/'_system/pylib'))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.colors as mcolors, matplotlib.patches as mpatches
from pycirclize import Circos

SIX = '/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/six_criteria_20261003/all_6849_programs_six_criteria.csv'
GATE = '/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/target_funnel_final_20261119/antigen_gate_matrix.csv'
TARGETS = ['CD71', 'CD119']
E = pd.read_csv(SIX, low_memory=False); E = E[E.protein.isin(TARGETS)].copy()
A = pd.read_csv(GATE).set_index('antigen')
tf = lambda v: str(v).strip().lower() in ('true', '1')
for g, col in [('MR', 'gate6'), ('Coloc', 'gate7'), ('Drug', 'gate8_expanded')]:
    E[g] = E.protein.map(lambda p: 'pass' if tf(A.loc[p, col]) else 'not_supported')
for k in ['C1', 'C2', 'C3', 'C4', 'C5']:
    E[k] = E[k].astype(str).str.strip().str.lower()
order = {'astro': 0, 'microglia_mhc2': 1, 'microglia_homeostatic': 2}
E['_b'] = E.combo.str.split('__').str[3].map(order).fillna(3); E['_s'] = E.combo
m = E.sort_values(['protein', '_b', '_s'], key=lambda s: s.map({t: i for i, t in enumerate(TARGETS)}) if s.name == 'protein' else s)
layers = ['C1', 'C2', 'C3', 'C4', 'C5', 'MR', 'Coloc', 'Drug']
names = ['C1 Antigen', 'C2 Dose', 'C3 Specificity', 'C4 PD brain', 'C5 Surface', 'MR', 'Coloc', 'Drug']
colors = ['#599CAF', '#D0AE67', '#C66C7E', '#879C66', '#9D88AE', '#4A9D8B', '#B07A4F', '#6D91AF']
(ROOT/'fig5/data').mkdir(exist_ok=True)
m[['combo', 'protein'] + layers].to_csv(ROOT/'fig5/data/evidence_wheel_final_targets.csv', index=False)

PT = 7.0; MM = 1/25.4; W = H = 85.5
plt.rcParams.update({'font.family': 'Arial', 'font.size': PT, 'pdf.fonttype': 42, 'savefig.bbox': None})
all_cmap = [mcolors.ListedColormap([mcolors.to_rgba(c, .12), mcolors.to_rgba(c, 1)]) for c in colors]
fig = plt.figure(figsize=(W*MM, H*MM))
ax = fig.add_axes([.03, .085, .94, .9], polar=True)
sector_group = {t: int(m.protein.eq(t).sum()) for t in TARGETS}
circo_group = Circos(sector_group, start=8, end=352, space=7, endspace=False)
NL = len(layers); step = 4.0
for sector in circo_group.sectors:
    z = m[m.protein.eq(sector.name)]
    start_radius = 72
    for i, (l, c) in enumerate(zip(layers, all_cmap)):
        radius = start_radius - i*step
        track_group = sector.add_track((radius - step*.9, radius))
        track_group.heatmap(z[l].eq('pass').astype(int).to_numpy().reshape(1, -1), vmin=0, vmax=1, cmap=c, rect_kws=dict(ec='white', lw=.2))
        for j, status in enumerate(z[l]):
            if status == 'unresolved': track_group.rect(j, j+1, fc='#F5F5F5', ec='#9EA8AE', lw=.15, hatch='//')
    x = np.arange(len(z)) + .5
    Sum = z[layers].eq('pass').sum(axis=1)
    track_sum = sector.add_track((76, 92))
    track_sum.bar(x, Sum.to_numpy(), vmin=0, vmax=NL, width=.76, color=['#B34F69' if v == NL else '#597B85' for v in Sum], ec='none')
    track_sum.line(x, Sum.to_numpy(), vmin=0, vmax=NL, color='#526770', lw=.35, ls='--')
    sector.text(sector.name, r=100, size=PT, weight='bold')
circo_group.plotfig(ax=ax)
# ring key in the centre, outer ring first (C1 is the outermost ring)
bb = ax.get_position(); cx = bb.x0 + bb.width/2; cy = bb.y0 + bb.height/2
dy = 2.5/H; y0 = cy + dy*(NL-1)/2; kx = cx - 10.6/W
for i, (name, c) in enumerate(zip(names, colors)):
    yy = y0 - i*dy
    fig.patches.append(mpatches.Rectangle((kx, yy - .8/H), 1.6/W, 1.6/H, transform=fig.transFigure, color=c, clip_on=False))
    fig.text(kx + 2.3/W, yy, name, va='center', ha='left', fontsize=PT)
handles = [mpatches.Patch(facecolor='#6B8F85', label='Met'), mpatches.Patch(facecolor='#EDF1EF', edgecolor='#BAC5C1', lw=.4, label='Not supported'),
           mpatches.Patch(facecolor='#F5F5F5', edgecolor='#9EA8AE', hatch='//', lw=.4, label='Unresolved'), mpatches.Patch(facecolor='#B34F69', label='All met')]
fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=PT, frameon=False, bbox_to_anchor=(.5, .005), handlelength=1, handleheight=1,
           handletextpad=.4, columnspacing=1.0, borderaxespad=0)
out = ROOT/'fig5/panels'; out.mkdir(exist_ok=True); stem = out/'5c_evidence_wheel_2x2'
fig.savefig(f'{stem}.pdf'); fig.savefig(f'{stem}.png', dpi=600); fig.savefig(f'{stem}_preview.png', dpi=150)
plt.close(fig)
json.dump({'id': '5c_evidence_wheel', 'block': '2x2', 'script': 'fig5/c_wheel.py', 'engine': 'python pycirclize 1.10.1 (original wheel code)',
           'sources': [SIX, GATE], 'n_programs': int(len(m)), 'rings': layers}, open(f'{stem}.manifest.json', 'w'), indent=1)
r = subprocess.run([os.environ['FIGV2_PY'], str(ROOT/'_system/check_panel.py'), f'{stem}.pdf', '--png', f'{stem}.png'], capture_output=True, text=True)
print(r.stdout.strip()[:1500], r.stderr.strip()[-800:])
