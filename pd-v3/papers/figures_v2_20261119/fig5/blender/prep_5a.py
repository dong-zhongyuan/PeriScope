# Build layout_5a.json for the Blender funnel (Figure 5a) from the antigen gate matrix: one disc per gate holding every
# antigen that reaches the gate (grey = passes, pink = removed here; CD71 and CD119 in their colours at the centre),
# sphere packing on a sunflower spiral; disc radius follows sqrt(count). Final disc = the two targets.
import json, math, sys
import pandas as pd
T = '/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/target_funnel_final_20261119/antigen_gate_matrix.csv'
A = pd.read_csv(T); tf = lambda s: s.astype(str).str.lower().isin(['true', '1'])
G = [('C1', 'gate1'), ('C2', 'gate2'), ('C3', 'gate3'), ('C4', 'gate4'), ('C5', 'gate5'), ('MR', 'gate6'), ('Coloc', 'gate7'), ('Direction', 'gate8a_direction'), ('Drug', 'gate8_expanded')]
alive = A.antigen.tolist(); tiers = []
for lab, col in G:
    ok = set(A.antigen[tf(A[col])]); keep = [a for a in alive if a in ok]; gone = [a for a in alive if a not in ok]
    tiers.append((lab, keep, gone)); alive = keep
assert sorted(alive) == ['CD119', 'CD71'], alive
TC = {'CD71': '#2E8B57', 'CD119': '#7B3F9E'}; KEEP, GONE = '#8C8C8C', '#E7A3BD'
W, H, cy, r = 174.0, 41.25, 18.6, .36
def sw(t): return .62 * 7 * .3528 * len(t) * 1.02          # approx. Arial width (mm) at 7 pt, for spacing only
def layout(c):
    xs, x = [], 1.0; items = tiers + [('Final', alive, [])]
    for lab, keep, gone in items:
        n = len(keep) + len(gone); R = c * math.sqrt(n) + r; w = max(2 * R, sw(lab) + 1.0); xs.append((x + w / 2, R, n)); x += w + 2.6
    return xs, x
c = .80
while layout(c)[1] > W - 1: c -= .01
xs, _ = layout(c); E, Lb = [], []
for (lab, keep, gone), (cx, R, n) in zip(tiers + [('Final', alive, [])], xs):
    order = [a for a in keep if a in TC] + sorted(a for a in keep if a not in TC) + sorted(gone)
    for k, a in enumerate(order):
        rr = c * math.sqrt(k + .5); th = k * 2.39996323
        col = TC.get(a, KEEP if a in keep else GONE) if a in keep else GONE
        E.append({'type': 'sphere', 'x': round(cx + rr * math.cos(th), 3), 'y': round(cy + rr * math.sin(th), 3), 'r': r if lab != 'Final' else .8, 'col': col})
    Lb.append({'t': lab, 'x': round(cx, 2), 'y': 34.4})
for (a, b) in zip(xs[:-1], xs[1:]):
    x0 = a[0] + max(a[1], .5) + .5; x1 = b[0] - max(b[1], .5) - .5
    if x1 - x0 > 1.2: E.append({'type': 'arrow', 'p': [[x0, cy], [x1, cy]], 'hl': .8, 'hw': .42, 'lw': .2})
L = {'W': W, 'H': H, 'dpi': 600, 'sun': .3, 'world': .9, 'elements': E, 'labels': Lb,
     'keys': {'keep': KEEP, 'gone': GONE, 'CD71': TC['CD71'], 'CD119': TC['CD119']}, 'spiral_c': c,
     'tiers': [{'gate': lab, 'n_in': len(k) + len(g), 'n_pass': len(k)} for lab, k, g in tiers]}
json.dump(L, open(sys.argv[1], 'w'), indent=1); print('tiers', [(t['gate'], t['n_in'], t['n_pass']) for t in L['tiers']], 'c', round(c, 2))
