"""assemble.py - build one figure page from a layout JSON by placing panel PDFs at 1:1 on the 4 x 5 unit grid.

Usage (after `. _system/env.sh`):  $FIGV2_PY _system/assemble.py _system/layout_fig1.json [--out-dir DIR]

Layout JSON: {"figure": 1, "letter_on_blank": false, "panels": [{"letter": "a", "id": "1a", "block": "1x3", "row": 1, "col": 1,
"file": "fig1/panels/1a_1x3.pdf", "status": "ready" | "blank-pending-ablation", ...}]}. block = rows x cols; row, col are
1-based grid cells (top-left of the block). Blocks must be allowed sizes, inside the 4 x 5 grid and must not overlap.
A panel whose status starts with "blank", or whose PDF does not exist yet, is drawn as an outlined blank block (no text;
no letter unless letter_on_blank is true). A panel PDF whose page is not exactly its block is refused (never scaled).
Writes Figure<n>.pdf (vector), Figure<n>.png (600 dpi), Figure<n>_preview.png (150 dpi), Figure<n>.assembly.json and the
page check (check_panel.py --page).
"""
import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figv2_pdf as fp  # noqa: E402

SYS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SYS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('layout')
    ap.add_argument('--out-dir')
    a = ap.parse_args()
    L = json.load(open(a.layout))
    n = L['figure']
    out_dir = a.out_dir or os.path.join(ROOT, f'fig{n}', 'page')
    os.makedirs(out_dir, exist_ok=True)
    occ = {}
    rows_used = 0
    placements, letters, blanks, report = [], [], [], []
    for p in L['panels']:
        b = p['block']
        if b not in fp.BLOCKS:
            raise SystemExit(f'panel {p["letter"]}: block {b} not allowed')
        r, c = (int(v) for v in b.split('x'))
        r0, c0 = int(p['row']), int(p['col'])
        if r0 < 1 or c0 < 1 or r0 + r - 1 > fp.ROWS or c0 + c - 1 > fp.COLS:
            raise SystemExit(f'panel {p["letter"]}: block {b} at row {r0}, col {c0} leaves the 4 x 5 grid')
        for i in range(r0, r0 + r):
            for j in range(c0, c0 + c):
                if (i, j) in occ:
                    raise SystemExit(f'panel {p["letter"]} overlaps panel {occ[(i, j)]} at cell ({i}, {j})')
                occ[(i, j)] = p['letter']
        rows_used = max(rows_used, r0 + r - 1)
        x, y = (c0 - 1) * (fp.U_MM + fp.G_MM), (r0 - 1) * (fp.U_MM + fp.G_MM)
        w, h = fp.block_size(b)
        f = os.path.join(ROOT, p['file']) if p.get('file') else None
        blank = str(p.get('status', '')).startswith('blank')
        state = 'blank (pending)' if blank else 'placed'
        if not blank and (f is None or not os.path.exists(f)):
            state = 'blank (panel not rendered yet)'
        if state == 'placed':
            pg = fp.inspect(f)[0]
            if abs(pg['w_mm'] - w) > 0.02 or abs(pg['h_mm'] - h) > 0.02:
                raise SystemExit(f'panel {p["letter"]}: {f} is {pg["w_mm"]:.2f} x {pg["h_mm"]:.2f} mm, block {b} is {w:.2f} x {h:.2f} mm')
            placements.append((f, x, y))
        else:
            blanks.append({'x': x, 'y': y, 'w': w, 'h': h})
        if state == 'placed' or L.get('letter_on_blank', False) or state.startswith('blank (panel not'):
            letters.append({'letter': p['letter'], 'x': x, 'y': y})
        report.append({'letter': p['letter'], 'id': p.get('id'), 'block': b, 'row': r0, 'col': c0, 'x_mm': x, 'y_mm': y,
                       'w_mm': w, 'h_mm': h, 'state': state, 'file': f})
    free = [(i, j) for i in range(1, rows_used + 1) for j in range(1, fp.COLS + 1) if (i, j) not in occ]
    W, H = fp.span_mm(fp.COLS), fp.span_mm(rows_used)
    spec = os.path.join(out_dir, f'.overlay_fig{n}.json')
    json.dump({'page_w': W, 'page_h': H, 'letters': letters, 'blanks': blanks}, open(spec, 'w'))
    ov = os.path.join(out_dir, f'.overlay_fig{n}.pdf')
    subprocess.run([os.environ['FIGV2_RSCRIPT'], os.path.join(SYS, 'overlay.R'), spec, ov], check=True)
    out = os.path.join(out_dir, f'Figure{n}.pdf')
    fp.merge(W, H, placements + [(ov, 0, 0)], out)
    os.remove(spec)
    os.remove(ov)
    for res, suffix in ((600, ''), (150, '_preview')):
        subprocess.run(['pdftoppm', '-r', str(res), '-png', '-singlefile', out, os.path.join(out_dir, f'Figure{n}{suffix}')], check=True)
    mans = [m for m in (f[:-4] + '.manifest.json' for f, _, _ in placements) if os.path.exists(m)]
    chk = subprocess.run([sys.executable, os.path.join(SYS, 'check_panel.py'), out, '--page', '--png', os.path.join(out_dir, f'Figure{n}.png'),
                          '--manifests', ','.join(mans)], capture_output=True, text=True)
    print(chk.stdout, end='')
    json.dump({'figure': n, 'page_mm': [W, H], 'unit_rows': rows_used, 'empty_cells': free, 'panels': report,
               'page_check_passed': chk.returncode == 0}, open(os.path.join(out_dir, f'Figure{n}.assembly.json'), 'w'), indent=1)
    print(f'Figure {n}: {len(placements)} panel(s) placed, {len(blanks)} blank, {len(free)} empty cell(s); {W:.2f} x {H:.2f} mm -> {out}')
    sys.exit(chk.returncode)


if __name__ == '__main__':
    main()
