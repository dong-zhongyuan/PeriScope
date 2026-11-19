"""check_panel.py - rule checker for PeriScope figures_v2 panels and pages.

Panel mode (default):  check_panel.py PANEL.pdf [--png PANEL.png] [--manifest PANEL.manifest.json] [--json REPORT.json]
Page mode:             check_panel.py PAGE.pdf --page [--png PAGE.png] [--json REPORT.json]

FAIL when: the PDF page is not an allowed block (panel) or not 174 mm x n rows (page); any font is not Arial or not
embedded; any text is not 7 pt (page mode also allows the 9-pt bold single-letter panel letters); any label looks like a
title or sentence (more than five words, sentence punctuation, caption prefix); the manifest records a FAIL in its grob
text audit or a bar chart with 1-3 bars; the PNG is not 600 dpi at the block size; the panel's letter zone (top-left
4.0 x 3.7 mm) is not empty. Exit status 1 on any FAIL. The report JSON is written next to the PDF unless --json is given.
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figv2_pdf as fp  # noqa: E402

PT, PT_LETTER = 7.0, 9.0
LETTER_ZONE = (4.0, 3.7)
SENTENCE = re.compile(r'[a-z]{2,}[.!?]( |$)|;')
CAPTION = re.compile(r':$|^(Note|Notes|n =|Source)\b')


def words(s):
    return sum(1 for t in s.split() if re.search(r'[A-Za-z].*[A-Za-z]', t))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pdf')
    ap.add_argument('--png')
    ap.add_argument('--manifest')
    ap.add_argument('--page', action='store_true')
    ap.add_argument('--manifests', help='page mode: comma-separated manifests of the placed panels')
    ap.add_argument('--json')
    a = ap.parse_args()
    fails, notes = [], []
    pages = fp.inspect(a.pdf)
    if len(pages) != 1:
        fails.append(f'{len(pages)} pages (one expected)')
    pg = pages[0]
    w, h = pg['w_mm'], pg['h_mm']
    if a.page:
        rows = [fp.span_mm(n) for n in range(1, fp.ROWS + 1)]
        if abs(w - fp.span_mm(fp.COLS)) > 0.02 or not any(abs(h - r) < 0.02 for r in rows):
            fails.append(f'page {w:.2f} x {h:.2f} mm is not 174 mm wide and an integer number of unit rows high')
        block = None
    else:
        block = fp.block_of_size(w, h)
        if block is None:
            fails.append(f'page {w:.3f} x {h:.3f} mm is not an allowed block')
    for f in pg['fonts']:
        if not f['family'].startswith('Arial'):
            fails.append(f'font {f["name"]} is not Arial')
        if not f['embedded']:
            fails.append(f'font {f["name"]} is not embedded')
    for lab in pg['labels']:
        t, s = lab['text'], lab['size']
        letter = a.page and abs(s - PT_LETTER) < 0.02 and re.fullmatch(r'[a-z]( +[a-z])*', t.strip()) and all('Bold' in x for x in lab['fonts'])
        if abs(s - PT) > 0.02 and not letter:
            fails.append(f'text "{t}" is {s:.2f} pt')
        if words(t) > 5:   # pypdf joins labels that share a baseline, so this is a note; the grob audit below is binding
            notes.append(f'PDF line "{t}" has {words(t)} words (may join neighbouring labels; see grob audit)')
        if SENTENCE.search(t) or CAPTION.search(t):
            fails.append(f'text "{t}" reads as a sentence or caption')
    man = None
    manifests = [a.manifest] if a.manifest else []
    if a.page and a.manifests:
        manifests = [m for m in a.manifests.split(',') if m]
    for mf in manifests:
        man = json.load(open(mf))
        if not a.page and man.get('block') != block:
            fails.append(f'manifest block {man.get("block")} differs from PDF block {block}')
        for row in man.get('text', []):
            if str(row.get('flag', '')).startswith('FAIL'):
                fails.append(f'grob text "{row["label"]}": {row["flag"]}')
        for b in man.get('bars', []):
            if 1 <= int(b['n_bars']) <= 3:
                fails.append(f'bar chart with {b["n_bars"]} bars ({b["what"]})')
        notes.append(f'{man.get("panel")} {man.get("block")}: {len(man.get("sources", []))} source file(s), {len(man.get("bars", []))} bar chart(s) recorded')
    if a.png:
        from PIL import Image
        im = Image.open(a.png).convert('L')
        ew, eh = round(w / 25.4 * 600), round(h / 25.4 * 600)
        if abs(im.size[0] - ew) > 1 or abs(im.size[1] - eh) > 1:
            fails.append(f'PNG {im.size} px is not 600 dpi at {w:.2f} x {h:.2f} mm ({ew} x {eh})')
        if not a.page:
            zx, zy = round(LETTER_ZONE[0] / 25.4 * 600), round(LETTER_ZONE[1] / 25.4 * 600)
            if im.crop((0, 0, zx, zy)).getextrema()[0] < 245:
                fails.append('letter zone (top-left 4.0 x 3.7 mm) is not empty')
    rep = {'pdf': os.path.abspath(a.pdf), 'mode': 'page' if a.page else 'panel', 'block': block, 'size_mm': [round(w, 3), round(h, 3)],
           'fonts': sorted({f['family'] for f in pg['fonts']}), 'n_labels': len(pg['labels']),
           'sizes_pt': sorted({lab['size'] for lab in pg['labels']}), 'images': pg['images'],
           'passed': not fails, 'fails': fails, 'notes': notes}
    out = a.json or re.sub(r'\.pdf$', '', a.pdf) + '.check.json'
    with open(out, 'w') as f:
        json.dump(rep, f, indent=1, ensure_ascii=False)
    status = 'PASS' if not fails else 'FAIL'
    print(f'{status} {os.path.basename(a.pdf)}: {rep["mode"]} {block or ""} {w:.2f} x {h:.2f} mm, fonts {rep["fonts"]}, '
          f'sizes {rep["sizes_pt"]}, {rep["n_labels"]} labels, {pg["images"]} raster image(s)')
    for x in fails:
        print('  FAIL', x)
    for x in notes:
        print('  note', x)
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
