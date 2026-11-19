"""figv2_pdf.py - PDF tools for the PeriScope figures_v2 system (crop to an exact block, inspect fonts and text, place
panels at 1:1). Used by render_panel.R (crop), check_panel.py (inspect) and assemble.py (merge). pypdf only."""
import math
import re
import sys

import pypdf
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject

MM = 72 / 25.4
U_MM, G_MM, COLS, ROWS = 41.25, 3.0, 4, 5
BLOCKS = ['1x1', '1x2', '1x3', '1x4', '2x2', '2x3', '2x4', '2x1', '3x1', '4x1', '3x2', '4x2']   # rows x cols


def span_mm(n):
    return n * U_MM + (n - 1) * G_MM


def block_size(b):
    r, c = (int(v) for v in b.split('x'))
    return span_mm(c), span_mm(r)


def block_of_size(w, h, tol=0.02):
    for b in BLOCKS:
        bw, bh = block_size(b)
        if abs(bw - w) < tol and abs(bh - h) < tol:
            return b
    return None


def crop(raw, final, w_mm, h_mm):
    """Set the page to exactly w x h mm, keeping the top-left corner (cairo truncates page sizes to whole points)."""
    w, h = float(w_mm) * MM, float(h_mm) * MM
    r = PdfReader(raw)
    assert len(r.pages) == 1, 'one page per panel'
    wr = PdfWriter()
    p = wr.add_page(r.pages[0])
    mb = p.mediabox
    dev_w, dev_h = float(mb.width), float(mb.height)
    assert dev_w + 1e-6 >= w and dev_h + 1e-6 >= h, f'device {dev_w}x{dev_h} pt smaller than block {w:.3f}x{h:.3f} pt'
    dy = float(mb.top) - h
    p.add_transformation(Transformation().translate(-float(mb.left), -dy))
    box = RectangleObject([0, 0, w, h])
    p.mediabox = box
    p.cropbox = box
    p.trimbox = box
    wr.add_metadata({'/Producer': 'PeriScope figures_v2 (R cairo_pdf + pypdf crop)'})
    with open(final, 'wb') as f:
        wr.write(f)


def _fonts(page):
    out = []
    res = page.get('/Resources')
    if res is None:
        return out
    res = res.get_object()
    seen = set()

    def walk(resources):
        resources = resources.get_object()
        fd = resources.get('/Font')
        if fd is not None:
            for _, ref in fd.get_object().items():
                f = ref.get_object()
                key = id(f)
                if key in seen:
                    continue
                seen.add(key)
                name = str(f.get('/BaseFont', ''))
                desc = f.get('/FontDescriptor')
                if desc is None and f.get('/DescendantFonts') is not None:
                    desc = f['/DescendantFonts'].get_object()[0].get_object().get('/FontDescriptor')
                desc = desc.get_object() if desc is not None else {}
                emb = any(k in desc for k in ('/FontFile', '/FontFile2', '/FontFile3'))
                out.append({'name': name.lstrip('/'), 'family': re.sub(r'^[A-Z]{6}\+', '', name.lstrip('/')), 'embedded': emb,
                            'subtype': str(f.get('/Subtype', ''))})
        xo = resources.get('/XObject')
        if xo is not None:
            for _, ref in xo.get_object().items():
                x = ref.get_object()
                if x.get('/Subtype') == '/Form' and x.get('/Resources') is not None:
                    walk(x['/Resources'])

    walk(res)
    return out


def _images(page):
    n = 0
    res = page.get('/Resources')
    if res is None:
        return 0

    def walk(resources):
        nonlocal n
        xo = resources.get_object().get('/XObject')
        if xo is None:
            return
        for _, ref in xo.get_object().items():
            x = ref.get_object()
            if x.get('/Subtype') == '/Image':
                n += 1
            elif x.get('/Subtype') == '/Form' and x.get('/Resources') is not None:
                walk(x['/Resources'])

    walk(res)
    return n


def inspect(path):
    """Page size (mm), fonts, image count, and every text label with its effective size (pt), font and origin."""
    r = PdfReader(path)
    pages = []
    for p in r.pages:
        mb = p.mediabox
        spans = []

        def visit(text, cm, tm, font_dict, font_size):
            if not text or not text.strip():
                return
            det_t = abs(tm[0] * tm[3] - tm[1] * tm[2])
            det_c = abs(cm[0] * cm[3] - cm[1] * cm[2])
            size = font_size * math.sqrt(det_t) * math.sqrt(det_c)
            x = tm[4] * cm[0] + tm[5] * cm[2] + cm[4]
            y = tm[4] * cm[1] + tm[5] * cm[3] + cm[5]
            name = str(font_dict.get('/BaseFont', '')) if font_dict else ''
            spans.append({'text': text.replace('\n', ''), 'size': round(size, 3), 'font': re.sub(r'^/?[A-Z]{6}\+', '', name.lstrip('/')),
                          'x': round(x, 2), 'y': round(y, 2)})

        p.extract_text(visitor_text=visit)
        labels = {}
        for s in spans:
            k = (s['x'], s['y'], s['size'])
            if k not in labels:
                labels[k] = {'text': '', 'size': s['size'], 'fonts': set(), 'x_mm': s['x'] / MM, 'y_mm': s['y'] / MM}
            labels[k]['text'] += s['text']
            labels[k]['fonts'].add(s['font'])
        lab = [dict(v, fonts=sorted(v['fonts'])) for v in labels.values()]
        pages.append({'w_mm': float(mb.width) / MM, 'h_mm': float(mb.height) / MM, 'origin': [float(mb.left), float(mb.bottom)],
                      'fonts': _fonts(p), 'images': _images(p), 'labels': lab})
    return pages


def merge(page_w_mm, page_h_mm, placements, out):
    """placements: list of (pdf_path, x_mm, y_mm) with x, y = top-left of the block on the page (mm from page top-left).
    Each panel page is translated only (no scaling); its MediaBox must already be the exact block."""
    W, H = page_w_mm * MM, page_h_mm * MM
    wr = PdfWriter()
    page = wr.add_blank_page(width=W, height=H)
    for path, x_mm, y_mm in placements:
        p = PdfReader(path).pages[0]
        mb = p.mediabox
        ph = float(mb.height)
        tx = x_mm * MM - float(mb.left)
        ty = H - y_mm * MM - ph - float(mb.bottom)
        page.merge_transformed_page(p, Transformation().translate(tx, ty))
    box = RectangleObject([0, 0, W, H])
    page.mediabox = box
    page.cropbox = box
    wr.add_metadata({'/Producer': 'PeriScope figures_v2 assembler (panels placed 1:1, pypdf)'})
    with open(out, 'wb') as f:
        wr.write(f)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'crop':
        crop(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif cmd == 'inspect':
        import json
        print(json.dumps(inspect(sys.argv[2]), indent=1, ensure_ascii=False))
    else:
        raise SystemExit('usage: figv2_pdf.py crop RAW FINAL W_MM H_MM | inspect PDF')
