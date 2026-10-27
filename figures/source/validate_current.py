"""Validate delivered files, native dimensions, source hashes and PPT placements."""
from pathlib import Path
from collections import Counter
import hashlib, json, zipfile, sys
import fitz
from PIL import Image
from lxml import etree
W = Path(__file__).resolve().parents[1]
ROOT = W.parents[2]
D = ROOT/'delivery'
manifest = json.loads((W/'panel_manifest.json').read_text())
layouts = json.loads((W/'layout.json').read_text())
assert len(manifest) == 49
report = {'panels': [], 'presentations': [], 'source_files': []}
for ident, row in sorted(manifest.items()):
    files = [D/'panels'/f'{ident}.{ext}' for ext in ['pdf', 'svg', 'png', 'tiff']]
    assert all(p.is_file() and p.stat().st_size for p in files), ident
    with fitz.open(files[0]) as doc:
        p = doc[0]
        assert abs(p.rect.width-row['width_in']*72)<.02, ident
        assert abs(p.rect.height-row['height_in']*72)<.02, ident
        fonts = Counter()
        outside = []
        for b in p.get_text('dict')['blocks']:
            for line in b.get('lines', []):
                for span in line['spans']:
                    if not span['text'].strip(): continue
                    fonts[round(span['size'], 1)] += 1
                    box = fitz.Rect(span['bbox'])
                    if not (p.rect+(-1,-1,1,1)).contains(box):
                        outside.append(span['text'])
        assert not outside, (ident, outside)
        if ident == '1d': assert '\nNA\n' not in p.get_text(), 'Invalid donor/stage facets'
        report['panels'].append(dict(id=ident, width_in=row['width_in'],
            height_in=row['height_in'], font_sizes=dict(fonts), outside_text=outside))
    for path in files[2:]:
        with Image.open(path) as im:
            assert im.size == (row['width_px'], row['height_px']), (ident, path.suffix)
            assert all(abs(x-600)<.1 for x in im.info['dpi']), ident
    svg = etree.parse(str(files[1]))
    for e in svg.iter():
        for key,value in e.attrib.items():
            if key.endswith('href'): assert value.startswith(('#','data:')), ident
ns = {'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
for layout in layouts:
    n = layout['figure']; ppt = D/'pptx'/f'Figure{n}.pptx'
    if '--skip-ppt' in sys.argv: continue
    with zipfile.ZipFile(ppt) as z:
        assert z.testzip() is None
        root = etree.fromstring(z.read('ppt/presentation.xml'))
        size = root.find('p:sldSz', ns)
        assert abs(int(size.get('cx'))/12700-layout['width'])<.1
        assert abs(int(size.get('cy'))/12700-layout['height'])<.1
        svg_count = len([s for s in z.namelist() if s.startswith('ppt/media/') and s.endswith('.svg')])
        assert svg_count == len(layout['all_panels'])-1, (n, svg_count)
        for panel in layout['all_panels']:
            ident = panel['id']
            if ident.endswith('a'): continue
            assert z.read(f'ppt/media/Figure{ident}.svg') == (D/'panels'/f'{ident}.svg').read_bytes(), ('Stale PPT panel', ident)
        slide = etree.fromstring(z.read('ppt/slides/slide1.xml'))
        native = len(slide.findall('.//p:sp', ns))
        assert native > 15, (n, native)
        report['presentations'].append(dict(figure=n, svg_panels=svg_count,
            native_shapes=native, sha256=hashlib.sha256(ppt.read_bytes()).hexdigest()))
    with fitz.open(D/'pdf'/f'Figure{n}.pdf') as doc:
        assert abs(doc[0].rect.width-layout['width'])<.02
        assert abs(doc[0].rect.height-layout['height'])<.02
    for panel in layout['all_panels']:
        r = manifest[panel['id']]
        assert abs(panel['width']-r['width_in']*72)<.02
        assert abs(panel['height']-r['height_in']*72)<.02
sources = json.loads((W/'data/source_manifest.json').read_text())
for name, expected in sources.items():
    path = ROOT/name
    assert path.exists(), name
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    assert actual == expected, ('Changed figure input', name)
    report['source_files'].append(dict(path=name, sha256=actual))
assert len(list((D/'panels').glob('*.pdf'))) == 49
assert len(list((D/'pptx').glob('Figure*.pptx'))) == 7
report['status'] = 'passed'
(W/'font_geometry_audit.json').write_text(json.dumps(report, indent=2))
print('PASS: 49 panels / 196 graphics files; '+('PPT inspection skipped' if '--skip-ppt' in sys.argv else '7 PPTs')+'; native placements; current source hashes.')
