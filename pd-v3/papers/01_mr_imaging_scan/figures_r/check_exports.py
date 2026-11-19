"""Check exact plotting units, native page dimensions, embedded fonts and text."""
import hashlib
import json
from pathlib import Path
import re
import sys
import pandas as pd
from PIL import Image


def main():
    folder = Path(sys.argv[1]).resolve()
    data = folder / 'plot_data'
    manifest = json.loads((data / 'MANIFEST.json').read_text())
    for name, expected in manifest['outputs'].items():
        assert hashlib.sha256((data / name).read_bytes()).hexdigest() == expected, name
    counts = {'A_landscape': 68175, 'B_matrix': 1170, 'C_black_substantia_nigra': 12,
              'D_tfrc_regions': 32, 'E_cross_modality': 24240, 'F_bilateral': 30300,
              'G_delta_se': 68175, 'H_palindromes': 269, 'S1_historical': 135}
    for name, n in counts.items():
        assert len(pd.read_csv(data / f'{name}.csv')) == n, name
    layout = json.loads((folder / 'qa/layout_checks.json').read_text())
    assert layout['panels'] == 8 and layout['plotting_areas'] == 22
    text = pd.concat([pd.read_csv(p) for p in (folder / 'qa').glob('text_*.csv')])
    assert text.pt.min() >= 6 - 1e-6
    exports = {}
    for name, mm in [('Figure1_expanded_MR', (174, 225)), ('historical_proteins_complete_panel', (174, 200))]:
        pdf = folder / f'{name}.pdf'
        b = pdf.read_bytes()
        font_names = sorted({v.decode() for v in re.findall(rb'/BaseFont\s*/([^\s/>]+)', b)})
        boxes = re.findall(rb'/MediaBox\s*\[([^]]+)\]', b)
        box = [float(x) for x in boxes[0].split()]
        page_mm = [(box[2]-box[0])*25.4/72, (box[3]-box[1])*25.4/72]
        assert all(abs(a-c) < .01 for a, c in zip(page_mm, mm))
        images = len(re.findall(rb'/Subtype\s*/Image', b))
        assert images == 0 and b.count(b'/FontFile2') >= 2 and all('Arial' in n for n in font_names)
        dimensions = Image.open(folder / f'{name}.png').size
        assert all(abs(actual - size/25.4*600) <= 1 for actual, size in zip(dimensions, mm))
        exports[name] = {'pdf_sha256': hashlib.sha256(b).hexdigest(), 'page_mm': page_mm,
                         'png_dimensions': dimensions, 'png_dpi': 600, 'fonts': font_names,
                         'embedded_true_type_streams': b.count(b'/FontFile2'), 'image_xobjects': images}
    report = {'passed': True, 'plotted_rows': counts, 'minimum_text_pt': float(text.pt.min()),
              'geometry': layout, 'exports': exports,
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'visual_review': 'Parent agent inspected full-page preview, individual A-H PNGs, and full 135-estimate supplement; corrected clipping, ticks, square scales, and legends'}
    (folder / 'qa/export_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS: plotted-unit counts, input hashes, layout, fonts, vector PDFs and 600-dpi exports')


if __name__ == '__main__':
    main()
