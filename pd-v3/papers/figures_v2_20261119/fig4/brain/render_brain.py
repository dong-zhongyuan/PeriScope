#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Brain sections for the two final targets, rendered on the 1 mm MNI152 template and written as
raw RGBA (+ a dims JSON) for the R panels to place with rasterGrob, exactly as the Blender art is.

CD71/TFRC acts on deep grey matter: the substantia nigra (Talairach Brodmann atlas) and the pallidum
(Harvard-Oxford subcortical), the structures carrying its MR estimates with colocalisation support.
CD119/IFNGR1's single imaging hit is one edge of the UK Biobank resting-state network matrix,
"rfMRI connectivity (ICA100 edge 868)" = BIG40 phenotype 3296 = UKB field 25753. UKB's own
edge_list_d100.txt gives edge 868 = good nodes 7 and 43, and rfMRI_GoodComponents_d100_v1.txt maps
those to group-ICA components 8 and 45 of rfMRI_ICA_d100.nii.gz.

Every section is cropped to the same field of view and the same aspect, so the four images are
directly comparable. Neurological convention: the left hemisphere is on the left; on the sagittal
section, anterior is on the right.
"""
import argparse
import json
import os
import warnings

warnings.filterwarnings('ignore')
import numpy as np
import nibabel as nib
from PIL import Image
from scipy import ndimage
from nilearn import datasets, image

DD = '/public/home/mengxl/dzy/.nilearn'
ICA = '/public/home/mengxl/dzy/pd-v3/data/ukb_group_ica/rfMRI_ICA_d100.nii.gz'
CD71, CD119 = '#2E8B57', '#7B3F9E'
FOV = 165.0               # one field of view for all four sections, so they share a scale
ASPECT = 148.0 / 116.0          # every section is cropped to this width:height ratio
DPI = 600
GREY_LO, GREY_HI = 0.22, 0.93

# stem, plane, coordinate (MNI mm), crop centre along the two in-plane axes, colour, what
SECTIONS = [
    ('sn',  'z', -9.0,  (0.0, -24.0), FOV, CD71,  'substantia nigra'),
    ('gp',  'z', -3.0,  (0.0, -24.0), FOV, CD71,  'pallidum'),
    ('dmn', 'x',  0.0,  (-16.0, 11.0), FOV, CD119, 'ICA d100 component 8'),
    ('amt', 'y', -12.0, (0.0,   6.0), FOV, CD119, 'ICA d100 component 45'),
]


def hex_rgb(h):
    h = h.lstrip('#')
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], float) / 255.0


def cut(vol, aff, plane, coord, centre, fov):
    """One section of a 1 mm volume on the template grid, cropped to FOV x FOV/ASPECT,
    returned as a 2-D array already oriented for display."""
    ax = 'xyz'.index(plane)
    k = int(round((coord - aff[ax, 3]) / aff[ax, ax]))
    sl = [slice(None)] * 3
    sl[ax] = k
    s = vol[tuple(sl)]                       # axes, in order, of the two remaining dims
    rest = [i for i in range(3) if i != ax]
    w, h = fov, fov / ASPECT
    # which remaining axis runs horizontally, and in which direction
    if plane == 'z':      # columns: x ascending (left hemisphere on the left); rows: y descending
        hax, vax, hflip, span = 0, 1, False, (w, h)
    elif plane == 'x':    # columns: y ascending (anterior to the right); rows: z descending
        hax, vax, hflip, span = 1, 2, False, (w, h)
    else:                 # coronal: columns x ascending; rows z descending
        hax, vax, hflip, span = 0, 2, False, (w, h)
    img = s if rest.index(hax) == 1 else s.T   # -> (vertical, horizontal)
    img = np.flipud(img)                        # vertical axis descends
    o_h, o_v = aff[hax, 3], aff[vax, 3]
    n_h, n_v = vol.shape[hax], vol.shape[vax]
    lo_h = int(round((centre[0] - span[0] / 2 - o_h) / aff[hax, hax]))
    hi_h = lo_h + int(round(span[0]))
    top_v = int(round((centre[1] + span[1] / 2 - o_v) / aff[vax, vax]))
    lo_r = n_v - 1 - top_v
    hi_r = lo_r + int(round(span[1]))
    pad = ((max(0, -lo_r), max(0, hi_r - img.shape[0])), (max(0, -lo_h), max(0, hi_h - img.shape[1])))
    if any(sum(p) for p in pad):
        img = np.pad(img, pad, constant_values=0)
        lo_r += pad[0][0]; hi_r += pad[0][0]; lo_h += pad[1][0]; hi_h += pad[1][0]
    out = img[lo_r:hi_r, lo_h:hi_h]
    return out[:, ::-1] if hflip else out


def upsample(a, size, mode):
    im = Image.fromarray((np.clip(a, 0, 1) * 255).astype(np.uint8), 'L')
    return np.asarray(im.resize(size, mode), float) / 255.0


def compose(anat, masks, size):
    v = anat.copy()
    hi = np.percentile(v[v > 0], 99.5) if (v > 0).any() else 1.0
    v = np.clip(v / hi, 0, 1)
    head = v > 0.04
    big = upsample(v, size, Image.LANCZOS)
    head_up = upsample(head.astype(float), size, Image.BILINEAR) > 0.5
    g = GREY_LO + (GREY_HI - GREY_LO) * np.clip(big, 0, 1)
    rgb = np.repeat(g[:, :, None], 3, axis=2)
    rgb[~head_up] = 1.0
    for m, col in masks:
        a = upsample(m.astype(float), size, Image.BILINEAR)
        c = hex_rgb(col)
        rgb = rgb * (1 - a[:, :, None]) + c[None, None, :] * a[:, :, None]
        core = ndimage.binary_erosion(a > 0.5, iterations=4)
        rim = (a > 0.5) & ~core
        rgb[rim] = c * 0.62
    out = np.concatenate([rgb, np.ones(rgb.shape[:2] + (1,))], axis=2)
    return (np.clip(out, 0, 1) * 255).astype(np.uint8)


def write(arr, stem, out):
    im = Image.fromarray(arr, 'RGBA')
    open(os.path.join(out, stem + '.rgba'), 'wb').write(im.tobytes())
    json.dump({'width': im.width, 'height': im.height}, open(os.path.join(out, stem + '.rgba.json'), 'w'))
    im.convert('RGB').save(os.path.join(out, stem + '.png'))
    print('%-10s %d x %d px' % (stem, im.width, im.height))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='.')
    ap.add_argument('--w-mm', type=float, default=41.10)
    ap.add_argument('--h-mm', type=float, default=32.20)
    ap.add_argument('--ica-z', type=float, default=8.0)   # both components are shown above this z
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    size = (int(round(a.w_mm / 25.4 * DPI)), int(round(a.h_mm / 25.4 * DPI)))

    tpl = datasets.load_mni152_template(resolution=1)
    aff, anat = tpl.affine, tpl.get_fdata()
    ho = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-1mm', data_dir=DD)
    ho_img = ho.maps if hasattr(ho.maps, 'get_fdata') else nib.load(ho.maps)
    tal = datasets.fetch_atlas_talairach('ba', data_dir=DD)
    tal_img = tal.maps if hasattr(tal.maps, 'get_fdata') else nib.load(tal.maps)
    rs = lambda im, i: image.resample_to_img(im, tpl, interpolation=i, force_resample=True, copy_header=True).get_fdata()
    hod, tald = rs(ho_img, 'nearest'), rs(tal_img, 'nearest')
    holab = [str(x) for x in ho.labels]
    tallab = [str(x) for x in tal.labels]
    X = (aff[0, 3] + np.arange(tpl.shape[0]))[:, None, None] * np.ones((1,) + tpl.shape[1:])

    ica = nib.load(ICA)
    comp = {c: rs(image.index_img(ica, c - 1), 'continuous') for c in (8, 45)}
    thr = {c: a.ica_z for c in comp}
    for c, v in comp.items():
        print('component %d  max %.1f  threshold %.1f  voxels %d (%d mm3)'
              % (c, v.max(), thr[c], int((v > thr[c]).sum()), int((v > thr[c]).sum())))

    MASKS = {
        'sn':  [((tald == tallab.index('Substania Nigra')) & (X < 0), CD71)],
        'gp':  [(hod == holab.index('Left Pallidum'), CD71), (hod == holab.index('Right Pallidum'), CD71)],
        'dmn': [(comp[8] > thr[8], CD119)],
        'amt': [(comp[45] > thr[45], CD119)],
    }
    meta = {'display_mm': [a.w_mm, a.h_mm], 'aspect': ASPECT, 'dpi': DPI,
            'ica_threshold_z': a.ica_z, 'sections': {}}
    for stem, plane, coord, centre, fov, col, what in SECTIONS:
        ms = [(cut(m.astype(float), aff, plane, coord, centre, fov), c) for m, c in MASKS[stem]]
        write(compose(cut(anat, aff, plane, coord, centre, fov), ms, size), stem, a.out)
        meta['sections'][stem] = {'plane': plane, 'coord': coord, 'centre': centre, 'fov_mm': fov, 'colour': col, 'what': what,
                                  'voxels': int(sum(m.sum() for m, _ in MASKS[stem]))}
    meta['sources'] = {'template': 'MNI152 1 mm (nilearn)',
                       'sn': 'Talairach Brodmann atlas, substantia nigra, left hemisphere',
                       'gp': 'Harvard-Oxford sub-maxprob-thr25-1mm, left and right pallidum',
                       'ica': ICA + ' (UK Biobank group-ICA d100); edge 868 = good nodes 7, 43 = components 8, 45'}
    json.dump(meta, open(os.path.join(a.out, 'brain.render.json'), 'w'), indent=1)
    print(json.dumps(meta['sections'], indent=1))


if __name__ == '__main__':
    main()
