# Generic Blender renderer for figure schematics (figures_v2). Reads a layout JSON (mm, y from the block top) and builds
# every element as a scene object: asset image planes (pd_product workflow_assets renders), shaded spheres and rounded
# cubes (Principled), flat fills, lines, arrows, small network glyphs, locks; a shadow catcher + sun give soft shadows.
# Renders with Cycles (CPU), orthographic, transparent film, exact block size at the layout dpi. Text is never drawn
# here: labels in the same JSON are placed as 7-pt vector text by R. Usage: blender.sh -b --python blender_render.py -- layout.json out.png
import bpy, json, math, sys
args = sys.argv[sys.argv.index('--') + 1:]; L = json.load(open(args[0])); OUT = args[1]
W, H = L['W'], L['H']; ASSETS = L.get('assets', '/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design/workflow_assets')
bpy.ops.wm.read_factory_settings(use_empty=True); sc = bpy.context.scene
def lin(hx, a=1.0):
    rgb = [int(hx[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    return (*[v / 12.92 if v <= .04045 else ((v + .055) / 1.055) ** 2.4 for v in rgb], a)
def mix_white(hx, a):
    rgb = [int(hx[i:i + 2], 16) for i in (1, 3, 5)]; return '#%02X%02X%02X' % tuple(round(255 - (255 - c) * a) for c in rgb)
MATS = {}
def emis(hx):
    k = 'E' + hx
    if k in MATS: return MATS[k]
    m = bpy.data.materials.new(k); m.use_nodes = True; n = m.node_tree.nodes; n.clear()
    e = n.new('ShaderNodeEmission'); e.inputs['Color'].default_value = lin(hx); o = n.new('ShaderNodeOutputMaterial')
    m.node_tree.links.new(e.outputs[0], o.inputs['Surface']); MATS[k] = m; return m
def shaded(hx, rough=.45):
    k = 'P' + hx
    if k in MATS: return MATS[k]
    m = bpy.data.materials.new(k); m.use_nodes = True; b = m.node_tree.nodes['Principled BSDF']
    b.inputs['Base Color'].default_value = lin(hx); b.inputs['Roughness'].default_value = rough
    if 'Specular IOR Level' in b.inputs: b.inputs['Specular IOR Level'].default_value = .35
    MATS[k] = m; return m
def P(x, y, z=0.0): return (x, H - y, z)
def link(ob): sc.collection.objects.link(ob); return ob
def poly(pts, z, hx, name='poly'):
    me = bpy.data.meshes.new(name); me.from_pydata([P(x, y, z) for x, y in pts], [], [list(range(len(pts)))]); me.update()
    ob = link(bpy.data.objects.new(name, me)); ob.data.materials.append(emis(hx)); return ob
def line(pts, z, hx, lw, closed=False, name='line'):
    cu = bpy.data.curves.new(name, 'CURVE'); cu.dimensions = '3D'; cu.bevel_depth = lw / 2; cu.bevel_resolution = 2
    sp = cu.splines.new('POLY'); sp.points.add(len(pts) - 1)
    for p, (x, y) in zip(sp.points, pts): p.co = (*P(x, y, z), 1)
    sp.use_cyclic_u = closed; ob = link(bpy.data.objects.new(name, cu)); ob.data.materials.append(emis(hx)); return ob
def dashed(pts, z, hx, lw, dash=.9):
    for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
        seg = math.hypot(x2 - x1, y2 - y1); n = max(1, int(seg / dash))
        for k in range(n):
            t0, t1 = k / n, (k + .55) / n; line([(x1 + (x2 - x1) * t0, y1 + (y2 - y1) * t0), (x1 + (x2 - x1) * t1, y1 + (y2 - y1) * t1)], z, hx, lw)
def disc(cx, cy, r, z, hx, n=48): return poly([(cx + r * math.cos(2 * math.pi * k / n), cy + r * math.sin(2 * math.pi * k / n)) for k in range(n)], z, hx, 'disc')
def sphere(cx, cy, r, hx, z=None):
    bpy.ops.mesh.primitive_uv_sphere_add(segments=32, ring_count=16, radius=r, location=P(cx, cy, (z if z is not None else 1.0) + r))
    ob = bpy.context.active_object; bpy.ops.object.shade_smooth(); ob.data.materials.append(shaded(hx)); return ob
def cube(cx, cy, s, hx, z=0.6):
    bpy.ops.mesh.primitive_cube_add(size=s, location=P(cx, cy, z + s / 2)); ob = bpy.context.active_object
    md = ob.modifiers.new('bev', 'BEVEL'); md.width = s * .18; md.segments = 3; ob.data.materials.append(shaded(hx, .5)); return ob
def image(d, z=2.0):
    x, y, s = d['x'], d['y'], d['s']
    me = bpy.data.meshes.new(d['id']); me.from_pydata([P(x, y + s, z), P(x + s, y + s, z), P(x + s, y, z), P(x, y, z)], [], [(0, 1, 2, 3)]); me.uv_layers.new()
    for lp, uv in zip(me.uv_layers.active.data, [(0, 0), (1, 0), (1, 1), (0, 1)]): lp.uv = uv
    ob = link(bpy.data.objects.new('img_' + d['id'], me)); m = bpy.data.materials.new('M_' + d['id'] + str(len(MATS))); m.use_nodes = True; n = m.node_tree.nodes; n.clear()
    tex = n.new('ShaderNodeTexImage'); tex.image = bpy.data.images.load(ASSETS + '/' + d['id'] + '.png', check_existing=True); tex.interpolation = 'Cubic'
    em = n.new('ShaderNodeEmission'); tr = n.new('ShaderNodeBsdfTransparent'); mx = n.new('ShaderNodeMixShader'); o = n.new('ShaderNodeOutputMaterial'); lk = m.node_tree.links
    lk.new(tex.outputs['Color'], em.inputs['Color']); lk.new(tex.outputs['Alpha'], mx.inputs[0]); lk.new(tr.outputs[0], mx.inputs[1]); lk.new(em.outputs[0], mx.inputs[2]); lk.new(mx.outputs[0], o.inputs['Surface'])
    ob.data.materials.append(m); return ob
def arrow(d, z=1.0):
    (x1, y1), (x2, y2) = d['p'][-2], d['p'][-1]; hx = d.get('col', '#595959' if d.get('dashed') else '#1A1A1A'); lw = d.get('lw', .26)
    ux, uy = x2 - x1, y2 - y1; Lh = math.hypot(ux, uy); ux, uy = ux / Lh, uy / Lh; hl, hw = d.get('hl', 1.0), d.get('hw', .5)
    def head(xe, ye, vx, vy): poly([(xe, ye), (xe - hl * vx - hw * vy, ye - hl * vy + hw * vx), (xe - hl * vx + hw * vy, ye - hl * vy - hw * vx)], z, hx, 'head')
    pts = [list(p) for p in d['p']]; pts[-1] = [x2 - hl * ux, y2 - hl * uy]
    if d.get('both'):
        (a1, b1), (a2, b2) = d['p'][0], d['p'][1]; vx, vy = a2 - a1, b2 - b1; lv = math.hypot(vx, vy); vx, vy = vx / lv, vy / lv; pts[0] = [a1 + hl * vx, b1 + hl * vy]; head(a1, b1, -vx, -vy)
    (dashed if d.get('dashed') else (lambda p, z_, c, w: line(p, z_, c, w)))(pts, z, hx, lw)
    if not d.get('nohead'): head(x2, y2, ux, uy)
def net(d):
    x, y, w, h = d['x'], d['y'], d['w'], d['h']; cy = y + h / 2; hl, hr = (h * .55, h) if d.get('widen') else (h, h * .55)
    pts = [(x, cy - hl / 2), (x + w, cy - hr / 2), (x + w, cy + hr / 2), (x, cy + hl / 2)]; poly(pts, 1.5, d['fill']); line(pts, 1.52, d['line'], .2, closed=True)
    nn = [2, 3, 4] if d.get('widen') else [4, 3, 2]; xs = [x + w * f for f in (.2, .5, .8)]; cols = []
    for k, m in enumerate(nn):
        span = (hl + (hr - hl) * (k + .5) / 3) * .66; cols.append([(xs[k], cy - span / 2 + span * i / (m - 1)) for i in range(m)])
    for k in range(2):
        for p in cols[k]:
            for q in cols[k + 1]: line([p, q], 1.55, mix_white(d['line'], .35), .07)
    for c in cols:
        for p in c: disc(p[0], p[1], .38, 1.6, d['line'])
def lock(d):
    lx, ly, s = d['x'], d['y'], d.get('s', 1.6)
    poly([(lx - s / 2, ly), (lx + s / 2, ly), (lx + s / 2, ly + s * .75), (lx - s / 2, ly + s * .75)], 2.2, '#1A1A1A', 'lock')
    line([(lx - s * .31 + i * s * .62 / 12, ly - s * .34 * math.sin(math.pi * i / 12)) for i in range(13)], 2.2, '#1A1A1A', s * .12)
def matrix(d):
    cs = d['s'] / d['n']; vals = d.get('vals', [.25, .8, .5, 1, .35, .6, .2, .9, .45, .7])
    for i in range(d['n']):
        for j in range(d.get('m', d['n'])):
            a = vals[(i * 3 + j * 7) % len(vals)]; x0, y0 = d['x'] + j * cs, d['y'] + i * cs
            poly([(x0 + .08, y0 + .08), (x0 + cs - .08, y0 + .08), (x0 + cs - .08, y0 + cs - .08), (x0 + .08, y0 + cs - .08)], 1.5, mix_white(d['col'], a), 'cell')
def rect(d):
    x, y, w, h, r = d['x'], d['y'], d['w'], d['h'], d.get('r', 0); n = 6; pts = []
    for cx, cy, a0 in [(x + w - r, y + r, -90), (x + w - r, y + h - r, 0), (x + r, y + h - r, 90), (x + r, y + r, 180)]:
        for k in range(n + 1): t = math.radians(a0 + 90 * k / n); pts.append((cx + r * math.cos(t), cy + r * math.sin(t)))
    poly(pts, d.get('z', 1.2), d['fill'])
    if d.get('line'): line(pts, d.get('z', 1.2) + .02, d['line'], d.get('lw', .18), closed=True)
for e in L['elements']:
    t = e['type']
    if t == 'image': image(e, e.get('z', 2.0))
    elif t == 'sphere': sphere(e['x'], e['y'], e['r'], e['col'], e.get('z'))
    elif t == 'cube': cube(e['x'], e['y'], e['s'], e['col'], e.get('z', .6))
    elif t == 'disc': disc(e['x'], e['y'], e['r'], e.get('z', 1.2), e['col'])
    elif t == 'line': (dashed if e.get('dashed') else (lambda p, z_, c, w: line(p, z_, c, w)))(e['p'], e.get('z', 1.0), e['col'], e.get('lw', .22))
    elif t == 'poly': poly(e['p'], e.get('z', 1.2), e['col'])
    elif t == 'arrow': arrow(e, e.get('z', 1.0))
    elif t == 'net': net(e)
    elif t == 'lock': lock(e)
    elif t == 'matrix': matrix(e)
    elif t == 'rect': rect(e)
    else: raise ValueError('unknown element ' + t)
me = bpy.data.meshes.new('ground'); me.from_pydata([(-50, -50, -.2), (W + 50, -50, -.2), (W + 50, H + 50, -.2), (-50, H + 50, -.2)], [], [(0, 1, 2, 3)])
gr = link(bpy.data.objects.new('shadow_catcher', me)); gr.is_shadow_catcher = True; gr.hide_render = not L.get('shadows', False)
sun = bpy.data.lights.new('Sun', 'SUN'); sun.energy = L.get('sun', .28); sun.angle = math.radians(L.get('sun_angle', 7)); so = link(bpy.data.objects.new('Sun', sun)); so.rotation_euler = (math.radians(-20), math.radians(16), 0)
wd = bpy.data.worlds.new('W'); wd.use_nodes = True; bg = wd.node_tree.nodes['Background']; bg.inputs['Color'].default_value = (1, 1, 1, 1); bg.inputs['Strength'].default_value = L.get('world', .9); sc.world = wd
cd = bpy.data.cameras.new('Cam'); cam = link(bpy.data.objects.new('Cam', cd)); cam.location = (W / 2, H / 2, 200); cd.type = 'ORTHO'; cd.ortho_scale = max(W, H); cd.clip_end = 1000; sc.camera = cam
sc.render.engine = 'CYCLES'; sc.cycles.device = 'CPU'; sc.cycles.samples = L.get('samples', 64); sc.cycles.use_denoising = True; sc.cycles.max_bounces = 4
sc.render.resolution_x = round(W / 25.4 * L.get('dpi', 600)); sc.render.resolution_y = round(H / 25.4 * L.get('dpi', 600)); sc.render.resolution_percentage = 100
sc.render.film_transparent = True; sc.render.image_settings.file_format = 'PNG'; sc.render.image_settings.color_mode = 'RGBA'
sc.view_settings.view_transform = 'Standard'; sc.view_settings.look = 'None'; sc.render.filepath = OUT; sc.render.threads_mode = 'FIXED'; sc.render.threads = 8
bpy.ops.render.render(write_still=True); bpy.ops.wm.save_as_mainfile(filepath=OUT.replace('.png', '.blend'), compress=True)
print('RENDER_DONE', OUT, sc.render.resolution_x, sc.render.resolution_y, len(L['elements']))
