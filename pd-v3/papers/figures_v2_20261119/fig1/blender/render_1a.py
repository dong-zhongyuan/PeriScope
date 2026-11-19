# Blender scene for Figure 1a (PeriScope architecture). Geometry, assets, arrows and shadows are built and rendered here
# (Cycles, orthographic, transparent film, exact block size at 600 dpi); text is added afterwards as 7-pt vector type in R
# (fig1/a_architecture.R) from the same layout file. Usage: blender.sh -b --python render_1a.py -- layout.json out.png
import bpy, json, math, sys
args = sys.argv[sys.argv.index('--') + 1:]; L = json.load(open(args[0])); OUT = args[1]
W, H = L['W'], L['H']
bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene
def lin(hx, a=1.0):
    rgb = [int(hx[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    return (*[v / 12.92 if v <= .04045 else ((v + .055) / 1.055) ** 2.4 for v in rgb], a)
MATS = {}
def emis(hx):
    if hx in MATS: return MATS[hx]
    m = bpy.data.materials.new('E' + hx); m.use_nodes = True; n = m.node_tree.nodes; n.clear()
    e = n.new('ShaderNodeEmission'); e.inputs['Color'].default_value = lin(hx); o = n.new('ShaderNodeOutputMaterial')
    m.node_tree.links.new(e.outputs[0], o.inputs['Surface']); MATS[hx] = m; return m
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
def disc(cx, cy, r, z, hx, n=48, name='disc'):
    return poly([(cx + r * math.cos(2 * math.pi * k / n), cy + r * math.sin(2 * math.pi * k / n)) for k in range(n)], z, hx, name)
def image(d, z=2.0):
    x, y, s = d['x'], d['y'], d['s']
    me = bpy.data.meshes.new(d['id']); me.from_pydata([P(x, y + s, z), P(x + s, y + s, z), P(x + s, y, z), P(x, y, z)], [], [(0, 1, 2, 3)]); me.uv_layers.new()
    for lp, uv in zip(me.uv_layers.active.data, [(0, 0), (1, 0), (1, 1), (0, 1)]): lp.uv = uv
    ob = link(bpy.data.objects.new('img_' + d['id'], me))
    m = bpy.data.materials.new('M_' + d['id']); m.use_nodes = True; n = m.node_tree.nodes; n.clear()
    tex = n.new('ShaderNodeTexImage'); tex.image = bpy.data.images.load(L['assets'] + '/' + d['id'] + '.png'); tex.interpolation = 'Cubic'
    em = n.new('ShaderNodeEmission'); tr = n.new('ShaderNodeBsdfTransparent'); mix = n.new('ShaderNodeMixShader'); o = n.new('ShaderNodeOutputMaterial')
    lk = m.node_tree.links; lk.new(tex.outputs['Color'], em.inputs['Color']); lk.new(tex.outputs['Alpha'], mix.inputs[0]); lk.new(tr.outputs[0], mix.inputs[1]); lk.new(em.outputs[0], mix.inputs[2]); lk.new(mix.outputs[0], o.inputs['Surface'])
    m.blend_method = 'BLEND' if hasattr(m, 'blend_method') else None; ob.data.materials.append(m); return ob
def mix_white(hx, a):
    rgb = [int(hx[i:i + 2], 16) for i in (1, 3, 5)]; return '#%02X%02X%02X' % tuple(round(255 - (255 - c) * a) for c in rgb)
def arrow(a, z=1.0, hx='#1A1A1A', lw=.26, dashed=False, both=False):
    (x1, y1), (x2, y2) = a; dx, dy = x2 - x1, y2 - y1; Lh = math.hypot(dx, dy); ux, uy = dx / Lh, dy / Lh; hl, hw = 1.0, .5
    def head(xe, ye, vx, vy):
        poly([(xe, ye), (xe - hl * vx - hw * vy, ye - hl * vy + hw * vx), (xe - hl * vx + hw * vy, ye - hl * vy - hw * vx)], z, hx, 'head')
    s0 = (x1 + (hl * ux if both else 0), y1 + (hl * uy if both else 0)); s1 = (x2 - hl * ux, y2 - hl * uy)
    if dashed:
        seg = math.hypot(s1[0] - s0[0], s1[1] - s0[1]); n = max(1, int(seg / .9))
        for k in range(n):
            t0, t1 = k / n, (k + .55) / n
            line([(s0[0] + (s1[0] - s0[0]) * t0, s0[1] + (s1[1] - s0[1]) * t0), (s0[0] + (s1[0] - s0[0]) * t1, s0[1] + (s1[1] - s0[1]) * t1)], z, hx, lw)
    else:
        line([s0, s1], z, hx, lw)
    head(x2, y2, ux, uy)
    if both: head(x1, y1, -ux, -uy)
for d in L['images']: image(d)
for d in L['nets']:
    x, y, w, h = d['x'], d['y'], d['w'], d['h']; cy = y + h / 2; hl, hr = (h * .55, h) if d['widen'] else (h, h * .55)
    pts = [(x, cy - hl / 2), (x + w, cy - hr / 2), (x + w, cy + hr / 2), (x, cy + hl / 2)]
    poly(pts, 1.5, d['fill'], d['id']); line(pts, 1.52, d['line'], .22, closed=True)
    nn = [2, 3, 4] if d['widen'] else [4, 3, 2]; xs = [x + w * f for f in (.2, .5, .8)]
    cols = []
    for k, m in enumerate(nn):
        span = (hl + (hr - hl) * (k + .5) / 3) * .66; cols.append([(xs[k], cy - span / 2 + span * i / (m - 1)) for i in range(m)])
    for k in range(2):
        for p in cols[k]:
            for q in cols[k + 1]: line([p, q], 1.55, mix_white(d['line'], .35), .07)
    for c in cols:
        for p in c: disc(p[0], p[1], .42, 1.6, d['line'])
for d in L['states']:
    for i in range(d['n']):
        disc(d['x'], d['y'] + i * d['d'] * 1.12 + d['d'] / 2, d['d'] / 2, 1.8, mix_white(d['col'], [.45, .9, .6, 1, .35][i % 5]))
    if d['lock']:
        lx, ly = d['x'] + 1.9, d['y'] - .3
        poly([(lx - .8, ly), (lx + .8, ly), (lx + .8, ly + 1.2), (lx - .8, ly + 1.2)], 1.9, '#1A1A1A', 'lock')
        line([(lx - .5 + i * 1.0 / 12, ly - .55 * math.sin(math.pi * i / 12)) for i in range(13)], 1.9, '#1A1A1A', .2)
M = L['matrix']; cs = M['s'] / M['n']; vals = [.25, .8, .5, 1, .35, .6, .2, .9, .45, .7]
for i in range(M['n']):
    for j in range(M['n']):
        a = vals[(i * 3 + j * 7) % len(vals)]
        poly([(M['x'] + j * cs + .1, M['y'] + i * cs + .1), (M['x'] + (j + 1) * cs - .1, M['y'] + i * cs + .1), (M['x'] + (j + 1) * cs - .1, M['y'] + (i + 1) * cs - .1), (M['x'] + j * cs + .1, M['y'] + (i + 1) * cs - .1)], 1.5, mix_white(M['col'], a), 'cell')
for a in L['arrows']: arrow(a['p'], dashed=a.get('dashed', False), both=a.get('both', False), hx='#595959' if a.get('dashed') else '#1A1A1A')
G = L['plot']; x0, y0, w, h = G['x'], G['y'], G['w'], G['h']
line([(x0, y0), (x0, y0 + h), (x0 + w, y0 + h)], 1.0, '#1A1A1A', .2)
for k, c in enumerate(G['cols']):
    for s in range(3):
        sh = (s - 1) * .06; f = (lambda t: 1 / (1 + math.exp(-9 * (t - .5 - sh)))) if k == 0 else (lambda t: 1 - 1 / (1 + math.exp(-7 * (t - .45 - sh))))
        pts = [(x0 + .8 + (w - 1.6) * t / 40, y0 + h - 1 - (h - 2.5) * (.15 + .7 * f(t / 40))) for t in range(41)]
        line(pts, 1.1 + s * .01, c if s == 1 else mix_white(c, .45), .34 if s == 1 else .16)
# shadow catcher, light, camera, render
me = bpy.data.meshes.new('ground'); me.from_pydata([(-50, -50, -.2), (W + 50, -50, -.2), (W + 50, H + 50, -.2), (-50, H + 50, -.2)], [], [(0, 1, 2, 3)])
gr = link(bpy.data.objects.new('shadow_catcher', me)); gr.is_shadow_catcher = True; gr.hide_render = True
sun = bpy.data.lights.new('Sun', 'SUN'); sun.energy = .45; sun.angle = math.radians(7); so = link(bpy.data.objects.new('Sun', sun)); so.rotation_euler = (math.radians(-20), math.radians(16), 0)
wd = bpy.data.worlds.new('W'); wd.use_nodes = True; wd.node_tree.nodes['Background'].inputs['Color'].default_value = (1, 1, 1, 1); wd.node_tree.nodes['Background'].inputs['Strength'].default_value = 1.6; sc.world = wd
cd = bpy.data.cameras.new('Cam'); cam = link(bpy.data.objects.new('Cam', cd)); cam.location = (W / 2, H / 2, 200); cd.type = 'ORTHO'; cd.ortho_scale = max(W, H); cd.clip_end = 1000; sc.camera = cam
sc.render.engine = 'CYCLES'; sc.cycles.device = 'CPU'; sc.cycles.samples = 64; sc.cycles.use_denoising = True; sc.cycles.max_bounces = 3
sc.render.resolution_x = round(W / 25.4 * L['dpi']); sc.render.resolution_y = round(H / 25.4 * L['dpi']); sc.render.resolution_percentage = 100
sc.render.film_transparent = True; sc.render.image_settings.file_format = 'PNG'; sc.render.image_settings.color_mode = 'RGBA'
sc.view_settings.view_transform = 'Standard'; sc.view_settings.look = 'None'; sc.render.filepath = OUT; sc.render.threads_mode = 'FIXED'; sc.render.threads = 8
bpy.ops.render.render(write_still=True)
bpy.ops.wm.save_as_mainfile(filepath=OUT.replace('.png', '.blend'), compress=True)
print('RENDER_DONE', OUT, sc.render.resolution_x, sc.render.resolution_y)
