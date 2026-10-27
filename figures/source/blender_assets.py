# -*- coding: utf-8 -*-
import bpy,math,random,pathlib,json
from mathutils import Vector
R=pathlib.Path(__file__).resolve().parents[1];O=R.parents[1]/'workflow_assets';random.seed(16)
def mat(name,c,rough=.72):
 m=bpy.data.materials.new(name);m.diffuse_color=(*c,1);m.use_nodes=True;p=m.node_tree.nodes.get('Principled BSDF');p.inputs['Base Color'].default_value=(*(v/12.92 if v<=.04045 else ((v+.055)/1.055)**2.4 for v in c),1);p.inputs['Roughness'].default_value=rough;return m
def sphere(name,pos,scale,m):
 bpy.ops.mesh.primitive_uv_sphere_add(segments=32,ring_count=20,location=pos);o=bpy.context.object;o.name=name;o.scale=scale;o.data.materials.append(m)
 for p in o.data.polygons:p.use_smooth=True
 return o
def curve(name,pts,r,m):
 c=bpy.data.curves.new(name,'CURVE');c.dimensions='3D';c.resolution_u=12;c.bevel_depth=r;c.bevel_resolution=3;s=c.splines.new('BEZIER');s.bezier_points.add(len(pts)-1)
 for p,co in zip(s.bezier_points,pts):p.co=co;p.handle_left_type='AUTO';p.handle_right_type='AUTO'
 o=bpy.data.objects.new(name,c);bpy.context.collection.objects.link(o);c.materials.append(m);return o
def antibody(name,p,s,m):
 x,y,z=p;curve(name+'_stem',[(x,y,z-s*.5),(x,y,z),(x-s*.36,y,z+s*.36)],s*.065,m);curve(name+'_arm',[(x,y,z),(x+s*.36,y,z+s*.36)],s*.065,m)
def cell(name,p,s,color,kind):
 x,y,z=p;sphere(name+'_membrane',p,(s,s*.85,s),color)
 nuc=mat(name+'_nucleus',(.29,.25,.52));sphere(name+'_nucleus',(x-s*.12,y-s*.76,z+s*.12),(s*.53,s*.3,s*.55),nuc)
 if kind=='mono':
  sphere(name+'_kidney_nucleus_lobe',(x+s*.28,y-s*.8,z-s*.16),(s*.38,s*.23,s*.37),nuc)
 for i in range(23):
  a=random.uniform(0,math.tau);r=random.uniform(.62,.86)*s;xx=x+r*math.cos(a);zz=z+r*math.sin(a)
  sphere(name+f'_granule_{i}',(xx,y-s*.62,zz),(s*.035,s*.035,s*.035),nuc)
 if kind=='dendritic':
  for i in range(10):
   a=i*math.tau/10;pts=[(x+s*.8*math.cos(a),y,z+s*.8*math.sin(a)),(x+s*1.17*math.cos(a+.08),y-.03,z+s*1.17*math.sin(a+.08)),(x+s*1.43*math.cos(a),y+.03,z+s*1.43*math.sin(a))];curve(name+f'_process_{i}',pts,s*.043,color)
for key in ['pig','antibody','binder','inhibitor']:
 bpy.ops.wm.read_factory_settings(use_empty=True);sc=bpy.context.scene;sc.name=key
 blue=mat('Blood_blue',(.42,.68,.83));violet=mat('Cortex_lilac',(.66,.58,.77));pink=mat('Tissue_pink',(.81,.52,.59));green=mat('Glia_sage',(.52,.69,.54));gold=mat('ADT_gold',(.88,.62,.25));dark=mat('Ink_navy',(.14,.22,.32));light=mat('Pale_blue',(.7,.84,.9))
 extent=4.1
 if key=='pig':
  sphere('Pig_head',(0,0,0),(.78,.46,.69),pink)
  sphere('Snout',(0,-.47,-.16),(.43,.24,.26),mat('Snout_pink',(.91,.63,.68)))
  for side in [-1,1]:
   sphere('Nostril_'+str(side),(side*.16,-.68,-.15),(.07,.035,.065),dark)
   sphere('Eye_'+str(side),(side*.4,-.4,.16),(.055,.04,.07),dark)
   o=sphere('Ear_'+str(side),(side*.62,-.02,.60),(.20,.12,.39),pink);o.rotation_euler[1]=side*.48
  extent=2.9
 elif key=='antibody':
  antibody('Antibody_Y',(0,0,0),1.8,blue)
  curve('Fc_light_chain',[(-.12,-.08,-.88),(-.12,-.08,-.05),(-.77,-.08,.70)],.055,gold)
  curve('Light_chain_right',[(.12,-.08,-.05),(.77,-.08,.70)],.055,gold);extent=2.8
 elif key=='binder':
  pts=[(-.65+1.3*i/35,-.05+.15*math.sin(i*.63),.28*math.cos(i*.63))for i in range(36)]
  curve('Engineered_binder_backbone',pts,.09,violet)
  for i in range(0,36,5):sphere('Binder_contact_'+str(i),pts[i],(.14,.13,.14),gold)
  extent=2.4
 elif key=='inhibitor':
  # Open binding cleft and a small molecule, all retained as separate objects.
  for i in range(24):
   a=math.radians(45+i*270/23);sphere('Protein_cleft_domain_'+str(i),(.62*math.cos(a),0,.62*math.sin(a)),(.26,.22,.26),blue)
  atom=[(.48,-.15,.48),(.73,-.15,.73),(1.03,-.15,.66),(.90,-.15,.39)]
  curve('Inhibitor_bonds',atom+[atom[0]],.025,gold)
  for i,pt in enumerate(atom):sphere('Inhibitor_atom_'+str(i),pt,(.10,.10,.10),violet if i%2 else gold)
  extent=3.0
 elif key=='pbmc':
  cell('T_lymphocyte',(-.82,0,.52),.54,blue,'lymph');cell('Monocyte',(.75,0,.35),.68,pink,'mono');cell('Dendritic_cell',(-.27,-.15,-.83),.52,green,'dendritic');extent=3.5
 elif key=='brain':
  import numpy as np
  text=(R/'source/BrainMesh_ICBM152.nv').read_text();vals=np.fromstring('\n'.join(l for l in text.splitlines()if not l.lstrip().startswith('#')),sep=' ');nv=int(vals[0]);v=vals[1:1+nv*3].reshape(nv,3);nf=int(vals[1+nv*3]);f=vals[2+nv*3:].reshape(nf,3).astype(int)-1
  v=v-(v.max(axis=0)+v.min(axis=0))/2;v=v*(2.6/(v.max(axis=0)-v.min(axis=0)).max());v=v[:,[1,0,2]]
  me=bpy.data.meshes.new('ICBM152_cortical_surface');me.from_pydata(v.tolist(),[],f.tolist());me.update();ob=bpy.data.objects.new('ICBM152_anatomical_brain_surface',me);bpy.context.collection.objects.link(ob);ob.data.materials.append(violet)
  for p in me.polygons:p.use_smooth=True
  ob['source']='BrainNet Viewer / ICBM152 cortical surface';ob['url']='https://github.com/mingruixia/BrainNet-Viewer/blob/master/Data/SurfTemplate/BrainMesh_ICBM152.nv';ob['license']='GPL-3.0; source mesh included';extent=3.15
 elif key=='cite_cell':
  cell('CITE_cell',(0,.08,0),.95,green,'lymph')
  for i in range(12):
   a=i*math.tau/12;x=1.04*math.cos(a);z=1.04*math.sin(a);antibody(f'Surface_ADT_{i}',(x,-.05,z),.34,gold if i%3 else pink)
  for i in range(5):
   pts=[(-.45+j*.11,-.93,-.24-i*.095+.045*math.sin(j*1.8+i))for j in range(8)];curve(f'Intracellular_RNA_{i}',pts,.015,light)
  extent=3.6
 elif key=='neuroglia':
  sphere('Microglia_soma',(-.7,0,.08),(.25,.22,.3),blue);sphere('Microglia_nucleus',(-.7,-.21,.08),(.095,.04,.12),dark)
  for i in range(7):
   a=i*math.tau/7
   base=Vector((-.7,0,.08));end=base+Vector((math.cos(a),0,math.sin(a)))*.73;mid=(base+end)/2+Vector((0,-.05,0));curve(f'Microglia_process_{i}',[base,mid,end],.025,blue)
   for j,off in enumerate([-.35,.35]):
    dst=end+Vector((math.cos(a+off),0,math.sin(a+off)))*.3;curve(f'Microglia_branch_{i}_{j}',[mid,end,dst],.012,blue)
  sphere('Astrocyte_soma',(.77,0,.27),(.2,.18,.2),pink)
  for i in range(12):
   a=i*math.tau/12;st=Vector((.77,0,.27));end=st+Vector((math.cos(a),0,math.sin(a)))*.64;curve(f'Astrocyte_process_{i}',[st,(st+end)/2+Vector((0,-.03,.04)),end],.022,pink)
   for j,delta in enumerate([-.65,.55]):
    ep=end+Vector((math.cos(a+delta),-.02,math.sin(a+delta)))*.22;curve(f'Astrocyte_endfoot_{i}_{j}',[end,ep],.012,pink);sphere(f'Endfoot_{i}_{j}',ep,(.035,.035,.035),pink)
  extent=3.6
 elif key=='spatial_section':
  # Tissue section on a capture grid: individual spots and polygon tissue mesh.
  verts=[];N=80
  for i in range(N):
   a=i*math.tau/N;r=1+.13*math.sin(3*a)+.08*math.sin(7*a);verts.append((1.12*r*math.cos(a),.06,r*.8*math.sin(a)))
  me=bpy.data.meshes.new('Section_mesh');me.from_pydata(verts,[],[list(range(N))]);me.update();ob=bpy.data.objects.new('Tissue_section',me);bpy.context.collection.objects.link(ob);ob.data.materials.append(mat('Histology_light',(.88,.74,.8)))
  for row in range(12):
   for col in range(15):
    x=-1.03+col*.145+(row%2)*.072;z=-.7+row*.13
    if (x/1.08)**2+(z/.77)**2<.96:sphere(f'Capture_spot_{row}_{col}',(x,-.045,z),(.028,.013,.028),[blue,pink,violet,gold,green][(row*3+col)%5])
  for j in range(12):
   a=j*math.tau/12;curve('Tissue_fiber_%02d'%j,[(.4*math.cos(a),.01,.3*math.sin(a)),(.75*math.cos(a+.13),.01,.5*math.sin(a+.13)),(1.03*math.cos(a+.2),.01,.62*math.sin(a+.2))],.012,violet)
  extent=3.1
 elif key=='protein':
  sphere('Protein_core',(0,0,0),(.45,.35,.55),blue)
  for i in range(24):
   a=i*2.399;v=1-2*(i+.5)/24;rad=math.sqrt(1-v*v);x=.5*rad*math.cos(a);y=.4*rad*math.sin(a);z=.62*v;sphere(f'Protein_domain_{i}',(x,y,z),(.22,.18,.22),blue if i%3 else violet)
  antibody('Binding_antibody',(.65,-.1,.7),.75,gold);extent=3.1
 elif key=='energy':
  verts=[];faces=[];N=45
  for j in range(N):
   for i in range(N):
    x=-1.4+i/(N-1)*2.8;y=-1.25+j/(N-1)*2.5;z=.24*(x*x+y*y)+.045*math.sin(x*3)*math.cos(y*2)-.4;verts.append((x,y,z))
  for j in range(N-1):
   for i in range(N-1):p=j*N+i;faces.append((p,p+1,p+1+N,p+N))
  me=bpy.data.meshes.new('Potential_surface_mesh');me.from_pydata(verts,[],faces);me.update();o=bpy.data.objects.new('Illustrative_potential_surface',me);bpy.context.collection.objects.link(o)
  for k in range(12):o.data.materials.append(mat(f'Energy_band_{k}',(.58+k*.016,.70+k*.009,.79+k*.008)))
  for f in me.polygons:f.material_index=min(11,int((f.center.z+.5)*13));f.use_smooth=True
  for j in range(0,N,5):curve('Potential_grid_%02d'%j,[verts[j*N+i]for i in range(N)],.009,dark)
  for i in range(0,N,5):curve('Potential_grid_v%02d'%i,[verts[j*N+i]for j in range(N)],.009,dark)
  pts=[]
  for i in range(20):t=i/19;x=1.15*(1-t);y=-.86*(1-t);z=.24*(x*x+y*y)+.045*math.sin(x*3)*math.cos(y*2)-.37;pts.append((x,y,z))
  curve('Descent_path',pts,.035,pink);sphere('Equilibrium',pts[-1],(.09,.09,.09),gold);extent=3.8
 world=bpy.data.worlds.new('World');sc.world=world;world.use_nodes=True;world.node_tree.nodes['Background'].inputs[0].default_value=(.95,.97,1,1);world.node_tree.nodes['Background'].inputs[1].default_value=.35
 for name,loc,energy,size in [('Key',(-3,-4,6),180,5),('Fill',(4,-1,3),100,4),('Rim',(0,3,4),150,3)]:
  d=bpy.data.lights.new(name,'AREA');d.energy=energy;d.shape='DISK';d.size=size;ob=bpy.data.objects.new(name,d);sc.collection.objects.link(ob);ob.location=loc;ob.rotation_euler=(-ob.location).to_track_quat('-Z','Y').to_euler()
 d=bpy.data.cameras.new('Camera');ob=bpy.data.objects.new('Camera',d);sc.collection.objects.link(ob);ob.location=(0,-8,3.0 if key!='energy'else 6.5);target=Vector((0,0,-.05));ob.rotation_euler=(target-ob.location).to_track_quat('-Z','Y').to_euler();d.type='ORTHO';d.ortho_scale=extent;sc.camera=ob
 sc.render.engine='BLENDER_EEVEE_NEXT';sc.render.resolution_x=900;sc.render.resolution_y=900;sc.render.resolution_percentage=100;sc.render.image_settings.file_format='PNG';sc.render.image_settings.color_mode='RGBA';sc.render.film_transparent=True;sc.render.filepath=str(O/f'{key}.png');sc.view_settings.view_transform='Standard';sc.view_settings.look='None';sc['purpose']='Stylized scientific workflow illustration; not an anatomical atlas or empirical result'
 bpy.context.preferences.filepaths.save_version=0;bpy.ops.wm.save_as_mainfile(filepath=str(O/'blender'/f'Asset_{key}.blend'),compress=True);bpy.ops.render.render(write_still=True)
 print('ASSET_READY',key,len(sc.objects),flush=True)
