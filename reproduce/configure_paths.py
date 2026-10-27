"""Prepare a relocated source copy; never edit the released checkout."""
from pathlib import Path
import argparse,re,shutil,sys
p=argparse.ArgumentParser();p.add_argument('--assets',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
root=Path(__file__).resolve().parents[1];out=a.output.expanduser().resolve();assets=a.assets.expanduser().resolve()
if out.exists():p.error('output must be a new directory')
if out==root or root in out.parents:p.error('output must be outside the released checkout')
out.mkdir(parents=True)
for name in ['src','scripts','_work','config','contracts','registry','results','figure_source_data','reproduce','figures']:
 shutil.copytree(root/name,out/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
for name in ['pyproject.toml','requirements.txt']:shutil.copy2(root/name,out/name)
# Scripts that concatenate W + '/diag/...' need the historical directory alias.
revision=out/'manuscript_revision_20261017';revision.mkdir()
(revision/'diag').symlink_to(Path('../results/nmi_revision_20261017'),target_is_directory=True)
(revision/'notes').symlink_to(Path('../_work/nmi_revision_20261017'),target_is_directory=True)
base='/public/home/mengxl/dzy/pd_product'
replacements=[
 (base+'/publication/analysis/pd_product',str(out)),
 (base+'/manuscript_revision_20261017/diag',str(out/'results/nmi_revision_20261017')),
 (base+'/manuscript_revision_20261017/notes',str(out/'_work/nmi_revision_20261017')),
 ('/public/home/mengxl/dzy/pd_product_assets',str(assets)),
 (base,str(out)),
 ('/public/home/mengxl/dzy/envs/pd_bbm/bin/python',sys.executable),
]
changed=0
for f in out.rglob('*'):
 if f.is_file() and f.suffix.lower() in {'.py','.r','.sh','.json','.yaml','.yml','.toml'}:
  old=f.read_text();new=old
  for x,y in replacements:new=new.replace(x,y)
  # Some entry points derive the assets root rather than storing its absolute path.
  new=re.sub(r"ROOT\.parent\s*/\s*(['\"])pd_product_assets\1",lambda _: 'Path('+repr(str(assets))+')',new)
  if new!=old:f.write_text(new);changed+=1
print(f'Prepared {out}; relocated {changed} source/config files. Configure remaining external data and R-library paths before full execution.')
