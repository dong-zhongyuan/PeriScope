"""Prepare a relocated source copy; never edit the released checkout."""
from pathlib import Path
import argparse,shutil,sys
p=argparse.ArgumentParser();p.add_argument('--assets',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
root=Path(__file__).resolve().parents[1];out=a.output.expanduser().resolve();assets=a.assets.expanduser().resolve()
if out.exists():p.error('output must be a new directory')
if out==root or root in out.parents:p.error('output must be outside the released checkout')
out.mkdir(parents=True)
for name in ['src','scripts','_work','config','contracts','registry']:
 shutil.copytree(root/name,out/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
for name in ['pyproject.toml','requirements.txt']:shutil.copy2(root/name,out/name)
replacements=[('/public/home/mengxl/dzy/pd_product_assets',str(assets)),('/public/home/mengxl/dzy/pd_product',str(out)),('/public/home/mengxl/dzy/envs/pd_bbm/bin/python',sys.executable)]
changed=0
for f in out.rglob('*'):
 if f.is_file() and f.suffix.lower() in {'.py','.r','.sh','.json','.yaml','.yml','.toml'}:
  old=f.read_text();new=old
  for x,y in replacements:new=new.replace(x,y)
  if new!=old:f.write_text(new);changed+=1
print(f'Prepared {out}; relocated {changed} source/config files. Configure remaining external data and R-library paths before full execution.')
