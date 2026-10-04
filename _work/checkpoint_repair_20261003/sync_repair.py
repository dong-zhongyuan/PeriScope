"""Follow the currently running repair and synchronize its completed artifacts."""
import hashlib,json,subprocess,sys,time
from pathlib import Path
P=Path('/Users/dawnmeng/Desktop/董仲元');W=P/'pd_product_figure_design/manuscript_review/checkpoint_repair_20261003';O=P/'delivery/benchmark_20261002';R=O/'seed_repair';R.mkdir(exist_ok=True)
REMOTE='/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002'
SSH=['ssh','-S','/tmp/pd_rescreen_20261002_ssh','-o','BatchMode=yes','-o','ConnectTimeout=10']
script="""import json
from pathlib import Path
r=Path('/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002/seed_repair')
out={}
for n in ['run_state.json','completion.json','finalization_error.json']:
 if (r/n).exists():out[n]=json.loads((r/n).read_text())
for s in [47,48,49]:
 p=r/f'models/seed_{s}/validation.json'
 if p.exists():
  v=json.loads(p.read_text())[-1];out[str(s)]={k:v[k] for k in ['epoch','step','seconds','checkpoint_eligible','selection_score']}
print(json.dumps(out))
"""
def save(obj):
 p=R/'local_progress.json.tmp';p.write_text(json.dumps(obj,indent=2));p.replace(R/'local_progress.json')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<22),b''):h.update(b)
 return h.hexdigest()
def main():
 failures=0
 for _ in range(240):
  try:
   r=subprocess.run(SSH+['a6000','/public/home/mengxl/dzy/envs/pd_bbm/bin/python -'],input=script,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30,check=True)
   state=json.loads(r.stdout);save({'status':'running','updated_unix':time.time(),'server':state});failures=0
  except Exception as e:
   failures+=1
   if failures>=3:raise RuntimeError('Could not read server state after three attempts; server training remains independent') from e
   time.sleep(20);continue
  if 'finalization_error.json' in state:raise RuntimeError(state['finalization_error.json'])
  if state.get('run_state.json',{}).get('status')=='failed':raise RuntimeError(state['run_state.json'])
  if state.get('completion.json',{}).get('status')=='complete':break
  time.sleep(45)
 else:raise TimeoutError('Three-hour local follow-up elapsed; inspect server status')
 subprocess.run([sys.executable,'-u',str(P/'pd_product_figure_design/manuscript_review/benchmark_full_20261003/sync_delivery.py')],check=True)
 subprocess.run(['rsync','-az','--exclude=last_model.pt','--exclude=launch.lock','-e',' '.join(SSH),'a6000:'+REMOTE+'/seed_repair/',str(R)+'/'],check=True)
 for name in ['checkpoint_overrides.json','补训结果说明.txt']:
  subprocess.run(['scp','-q','-o','ControlPath=/tmp/pd_rescreen_20261002_ssh','-o','BatchMode=yes','a6000:'+REMOTE+'/'+name,str(O/name)],check=True)
 checks=json.loads((R/'checkpoint_validation.json').read_text())
 for x in checks['records']:assert sha(R/f"models/seed_{x['seed']}/model.pt")==x['model_sha256']
 save({'status':'complete','updated_unix':time.time(),'checks':checks,'delivery':str(O)})
 (W/'完成状态.txt').write_text('47/48/49补训、检查点验证、四轴评估及十种子汇总均已完成，服务器与本地文件哈希复核通过。\n'+str(O/'补训结果说明.txt')+'\n')
 print('REPAIR_LOCAL_DELIVERY_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save({'status':'needs_attention','error':str(e),'time':time.time()});raise
