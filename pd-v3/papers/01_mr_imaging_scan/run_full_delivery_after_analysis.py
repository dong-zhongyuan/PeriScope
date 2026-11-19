"""Continue the already specified full-result audits and rendering on A6000.

Never restarts acquisition/MR, changes estimates, overwrites existing results,
or treats automated rendering as visual acceptance.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import fcntl,hashlib,json,os,subprocess,sys,threading,time,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/'papers/01_mr_imaging_scan'
OUT=HERE/'results/full_delivery_20261117'
BASE=HERE/'results/full_atlas_20261116/full'
PRECISION=HERE/'results/precision_overlap_20261116/full'
STATE_LOCK=threading.Lock()

def now():return datetime.now(timezone.utc).isoformat()
def sha(p):
    with p.open('rb') as s:return hashlib.file_digest(s,'sha256').hexdigest()
def save(path,value):
    tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)
def require(ok,message):
    if not ok:raise RuntimeError(message)

def live_analysis(pid):
    try:
        cmd=(Path('/proc')/str(pid)/'cmdline').read_bytes()
        stat=(Path('/proc')/str(pid)/'stat').read_text().split(') ',1)[1].split()[0]
        return b'run_full_atlas_after_acquisition.py' in cmd and stat!='Z'
    except FileNotFoundError:return False

def verified_stage(name,path):
    report=json.loads(path.read_text())
    require(report.get('complete') is True and report.get('mode')=='full','Incomplete stage output: '+name)
    if name=='complete_result_audit':
        require(report.get('formal_full_atlas_verified') is True,'Not a formal full result verification')
        require(all(report['policies'][p]['rows']==11655183 for p in ('drop','infer')),'Wrong numerical verification scope')
    elif name=='precision':
        require(report.get('analysed_pairs_per_policy')==11655183,'Wrong precision family')
        require(all(report['policies'][p]['rows']==11655183 for p in ('drop','infer')),'Wrong precision rows')
    elif name=='precision_verification':require(report.get('formal_full_precision_verified') is True,'Not full precision verification')
    elif name=='presentation':require(report.get('formal_full_atlas_result') is True,'Not full presentation')
    elif name=='presentation_verification':
        require(report.get('formal_full_presentation_verified') is True,'Not full presentation verification')
        require(all(report['policies'][p]['planned_rows']==11655183 for p in ('drop','infer')),'Wrong presentation verification rows')
    elif name=='render':
        require(report.get('formal_full_atlas_result') is True,'Not full rendering')
        require(report.get('visual_review_complete') is False,'Unexpected automatic visual acceptance')
    else:raise RuntimeError('Unknown stage')
    return dict(report=str(path.relative_to(ROOT)),report_sha256=sha(path),mode='full')

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'.runner.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        require(not (OUT/'RUN.json').exists(),'Prior delivery run exists; inspect before restart')
        plan_path=OUT/'PLAN.json';plan=json.loads(plan_path.read_text())
        def check_code():
            for relative,expected in plan['source_code_sha256'].items():require(sha(ROOT/relative)==expected,'Frozen delivery code changed: '+relative)
            for relative,expected in plan['implementation_evidence_sha256'].items():require(sha(ROOT/relative)==expected,'Implementation evidence changed: '+relative)
        check_code()
        state=dict(pid=os.getpid(),status='waiting_for_complete_MR',started_at_utc=now(),
          plan_sha256=sha(plan_path),steps={},complete=False,computational_steps_complete=False,
          visual_review_complete=False,new_MR_branch=False,source_code_sha256=plan['source_code_sha256'])
        save(OUT/'RUN.json',state)
        def update(**kwargs):
            with STATE_LOCK:state.update(kwargs);save(OUT/'RUN.json',state)
        env=os.environ.copy();env.update(LC_ALL='C',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
          PDV3_RSCRIPT='/public/home/mengxl/dzy/pd-v3_runtime/bin/Rscript')
        stages={
          'complete_result_audit':('verify_complete_atlas.py',BASE/'COMPLETE_RESULT_VERIFICATION.json'),
          'precision':('audit_mr_precision.py',PRECISION/'AUDIT.json'),
          'precision_verification':('verify_complete_precision.py',PRECISION/'COMPLETE_PRECISION_VERIFICATION.json'),
          'presentation':('s19_prepare_atlas_presentation.py',BASE/'presentation/MANIFEST.json'),
          'presentation_verification':('verify_complete_presentation.py',BASE/'presentation/COMPLETE_PRESENTATION_VERIFICATION.json'),
          'render':('s20_render_atlas_overview.py',BASE/'overview_figure/MANIFEST.json')}
        def run_stage(name):
            check_code();script,report=stages[name]
            require(not report.exists(),'Stage output already exists; no automatic overwrite: '+str(report))
            cmd=[sys.executable,'-u',str(HERE/script),'--mode','full']
            logpath=OUT/(name+'.log')
            with logpath.open('xb') as log:
                process=subprocess.Popen(cmd,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
                with STATE_LOCK:
                    state['steps'][name]=dict(status='running',pid=process.pid,command=cmd,log=str(logpath.relative_to(ROOT)),started_at_utc=now())
                    save(OUT/'RUN.json',state)
                code=process.wait()
            with STATE_LOCK:
                state['steps'][name].update(exit_code=code,ended_at_utc=now(),status='verifying_output' if code==0 else 'failed')
                save(OUT/'RUN.json',state)
            require(code==0,f'{name} exited {code}; preserved '+str(logpath))
            evidence=verified_stage(name,report);check_code()
            with STATE_LOCK:
                state['steps'][name].update(status='complete',**evidence);save(OUT/'RUN.json',state)
            print(name+' complete',flush=True)
        def branch(names):
            for name in names:run_stage(name)
        try:
            for name,(_,report) in stages.items():require(not report.exists(),'Existing stage requires review: '+name)
            while True:
                a=json.loads((ROOT/'data/big40_atlas/ANALYSIS_RUN.json').read_text())
                if a.get('status')=='complete':
                    result=json.loads((BASE/'RUN.json').read_text())
                    require(result.get('complete') is True and result.get('status')=='complete'
                      and result.get('planned_pairs')==11655183 and result.get('completed_traits')==3963,
                      'Analysis completion receipt has wrong scope')
                    update(analysis_run_sha256=sha(ROOT/'data/big40_atlas/ANALYSIS_RUN.json'),
                      full_run_sha256=sha(BASE/'RUN.json'),status='running_complete_result_audit')
                    break
                require(a.get('status') in ('waiting_for_complete_acquisition','running_statistics'),
                  'MR requires attention: '+str(a.get('status')))
                require(live_analysis(a['pid']),'Authoritative MR process is missing; do not restart automatically')
                update(last_live_analysis_check_utc=now(),analysis_pid=a['pid'],analysis_status=a['status'])
                time.sleep(45)
            run_stage('complete_result_audit')
            update(status='running_original_followup_steps')
            errors=[]
            with ThreadPoolExecutor(max_workers=2) as pool:
                jobs=[pool.submit(branch,['precision','precision_verification']),
                      pool.submit(branch,['presentation','presentation_verification','render'])]
                for job in jobs:
                    try:job.result()
                    except Exception as error:errors.append(repr(error))
            require(not errors,'Followup step failed: '+'; '.join(errors))
            update(status='awaiting_visual_review',computational_steps_complete=True,ended_at_utc=now(),
              remaining_scope='Inspect full rendered figures and scientific results; candidate regional colocalization/reverse MR and actual overlap questions remain. Papers 2/3 are not completed by this runner.')
        except Exception as error:
            update(status='needs_attention',error=repr(error),ended_at_utc=now());traceback.print_exc();raise

if __name__=='__main__':main()
