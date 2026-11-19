"""Continue the existing verified cache batches; hand off before final acquisition."""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import precompute_atlas_cache as pc

BASE = pc.RESULTS / 'cache_continuation_nativefix'
VERIFIER = HERE / 'verify_staged_atlas_snapshot.py'


def utc():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    return json.loads(path.read_text())


def write_new(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


def alive(pid):
    result = subprocess.run(['ps', '-p', str(pid), '-o', 'pid=,stat=,args='],
                            capture_output=True, text=True)
    fields = result.stdout.strip().split()
    return len(fields) >= 2 and int(fields[0]) == pid and not fields[1].startswith('Z')


def decision(progress, acquisition, analysis, previous_count, cached, minimum_new, free, full_exists):
    """No writes or new jobs are permitted after the original full-scan handoff."""
    if (full_exists or acquisition == 'acquisition_complete'
            or analysis in ('running_statistics', 'complete')):
        return 'handoff_to_full_scan'
    if acquisition != 'running' or analysis != 'waiting_for_complete_acquisition':
        raise RuntimeError('Main acquisition/statistics requires attention')
    if progress['requested'] != 3935 or not 0 <= progress['completed'] <= 3935:
        raise RuntimeError('Acquisition plan dimensions changed')
    if progress['complete'] or 3935 - progress['completed'] <= 128:
        return 'handoff_before_last_128'
    if progress['completed'] < previous_count:
        return 'handoff_acquisition_attempt_changed'
    if free < 2 * 2**30:
        raise RuntimeError('Less than 2 GiB free; retain existing caches and stop')
    # The 28 completed extension traits are guaranteed. Do not assume the eight
    # early tail imports remain outside the main counter; this lower bound can
    # delay a batch but cannot over-count its ready inputs.
    return 'batch_ready' if progress['completed'] + 28 - cached >= minimum_new else 'waiting_for_batch'


def checked_snapshot(folder):
    manifest = read(folder / 'DELIVERABLES.json')
    if manifest.get('complete') is not True:
        raise RuntimeError('Baseline snapshot has no complete fixed manifest')
    for item in manifest['artifacts']:
        if pc.atlas.sha(ROOT / item['path']) != item['sha256']:
            raise RuntimeError(f'Baseline artifact changed: {item["path"]}')
    proof = read(folder / 'VERIFICATION.json')
    run = read(folder / 'SNAPSHOT_RUN.json')
    if (proof.get('complete') is not True or run.get('status') != 'snapshot_complete'
            or run.get('mode') != 'full' or run.get('planned_traits') != 3963
            or proof['snapshot_traits'] != run['completed_traits']
            or run['completed_traits'] != run['imported_traits']):
        raise RuntimeError('Baseline snapshot scope or completion differs')
    for live, frozen in [('RUN.json', 'SNAPSHOT_RUN.json'), ('CHUNKS.json', 'SNAPSHOT_CHUNKS.json')]:
        if (pc.CACHE_FULL / live).read_bytes() != (folder / frozen).read_bytes():
            raise RuntimeError('Current cache is not the specified verified baseline')
    if alive(run['pid']):
        raise RuntimeError('Previous precompute process is still present')
    return proof, run


def reuse_summary(previous, current):
    old_records, records = read(previous / 'SNAPSHOT_CHUNKS.json'), read(current / 'SNAPSHOT_CHUNKS.json')
    old_files, files = read(previous / 'CACHE_FILES.json'), read(current / 'CACHE_FILES.json')
    record_map, file_map = {x['index']: x for x in records}, {x['name']: x for x in files}
    if len(record_map) != len(records) or len(file_map) != len(files):
        raise RuntimeError('New snapshot contains duplicate identities')
    if any(record_map.get(x['index']) != x for x in old_records):
        raise RuntimeError('A previously verified cache record changed')
    if any(file_map.get(x['name']) != x for x in old_files):
        raise RuntimeError('A previously verified cache file identity changed')
    return dict(existing_traits_with_identical_records=len(old_records),
                existing_files_with_identical_sha_size_and_inode=len(old_files),
                new_traits=len(records) - len(old_records), new_files=len(files) - len(old_files),
                comparison_basis='Current full verifier freshly read and hashed every payload; compare its fixed file records to the previous verified snapshot.')


def run_child(command, log_path, stage, state):
    with log_path.open('x') as stream:
        child = subprocess.Popen(command, cwd=ROOT, stdout=stream,
                                 stderr=subprocess.STDOUT, start_new_session=True)
        state.update(status=stage, child_pid=child.pid, child_command=command, updated_at_utc=utc())
        pc.atlas.save_json(BASE / 'RUN.json', state)
        while True:
            try:
                code = child.wait(timeout=30)
                break
            except subprocess.TimeoutExpired:
                state['updated_at_utc'] = utc()
                pc.atlas.save_json(BASE / 'RUN.json', state)
        receipt = dict(pid=child.pid, command=command, actual_exit_code=code,
                       observed_at_utc=utc(), log_sha256=pc.atlas.sha(log_path))
        write_new(log_path.with_suffix('.exit.json'), receipt)
        if code:
            raise RuntimeError(f'{stage} exited {code}; no automatic restart')
        return receipt


def finish_snapshot(previous, current, attempt, exits, code):
    proof, run = read(current / 'VERIFICATION.json'), read(current / 'SNAPSHOT_RUN.json')
    if proof.get('complete') is not True or run.get('status') != 'snapshot_complete':
        raise RuntimeError('Current payload verification is incomplete')
    if proof['snapshot_traits'] != run['completed_traits'] or run['imported_traits'] != run['completed_traits']:
        raise RuntimeError('Current snapshot scope differs')
    reuse = reuse_summary(previous, current)
    if reuse['new_traits'] <= 0 or reuse['new_files'] != 3 * reuse['new_traits']:
        raise RuntimeError('Continuation did not add complete trait caches')
    reuse.update(observed_at_utc=utc(), previous_snapshot=str(previous.relative_to(ROOT)),
                 previous_snapshot_manifest_sha256=pc.atlas.sha(previous / 'DELIVERABLES.json'),
                 reused_cache_contents_unchanged=True, old_snapshots_not_replaced=True)
    write_new(current / 'REUSE_VERIFICATION.json', reuse)
    receipt = dict(complete=True, observed_at_utc=utc(), supervisor_pid=os.getpid(),
        attempt=str(attempt.relative_to(ROOT)), children=exits, supervisor_code_sha256=code,
        full_verifier_report_sha256=pc.atlas.sha(current / 'VERIFICATION.json'),
        precompute_exit_observed_by_supervisor=True, formal_full_scan_complete=False,
        clarification='The payload verifier does not itself wait for precompute; this supervisor directly waited for both children and observed exit 0.')
    write_new(current / 'CONTINUATION_EXECUTION.json', receipt)
    (current / 'README.md').write_text(
        '# 同一主扫描的自动接续缓存快照\n\n'
        f'累计{proof["snapshot_traits"]:,}项表型，新增{reuse["new_traits"]:,}项。'
        f'两策略各{proof["actual_rows_read_per_policy"]["drop"]:,}行均已完整读取核验，'
        f'{proof["cache_files_checked"]:,}个文件及主目录链接通过。\n\n'
        f'此前{reuse["existing_traits_with_identical_records"]:,}项记录与'
        f'{reuse["existing_files_with_identical_sha_size_and_inode"]:,}个文件的SHA、大小、inode保持。'
        '原快照保留；两子进程均由监督程序实际等待并观察退出0。\n\n'
        '这是原主扫描的未校正缓存，不是正式全量MR完成，不产生最终q值或候选。'
        '代码、日志、逐文件记录和退出证据见本目录清单及CONTINUATION_EXECUTION.json。\n', encoding='utf-8')
    paths = [p for p in current.iterdir() if p.is_file()]
    paths += [p for p in attempt.iterdir() if p.is_file()]
    paths += [ROOT / relative for relative in code]
    artifacts = [dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, sha256=pc.atlas.sha(p))
                 for p in sorted(set(paths))]
    write_new(current / 'DELIVERABLES.json', dict(complete=True, created_at_utc=utc(),
        snapshot_traits=proof['snapshot_traits'], artifact_count=len(artifacts), artifacts=artifacts,
        formal_full_scan_complete=False, final_corrections_computed=False, candidates_selected=False))
    return proof


def run(baseline, minimum_new=128, check_only=False):
    if minimum_new < 100:
        raise ValueError('Use a substantial batch of at least 100 newly ready traits')
    baseline = baseline.resolve()
    proof, _ = checked_snapshot(baseline)
    _, _, _, _, frozen_science, frozen_inputs = pc.context()
    code = {str(p.relative_to(ROOT)): pc.atlas.sha(p) for p in
            [Path(__file__), Path(pc.__file__), VERIFIER]}
    def check_code():
        for relative, expected in code.items():
            if pc.atlas.sha(ROOT / relative) != expected:
                raise RuntimeError(f'Continuation code changed: {relative}')
        _, _, _, _, current_code, current_inputs = pc.context()
        if current_code != frozen_science or current_inputs != frozen_inputs:
            raise RuntimeError('Frozen analysis inputs or implementation changed')
    progress = read(pc.DATA / 'extraction_progress.json')
    previous_count = progress['completed']
    acquisition, analysis = read(pc.DATA / 'RUN.json'), read(pc.DATA / 'ANALYSIS_RUN.json')
    state = dict(pid=os.getpid(), status='initializing', started_at_utc=utc(),
        minimum_new_traits=minimum_new, latest_verified_snapshot=str(baseline.relative_to(ROOT)),
        latest_verified_traits=proof['snapshot_traits'], source_code_sha256=code,
        original_scientific_code_sha256=frozen_science, original_input_sha256=frozen_inputs,
        no_final_corrections=True, no_candidate_selection=True, main_scan_complete=False,
        completed_batches=[])
    action = decision(progress, acquisition['status'], analysis['status'], previous_count,
                      proof['snapshot_traits'], minimum_new, shutil.disk_usage(ROOT).free,
                      (pc.RESULTS / 'full/RUN.json').exists())
    if check_only:
        print(json.dumps(dict(action=action, verified_traits=proof['snapshot_traits'],
            main_queue_completed=previous_count, minimum_new=minimum_new, no_outputs_written=True), indent=2))
        return
    BASE.mkdir(exist_ok=True)
    with (BASE / '.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (BASE / 'RUN.json').exists():
            raise RuntimeError('A continuation run record already exists; inspect it before restarting')
        pc.atlas.save_json(BASE / 'RUN.json', state)
        try:
            while True:
                check_code()
                progress = read(pc.DATA / 'extraction_progress.json')
                acquisition, analysis = read(pc.DATA / 'RUN.json'), read(pc.DATA / 'ANALYSIS_RUN.json')
                action = decision(progress, acquisition['status'], analysis['status'], previous_count,
                    proof['snapshot_traits'], minimum_new, shutil.disk_usage(ROOT).free,
                    (pc.RESULTS / 'full/RUN.json').exists())
                state.update(status=action, main_queue_completed=progress['completed'],
                             updated_at_utc=utc(), child_pid=None)
                pc.atlas.save_json(BASE / 'RUN.json', state)
                if action.startswith('handoff_'):
                    state['ended_at_utc'] = utc()
                    pc.atlas.save_json(BASE / 'RUN.json', state)
                    print(json.dumps(state, indent=2), flush=True)
                    return
                if not alive(acquisition['pid']) or not alive(analysis['pid']):
                    raise RuntimeError('Main acquisition/statistics process is missing')
                pc.import_guard()
                previous_count = progress['completed']
                if action == 'waiting_for_batch':
                    time.sleep(45)
                    continue
                checked_snapshot(baseline)
                stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
                attempt = BASE / 'attempts' / stamp
                attempt.mkdir(parents=True, exist_ok=False)
                write_new(attempt / 'LAUNCH.json', dict(observed_at_utc=utc(), supervisor_pid=os.getpid(),
                    baseline=str(baseline.relative_to(ROOT)), baseline_manifest_sha256=pc.atlas.sha(baseline / 'DELIVERABLES.json'),
                    main_completed=progress['completed'], minimum_new_traits=minimum_new, code_sha256=code))
                exits = {'precompute': run_child([sys.executable, '-u', str(Path(pc.__file__)), '--mode', 'full', '--import-main'],
                    attempt / 'precompute.log', 'running_precompute', state)}
                exits['verification'] = run_child([sys.executable, '-u', str(VERIFIER)],
                    attempt / 'verification.log', 'running_verification', state)
                check_code()
                current_run = read(pc.CACHE_FULL / 'RUN.json')
                finished = datetime.fromisoformat(current_run['finished_at_utc']).strftime('%Y%m%dT%H%M%SZ')
                current = pc.CACHE_FULL / 'snapshots' / finished
                proof = finish_snapshot(baseline, current, attempt, exits, code)
                baseline = current
                state['completed_batches'].append(dict(snapshot=str(current.relative_to(ROOT)),
                    traits=proof['snapshot_traits'], verified_at_utc=proof['verified_at_utc']))
                state.update(latest_verified_snapshot=str(current.relative_to(ROOT)),
                             latest_verified_traits=proof['snapshot_traits'])
                pc.atlas.save_json(BASE / 'RUN.json', state)
                print(json.dumps(state['completed_batches'][-1]), flush=True)
        except Exception as error:
            state.update(status='stopped_with_reason', error=str(error), ended_at_utc=utc())
            pc.atlas.save_json(BASE / 'RUN.json', state)
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--minimum-new', type=int, default=128)
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    run(args.baseline, args.minimum_new, args.check_only)
