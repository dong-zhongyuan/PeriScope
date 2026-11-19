"""Operational failures the cache continuation must reject before more work."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('cache_continuation', ROOT/'papers/01_mr_imaging_scan/continue_atlas_cache.py')
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)


class CacheContinuationTests(unittest.TestCase):
    def decision(self, completed=2000, previous=1990, cached=1800, free=4*2**30,
                 full_exists=False, acquisition='running', analysis='waiting_for_complete_acquisition'):
        return C.decision(dict(requested=3935, completed=completed, complete=False),
            acquisition, analysis, previous, cached, 128, free, full_exists)

    def test_complete_scope_handoff_and_last_128(self):
        self.assertEqual(self.decision(full_exists=True), 'handoff_to_full_scan')
        self.assertEqual(self.decision(analysis='running_statistics'), 'handoff_to_full_scan')
        self.assertEqual(self.decision(completed=3807, previous=3806), 'handoff_before_last_128')
        self.assertEqual(self.decision(completed=3806, previous=3805), 'batch_ready')

    def test_waits_for_a_substantial_batch_without_counting_tail_twice(self):
        self.assertEqual(self.decision(completed=1713, previous=1712, cached=1725), 'waiting_for_batch')
        self.assertEqual(self.decision(completed=1825, previous=1824, cached=1725), 'batch_ready')

    def test_attempt_reset_and_terminal_failure_do_not_restart_jobs(self):
        self.assertEqual(self.decision(completed=17, previous=3800), 'handoff_acquisition_attempt_changed')
        with self.assertRaises(RuntimeError):
            self.decision(acquisition='needs_attention')

    def test_low_space_and_wrong_plan_refused(self):
        with self.assertRaises(RuntimeError):
            self.decision(free=2*2**30-1)
        with self.assertRaises(RuntimeError):
            C.decision(dict(requested=45, completed=45, complete=False), 'running',
                'waiting_for_complete_acquisition', 0, 0, 128, 4*2**30, False)

    def test_changed_existing_record_or_payload_identity_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            old, new = Path(temp)/'old', Path(temp)/'new'
            old.mkdir(); new.mkdir()
            record = dict(index=0, identity='old', imported=True)
            file = dict(name='0000_drop.csv.gz', sha256='old', shared_inode=1, size_bytes=10)
            for folder in [old,new]:
                (folder/'SNAPSHOT_CHUNKS.json').write_text(json.dumps([record]))
                (folder/'CACHE_FILES.json').write_text(json.dumps([file]))
            self.assertEqual(C.reuse_summary(old,new)['existing_traits_with_identical_records'],1)
            (new/'SNAPSHOT_CHUNKS.json').write_text(json.dumps([{**record,'identity':'different'}]))
            with self.assertRaises(RuntimeError): C.reuse_summary(old,new)
            (new/'SNAPSHOT_CHUNKS.json').write_text(json.dumps([record]))
            (new/'CACHE_FILES.json').write_text(json.dumps([{**file,'shared_inode':2}]))
            with self.assertRaises(RuntimeError): C.reuse_summary(old,new)

    def test_actual_failed_child_exit_is_saved_and_raises(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(C,'BASE',Path(temp)):
            log=Path(temp)/'failure.log'
            with self.assertRaisesRegex(RuntimeError,'exited 7'):
                C.run_child([sys.executable,'-c','raise SystemExit(7)'],log,'test_child',{})
            receipt=json.loads(log.with_suffix('.exit.json').read_text())
            self.assertEqual(receipt['actual_exit_code'],7)
            self.assertFalse(C.alive(receipt['pid']))

    def test_actual_success_is_observed_before_receipt(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(C,'BASE',Path(temp)):
            log=Path(temp)/'success.log'
            result=C.run_child([sys.executable,'-c','print("operational child test")'],log,'test_child',{})
            self.assertEqual(result['actual_exit_code'],0)
            self.assertEqual(result['log_sha256'],C.pc.atlas.sha(log))
            self.assertFalse(C.alive(result['pid']))


if __name__ == '__main__':
    unittest.main()
