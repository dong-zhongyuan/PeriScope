import importlib.util
from pathlib import Path
import sys
import unittest
import io
import gzip
import tarfile
import tempfile

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[1] / 'papers/01_mr_imaging_scan'
sys.path.insert(0, str(PAPER))
spec = importlib.util.spec_from_file_location('expanded_cli', PAPER / 's11_run_expanded.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
qspec = importlib.util.spec_from_file_location('qsm_parser', PAPER / 's13_extract_qsm_panel.py')
qmodule = importlib.util.module_from_spec(qspec)
qspec.loader.exec_module(qmodule)


class FrequencyJoinTests(unittest.TestCase):
    def test_exact_and_reversed_frequency_orientation(self):
        y = pd.DataFrame([dict(build_y='GRCh37', chrom_y='01', pos_y=100, ea_y='A', oa_y='T'),
                          dict(build_y='GRCh37', chrom_y=1, pos_y=100, ea_y='T', oa_y='A')])
        af = pd.DataFrame([dict(chr=1, pos37=100, ea_y='A', oa_y='T', eaf_y=.2, source_y='source')])
        out = module.attach_frequencies(y, af)
        np.testing.assert_allclose(out.eaf_y, [.2, .8])

    def test_no_cross_build_or_complement_guess(self):
        y = pd.DataFrame([dict(build_y='GRCh38', chrom_y='1', pos_y=100, ea_y='A', oa_y='G'),
                          dict(build_y='GRCh37', chrom_y=1, pos_y=100, ea_y='T', oa_y='C')])
        af = pd.DataFrame([dict(chr=1, pos37=100, ea_y='A', oa_y='G', eaf_y=.2, source_y='source')])
        self.assertTrue(module.attach_frequencies(y, af).eaf_y.isna().all())

    def test_ambiguous_frequency_rows_remain_missing(self):
        y = pd.DataFrame([dict(build_y='GRCh37', chrom_y=1, pos_y=100, ea_y='A', oa_y='G')])
        af = pd.DataFrame([dict(chr=1, pos37=100, ea_y='A', oa_y='G', eaf_y=f, source_y='source') for f in [.2, .3]])
        out = module.attach_frequencies(y, af)
        self.assertTrue(out.eaf_y.isna().all())
        self.assertEqual(out.iloc[0].eaf_attachment_status, 'ambiguous_frequency_rows')


class ArchiveExtractionTests(unittest.TestCase):
    def test_sequential_tar_members_and_x_metadata(self):
        payload = b'chr rsid pos a1 a2 beta se pval(-log10)\n01 rs1 100 A G .1 .01 3\n0X rs2 200 C T -.2 .03 4\n01 rs3 300 A T .01 .01 1\n'
        compressed = io.BytesIO()
        with tarfile.open(fileobj=compressed, mode='w:gz') as archive:
            for name in ['0001.txt', '0002.txt']:
                info = tarfile.TarInfo(name)
                info.size = len(payload)
                archive.addfile(info, io.BytesIO(payload))
        compressed.seek(0)
        with tempfile.TemporaryDirectory(dir=PAPER.parents[1]) as directory:
            with gzip.GzipFile(fileobj=compressed) as uncompressed:
                with tarfile.open(fileobj=uncompressed, mode='r|') as archive:
                    for member in archive:
                        meta = dict(member_basename=member.name, phenotype=member.name, trait='QSM test',
                                    family='QSM', source_y='fixture')
                        path = Path(directory) / (member.name + '.csv')
                        result = qmodule.parse_member(archive.extractfile(member), {('1', 100), ('X', 200)}, meta, path)
                        data = pd.read_csv(path)
                        self.assertEqual(result['source_rows'], 3)
                        self.assertEqual(len(data), 2)
                        self.assertTrue(data.n_y.isna().all())
                        self.assertEqual(data.ea_y.tolist(), ['G', 'T'])
                while uncompressed.read(1024):
                    pass


if __name__ == '__main__':
    unittest.main()
