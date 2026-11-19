import gzip
import io
import unittest
from unittest.mock import patch
from pdv3.native_imaging import stream_oxford


class NativeTests(unittest.TestCase):
    def response(self,data):
        class Response:
            status_code=200
            headers={'Content-Length':str(len(data))}
            raw=io.BytesIO(data)
            def raise_for_status(self):pass
            def __enter__(self):return self
            def __exit__(self,*args):pass
        return Response()

    def test_full_stream_keeps_multialleles_and_missing_x_n(self):
        raw=b'chr rsid pos a1 a2 beta se pval(-log10)\n01 rs1 42 C T 0.2 0.1 1\n01 rs2 42 C A -0.1 0.2 0.1\n0X rs3 50 G A 0.1 0.1 0.5\n'
        data=gzip.compress(raw)
        with patch('pdv3.native_imaging.requests.get',return_value=self.response(data)):
            rows,audit=stream_oxford('https://source', [('1',42),('23',50)],
                                    {'phenotype':'test','n_discovery':'100'},expected_rows=3)
        self.assertEqual(len(rows),3)
        self.assertEqual(rows[-1]['n_y'],'')
        self.assertEqual(rows[0]['n_y'],'100')
        self.assertTrue(audit['gzip_crc_eof_verified'])
        self.assertEqual(audit['source_size_bytes'],len(data))

    def test_truncated_source_cannot_be_complete(self):
        raw=b'chr rsid pos a1 a2 beta se pval(-log10)\n01 rs1 42 C T 0.2 0.1 1\n'
        data=gzip.compress(raw)[:-5]
        with patch('pdv3.native_imaging.requests.get',return_value=self.response(data)):
            with self.assertRaises(EOFError):stream_oxford('https://source', [('1',42)],
                                        {'phenotype':'test','n_discovery':'100'},expected_rows=1)


if __name__=='__main__':unittest.main()
