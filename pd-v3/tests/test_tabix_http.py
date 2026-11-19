"""Scientific input integrity: point boundaries, multiple alleles, CRC and HTTP identity."""
import gzip
import struct
import unittest
import zlib

from pdv3.tabix_http import (bgzf_blocks, parse_tbi, point_chunks, merge_intervals,
                            split_ranges, HTTPStore, RemotePointTable)


def block(raw):
    comp = zlib.compressobj(wbits=-15)
    payload = comp.compress(raw) + comp.flush()
    size = 18 + len(payload) + 8
    return (b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00' +
            struct.pack('<H', size-1) + payload + struct.pack('<II', zlib.crc32(raw), len(raw)))


class PointTests(unittest.TestCase):
    def test_bgzf_crc_rejects_corruption(self):
        valid=block(b'abc\ndef\n')
        self.assertEqual(list(bgzf_blocks(valid)),[(0,b'abc\ndef\n')])
        bad=valid[:-8]+bytes([valid[-8]^1])+valid[-7:]
        with self.assertRaises(zlib.error):list(bgzf_blocks(bad))

    def test_bgzf_final_partial_is_not_decoded(self):
        valid=block(b'abc\n')
        self.assertEqual(list(bgzf_blocks(valid+valid[:20])),[(0,b'abc\n')])

    def test_16kb_one_based_boundary(self):
        refs={'1':({4681:[(5,20)],4682:[(20,40)]},[5,20])}
        self.assertEqual(point_chunks(refs,'1',16384),[(5,20)])
        self.assertEqual(point_chunks(refs,'1',16385),[(20,40)])
        self.assertEqual(point_chunks(refs,'23',100),[])
        with self.assertRaises(ValueError):point_chunks(refs,'1',0)

    def test_tbi_format_and_truncation(self):
        raw=b'TBI\x01'+struct.pack('<8i',1,0,1,2,0,35,1,2)+b'1\0'
        raw+=struct.pack('<iIiQQiQ',1,4681,1,5,20,1,5)
        self.assertEqual(point_chunks(parse_tbi(gzip.compress(raw)),'1',3),[(5,20)])
        bad=raw[:8]+struct.pack('<i',65536)+raw[12:]
        with self.assertRaises(ValueError):parse_tbi(gzip.compress(bad))
        with self.assertRaises(struct.error):parse_tbi(gzip.compress(raw[:-5]))

    def test_http_rejects_wrong_range_and_short_body(self):
        info={'content_type':'application/octet-stream','content_range':'bytes 10-12/100'}
        self.assertEqual(split_ranges(b'abc',info,[(10,12)])[1],100)
        with self.assertRaises(ValueError):split_ranges(b'ab',info,[(10,12)])
        with self.assertRaises(ValueError):split_ranges(b'abc',info,[(11,13)])
        with self.assertRaises(ValueError):HTTPStore.check_identity({'etag':'new'},{'etag':'old'})

    def test_multipart_missing_or_duplicate_part_rejected(self):
        part=b'--test\r\nContent-Type: application/octet-stream\r\nContent-Range: bytes 10-12/100\r\n\r\nabc\r\n'
        info={'content_type':'multipart/byteranges; boundary=test'}
        good=part+b'--test--\r\n'
        self.assertEqual(split_ranges(good,info,[(10,12)])[0],{(10,12):b'abc'})
        with self.assertRaises(ValueError):split_ranges(good,info,[(10,12),(20,22)])
        with self.assertRaises(ValueError):split_ranges(part+part+b'--test--\r\n',info,[(10,12)])

    def test_record_crossing_blocks_and_all_alleles(self):
        header=b'chromosome\tbase_pair_location\teffect_allele\tbeta\n'
        record1=b'1\t16384\tA\t0.1\n'
        record2=b'1\t16384\tG\t0.2\n'
        record3=b'1\t16385\tT\t0.3\n'
        first=block(header+record1[:7]); second=block(record1[7:]+record2+record3)
        raw=first+second
        end=(len(first)<<16)+len(record1[7:]+record2+record3)
        table=RemotePointTable.__new__(RemotePointTable)
        table.url='https://test/table';table.total=len(raw);table.identity={}
        table.index_bytes=b'test';table.header=header.decode().strip().split('\t')
        table.refs={'1':({4681:[(len(header),end)],4682:[(len(header),end)]},[len(header),len(header)])}
        class Store:
            def get(self,url,ranges,identity):
                self_test.assertEqual(len(ranges),1)
                a,b=ranges[0]
                return raw[a:b+1],{'content_type':'application/octet-stream','content_range':f'bytes {a}-{b}/{len(raw)}'}
        self_test=self;table.store=Store()
        rows,info=table.fetch_many([('1',16384)])
        self.assertEqual([r['effect_allele'] for r in rows],['A','G'])
        self.assertEqual(info['returned_rows'],2)


if __name__=='__main__':unittest.main()
