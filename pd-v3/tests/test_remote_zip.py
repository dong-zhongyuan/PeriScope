import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch
import zipfile
from pdv3.remote_zip import RangeSource, RangeReader, nested_reader


class MemorySource:
    def __init__(self, data):
        self.data, self.size = data, len(data)

    def get(self, start, length):
        return self.data[start:start+length]


def archive_bytes(compression=zipfile.ZIP_STORED):
    child, outer = io.BytesIO(), io.BytesIO()
    with zipfile.ZipFile(child, 'w', compression=zipfile.ZIP_STORED) as z:
        z.writestr('image.bin', b'uniquepixeldata' * 1000)
    with zipfile.ZipFile(outer, 'w') as z:
        z.writestr('prefix', b'other data')
        z.writestr('child.zip', child.getvalue(), compress_type=compression)
    return outer.getvalue()


class TestRemoteZip(unittest.TestCase):
    def test_nested_offset_and_crc(self):
        reader = RangeReader(MemorySource(archive_bytes()))
        with zipfile.ZipFile(reader) as outer:
            nested = nested_reader(outer, reader, outer.getinfo('child.zip'))
            with zipfile.ZipFile(nested) as inner:
                self.assertEqual(inner.read('image.bin'), b'uniquepixeldata' * 1000)

    def test_corrupt_inner_payload_is_rejected(self):
        data = bytearray(archive_bytes())
        data[data.index(b'uniquepixeldata')] ^= 1
        reader = RangeReader(MemorySource(data))
        with zipfile.ZipFile(reader) as outer:
            nested = nested_reader(outer, reader, outer.getinfo('child.zip'))
            with zipfile.ZipFile(nested) as inner:
                with self.assertRaises(zipfile.BadZipFile):
                    inner.read('image.bin')

    def test_compressed_child_and_out_of_range_seek_rejected(self):
        reader = RangeReader(MemorySource(archive_bytes(zipfile.ZIP_DEFLATED)))
        with zipfile.ZipFile(reader) as outer:
            with self.assertRaises(ValueError):
                nested_reader(outer, reader, outer.getinfo('child.zip'))
        with self.assertRaises(ValueError): reader.seek(-1)
        with self.assertRaises(ValueError): reader.seek(1, 2)

    @staticmethod
    def response(data=b'abc', content_range='bytes 0-2/10', etag='version1'):
        r=Mock(); r.status_code=206
        r.headers={'Content-Range':content_range,'ETag':etag}
        buf=io.BytesIO(data)
        r.raw=Mock(); r.raw.read=Mock(side_effect=lambda n,decode_content=False:buf.read(n))
        r.raise_for_status=Mock()
        r.__enter__=Mock(return_value=r); r.__exit__=Mock(return_value=False)
        return r

    def test_wrong_http_range_and_short_payload_rejected(self):
        with tempfile.TemporaryDirectory() as tmp, patch('pdv3.remote_zip.time.sleep'):
            source=RangeSource('https://example.org/archive',10,tmp)
            source.session.get=Mock(return_value=self.response(content_range='bytes 1-3/10'))
            with self.assertRaises(ValueError): source.get(0,3)
            source.session.get=Mock(return_value=self.response(data=b'ab'))
            with self.assertRaises(ValueError): source.get(0,3)

    def test_missing_validators_require_pinned_metadata_and_cache_is_checked(self):
        with tempfile.TemporaryDirectory() as tmp, patch('pdv3.remote_zip.time.sleep'):
            source=RangeSource('https://example.org/archive',10,tmp)
            source.session.get=Mock(return_value=self.response(etag=None))
            with self.assertRaises(ValueError): source.get(0,3)
            source=RangeSource('https://example.org/archive',10,tmp,published_identity={'version':'1','checksum':'md5:example'})
            source.session.get=Mock(return_value=self.response(etag=None))
            self.assertEqual(source.get(0,3),b'abc')
            self.assertEqual(source.get(0,3),b'abc')
            self.assertEqual(source.session.get.call_count,1)
            next(Path(tmp).glob('*.bin')).write_bytes(b'bad')
            with self.assertRaises(ValueError): source.get(0,3)

    def test_http_identity_change_rejected(self):
        with tempfile.TemporaryDirectory() as tmp, patch('pdv3.remote_zip.time.sleep'):
            source=RangeSource('https://example.org/archive',10,tmp)
            source.session.get=Mock(return_value=self.response())
            source.get(0,3)
            source.session.get=Mock(return_value=self.response(content_range='bytes 3-5/10',etag='version2'))
            with self.assertRaises(ValueError): source.get(3,3)


if __name__ == '__main__':
    unittest.main()
