"""Retrieve one public ZIP member by HTTP range, verify ZIP CRC and provenance."""
from pathlib import Path
import urllib.request,struct,zlib,json,hashlib,time
from concurrent.futures import ThreadPoolExecutor
dest=Path('/public/home/mengxl/dzy/pd_product_assets/raw/ma_published_reference_20261003');dest.mkdir(parents=True,exist_ok=True)
url='https://zenodo.org/api/records/19454950/files/data.zip/content'
op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:18792','http':'http://127.0.0.1:18792'}))
offset=911856379;compressed_size=1728902953;expected_size=1770471811
def request(start,end):
    opener=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:18792','http':'http://127.0.0.1:18792'}))
    r=opener.open(urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}'}),timeout=60)
    assert r.status==206 and r.headers.get('Content-Range','').startswith(f'bytes {start}-{end}/'),r.headers
    return r
with request(offset,offset+511) as r:header=r.read()
h=struct.unpack_from('<4s5H3I2H',header);assert h[0]==b'PK\x03\x04' and h[3]==8
name=header[30:30+h[-2]].decode();assert name=='data/processed_sn_obj.rds'
start=offset+30+h[-2]+h[-1];end=start+compressed_size-1
out=dest/'processed_sn_obj.rds';partial=dest/'processed_sn_obj.rds.partial'
if out.exists() and (dest/'download_provenance.json').exists():print('REFERENCE_ALREADY_DOWNLOADED');raise SystemExit(0)
dec=zlib.decompressobj(-15);crc=0;sha=hashlib.sha256();written=0;position=start
cache=dest/'compressed_member_parts';cache.mkdir(exist_ok=True)
def fetch_block(bounds):
    a,b=bounds
    path=cache/f'{a}-{b}.bin'
    part=path.with_suffix('.partial')
    if path.exists() and path.stat().st_size==b-a+1:return path
    for attempt in range(8):
        try:
            have=part.stat().st_size if part.exists() else 0
            if have<b-a+1:
                with request(a+have,b) as r,part.open('ab') as f:
                    while True:
                        data=r.read(min(1024*1024,b-a+1-have))
                        if not data:break
                        f.write(data);have+=len(data)
                        if have==b-a+1:break
            assert part.stat().st_size==b-a+1,(a,b,part.stat().st_size)
            part.rename(path);print('downloaded_range',a,b,flush=True);return path
        except Exception as exc:
            print('range_retry',a,attempt+1,repr(exc),flush=True)
            if attempt==7:raise
            time.sleep(2)
blocks=[(i,min(i+8*1024*1024-1,end)) for i in range(start,end+1,8*1024*1024)]
with partial.open('wb') as f:
    with ThreadPoolExecutor(max_workers=2) as pool:
        for part in pool.map(fetch_block,blocks):
            b=part.read_bytes()
            position+=len(b);data=dec.decompress(b);f.write(data);crc=zlib.crc32(data,crc);sha.update(data);written+=len(data)
            print(json.dumps(dict(compressed_bytes=position-start,total_compressed=compressed_size,output_bytes=written)),flush=True)
    data=dec.flush();f.write(data);crc=zlib.crc32(data,crc);sha.update(data);written+=len(data)
assert dec.eof and written==expected_size,(dec.eof,written)
assert crc&0xffffffff==h[6],(crc,h[6])
partial.rename(out)
(dest/'download_provenance.json').write_text(json.dumps(dict(
  source_record='https://doi.org/10.5281/zenodo.19454950',source_url=url,
  source_repository='https://github.com/dalhoomist/T-cell_and_glial_pathology_in_PD',
  member=name,member_header_offset=offset,compressed_range=[start,end],compressed_size=compressed_size,
  uncompressed_size=written,zip_crc32=f'{crc&0xffffffff:08x}',crc_verified=True,sha256=sha.hexdigest()),indent=2))
print('MA_REFERENCE_READY',written,sha.hexdigest(),flush=True)
for p in cache.iterdir():p.unlink()
cache.rmdir()
