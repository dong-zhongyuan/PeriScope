"""Reproduce the public CD40 source and pinned coloc function downloads."""
import hashlib
import json
import os
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/scallop_cd40'
URL='https://zenodo.org/records/2615265/files/CD40.txt.gz?download=1'
SIZE=347619368
MD5='5b909d2ce35f144c74fd73271a3a0287'


def hashes(path):
    md5=hashlib.md5();sha=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):md5.update(b);sha.update(b)
    return md5.hexdigest(),sha.hexdigest()


def main():
    target=DATA/'download/CD40.txt.gz';target.parent.mkdir(parents=True,exist_ok=True)
    if not target.exists():
        record=requests.get('https://zenodo.org/api/records/2615265',timeout=30)
        record.raise_for_status();obj=record.json()
        source=next(r for r in obj['files'] if r['key']==target.name)
        assert source['size']==SIZE and source['checksum']=='md5:'+MD5
        (DATA/'ZENODO_RECORD.json').write_text(json.dumps(obj,indent=2)+'\n')
        part=target.with_suffix(target.suffix+'.part')
        with requests.get(URL,stream=True,timeout=(20,90),headers={'Accept-Encoding':'identity'}) as r:
            r.raise_for_status()
            if r.status_code!=200:raise ValueError('Expected full response')
            with part.open('wb') as f:
                for b in r.iter_content(1<<20):f.write(b)
                f.flush();os.fsync(f.fileno())
        assert part.stat().st_size==SIZE and hashes(part)[0]==MD5
        part.replace(target)
    md5,sha=hashes(target)
    assert target.stat().st_size==SIZE and md5==MD5
    manifest=DATA/'DOWNLOAD_MANIFEST.json'
    if manifest.exists():assert json.loads(manifest.read_text())['sha256']==sha
    else:
        manifest.write_text(json.dumps(dict(url=URL,size_bytes=SIZE,published_md5=MD5,verified_md5=md5,
             sha256=sha,source_record='https://zenodo.org/records/2615265',gzip_crc='verified by subsequent full regional stream'),indent=2)+'\n')
    ref=json.loads((DATA/'COLOC_REFERENCE.json').read_text())
    for info in ref['files'].values():
        path=ROOT/info['path'];path.parent.mkdir(parents=True,exist_ok=True)
        if not path.exists():
            r=requests.get(info['url'],timeout=30);r.raise_for_status()
            assert hashlib.sha256(r.content).hexdigest()==info['sha256'];path.write_bytes(r.content)
        assert hashlib.sha256(path.read_bytes()).hexdigest()==info['sha256']
    print('Public CD40 source and pinned coloc reference verified')


if __name__=='__main__':main()
