"""Acquire all five predefined AGES assays with published MD5 and full CRC checks."""
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from s07_acquire_big40 import acquire
from pdv3.tabix_http import save_json

DATA=ROOT/'data/ages_public'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def one(record,instruments):
    accession=record['accession'];gene=record['gene']
    lead=instruments[instruments.gene.eq(gene)].iloc[0]
    chrom=int(lead.chrom);center=int(lead.pos37);low,high=center-1000000,center+1000000
    folder=DATA/accession
    md5s={v[1]:v[0] for l in (folder/'md5sum.txt').read_text().splitlines() if len(v:=l.split())==2}
    expected=md5s[record['file_name']]
    metadata=(folder/f'{accession}_buildGRCh37.tsv-meta.yaml').read_text()
    assert 'genome_assembly: GRCh37' in metadata and f'protein {gene}' in metadata
    n=int(re.search(r'sample_size: (\d+)',metadata).group(1))
    original_md5=re.search(r'data_file_md5sum: (\w+)',metadata).group(1)
    assert expected==original_md5 and record['file_name']==accession+'_buildGRCh37.tsv'
    output=DATA/'extracted'/f'{accession}.region.csv.gz';output.parent.mkdir(exist_ok=True)
    def validate(path):
        h=hashlib.md5()
        with path.open('rb') as f:
            for b in iter(lambda:f.read(1<<20),b''):h.update(b)
        if h.hexdigest()!=expected:raise ValueError(f'Published MD5 mismatch {accession}')
        rows=[];total=0
        with path.open('rt') as f:
            header=next(f).strip().split('\t')
            assert header=='variant_id p_value chromosome base_pair_location effect_allele other_allele effect_allele_frequency beta standard_error'.split()
            for line in f:
                v=line.strip().split('\t');total+=1
                if len(v)!=len(header):raise ValueError('Malformed source row')
                if v[2]!=str(chrom) or not low<=int(v[3])<=high:continue
                r=dict(zip(header,v));ea,oa=v[4:6]
                r.update(gene=gene,assay=record['assay'],accession=accession,build='GRCh37',
                         pos37=int(v[3]),chrom=chrom,ea_x=ea,oa_x=oa,eaf_x=float(v[6]),
                         beta_x=float(v[7]),se_x=float(v[8]),p_x=float(v[1]),n_x=n,source_row=total,
                         variant_key=f'{chrom}:{v[3]}:'+ '_'.join(sorted([ea,oa])) if len(ea)==len(oa)==1 and ea in 'ACGT' and oa in 'ACGT' and ea!=oa else None)
                if not np.isfinite([r['beta_x'],r['se_x'],r['eaf_x'],r['p_x']]).all() or r['se_x']<=0 or not 0<=r['p_x']<=1 or not 0<=r['eaf_x']<=1:
                    raise ValueError('Invalid regional statistics')
                rows.append(r)
        frame=pd.DataFrame(rows)
        assert len(frame)>0
        frame.to_csv(output,index=False,compression={'method':'gzip','mtime':0})
        return dict(published_md5_verified=True,md5=expected,uncompressed_text_eof_verified=True,source_rows=total,
                    regional_rows=len(frame),output=str(output.relative_to(ROOT)),output_sha256=sha(output),
                    region=dict(chrom=chrom,start=low,end=high,center=center,lead=lead.variant),sample_size=n,
                    original_source_md5=original_md5,
                    build_evidence='Original unformatted source with GRCh37 metadata, exact published original-source MD5 verified',
                    formatted_sidecar_build_conflict='Formatted files excluded after build metadata conflict and failed original checksum reconstruction; original source used directly',
                    allele_or_p_value_filter=None)
    path=DATA/'download'/record['file_name']
    receipt=acquire(record['url'],path,validator=validate)
    # A valid cached full file may need its derived region rebuilt.
    if not output.exists() or sha(output)!=receipt['validation']['output_sha256']:
        receipt['validation']=validate(path)
    result=dict(**record,complete=True,source_file=str(path.relative_to(ROOT)),download=receipt,
                metadata_sha256=sha(folder/f'{accession}_buildGRCh37.tsv-meta.yaml'),source_code_sha256=sha(Path(__file__)))
    save_json(folder/'ACQUISITION.json',result)
    print(accession,'COMPLETE',receipt['validation']['regional_rows'],flush=True)
    return result


def main():
    records=json.loads((DATA/'SOURCE_CATALOG.json').read_text())
    assert len(records)==5
    inst=pd.read_csv(ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv')
    with ThreadPoolExecutor(max_workers=3) as ex:results=list(ex.map(lambda r:one(r,inst),records))
    save_json(DATA/'MANIFEST.json',dict(complete=True,assays=5,selection='All AGES assays for three historical proteins, identified before regional analysis',datasets=results))


if __name__=='__main__':main()
