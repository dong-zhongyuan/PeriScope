"""Extract public EUR reference genotypes, retaining ALT dosage and SNP identity.

This is an external LD reference, not cohort-matched LD. No association or MAF
filter is applied at extraction. Remote BGZF/VCF access uses pinned pysam.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import requests

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'.deps/python'))
sys.path.insert(0,str(ROOT))
import pysam
from pdv3.tabix_http import save_json

BASE='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
PANEL='integrated_call_samples_v3.20130502.ALL.panel'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene',choices=['CD22','CD40','TFRC'],required=True)
    args=parser.parse_args()
    inst=pd.read_csv(ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv')
    lead=inst[inst.gene.eq(args.gene)]
    assert len(lead)==1
    lead=lead.iloc[0];chrom=str(int(lead.chrom));center=int(lead.pos37)
    start,end=center-1000000,center+1000000
    out=ROOT/'data/ld_1000g';out.mkdir(exist_ok=True)
    ex=out/'extracted'/args.gene;ex.mkdir(parents=True,exist_ok=True)
    manifest=ex/'MANIFEST.json'
    if manifest.exists():
        m=json.loads(manifest.read_text())
        if m['complete'] and all(sha(ROOT/p)==v for p,v in m['outputs'].items()):
            print('Verified cached',args.gene,flush=True);return
        raise ValueError('Existing reference extraction failed hash validation')
    r=requests.get(BASE,timeout=40);r.raise_for_status()
    names=re.findall(r'href="([^"]+)"',r.text)
    files=[n for n in names if n.startswith(f'ALL.chr{chrom}.') and n.endswith('.genotypes.vcf.gz')]
    assert len(files)==1,files
    name=files[0];url=BASE+name
    for n in [name+'.tbi',PANEL]:
        p=out/n
        if not p.exists():
            r=requests.get(BASE+n,timeout=40);r.raise_for_status();p.write_bytes(r.content)
    panel=pd.read_csv(out/PANEL,sep='\t',usecols=['sample','pop','super_pop'])
    eur=panel[panel.super_pop.eq('EUR')];assert len(eur)==503 and eur['sample'].is_unique
    head=requests.head(url,timeout=40);head.raise_for_status()
    identity={k:head.headers.get(k) for k in ['ETag','Last-Modified','Content-Length']}
    counts=Counter();rows=[];dosages=[]
    with pysam.VariantFile(url,index_filename=str(out/(name+'.tbi'))) as source:
        source.subset_samples(eur['sample'].tolist())
        samples=list(source.header.samples)
        assert set(samples)==set(eur['sample'])
        for variant in source.fetch(chrom,start-1,end):
            if not start<=variant.pos<=end:continue
            counts['source_region_variants']+=1
            if len(variant.ref)!=1 or not variant.alts or len(variant.alts)!=1 or len(variant.alts[0])!=1 or any(a not in 'ACGT' for a in [variant.ref,variant.alts[0]]):
                counts['not_biallelic_snv']+=1;continue
            dosage=[]
            for sample in samples:
                gt=variant.samples[sample]['GT']
                if len(gt)!=2 or any(a is None for a in gt):dosage.append(-1)
                elif any(a not in [0,1] for a in gt):raise ValueError('Non-biallelic genotype in biallelic record')
                else:dosage.append(sum(gt))
            row=dict(chromosome=chrom,pos37=variant.pos,ref=variant.ref,alt=variant.alts[0],source_variant_id=variant.id,
                     variant_key=f'{chrom}:{variant.pos}:'+ '_'.join(sorted([variant.ref,variant.alts[0]])),
                     genotype_coding='ALT dosage',missing_calls=dosage.count(-1))
            rows.append(row);dosages.append(dosage)
    after=requests.head(url,timeout=40);after.raise_for_status()
    assert identity=={k:after.headers.get(k) for k in identity},'Reference HTTP identity changed during extraction'
    variants=pd.DataFrame(rows)
    if variants.variant_key.duplicated().any():raise ValueError('Duplicate reference SNP identity')
    g=np.asarray(dosages,dtype=np.int8)
    assert g.shape==(len(rows),503)
    variants['alt_frequency']=np.where((g>=0).sum(axis=1)>0,np.maximum(g,0).sum(axis=1)/(2*np.maximum((g>=0).sum(axis=1),1)),np.nan)
    np.savez_compressed(ex/'alt_dosages.npz',dosages=g,samples=np.asarray(samples),variant_keys=variants.variant_key.to_numpy(dtype=str))
    variants.to_csv(ex/'variants.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    eur.set_index('sample').loc[samples].reset_index().to_csv(ex/'sample_population.csv',index=False)
    outputs={str(p.relative_to(ROOT)):sha(p) for p in ex.iterdir() if p.is_file() and p.name!='MANIFEST.json'}
    save_json(manifest,dict(complete=True,gene=args.gene,chromosome=chrom,start=start,end=end,center=center,
         historical_lead=lead.variant,genome_assembly='GRCh37',url=url,http_identity=identity,
         index_sha256=sha(out/(name+'.tbi')),panel_sha256=sha(out/PANEL),
         sample_count=503,population_counts=eur['pop'].value_counts().to_dict(),counts=dict(counts),
         retained_biallelic_snvs=len(variants),missing_calls=int((g<0).sum()),
         maf_filter=None,association_filter=None,cohort_matched=False,
         whole_remote_file_checksum_verified=False,local_output_integrity='SHA256; HTSlib BGZF decompression and exact interval/allele checks',
         pysam=pysam.__version__,htslib=pysam.__samtools_version__,outputs=outputs,source_code_sha256=sha(Path(__file__))))
    print(json.dumps(dict(gene=args.gene,variants=len(variants),samples=503,counts=dict(counts)),indent=2),flush=True)


if __name__=='__main__':main()
