"""Stream an official Oxford gzip to EOF and extract supplied GRCh37 positions."""
import collections
import csv
import gzip
import hashlib
import io
import math

import requests

from .tabix_http import digest

HEADER=['chr','rsid','pos','a1','a2','beta','se','pval(-log10)']
FIELDS=['chromosome','base_pair_location','effect_allele','other_allele','beta',
        'standard_error','p_value','neg_log10_p_value','variant_id','rsid',
        'phenotype','source_y','build_y','n_y','n_y_status']


class HashReader:
    def __init__(self,raw):
        self.raw=raw;self.h=hashlib.sha256();self.n=0
    def read(self,size=-1):
        data=self.raw.read(size);self.h.update(data);self.n+=len(data);return data


def stream_oxford(url,points,meta,expected_rows=17103079):
    wanted={(str(c).lstrip('0').replace('X','23'),int(p)) for c,p in points}
    with requests.get(url,headers={'Accept-Encoding':'identity','User-Agent':'pd-v3-full-catalog/1.0'},
                      stream=True,timeout=(20,90)) as response:
        response.raise_for_status()
        if response.status_code!=200:raise ValueError('Expected a complete native source response')
        expected=int(response.headers['Content-Length'])
        raw=HashReader(response.raw)
        rows=[];total=0;chromosomes=collections.Counter()
        with gzip.GzipFile(fileobj=raw) as zipped, io.TextIOWrapper(zipped,encoding='ascii') as source:
            if source.readline().split()!=HEADER:raise ValueError('Unexpected Oxford columns')
            for line in source:
                v=line.split()
                if len(v)!=8:raise ValueError('Invalid native record width')
                c,rs,p,oa,ea,b,s,lp=v
                c=c.lstrip('0').replace('X','23');p=int(p)
                b,s,lp=float(b),float(s),float(lp)
                if not (c in {str(i) for i in range(1,24)} and p>0 and math.isfinite(b)
                        and math.isfinite(s) and s>0 and math.isfinite(lp) and lp>=0 and oa!=ea):
                    raise ValueError('Invalid native coordinate/estimate')
                total+=1;chromosomes[c]+=1
                if (c,p) in wanted:
                    rows.append(dict(zip(FIELDS,[c,p,ea,oa,b,s,10**(-lp),lp,rs,rs,
                                meta['phenotype'],url,'GRCh37',meta['n_discovery'] if c!='23' else '',
                                'catalog_autosomal_exact' if c!='23' else 'X_not_assigned'])))
        # gzip has consumed every member and verified CRC at EOF.
        if raw.n!=expected or total!=expected_rows:
            raise ValueError(f'Native source incomplete: bytes={raw.n}/{expected}; rows={total}/{expected_rows}')
        receipt={'url':url,'source_size_bytes':raw.n,'source_sha256':raw.h.hexdigest(),
                 'etag':response.headers.get('ETag'),'last_modified':response.headers.get('Last-Modified'),
                 'source_rows':total,'chromosome_rows':dict(chromosomes),'gzip_crc_eof_verified':True,
                 'source_retention':'streamed with SHA256; full compressed payload not retained',
                 'all_records_numerically_validated':True}
    return rows,receipt
