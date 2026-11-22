#!/usr/bin/env python
"""Stream the Oxford QSM/T2* discovery tarball once and keep the TFRC window (chr3 195.3-196.3 Mb,
GRCh37) for nigral QSM right (0016) and nigral T2* left/right (0033, 0034). Header of every member:
chr rsid pos a1 a2 beta se pval(-log10); a2 is the effect allele."""
import gzip, os, sys, tarfile
TAR = '/public/home/mengxl/dzy/pd-v3/data/oxford_qsm/download/QSMandT2s_IDPs.tar.gz'
OUT = sys.argv[1]
WANT = {'0015.txt', '0016.txt', '0033.txt', '0034.txt'}
os.makedirs(OUT, exist_ok=True)
done = set()
with tarfile.open(TAR, 'r|gz') as tf:
    for m in tf:
        name = os.path.basename(m.name)
        if name not in WANT:
            continue
        fh = tf.extractfile(m)
        out = gzip.open(os.path.join(OUT, name.replace('.txt', '_chr3_tfrc.tsv.gz')), 'wt')
        hdr = fh.readline().decode().split()
        out.write('\t'.join(hdr) + '\n')
        n = 0
        for raw in fh:
            f = raw.split()
            if f[0] in (b'3', b'03') and 195300000 <= int(f[2]) <= 196300000:
                out.write(raw.decode().replace(' ', '\t'))
                n += 1
            elif n and f[0] not in (b'3', b'03'):
                break
        out.close()
        done.add(name)
        print(name, n, 'rows', flush=True)
        if done == WANT:
            break
