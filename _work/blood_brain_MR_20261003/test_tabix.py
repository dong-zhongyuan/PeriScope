import pysam,urllib.request,json
u='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90003001-GCST90004000/GCST90003867/harmonised/GCST90003867.h.tsv.gz'
print('OPEN',flush=True)
t=pysam.TabixFile(u)
print('HEADER',list(t.header),flush=True)
print('CONTIG',t.contigs[:5],flush=True)
print('ROW',list(t.fetch('1',154453787,154453788)),flush=True)
