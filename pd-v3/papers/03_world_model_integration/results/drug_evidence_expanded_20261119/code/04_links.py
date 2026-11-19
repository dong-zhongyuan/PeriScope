#!/usr/bin/env python
"""Candidate ligand->receptor and complex-partner links for every antigen gene, with drugs on the linked gene."""
import sys; sys.path.insert(0,'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119/code')
from common import *
import collections
R=O+'/raw/'; I=O+'/interim/'
AG=antigens(); ag_genes=sorted({g for a in AG for g in a['genes']}); by_gene=collections.defaultdict(list)
for a in AG:
    for g in a['genes']: by_gene[g].append(a['target'])
gi=pd.read_csv(R+'cellphonedb/gene_input.csv').fillna(''); u2g={}
for r in gi.itertuples():
    if r.hgnc_symbol: u2g.setdefault(r.uniprot,r.hgnc_symbol)
cx=pd.read_csv(R+'cellphonedb/complex_input.csv',dtype=str).fillna('')
C={}
for r in cx.to_dict('records'):
    C[r['complex_name']]=[u2g.get(r[f'uniprot_{i}'],r[f'uniprot_{i}']) for i in range(1,6) if r[f'uniprot_{i}']]
def members(p): return C.get(p) or ([u2g[p]] if p in u2g else [])
L=collections.defaultdict(lambda:dict(src=set(),refs=set(),ty=set()))
ip=pd.read_csv(R+'cellphonedb/interaction_input.csv',dtype=str).fillna('')
for r in ip[ip.directionality.eq('Ligand-Receptor')].to_dict('records'):
    for rg in members(r['partner_b']):
        if rg in by_gene:
            for lg in members(r['partner_a']):
                x=L[(rg,lg)]; x['src'].add('CellPhoneDB'); x['refs'].add(r['source'][:70])
                if r['modulatory_effect']: x['ty'].add('CPDB:'+r['modulatory_effect'])
g=pd.read_csv(R+'gtopdb/interactions.csv',skiprows=1,low_memory=False,dtype=str).fillna('')
e=g[g.Endogenous.str.lower().isin(['true','t'])&g['Target Species'].str.lower().eq('human')&g['Ligand Gene Symbol'].ne('')]
for r in e.to_dict('records'):
    for rg in r['Target Gene Symbol'].split('|'):
        if rg in by_gene:
            for lg in re.split(r'[|;]',r['Ligand Gene Symbol']):
                lg=lg.strip().upper()
                x=L[(rg,lg)]; x['src'].add('GtoPdb_endogenous'); x['ty'].add(('GTOP:'+r['Type']+'/'+r['Action']).upper())
                if r['PubMed ID']: x['refs'].add('PMID:'+r['PubMed ID'][:40])
o=pd.read_csv(R+'omnipath/interactions.tsv',sep='\t',low_memory=False,dtype=str).fillna('')
ic=pd.read_csv(R+'omnipath/intercell_generic.tsv',sep='\t',low_memory=False,dtype=str).fillna('')
lig=set(ic[ic.category.eq('ligand')].genesymbol); rec=set(ic[ic.category.eq('receptor')].genesymbol)
o2=o[(o.ligrecextra.str.lower().eq('true')|o.sources.str.contains('CellChatDB|CellPhoneDB|Guide2Pharma|ICELLNET|CellTalkDB|connectomeDB|Ramilowski',regex=True))&o.target_genesymbol.isin(by_gene)&o.source_genesymbol.isin(lig)&o.target_genesymbol.isin(rec)]
for r in o2.to_dict('records'):
    x=L[(r['target_genesymbol'],r['source_genesymbol'])]; x['src'].add('OmniPath'); s=[y for y in r['sources'].split(';')]; x['ty'].add(f'OP_n={len(s)}')
    if r['is_inhibition']=='True' and r['is_stimulation']!='True': x['ty'].add('OP:inhibition')
    if r['references']: x['refs'].update(y for y in r['references'].split(';')[:3])
rows=[dict(receptor_gene=k[0],ligand_gene=k[1],antigens=';'.join(by_gene[k[0]]),link_sources=';'.join(sorted(v['src'])),ligand_type=';'.join(sorted(v['ty']))[:120],refs=';'.join(sorted(v['refs']))[:300]) for k,v in L.items() if k[0]!=k[1]]
LL=pd.DataFrame(rows); LL.to_csv(I+'ligand_link_candidates.csv',index=False); print('ligand links',len(LL),LL.receptor_gene.nunique(),LL.ligand_gene.nunique())
# complexes
P=collections.defaultdict(lambda:dict(src=set(),cx=set()))
for n,m in C.items():
    for g0 in m:
        if g0 in by_gene:
            for h in m:
                if h!=g0: x=P[(g0,h)]; x['src'].add('CellPhoneDB_complex'); x['cx'].add(n[:40])
cp=pd.read_csv(R+'complexportal/9606.tsv',sep='\t',dtype=str).fillna('')
idn=cp.columns[4]; acc2=gi.set_index('uniprot').hgnc_symbol.to_dict()
for r in cp.to_dict('records'):
    acc=[re.sub(r'[-_].*$','',x.split('(')[0]) for x in r[idn].split('|')]; gs=sorted({acc2[a] for a in acc if acc2.get(a)})
    if len(gs)>12: continue
    for g0 in gs:
        if g0 in by_gene:
            for h in gs:
                if h!=g0: x=P[(g0,h)]; x['src'].add('ComplexPortal'); x['cx'].add(r['#Complex ac']+':'+r['Recommended name'][:40])
EV=pd.read_csv(I+'evidence_long.csv.gz',low_memory=False,dtype={'target_genes':str}).fillna({'target_genes':'','action':'','target_type':''})
cxr=EV[EV.target_type.isin(['protein complex','protein complex group','protein-protein interaction','chimeric protein'])]
for r in cxr[['target_genes','target_name','source']].drop_duplicates().itertuples():
    gs=norm_list(r.target_genes)
    if len(gs)>12: continue
    for g0 in gs:
        if g0 in by_gene:
            for h in gs:
                if h!=g0: x=P[(g0,h)]; x['src'].add('drug_target_complex:'+r.source.split('_')[0]); x['cx'].add(str(r.target_name)[:40])
PP=pd.DataFrame([dict(antigen_gene=k[0],partner_gene=k[1],antigens=';'.join(by_gene[k[0]]),link_sources=';'.join(sorted(v['src'])),complexes=';'.join(sorted(v['cx']))[:200]) for k,v in P.items()])
PP.to_csv(I+'complex_link_candidates.csv',index=False); print('complex links',len(PP),PP.antigen_gene.nunique(),PP.partner_gene.nunique())
extra=sorted((set(LL.ligand_gene)|set(PP.partner_gene))-set(ag_genes)); extra=[x for x in extra if re.match(r'^[A-Z0-9][A-Z0-9-]+$',x)]
open(O+'/code/genes_linked.txt','w').write('\n'.join(extra)+'\n'); print('extra genes',len(extra))
