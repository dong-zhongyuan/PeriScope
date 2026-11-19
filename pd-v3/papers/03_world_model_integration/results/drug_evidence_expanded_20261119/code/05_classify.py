#!/usr/bin/env python
"""Classify every antigen-drug-route record, then write drug_records.csv.gz and antigen_drug_summary.csv."""
import sys; sys.path.insert(0,'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119/code')
from common import *
import collections, os
I=O+'/interim/'; C=O+'/curation/'
EV=pd.read_csv(I+'evidence_long.csv.gz',low_memory=False).fillna({'target_genes':'','action':'','target_type':'','mechanism':'','refs':'','link':'','extra':'','act_type':'','target_name':''})
DM=pd.read_csv(I+'drug_master.csv.gz',low_memory=False).set_index('drug_key')
V=pd.read_csv(S6+'verified_drug_target_mechanisms.csv'); AM=pd.read_csv(S6+'all_drug_target_mechanisms.csv'); ANN=pd.read_csv(S6+'all_antigen_annotations.csv').set_index('protein')
LP=pd.read_csv(C+'ligand_pairs.csv').fillna(''); CP=pd.read_csv(C+'complex_pairs.csv').fillna('')
def opt(n,cols):
    p=C+n
    return pd.read_csv(p,dtype=str).fillna('') if os.path.exists(p) else pd.DataFrame(columns=cols)
MAN=opt('manual_records.csv',['antigen','acted_gene','drug','chembl_id','route_class','direction','level','max_phase','link_quality','mechanism','evidence','note'])
ABE=opt('antibody_effect.csv',['drug','antigens','routes','decision','excerpt'])
EXC=opt('route_exceptions.csv',['receptor_gene','linked_gene','drug','effect','reason'])
LAB={}
if os.path.exists(O+'/raw/openfda/labels.jsonl'):
    for l in open(O+'/raw/openfda/labels.jsonl'):
        o=json.loads(l); n=re.sub('[^A-Z0-9]','',str(o['name']).upper())
        nm=lambda h:[re.sub('[^A-Z0-9]','',str(g).upper()) for g in (h.get('generic') or [])+(h.get('brand') or [])]
        hit=[h for h in o['hits'] if any(n and x.startswith(n) for x in nm(h))] or [h for h in o['hits'] if any(n and n[:12] in x for x in nm(h))]
        if hit: LAB[o['drug_key']]=hit[0]
COMPLEX={'protein complex','protein complex group','protein-protein interaction','chimeric protein'}
FAMILY={'protein family','selectivity group'}
NUC={'protein nucleic-acid complex','nucleic-acid'}
CD3G={'CD3E','CD3D','CD3G','CD247'}
s2k=EV[EV.source.eq('ChEMBL_REST_20261003')].drop_duplicates('src_drug').set_index('src_drug').drug_key.to_dict()
ver=collections.defaultdict(set)
for r in V.itertuples(): ver[r.protein].add(s2k.get(r.molecule_chembl_id,r.molecule_chembl_id))
SPECIAL={}
for r in AM.itertuples():
    k=s2k.get(r.molecule_chembl_id,r.molecule_chembl_id); m=r.match_status
    if m in('pan_PTPRC_drug_isoform_specificity_unresolved','Fc_receptor_engagement_not_independent_CD16_intervention','CD44v6_epitope_specificity_unresolved','CR1_functional_analogue_not_direct_CR1_targeting','HLA_DR_allele_epitope_requires_resolution'): SPECIAL[(r.protein,k)]=('family_member',m)
    elif m=='incorrect_source_target_annotation': SPECIAL[(r.protein,k)]=('drop',m)
    elif m=='ligand_trap_not_direct_TREM1_binding': SPECIAL[(r.protein,k)]=('ligand_route',m)
names=collections.defaultdict(set)
for k,n in DM.names.dropna().items():
    for x in n.split('|'): names[re.sub('[^A-Z0-9]','',x.upper())].add(k)
def res(chembl,name):
    if chembl and chembl in DM.index: return chembl
    if chembl:
        h=EV[EV.src_drug.eq(chembl)].drug_key
        if len(h): return h.iloc[0]
    h=names.get(re.sub('[^A-Z0-9]','',str(name).upper()),set())
    return next(iter(h)) if len(h)==1 else 'MANUAL:'+re.sub('[^A-Z0-9]','',str(name).upper())
EV['genes']=EV.target_genes.str.split(';')
d=[action_direction(a) for a in EV.action]; EV['dirn']=[x[0] for x in d]; EV['sil']=[x[1] for x in d]
def tier(r):
    s=r['source']
    if s.startswith('Open') or s.startswith('ChEMBL'): return 2 if r['primary'] else 4
    if s.startswith('Gto') or s.startswith('Drug'): return 3 if r['primary'] else 4
    return 5
recs=EV.to_dict('records')
idx=collections.defaultdict(list)
for i,r in enumerate(recs):
    for g in r['genes']:
        if g: idx[g].append(i)
cd3drug={r['drug_key'] for g in CD3G for r in (recs[i] for i in idx[g]) if r['source'].startswith(('Open','ChEMBL'))}
lp=collections.defaultdict(list); cp=collections.defaultdict(list)
for r in LP.itertuples(): lp[r.receptor_gene].append((r.ligand_gene,r.quality,r.note))
for r in CP.itertuples(): cp[r.antigen_gene].append((r.partner_gene,r.quality,r.note))
cpd={(r.antigen_gene,r.partner_gene):r.quality for r in CP.itertuples()}
exc={}
for r in EXC.itertuples(): exc[(r.receptor_gene,r.linked_gene,r.drug.upper())]=(r.effect,r.reason)
abe={}
for r in ABE.itertuples():
    for a_ in r.antigens.split(';'):
        for rt in r.routes.split(';'): abe[(a_,res('',r.drug),rt)]=r
AG=antigens(); out=collections.defaultdict(list)
UPLIG={'EXOGENOUS PROTEIN','EXOGENOUS GENE'}
for a in AG:
    p=a['target']; G=set(a['genes'])
    def put(k,r,route,gene,link='',q='',note='',dirn=None):
        out[(p,k,route)].append(dict(r,_gene=gene,_link=link,_q=q,_note=note,_d=r['dirn'] if dirn is None else dirn))
    for g in G:
        for i in idx[g]:
            r=recs[i]; k=r['drug_key']; tt=r['target_type']; S=set(r['genes'])
            sp=SPECIAL.get((p,k))
            if sp and sp[0]=='drop': continue
            if sp and sp[0]=='family_member': put(k,r,'family_member',g,note=sp[1]); continue
            if sp and sp[0]=='ligand_route': put(k,r,'ligand_route',g,link='TREM1 ligand',q='canonical',note=sp[1],dirn='down'); continue
            if p in('CD45RA','CD45RB','CD45RO'): put(k,r,'family_member',g,note='pan_PTPRC_drug_isoform_specificity_unresolved'); continue
            if g.startswith('FCGR') and str(DM.drug_type.get(k,'')).startswith('Antibody'): put(k,r,'family_member',g,note='Fc_receptor_engagement_not_independent_intervention_on_this_antigen'); continue
            if r['sil'] and (tt in NUC or tt=='single protein'): put(k,r,'gene_silencing',g); continue
            if tt in NUC: continue
            if tt in('single protein','ligand_target'): put(k,r,'direct',g)
            elif tt in COMPLEX:
                oth=S-{g}
                if S<=G or k in ver[p]: put(k,r,'direct',g)
                else:
                    ok=[o for o in oth if (g,o) in cpd]
                    if ok:
                        q='obligate' if all(cpd.get((g,o))=='obligate' for o in oth) else 'partial'
                        put(k,r,'complex_partner',g,link=';'.join(sorted(ok)),q=q,note='drug target is the complex containing the antigen; contacted subunit not resolved')
                    else: put(k,r,'family_member',g,note='complex member; contacted subunit not resolved')
            elif tt in FAMILY: put(k,r,'family_member',g,note='pan-family mechanism record')
        for L,q,nt in lp.get(g,[]):
            for i in idx.get(L,[]):
                r=recs[i]; k=r['drug_key']; tt=r['target_type']
                if not(tt in('single protein','ligand_target') or tt in COMPLEX): continue
                if set(r['genes'])&G and tt in COMPLEX: continue
                dn=r['dirn']; ac=r['action'].upper()
                if dn=='up' and not any(x in ac for x in UPLIG): dn='unknown'
                if q=='ambiguous': dn='unknown'
                ex=exc.get((g,L,str(DM.name.get(k,'')).upper()))
                if ex:
                    if ex[0]=='exclude': continue
                    if ex[0]=='unknown': dn='unknown'
                put(k,r,'ligand_route',L,link=f'{g}<-{L}',q=q,note=(ex[1] if ex else nt),dirn=dn)
        for P,q,nt in cp.get(g,[]):
            for i in idx.get(P,[]):
                r=recs[i]; k=r['drug_key']; tt=r['target_type']
                if tt!='single protein': continue
                ex=exc.get((g,P,str(DM.name.get(k,'')).upper()))
                if ex and ex[0]=='exclude': continue
                put(k,r,'complex_partner',P,link=f'{g}+{P}',q=q,note=(ex[1] if ex else nt),dirn='unknown' if ex else None)
MR=[]
for r in MAN.to_dict('records'):
    k=res(r['chembl_id'],r['drug'])
    if k not in DM.index:
        DM.loc[k,['name','drug_type']]=[r['drug'],r.get('drug_type','') or 'manual']
        DM.loc[k,['gtop_approved','dgidb_approved','drugcentral','withdrawn']]=[False,False,False,False]
    out[(r['antigen'],k,r['route_class'])].append(dict(source='Manual_curation',drug_key=k,src_drug=r['chembl_id'] or r['drug'],target_genes=r['acted_gene'],target_type='manual',target_name='',action=r['mechanism'][:60],mechanism=r['mechanism'],primary=True,pact=float(r['pact']) if r.get('pact') else np.nan,act_type=r.get('act_type',''),refs=r['evidence'],link=r['evidence'].split(';')[0],extra='',genes=[r['acted_gene']],dirn=r['direction'],sil=False,_gene=r['acted_gene'],_link=r['link'],_q=r['link_quality'],_note=r['note'],_d=r['direction'],_man=r))
LV={'approved':3,'clinical':2,'preclinical_tool':1,'no_stage':0}
T_LAB={0:'manual',1:'label',2:'mechanism_record',3:'curated_primary',4:'activity_only',5:'dgidb_only'}
agm={a['target']:a for a in AG}; R=[]
for (p,k,route),rows in out.items():
    a=agm[p]; G=set(a['genes']); m=DM.loc[k] if k in DM.index else None
    if isinstance(m,pd.DataFrame): m=m.iloc[0]
    dt=str(m.drug_type) if m is not None else ''; nm=str(m['name']) if m is not None else k
    for r in rows: r['_t']=0 if r['source']=='Manual_curation' else tier(r)
    tiers=sorted({r['_t'] for r in rows}); dirn='unknown'; basis=''; conflict=False
    ab=abe.get((p,k,route)); manual_only=any(r['_t']==0 for r in rows)
    if ab is not None and not manual_only:
        dirn=ab.decision; basis='label_or_cited_effect: '+ab.excerpt[:90]; tb=1
    else:
        tb=min(tiers)
        for t in tiers:
            ds={r['_d'] for r in rows if r['_t']==t}-{'unknown',''}
            if len(ds)==1: dirn=next(iter(ds)); basis=';'.join(sorted({r['source'].split('_')[0]+':'+str(r['action'])[:30] for r in rows if r['_t']==t and r['_d']==dirn}))[:120]; tb=t; break
            if len(ds)>1: conflict=True; basis='conflict:'+';'.join(sorted(ds)); tb=t; break
    txt=' '.join(str(r['mechanism'])+' '+str(r.get('extra','')) for r in rows)
    kill=False
    is_cd3=bool(G&CD3G)
    if ab is None and not manual_only and route in('direct','complex_partner') and (dirn!='up' or tb>1):
        if dt=='Antibody drug conjugate': dirn='depletion'; basis='record drug type: antibody-drug conjugate'
        elif dt=='Cell' and not is_cd3 and k in LAB and depl_hit(LAB[k].get('moa','')): dirn='depletion'; basis='label: cell therapy directed at the antigen'
        elif depl_hit(txt) and tb>=2:
            if dirn in('up','down') and re.search(r'\bblocks?\b|neutrali[sz]',txt,re.I): kill=True
            else: dirn='depletion'; basis='mechanism text: '+depl_hit(txt).group(0)
        elif k in cd3drug and not is_cd3 and route=='direct' and dt.startswith('Antibody'): dirn='depletion'; basis='mechanism records on CD3 and this antigen (T-cell engager)'
        elif k in cd3drug and is_cd3 and dt.startswith('Antibody') and len({g for r in rows for g in r['genes']})<=4 and dirn=='unknown': basis='CD3 arm record; no direction assigned'
    if route=='direct' and dirn in('up','down') and k in LAB and depl_hit(LAB[k].get('moa','')): kill=True
    if m is not None:
        os_,cp_=m.ot_stage,m.chembl_phase; mp=np.nanmax([os_,cp_]) if not(np.isnan(os_) and np.isnan(cp_)) else np.nan
        otap=str(m.ot_stage_label)=='APPROVAL'
        ap=(cp_>=4) or bool(m.gtop_approved) or (k in LAB) or (otap and (bool(m.drugcentral) or bool(m.dgidb_approved)))
        apn=''
        if ap: mp=4.0
        elif mp>=4: mp=3.0; apn='Open Targets lists approval but no second source confirmed it; treated as clinical, phase capped at 3'
        wd=bool(m.withdrawn); wdd=str(m.withdrawn_detail) if pd.notna(m.withdrawn_detail) else ''
    else: mp=np.nan; ap=False; wd=False; wdd=''; apn=''
    man=[r['_man'] for r in rows if '_man' in r]
    lvl='approved' if ap else ('clinical' if mp>=.5 else None)
    if man and man[0].get('level'):
        lvl=man[0]['level']; 
        if man[0].get('max_phase'): mp=float(man[0]['max_phase'])
    pc=[r['pact'] for r in rows if pd.notna(r['pact'])]; pb=max(pc) if pc else np.nan
    pty=';'.join(sorted({str(r['act_type']) for r in rows if pd.notna(r['pact']) and r['act_type']}))[:30]
    if lvl is None: lvl='preclinical_tool' if (dirn in('up','down') and pd.notna(pb) and tb<=4) else 'no_stage'
    sup=T_LAB[tb]
    fam=(p,k,'family_member') in out
    if route=='direct' and sup=='activity_only' and fam and pd.notna(pb) and pb>=7: sup='family_with_measured_activity'
    qs={r['_q'] for r in rows if r['_q']}
    if route=='ligand_route': q='canonical' if 'canonical' in qs else ('partial' if 'partial' in qs else next(iter(qs),''))
    else: q='partial' if 'partial' in qs else ('obligate' if 'obligate' in qs else next(iter(qs),''))
    lk=';'.join(sorted({r['_link'] for r in rows if r['_link']}))[:100]
    note=';'.join(sorted({str(r['_note']) for r in rows if r['_note']}))[:220]
    srcs=sorted({r['source'] for r in rows})
    lab_url=('https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid='+str(LAB[k]['set_id'])) if k in LAB and LAB[k].get('set_id') else ''
    if ab is not None and lab_url: srcs=srcs+['openFDA_label']
    refs=list(dict.fromkeys(x for r in rows for x in str(r['refs']).split(';') if x))[:8]; lks=list(dict.fromkeys(([lab_url] if ab is not None and lab_url else [])+[r['link'] for r in rows if r['link']]))[:4]
    acts=';'.join(sorted({str(r['action'])[:36] for r in rows if r['action']}))[:160]; mech=';'.join(sorted({str(r['mechanism'])[:70] for r in rows if r['mechanism'] and not str(r['mechanism']).startswith('sources=')}))[:220]
    inv=k in ver[p]
    if route=='family_member': why='family_or_subunit_unresolved'
    elif conflict: why='direction_conflict'
    elif dirn=='unknown': why='direction_unknown'
    elif dirn=='depletion': why='depletion_is_not_up_or_down'
    elif lvl=='preclinical_tool': why='preclinical_tool_only'
    elif lvl=='no_stage': why='no_clinical_stage'
    elif sup in('activity_only','dgidb_only'): why='no_mechanism_grade_record'
    elif route in('ligand_route','complex_partner') and q not in('canonical','obligate'): why='partial_or_ambiguous_link'
    else: why=''
    tier_u='' if why else ('A_direct' if route in('direct','gene_silencing') else 'B_ligand_or_partner')
    if why=='partial_or_ambiguous_link' and q=='partial': tier_u='C_partial_arm'
    R.append(dict(antigen=p,antigen_genes=';'.join(a['genes']),identity_mapping=a['mapping'],route_class=route,acted_gene=';'.join(sorted({r['_gene'] for r in rows})),link=lk,link_quality=q,drug=nm,drug_key=k,chembl_id=k if k.startswith('CHEMBL') else '',drug_type=dt,direction=dirn,direction_basis=basis,direction_conflict=conflict,also_cell_killing=kill,action_records=acts,mechanism_text=mech,level=lvl,max_phase=mp,approval_note=apn,label_set_id=(LAB[k]['set_id'] if k in LAB else ''),withdrawn=wd,withdrawn_detail=wdd,support=sup,sources=';'.join(srcs),n_sources=len(srcs),best_pX=pb,potency_type=pty,in_verified_list=inv,usable=not why,usable_tier=tier_u,not_usable_reason=why,curation_note=note,references=';'.join(refs)[:500],evidence_links=';'.join(lks)[:300]))
D=pd.DataFrame(R)
# drop family_member rows where a direct record exists for the same antigen-drug
dk=set(zip(D[D.route_class.eq('direct')].antigen,D[D.route_class.eq('direct')].drug_key))
D=D[~(D.route_class.eq('family_member')&pd.Series([x in dk for x in zip(D.antigen,D.drug_key)],index=D.index))]
D['_l']=D.level.map(LV); D=D.sort_values(['antigen','route_class','_l','max_phase'],ascending=[True,True,False,False]).drop(columns='_l')
D.to_csv(O+'/drug_records.csv.gz',index=False,compression='gzip')
# verified baseline directions
base={}
for p in agm:
    x=D[(D.antigen.eq(p))&D.in_verified_list&D.route_class.eq('direct')]
    base[p]=set(x[x.direction.isin(['up','down'])&x.level.isin(['approved','clinical'])].direction)
ROUTES=['direct','ligand_route','complex_partner','gene_silencing','family_member']; DIRS=['up','down','depletion','unknown']
S=[]
for a in AG:
    p=a['target']; x=D[D.antigen.eq(p)]; o=dict(antigen=p,genes=';'.join(a['genes']),identity_mapping=a['mapping'],n_records=len(x))
    for rt in ROUTES:
        for dr in DIRS:
            y=x[x.route_class.eq(rt)&x.direction.eq(dr)]; c=y[y.level.isin(['approved','clinical'])]
            o[f'{rt}_{dr}_n']=c.drug_key.nunique(); o[f'{rt}_{dr}_best']='approved' if (c.level=='approved').any() else (f'phase_{c.max_phase.max():g}' if len(c) and c.max_phase.notna().any() else ('clinical' if len(c) else ''))
            o[f'{rt}_{dr}_examples']='; '.join(c.drug.drop_duplicates().head(4).tolist())
            o[f'{rt}_{dr}_preclinical_n']=y[y.level.eq('preclinical_tool')].drug_key.nunique()
    for t,l in (('A_direct','direct'),('B_ligand_or_partner','ligand_or_partner'),('C_partial_arm','partial_arm')):
        y=x[x.usable_tier.eq(t)]
        for dr in('up','down'):
            z=y[y.direction.eq(dr)]; o[f'usable_{l}_{dr}_n']=z.drug_key.nunique(); o[f'usable_{l}_{dr}_examples']='; '.join(z.drug.drop_duplicates().head(4).tolist())
    U=x[x.usable]; uu=set(U.direction); ua=set(U[U.usable_tier.eq('A_direct')].direction); vb=base[p]
    o.update(usable_drugs_not_in_verified_list=int(U[~U.in_verified_list].drug_key.nunique()),usable_drugs_not_in_verified_list_direct_only=int(U[(~U.in_verified_list)&U.usable_tier.eq('A_direct')].drug_key.nunique()),verified_list_drugs=len(ver[p]),verified_directions=';'.join(sorted(vb)),usable_directions=';'.join(sorted(uu)),usable_directions_direct_only=';'.join(sorted(ua)),gain_new_direction=';'.join(sorted(uu-vb)),gain_new_direction_direct_only=';'.join(sorted(ua-vb)),gain_first_usable=bool(uu) and not vb,gain_first_usable_direct_only=bool(ua) and not vb,any_depletion_drug=bool(((x.direction=='depletion')&x.level.isin(['approved','clinical'])).any()))
    S.append(o)
S=pd.DataFrame(S); S.to_csv(O+'/antigen_drug_summary.csv',index=False)
gs=dict(antigens_with_new_usable_drug_any_route=int((S.usable_drugs_not_in_verified_list>0).sum()),antigens_with_new_usable_drug_direct_only=int((S.usable_drugs_not_in_verified_list_direct_only>0).sum()),n_antigens=len(S),n_records=len(D),n_drugs=int(D.drug_key.nunique()),verified_antigens_with_any_drug=int((S.verified_list_drugs>0).sum()),verified_antigens_with_up_or_down=int((S.verified_directions!='').sum()),usable_antigens_all=int((S.usable_directions!='').sum()),usable_antigens_direct_only=int((S.usable_directions_direct_only!='').sum()),gain_first_usable_all=int(S.gain_first_usable.sum()),gain_first_usable_direct_only=int(S.gain_first_usable_direct_only.sum()),gain_new_direction_all=int((S.gain_new_direction!='').sum()),gain_new_direction_direct_only=int((S.gain_new_direction_direct_only!='').sum()))
json.dump(gs,open(O+'/interim/gain_summary.json','w'),indent=1); print(json.dumps(gs)); print(D.route_class.value_counts().to_dict()); print(D.direction.value_counts().to_dict()); print(D.level.value_counts().to_dict()); print(D.not_usable_reason.value_counts().to_dict())
