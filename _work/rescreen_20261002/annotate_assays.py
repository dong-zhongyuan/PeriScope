"""Measured CITE channels and gene-level public Olink coverage for the complete panel."""
from pathlib import Path
import json,hashlib,time,urllib.request
import pandas as pd
from input_utils import canonical_symbols

O=Path('/public/home/mengxl/dzy/pd_product_assets/results/rescreen_20261002');REF=Path(__file__).parent/'reference'
url='https://biobank.ndph.ox.ac.uk/showcase/codown.cgi?id=143';path=REF/'olink_ukb_assays.tsv'
if not path.exists():urllib.request.urlretrieve(url,path)
table=pd.read_csv(path,sep='\t');assert len(table)==2923 and {'coding','meaning'}<=set(table.columns)
table['gene']=canonical_symbols(table.meaning.str.split(';').str[0]);rows=[]
for r in json.loads((O/'candidate_registry.json').read_text()):
    if r['is_control']:continue
    match=table[table.gene.isin(r['genes'])]
    rows.append(dict(target=r['target'],genes=';'.join(r['genes']),mapping=r['mapping'],CITE_measured=True,
        CITE_channels=';'.join(r['channel_names']),Olink_covered_genes=';'.join(sorted(set(match.gene))),
        Olink_assay_ids=';'.join(match.coding.astype(str)),Olink_gene_level_coverage=bool(len(match)),
        molecular_form_note='RNA gene mapping does not resolve CD45 isoforms or intact heteromeric complexes' if r['target'] in ['CD45RA','CD45RB','CD45RO'] or r['mapping']=='defined_complex' else '',
        participant_level_plasma_measurements_used=False))
pd.DataFrame(rows).to_csv(O/'all_target_assay_coverage.csv',index=False)
(O/'assay_annotation_source.json').write_text(json.dumps(dict(url=url,catalogue='https://biobank.ctsu.ox.ac.uk/ukb/coding.cgi?id=143',
    checked=time.strftime('%Y-%m-%d %H:%M:%S'),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),n_assays=len(table),
    scope='panel inclusion only; no new participant-level plasma measurements or plasma accuracy claims',
    soma_scan='Full current panel not retrieved for this rerun; old seven-target annotations are not extended to new candidates'),indent=2))
print('ASSAY COVERAGE',len(rows),sum(r['Olink_gene_level_coverage'] for r in rows),flush=True)
