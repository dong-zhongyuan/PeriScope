"""Package reproducible inputs, code and provenance after completed analyses."""
from pathlib import Path
import hashlib,json,subprocess,tarfile,shutil,os,platform
import importlib.metadata

W=Path(__file__).resolve().parent
A=Path('/public/home/mengxl/dzy/pd_product_assets')
P=A/'results/optimization_20261003'
S=Path(os.environ.get('SPATIAL_OUTPUT_DIR',str(P/'spatial_subtype_20261003')))
summary=json.loads((S/'completed.json').read_text())
assert summary['status']=='complete'
assert all(x['returncode']==0 for x in json.loads((S/'run_status.json').read_text())['completed'])
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
inputs=[A/'processed/gse178265/v0.1/GSE178265_sn_annotated.h5ad',
        A/'processed/gse253975/v0.1/GSE253975_geomx.h5ad',
        A/'interim/v0.1/purification/c178_brain_astro.csv',
        A/'interim/v0.1/purification/c178_brain_microglia.csv',
        A/'interim/v0.1/splits/GSE178265_sn_donor_split_v1.json',
        A/'raw/citeseq_hao/GSE164378_sc.meta.data_3P.csv.gz']
inputs += list((A/'raw/citeseq_hao/raw_3p').glob('GSM5008737_RNA_3P-*.gz'))
inputs += [P/'programs.json',P/'target_specific/programs.json',
           P/'six_criteria_20261003/all_6849_programs_six_criteria.csv',
           W.parent/'optimization_20261003/reference/spatial_coordinates.csv',
           W.parent/'optimization_20261003/spatial_analysis.py',
           A/'interim/rescreen_20261002/brain_train.npz']
inputs += [A/'raw/ma_published_reference_20261003/processed_sn_obj.rds', A/'raw/ma_published_reference_20261003/download_provenance.json', W/'ma_reference_sources/GEO_author_donor_crosswalk.csv', W.parent/'optimization_20261003/reference/hgnc_complete_set.txt']
records=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in inputs]
commit=subprocess.check_output(['git','--git-dir='+str(W/'spacexr/.git'),'rev-parse','HEAD'],text=True).strip()
(S/'source_provenance.json').write_text(json.dumps(dict(
    spacexr_repository='https://github.com/dmcable/spacexr',spacexr_commit=commit,
    biological_data='Real raw counts; synthetic known mixtures used only to assess state recovery',
    sources=records),indent=2))
shutil.copy2(W/'R_sessionInfo.txt',S/'R_sessionInfo.txt')
shutil.copy2(W/'复现说明.txt',S/'复现说明.txt')
packages={k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','h5py','anndata']}
(S/'Python_environment.json').write_text(json.dumps(dict(python=platform.python_version(),packages=packages),indent=2))
(S/'requirements.txt').write_text(''.join(k+'=='+v+'\n' for k,v in packages.items()))
with tarfile.open(S/'reproduction_source.tar.gz','w:gz') as tar:
    for p in sorted(W.iterdir()):
        if p.suffix in ['.py','.R'] and not p.name.startswith(('CSIDE','RCTD')):
            tar.add(p,arcname=p.name)
    tar.add(W.parent/'optimization_20261003/input_utils.py',arcname='reference_dependency/input_utils.py')
    tar.add(W.parent/'optimization_20261003/reference/hgnc_complete_set.txt',arcname='reference_dependency/hgnc_complete_set.txt')
    tar.add(W/'ma_reference_sources/GEO_author_donor_crosswalk.csv',arcname='ma_reference_sources/GEO_author_donor_crosswalk.csv')
    tar.add(W.parent/'optimization_20261003/screen_six_criteria.py',arcname='pipeline_integration/screen_six_criteria.py')
    tar.add(W.parent/'optimization_20261003/reapply_completed_spatial.py',arcname='pipeline_integration/reapply_completed_spatial.py')
    for name in ['DESCRIPTION','LICENSE']:
        p=W/'spacexr'/name
        if p.exists():tar.add(p,arcname='spacexr_version/'+name)
manifest=[dict(path=str(p.relative_to(S)),bytes=p.stat().st_size,sha256=sha(p))
          for p in sorted(S.rglob('*')) if p.is_file() and p.name!='output_manifest.json']
(S/'output_manifest.json').write_text(json.dumps(manifest,indent=2))
status=json.loads((P/'delivery_status.json').read_text()) if (P/'delivery_status.json').exists() else {}
status['spatial_subtype_directory']=S.name
status['spatial_subtype_result']=summary
status['spatial_seventh_criterion_adopted']=summary['spatial_layer_adoptable']
status['spatial_supported_nomination_directory']='seven_criteria_20261003'
if not summary['spatial_layer_adoptable']:
    status['final_nomination_directory']='six_criteria_20261003'
    status['nomination_scope']='Six-criterion primary screen plus separate seven-criterion spatial-supported tier; unresolved spatial evidence is retained.'
if summary['spatial_layer_adoptable']:
    status['final_nomination_directory']='seven_criteria_20261003'
    status['nomination_scope']='Unchanged C1-C6 plus matched-state spatial C7; unresolved spatial coverage remains distinct from tested non-support.'
    status['result']=summary
(P/'delivery_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
print('SPATIAL_DELIVERY_PACKAGED',len(manifest),commit)
