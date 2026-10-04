from pathlib import Path
import json,os,hashlib
import pandas as pd,numpy as np
P=Path(os.environ.get('PERISCOPE_RESULT_ROOT','/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003'))
S=Path(os.environ.get('SPATIAL_OUTPUT_DIR',str(P/'spatial_subtype_author_20261003')));O=P/'seven_criteria_20261003'
e=pd.read_csv(O/'all_6849_programs_seven_criteria.csv');old=pd.read_csv(P/'six_criteria_20261003/all_6849_programs_six_criteria.csv')
a=old.set_index('combo').sort_index();b=e.set_index('combo').sort_index()
assert len(e)==6849 and e.protein.nunique()==213
for c in a:pd.testing.assert_series_equal(a[c],b[c],check_names=False)
assert e.combo.is_unique
assert (e.seven_criteria_pass==(e.six_criteria_pass&e.C7.eq('pass'))).all()
leaf=pd.read_csv(S/'all_leaf_program_spatial_tests.csv')
assert not leaf.duplicated(['combo','spatial_subtype']).any()
assert len(leaf)==len(e)+int(e.brain_state.eq('astro').sum())
assert (leaf.loc[leaf.leaf_C7.eq('pass'),'effect']*leaf.loc[leaf.leaf_C7.eq('pass'),'matched_brain_effect']>0).all()
for combo,g in leaf.groupby('combo'):
 expect='pass' if g.leaf_C7.eq('pass').any() else 'unresolved' if g.leaf_C7.eq('unresolved').any() else 'not_supported'
 assert b.loc[combo,'C7']==expect
reserved=set((S/'inputs/deconvolution_genes.txt').read_text().splitlines());evals=set((S/'inputs/evaluation_genes.txt').read_text().splitlines())
assert not reserved&evals
roles=pd.read_csv(S/'donor_roles.csv');assert not roles.loc[roles.role.isin(['reference','QC_holdout']),'spatial_overlap'].any()
counts={}
for p in sorted((S/'samples').iterdir()):
 assert (p/'completed.txt').exists()
 w=pd.read_csv(p/'cell_state_RNA_fractions.csv',index_col=0)
 assert np.isfinite(w.to_numpy()).all() and w.min().min()>=0 and np.allclose(w.sum(1),1)
 assert not set((p/'actual_deconvolution_genes.txt').read_text().splitlines())&evals
 counts[p.name]=len(w)
report=dict(status='passed',original_six_criteria_all_columns_unchanged=True,n_programs=len(e),n_antigens=e.protein.nunique(),n_leaf_tests=len(leaf),source_donor_overlap=False,composition_evaluation_gene_overlap=False,n_spots_by_sample=counts,classification_verified=True)
(S/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
