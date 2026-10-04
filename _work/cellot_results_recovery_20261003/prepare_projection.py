"""Reconstruct the original PCA using its verified two-thread numerical setting."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='2'
import json
from pathlib import Path
import numpy as np
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
B=Path('/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002')
R=B/'cellot_recomputed_20261003';P=R/'projection';P.mkdir(parents=True,exist_ok=True)
axes=[bs+'__'+rs for bs in ['cDC','classical_mono','nonclassical_mono','pDC'] for rs in ['astro','microglia_mhc2']]
for axis in axes:
 dest=P/f'{axis}.npz';assert not dest.exists()
 cache=dict(np.load(B/'cache'/f'{axis}.npz'))
 with threadpool_limits(limits=1 if axis.startswith('pDC__') else 2, user_api='blas'):
  p=PCA(64,svd_solver='randomized',random_state=42).fit(np.concatenate([cache[f'train_{k}{c}'] for k in ['x','y'] for c in (0,1)]))
 np.savez_compressed(dest,pca_mean=p.mean_,pca_components=p.components_)
 print('PCA_READY',axis,flush=True)
(R/'projection_reconstruction.json').write_text(json.dumps(dict(numerical_threads=dict(original_six_axes=2,pDC_extension=1),solver='randomized',n_components=64,random_state=42,confirmed_original_log={'classical_mono__astro_seed43_condition0_step1000':0.08640657587276501,'classical_mono__astro_seed43_condition1_step1000':0.08771100400034615}),indent=2))
