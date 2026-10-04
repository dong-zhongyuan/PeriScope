"""Record the numerical environment and input identities for CellOT reconstruction."""
import sys,json,platform,hashlib
from pathlib import Path
import numpy,scipy,sklearn,torch
from threadpoolctl import threadpool_info
W=Path(__file__).resolve().parent
B=Path('/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002');R=B/'cellot_recomputed_20261003';P=B.parent/'pig_external_validation_20261003'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(4*1024*1024),b''):h.update(v)
 return h.hexdigest()
files=list((B/'cache').glob('*.npz'))+list((R/'projection').glob('*.npz'))+[P/'adapted_pig_blood_reference_inputs.npz',P/'human_reference_cells.npz',P/'inputs/brain_one2one_sum_provided_transcript_TPM.csv']
obj=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,torch=torch.__version__,cuda=torch.version.cuda,blas_libraries=threadpool_info(),training_function_sha256=sha(W/'original_cellot_function.py'),inputs={str(p):sha(p) for p in files},numerical_settings={'PCA_original_six_axes_threads':2,'PCA_pDC_extension_threads':1,'CPU_worker_torch_threads':1,'CPU_worker_BLAS_threads':1,'training_seeds':[42,43,44],'PCA_seed':42,'seed42_original_six_axes_device':'cpu','seed42_original_six_axes_torch_threads':2,'seed42_original_six_axes_BLAS_threads':2,'residual_BLAS_threads':1,'pig_BLAS_threads':2,'pig_device':'cuda:1'})
(R/'recovery_environment.json').write_text(json.dumps(obj,indent=2));print('ENVIRONMENT_AND_INPUT_HASHES_SAVED',len(files))
