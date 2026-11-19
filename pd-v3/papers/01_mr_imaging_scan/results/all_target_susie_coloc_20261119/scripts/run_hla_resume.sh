# HLA-DRA group 8b6b5e72ab13 (MHC, 19803 SNPs): same s03_susie.R and same settings, one process per pair, 6 BLAS threads each.
O=/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/all_target_susie_coloc_20261119; cd $O/scripts
export OPENBLAS_NUM_THREADS=6 OMP_NUM_THREADS=6 MKL_NUM_THREADS=6
i=$1
SUSIE_L=5 SUSIE_MAXIT=2000 SUSIE_PILOT_N=0 nice -n 5 /public/home/mengxl/dzy/envs/rcoloc/bin/Rscript s03_susie.R $O/work $O/work/results 1 $O/work/resume_hla/g.txt $O/work/resume_hla/p$i.txt
echo FIN_HLA_$i
