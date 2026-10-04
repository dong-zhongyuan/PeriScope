main <- function() {
args <- commandArgs(trailingOnly=TRUE)
sample_id <- if(length(args)) args[1] else 'QC'
work <- Sys.getenv('SPATIAL_WORK_DIR','/public/home/mengxl/dzy/pd_product/_work/spatial_subtype_20261003')
default_out <- '/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/spatial_subtype_20261003'
out <- Sys.getenv('SPATIAL_OUTPUT_DIR',default_out)
inp <- Sys.getenv('SPATIAL_INPUT_DIR',file.path(out,'inputs')); dest <- file.path(out,'samples',sample_id)
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
# Per-sample lock permits independent scheduling without duplicate fits.
if(file.exists(file.path(dest,'completed.txt')))return(invisible(NULL))
lock<-file.path(dest,'fit.lock')
while(!dir.create(lock,showWarnings=FALSE)) {
 if(file.exists(file.path(dest,'completed.txt')))return(invisible(NULL))
 Sys.sleep(2)
}
on.exit(unlink(lock,recursive=TRUE),add=TRUE)
writeLines(as.character(Sys.getpid()),file.path(lock,'pid'))
.libPaths(c(file.path(work,'R_library'),.libPaths()))
Sys.setenv(OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,R_LIBS_USER=file.path(work,'R_library'))
library(Matrix); library(spacexr)
set.seed(20261003)
genes <- readLines(file.path(inp,'genes.txt'))
deconv <- readLines(Sys.getenv('SPATIAL_DECONV_GENE_FILE',file.path(inp,'deconvolution_genes.txt')))
evaluation <- readLines(file.path(inp,'evaluation_genes.txt'))
profiles <- as.matrix(read.csv(file.path(inp,'reference_profiles.csv.gz'),row.names=1,check.names=FALSE))
read_matrix <- function(prefix,meta) {
  con <- gzfile(file.path(inp,paste0(prefix,'_counts.mtx.gz')),'rt')
  on.exit(close(con)); x <- as(readMM(con),'CsparseMatrix')
  rownames(x) <- genes;colnames(x) <- rownames(meta);x
}
rm <- read.csv(file.path(inp,'reference_metadata.csv'),row.names=1,check.names=FALSE)
rc <- read_matrix('reference',rm)
ct <- factor(rm$state);names(ct) <- rownames(rm)
nu <- rm$nUMI;names(nu)<-rownames(rm)
reference <- Reference(rc,ct,nu)
if(sample_id=='QC') {
  md <- read.csv(file.path(inp,'qc_metadata.csv'),row.names=1)
  counts <- read_matrix('qc',md)
} else {
  md <- read.csv(file.path(inp,'spatial_metadata.csv'),row.names=1)
  counts <- read_matrix('spatial',md)
  sel <- md$sample_gsm==sample_id;md<-md[sel,,drop=FALSE];counts<-counts[,rownames(md),drop=FALSE]
}
nUMI <- md$nUMI;names(nUMI)<-rownames(md)
coords<-md[,c('x','y')]
if(!file.exists(file.path(dest,'rctd.rds'))) {
  puck <- SpatialRNA(coords,counts[deconv,,drop=FALSE],nUMI)
  obj <- create.RCTD(puck,reference,max_cores=as.integer(Sys.getenv('SPATIAL_MAX_CORES','8')),cell_type_profiles=profiles,keep_reference=FALSE,UMI_min=100)
  diagnostic <- sample_id=='QC' && Sys.getenv('SPATIAL_QC_DIAGNOSTIC')=='TRUE'
  if(!diagnostic) stopifnot(length(intersect(obj@internal_vars$gene_list_reg,evaluation))==0)
  obj <- run.RCTD(obj,doublet_mode='full')
  saveRDS(obj,file.path(dest,'rctd.rds'))
} else obj <- readRDS(file.path(dest,'rctd.rds'))
w <- as.matrix(obj@results$weights);w <- w/rowSums(w)
write.csv(w,file.path(dest,'cell_state_RNA_fractions.csv'))
writeLines(obj@internal_vars$gene_list_reg,file.path(dest,'actual_deconvolution_genes.txt'))
write.csv(data.frame(state=colnames(w),equivalent_spots=colSums(w),n_spots=nrow(w)),file.path(dest,'state_support.csv'),row.names=FALSE)
if(sample_id=='QC') {
  truth <- read.csv(file.path(inp,'qc_truth.csv'),row.names=1,check.names=FALSE)
  qc <- do.call(rbind,lapply(colnames(w),function(st) {
    x<-truth[rownames(w),st];y<-w[,st]
    data.frame(state=st,n=nrow(w),spearman=cor(x,y,method='spearman'),RMSE=sqrt(mean((x-y)^2)),mean_true=mean(x),mean_estimated=mean(y))
  }))
  qc$quality_pass<-with(qc,n>=100&is.finite(spearman)&spearman>=.5&RMSE<=.15)
  write.csv(qc,file.path(out,'state_identifiability_QC.csv'),row.names=FALSE)
  writeLines('complete',file.path(dest,'completed.txt'));return(invisible(NULL))
}

# Infer expression from real spatial counts, conditional on composition.
# The deconvolution input never contained the evaluation genes restored here.
barcodes<-rownames(w)
obj@originalSpatialRNA <- SpatialRNA(coords[barcodes,,drop=FALSE],counts[evaluation,barcodes,drop=FALSE],nUMI[barcodes])
state_keep<-colnames(w)[colSums(w)>=25]
writeLines(state_keep,file.path(dest,'CSIDE_initial_supported_states.txt'))
if(length(state_keep)<2) {
  writeLines('insufficient_mixture_components',file.path(dest,'unavailable.txt'))
  writeLines('complete_with_unavailable_expression',file.path(dest,'completed.txt'));return(invisible(NULL))
}
if(!file.exists(file.path(dest,'cside.rds'))) {
  # Let the package apply its composition-coverage filter to a fixed point.
  # Passing a preselected list makes expected low coverage a fatal API error.
  obj <- tryCatch(run.CSIDE.intercept(obj,cell_types=NULL,cell_type_threshold=25,
    gene_threshold=1e-5,doublet_mode=FALSE,weight_threshold=.8,sigma_gene=TRUE,normalize_expr=FALSE),
    error=function(e) {
      if(grepl('choose_cell_types: length(cell_types) is 0',conditionMessage(e),fixed=TRUE)) {
        writeLines(conditionMessage(e),file.path(dest,'unavailable.txt'));return(NULL)
      }
      stop(e)
    })
  if(is.null(obj)) {
    writeLines('complete_with_unavailable_expression',file.path(dest,'completed.txt'));return(invisible(NULL))
  }
  saveRDS(obj,file.path(dest,'cside.rds'))
} else obj <- readRDS(file.path(dest,'cside.rds'))
fit<-obj@de_results$gene_fits
expr<-fit$intercept_val
writeLines(colnames(expr),file.path(dest,'CSIDE_included_states.txt'))
se<-fit$s_mat;colnames(se)<-colnames(expr)
ok<-fit$con_mat & !fit$error_mat & is.finite(expr)&is.finite(se)
write.csv(expr,file.path(dest,'cell_state_log_expression.csv'))
write.csv(se,file.path(dest,'cell_state_log_expression_SE.csv'))
write.csv(ok,file.path(dest,'cell_state_gene_converged.csv'))
write.csv(data.frame(state=colnames(expr),n_genes=nrow(expr),converged=colSums(ok)),file.path(dest,'fit_summary.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
writeLines('complete',file.path(dest,'completed.txt'))
cat('SPATIAL_SAMPLE_COMPLETE',sample_id,'\n')

}
main()
