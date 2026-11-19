args <- commandArgs(trailingOnly=TRUE);ROOT <- normalizePath(args[1])
.libPaths(c(file.path(ROOT,'.deps/R'),.libPaths()))
suppressPackageStartupMessages({library(coloc);library(susieR);library(data.table);library(jsonlite)})
OUT <- file.path(ROOT,'papers/01_mr_imaging_scan/results/public_cd40_ld_20261116')
meta <- fromJSON(file.path(OUT,'INPUT_AUDIT.json'))
stopifnot(meta$complete)
d <- read.csv(gzfile(file.path(OUT,'model_inputs.csv.gz')),check.names=FALSE)
genotypes <- read.csv(gzfile(file.path(OUT,'alt_dosages.csv.gz')),check.names=FALSE)
stopifnot(identical(d$variant_key,genotypes$variant_key),nrow(d)>100)
G <- as.matrix(genotypes[,-1]);R <- cor(t(G))
dimnames(R) <- list(d$variant_key,d$variant_key)
stopifnot(all(is.finite(R)),max(abs(R-t(R)))<1e-12,max(abs(diag(R)-1))<1e-12)
checks <- read.csv(file.path(OUT,'python_ld_checks.csv'))
lead <- which(d$source_variant=='rs4810485');stopifnot(length(lead)==1)
stopifnot(max(abs(R[,lead]-checks$signed_r_to_rs4810485))<1e-10,
          max(abs(rowMeans(G)/2-checks$alt_frequency))<1e-12)
saveRDS(R,file.path(OUT,'signed_alt_ld.rds'))
make_data <- function(b,se,n)list(beta=b,varbeta=se^2,snp=d$variant_key,position=d$pos37,
                                type='quant',sdY=1,N=n,LD=R)
datasets <- list(protein=make_data(d$beta_x_alt,d$se_x,meta$n_x_model),
                 left=make_data(d$beta_y_left,d$se_y_left,meta$n_y_model),
                 right=make_data(d$beta_y_right,d$se_y_right,meta$n_y_model))
fits <- list();diagnostics <- list();warnings <- list()
for(name in names(datasets)) {
 cat('Running',name,'with',nrow(d),'SNPs\n');flush.console()
 ds <- datasets[[name]];w <- character()
 fit <- withCallingHandlers(runsusie(ds,L=10,maxit=1000,repeat_until_convergence=FALSE,
                                    estimate_residual_variance=FALSE,coverage=.95,min_abs_corr=.5),
                            warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')})
 fits[[name]] <- fit;warnings[[name]] <- w
 z <- ds$beta/sqrt(ds$varbeta)
 s <- estimate_s_rss(z=z,R=R,n=ds$N)
 diagnostics[[name]] <- list(converged=fit$converged,iterations=fit$niter,
     n_credible_sets=length(fit$sets$cs),ld_regularization_diagnostic=s,
     min_p=2*pnorm(-max(abs(z))),warnings=w)
 saveRDS(fit,file.path(OUT,paste0('susie_',name,'.rds')))
 write.csv(data.frame(variant_key=d$variant_key,pip=fit$pip),file.path(OUT,paste0('pip_',name,'.csv')),row.names=FALSE)
}
summaries <- list();baseline <- list()
for(name in c('left','right')) {
 cat('Colocalization',name,'\n');flush.console()
 base <- coloc.abf(datasets$protein,datasets[[name]],p1=1e-4,p2=1e-4,p12=1e-5)
 baseline[[name]] <- as.list(base$summary)
 if(length(fits$protein$sets$cs)==0 || length(fits[[name]]$sets$cs)==0) {
  summaries[[name]] <- list(status='not_estimable_no_credible_set',posterior=NULL)
 } else {
  res <- coloc.susie(fits$protein,fits[[name]],p1=1e-4,p2=1e-4,p12=1e-5)
  saveRDS(res,file.path(OUT,paste0('coloc_susie_',name,'.rds')))
  fwrite(res$summary,file.path(OUT,paste0('coloc_susie_',name,'.csv')))
  summaries[[name]] <- list(status='computed',summary=as.data.frame(res$summary))
 }
}
write_json(list(complete=TRUE,model_snps=nrow(d),reference_samples=ncol(G),
    reference='1000 Genomes phase 3 EUR, external and not cohort-matched',
    L=10,coverage=.95,min_abs_corr=.5,estimate_residual_variance=FALSE,
    diagnostics=diagnostics,coloc_susie=summaries,same_input_abf=baseline,
    independent_ld_check='All SNP correlations to historical lead and allele frequencies agree with Python',
    versions=list(coloc=as.character(packageVersion('coloc')),susieR=as.character(packageVersion('susieR')))),
    file.path(OUT,'RESULT.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=16)
writeLines(capture.output(sessionInfo()),file.path(OUT,'R_sessionInfo.txt'))
