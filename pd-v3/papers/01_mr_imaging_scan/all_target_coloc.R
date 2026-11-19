# Original pinned single-causal-variant ABF recipe, one planned region at a time.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
ROOT <- normalizePath(args[1]); OUT <- normalizePath(args[2])
source(file.path(ROOT,'reference/coloc/claudia.R'))
rec <- read.csv(file.path(OUT,'region.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(rec)==1)
d <- read.csv(gzfile(file.path(OUT,'shared_snps.csv.gz')),check.names=FALSE)
stopifnot(!anyDuplicated(d$variant_key),all(d$standard_error_x>0),all(d$standard_error_y>0))
priors <- sort(unique(c(10^seq(-6,log10(3e-5),length.out=31),1e-5)))
results <- list(); diagnostics <- list()
for(window in c(500000,1000000)) {
 v <- d[abs(d$pos37-rec$center)<=window,]
 ns <- nrow(v); status <- 'defined'; reason <- ''; sdY_x <- NA_real_
 maf <- pmin(v$effect_allele_frequency_x,1-v$effect_allele_frequency_x)
 if(ns==0) {status <- 'undefined_no_shared_snps';reason <- 'No exact biallelic SNP shared in this window'}
 else {
  stopifnot(all(is.finite(maf)),length(unique(v$n_x))==1,all(v$n_x>0))
  estimate <- tryCatch(sdY.est(v$standard_error_x^2,maf,unique(v$n_x)),error=function(e)e)
  if(inherits(estimate,'error')) {status <- 'undefined_protein_scale';reason <- conditionMessage(estimate)}
  else if(!is.finite(estimate)||estimate<=0) {status <- 'undefined_protein_scale';reason <- 'Nonpositive or nonfinite estimated protein SD'}
  else sdY_x <- estimate
 }
 p1 <- if(ns*1e-4>=1) 1/(ns+1) else 1e-4
 if(status=='defined') {
  l1 <- approx.bf.estimates(v$beta_x/v$standard_error_x,v$standard_error_x^2,type='quant',sdY=sdY_x)$lABF
  l2 <- approx.bf.estimates(v$beta_y/v$standard_error_y,v$standard_error_y^2,type='quant',sdY=1)$lABF
  stopifnot(all(is.finite(l1)),all(is.finite(l2)))
 } else {l1 <- rep(NA_real_,ns);l2 <- rep(NA_real_,ns)}
 for(p12 in priors) {
  effective_p12 <- if(ns*p12>=1) 1/(ns+1) else p12
  pp <- setNames(rep(NA_real_,5),paste0('PP.H',0:4,'.abf'))
  if(status=='defined') {
   pp <- combine.abf(l1,l2,p1,p1,effective_p12,quiet=TRUE)
   stopifnot(all(is.finite(pp)),abs(sum(pp)-1)<1e-10)
  }
  results[[length(results)+1]] <- data.frame(region_id=rec$region_id,accession=rec$accession,gene=rec$gene,
   phenotype=rec$phenotype,window_bp=window,nsnps=ns,p1_requested=1e-4,p2_requested=1e-4,p12_requested=p12,
   p1=p1,p2=p1,p12=effective_p12,sdY_protein_estimated=sdY_x,sdY_imaging=1,
   protein_N_metadata=rec$protein_N_metadata,effect_prior_sd_fraction=.15,t(pp),
   status=status,undefined_reason=reason,check.names=FALSE)
 }
 safe_min <- function(x) if(length(x)) min(x) else NA_real_
 safe_max <- function(x) if(length(x)) max(x) else NA_real_
 diagnostics[[length(diagnostics)+1]] <- data.frame(region_id=rec$region_id,accession=rec$accession,gene=rec$gene,
  window_bp=window,nsnps=ns,sdY_protein_estimated=sdY_x,maf_min=safe_min(maf),maf_max=safe_max(maf),
  protein_N_metadata=rec$protein_N_metadata,min_p_protein=safe_min(v$p_value_x),min_p_imaging=safe_min(v$p_value_y),
  min_se_protein=safe_min(v$standard_error_x),max_se_protein=safe_max(v$standard_error_x),status=status,undefined_reason=reason)
 v$log_abf_protein <- l1;v$log_abf_imaging <- l2
 v$snp_pp_H4_conditional <- if(status=='defined') exp(l1+l2-logsum(l1+l2)) else rep(NA_real_,ns)
 v$snp_pp_protein_conditional <- if(status=='defined') exp(l1-logsum(l1)) else rep(NA_real_,ns)
 v$snp_pp_imaging_conditional <- if(status=='defined') exp(l2-logsum(l2)) else rep(NA_real_,ns)
 con <- gzfile(file.path(OUT,paste0('window_',format(window,scientific=FALSE,trim=TRUE),'_snp_posteriors.csv.gz')),'wt')
 write.csv(v,con,row.names=FALSE);close(con)
}
write.csv(do.call(rbind,results),file.path(OUT,'coloc_prior_sensitivity.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(OUT,'coloc_scale_diagnostics.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'R_sessionInfo.txt'))
