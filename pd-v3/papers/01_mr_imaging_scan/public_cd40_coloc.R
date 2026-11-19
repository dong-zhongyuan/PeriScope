# Official, pinned coloc ABF functions; the complete package is not installed.
# Inputs have already been identity-matched. All SNPs are passed without a P
# filter. Known unit-variance phenotypes avoid estimating sdY from imputed N.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
ROOT <- normalizePath(args[1])
OUT <- file.path(ROOT,'papers/01_mr_imaging_scan/results/public_cd40_20261116')
source(file.path(ROOT,'reference/coloc/claudia.R'))
priors <- sort(unique(c(10^seq(-6,log10(3e-5),length.out=31),1e-5)))
results <- list()
for(idp in c('0017','0018')) {
  d <- read.csv(gzfile(file.path(OUT,paste0(idp,'_shared_snps.csv.gz'))),check.names=FALSE)
  stopifnot(!anyDuplicated(d$variant_key),all(d$se_x>0),all(d$se_y>0),nrow(d)>100)
  for(window in c(500000,1000000)) {
    v <- d[abs(d$pos37-44747947)<=window,]
    stopifnot(nrow(v)>100)
    l1 <- approx.bf.estimates(v$beta_x/v$se_x,v$se_x^2,type='quant',sdY=1)$lABF
    l2 <- approx.bf.estimates(v$beta_y/v$se_y,v$se_y^2,type='quant',sdY=1)$lABF
    # Apply coloc's prior adjustment to exactly the common SNP set supplied.
    p1 <- if(nrow(v)*1e-4>=1) 1/(nrow(v)+1) else 1e-4
    p2 <- p1
    for(p12 in priors) {
      effective_p12 <- if(nrow(v)*p12>=1) 1/(nrow(v)+1) else p12
      pp <- combine.abf(l1,l2,p1,p2,effective_p12,quiet=TRUE)
      results[[length(results)+1]] <- data.frame(idp=idp,window_bp=window,nsnps=nrow(v),
         p1_requested=1e-4,p2_requested=1e-4,p12_requested=p12,p1=p1,p2=p2,p12=effective_p12,
         sdY_protein=1,sdY_imaging=1,effect_prior_sd=.15,t(pp),check.names=FALSE)
    }
    v$log_abf_protein <- l1; v$log_abf_imaging <- l2
    v$snp_pp_H4_conditional <- exp(l1+l2-logsum(l1+l2))
    v$snp_pp_protein_conditional <- exp(l1-logsum(l1))
    v$snp_pp_imaging_conditional <- exp(l2-logsum(l2))
    con <- gzfile(file.path(OUT,paste0(idp,'_window_',format(window,scientific=FALSE,trim=TRUE),'_snp_posteriors.csv.gz')),'wt')
    write.csv(v,con,row.names=FALSE);close(con)
  }
}
write.csv(do.call(rbind,results),file.path(OUT,'coloc_prior_sensitivity.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'R_sessionInfo.txt'))
