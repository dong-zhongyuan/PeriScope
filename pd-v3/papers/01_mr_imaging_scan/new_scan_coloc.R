# Pinned official coloc ABFs. No decision threshold and no unit-variance assumption for AGES protein.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
ROOT <- normalizePath(args[1]); OUT <- file.path(ROOT,'papers/01_mr_imaging_scan/results/new_scan_regions_20261117')
source(file.path(ROOT,'reference/coloc/claudia.R'))
batch <- read.csv(file.path(OUT,'batch_regions.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(batch)>0,!anyDuplicated(batch$accession))
priors <- sort(unique(c(10^seq(-6,log10(3e-5),length.out=31),1e-5)))
results <- list(); diagnostics <- list()
for(i in seq_len(nrow(batch))) {
  rec <- batch[i,]
  accession <- rec$accession; folder <- file.path(OUT,accession)
  d <- read.csv(gzfile(file.path(folder,'shared_snps.csv.gz')),check.names=FALSE)
  stopifnot(!anyDuplicated(d$variant_key),all(d$standard_error_x>0),all(d$standard_error_y>0))
  for(window in c(500000,1000000)) {
    v <- d[abs(d$pos37-rec$center)<=window,]
    stopifnot(nrow(v)>0,all(is.finite(v$effect_allele_frequency_x)),length(unique(v$n_x))==1)
    maf <- pmin(v$effect_allele_frequency_x,1-v$effect_allele_frequency_x)
    sdY_x <- sdY.est(v$standard_error_x^2,maf,unique(v$n_x))
    stopifnot(is.finite(sdY_x),sdY_x>0)
    l1 <- approx.bf.estimates(v$beta_x/v$standard_error_x,v$standard_error_x^2,type='quant',sdY=sdY_x)$lABF
    l2 <- approx.bf.estimates(v$beta_y/v$standard_error_y,v$standard_error_y^2,type='quant',sdY=1)$lABF
    stopifnot(all(is.finite(l1)),all(is.finite(l2)))
    p1 <- if(nrow(v)*1e-4>=1) 1/(nrow(v)+1) else 1e-4; p2 <- p1
    for(p12 in priors) {
      effective_p12 <- if(nrow(v)*p12>=1) 1/(nrow(v)+1) else p12
      pp <- combine.abf(l1,l2,p1,p2,effective_p12,quiet=TRUE)
      stopifnot(all(is.finite(pp)),abs(sum(pp)-1)<1e-10)
      results[[length(results)+1]] <- data.frame(accession=accession,gene=rec$gene,phenotype=rec$phenotype,window_bp=window,nsnps=nrow(v),p1_requested=1e-4,p2_requested=1e-4,p12_requested=p12,p1=p1,p2=p2,p12=effective_p12,sdY_protein_estimated=sdY_x,sdY_imaging=1,protein_N_metadata=unique(v$n_x),effect_prior_sd_fraction=.15,t(pp),check.names=FALSE)
    }
    diagnostics[[length(diagnostics)+1]] <- data.frame(accession=accession,gene=rec$gene,window_bp=window,nsnps=nrow(v),sdY_protein_estimated=sdY_x,maf_min=min(maf),maf_max=max(maf),protein_N_metadata=unique(v$n_x),min_p_protein=min(v$p_value_x),min_p_imaging=min(v$p_value_y),min_se_protein=min(v$standard_error_x),max_se_protein=max(v$standard_error_x))
    v$log_abf_protein <- l1;v$log_abf_imaging <- l2
    v$snp_pp_H4_conditional <- exp(l1+l2-logsum(l1+l2));v$snp_pp_protein_conditional <- exp(l1-logsum(l1));v$snp_pp_imaging_conditional <- exp(l2-logsum(l2))
    con <- gzfile(file.path(folder,paste0('window_',format(window,scientific=FALSE,trim=TRUE),'_snp_posteriors.csv.gz')),'wt');write.csv(v,con,row.names=FALSE);close(con)
  }
}
write.csv(do.call(rbind,results),file.path(OUT,'coloc_prior_sensitivity.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(OUT,'coloc_scale_diagnostics.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'R_sessionInfo.txt'))
