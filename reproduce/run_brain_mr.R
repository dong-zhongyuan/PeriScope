# Single cis-instrument blood protein -> continuous brain phenotype MR.
args <- commandArgs(trailingOnly=TRUE)
o <- normalizePath(args[1], mustWork=TRUE)
e <- read.csv(file.path(o,'combined_exposure_instruments.csv'),check.names=FALSE)
y <- read.csv(file.path(o,'brain_sentinel_statistics.csv'),check.names=FALSE)
x <- merge(e,y,by.x='assay',by.y='exposure_assay',all.x=FALSE,suffixes=c('_exposure','_outcome'))
comp <- function(a) chartr('ACGT','TGCA',a)
x$harmonization <- 'allele_mismatch';x$sign <- NA_real_
for(i in seq_len(nrow(x))) {
 a<-x$effect_allele_exposure[i]; b<-x$other_allele_exposure[i]; c<-x$effect_allele_outcome[i];d<-x$other_allele_outcome[i]
 if(!all(nchar(c(a,b,c,d))==1) || !all(c(a,b,c,d)%in%c('A','C','G','T'))) {x$harmonization[i]<-'indel_or_nonSNV';next}
 if(paste(sort(c(a,b)),collapse='')%in%c('AT','CG')) {x$harmonization[i]<-'palindrome_without_outcome_frequency';next}
 if((a==c&&b==d)||(comp(a)==c&&comp(b)==d)) {x$sign[i]<-1;x$harmonization[i]<-'aligned'}
 if((a==d&&b==c)||(comp(a)==d&&comp(b)==c)) {x$sign[i]<--1;x$harmonization[i]<-'swapped'}
}
x$eligible <- is.finite(x$sign)&is.finite(x$beta_exposure)&is.finite(x$se_exposure)&is.finite(x$beta_outcome)&is.finite(x$se_outcome)&x$F_exposure>=10&x$beta_exposure!=0&x$se_outcome>0
# Cohort priority uses eligible instruments, never the outcome effect or P value.
ukb_genes <- unique(x$gene[x$source=='Sun2023_ST9_discovery' & x$eligible])
a <- e[e$source=='Gudjonsson2022_AGES_serum' & !e$gene%in%ukb_genes,]
a <- a[order(-a$F_exposure,a$assay),]
ages_primary <- a$assay[!duplicated(a$gene)]
x$selected <- x$eligible & (x$source=='Sun2023_ST9_discovery' | x$assay%in%ages_primary)
x$analysis_set <- ifelse(x$source=='Sun2023_ST9_discovery','UKB_primary',
 ifelse(x$gene%in%ukb_genes,'AGES_not_run_UKB_covered',
 ifelse(x$assay%in%ages_primary,'AGES_gap_fill','AGES_not_run_alternate_assay')))
x$selection_reason <- ifelse(x$selected,'selected_by_cohort_priority',
 ifelse(!x$eligible,x$harmonization,x$analysis_set))
x$beta_MR <- ifelse(x$selected,x$sign*x$beta_outcome/x$beta_exposure,NA_real_)
x$se_MR <- ifelse(x$selected,x$se_outcome/abs(x$beta_exposure),NA_real_)
x$se_second_order <- ifelse(x$selected,sqrt(x$se_outcome^2/x$beta_exposure^2+x$beta_outcome^2*x$se_exposure^2/x$beta_exposure^4),NA_real_)
x$p_MR <- 2*pnorm(-abs(x$beta_MR/x$se_MR));x$CI_low<-x$beta_MR-1.96*x$se_MR;x$CI_high<-x$beta_MR+1.96*x$se_MR
# One combined hypothesis family across UKB primary and AGES gap-fill.
x$q_BH <- NA_real_
ii <- which(x$selected)
x$q_BH[ii] <- p.adjust(x$p_MR[ii],method='BH')
x<-x[order(x$analysis_set,x$p_MR),]
write.csv(x[x$selected,],file.path(o,'all_brain_MR_estimates.csv'),row.names=FALSE,na='')
write.csv(x[x$selected & x$analysis_set=='UKB_primary',],file.path(o,'UKB_primary_MR.csv'),row.names=FALSE,na='')
write.csv(x[x$selected & x$analysis_set=='AGES_gap_fill',],file.path(o,'AGES_gap_fill_MR.csv'),row.names=FALSE,na='')
write.csv(x[,c('assay','gene','source','outcome_accession','eligible','selected','analysis_set','selection_reason')],file.path(o,'cohort_selection_audit.csv'),row.names=FALSE,na='')
write.csv(x[,c('assay','gene','outcome_accession','harmonization','eligible','selected','analysis_set')],file.path(o,'harmonization_audit.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(o,'R_sessionInfo.txt'))
for(g in unique(x$analysis_set)) {z<-x[x$analysis_set==g&x$selected,];cat(g,nrow(z),'eligible;',sum(z$p_MR<.05),'nominal;',sum(z$q_BH<.05),'FDR\n')}
