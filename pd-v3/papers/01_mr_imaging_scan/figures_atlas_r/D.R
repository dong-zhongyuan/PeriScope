make_D <- function() {
  d<-read_table('D_phenotype_coverage')
  stopifnot(nrow(d)==cfg$phenotypes*2,all(d$n_planned==cfg$assays))
  d$family_label<-factor(d$family,levels=cats$family,
    labels=paste0(cats$label,' (',cats$n_phenotypes,')'))
  # stat_ecdf uses every phenotype, including tied/zero fractions. The full
  # per-phenotype table remains available; tied curves are intentionally equal.
  ggplot(d,aes(estimated_fraction,colour=policy,linetype=policy))+
    stat_ecdf(geom='step',pad=TRUE,linewidth=.4/.pt)+
    facet_wrap(~family_label,ncol=if(nrow(cats)>12)4 else if(nrow(cats)>4)3 else 1)+
    scale_colour_manual(values=PCOL,labels=POLICY,name=NULL)+
    scale_linetype_manual(values=PLINE,labels=POLICY,name=NULL)+
    scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','.5','1'),expand=c(0,0))+
    scale_y_continuous(limits=c(0,1),breaks=c(0,1),expand=c(0,0))+
    labs(x='Estimated fraction of all 2,941 assays',y='Cumulative phenotype fraction')+
    theme(legend.position='bottom',panel.spacing=unit(2,'mm'),strip.text=element_text(size=6,face='bold'))
}
