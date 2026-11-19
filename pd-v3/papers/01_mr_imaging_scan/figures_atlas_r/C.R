make_C <- function() {
  d<-read_table('C_instrument_strength')
  stopifnot(!anyDuplicated(d$exposure_id),all(d$F_x>0),all(is.finite(d$absolute_best_z)))
  counts<-table(factor(d$availability,levels=c('Primary','EAF-only')))
  ggplot(d,aes(F_x,absolute_best_z,colour=availability,shape=availability))+
    geom_point(size=.55,alpha=.6,stroke=.2)+
    scale_x_log10(breaks=10^(0:8),labels=function(x)format(x,scientific=TRUE,trim=TRUE))+
    scale_colour_manual(values=c(Primary=PCOL[['drop']],`EAF-only`=PCOL[['infer']]),
      breaks=names(counts),labels=paste0(names(counts),' (',counts,')'),name=NULL)+
    scale_shape_manual(values=c(Primary=16,`EAF-only`=1),
      breaks=names(counts),labels=paste0(names(counts),' (',counts,')'),name=NULL)+
    labs(x='Exposure instrument F',y='Absolute best-pair z')+
    guides(colour=guide_legend(nrow=2,override.aes=list(alpha=1,size=1,shape=c(16,1))),shape='none')+
    theme(legend.position='bottom',legend.key.width=unit(3,'mm'))
}
