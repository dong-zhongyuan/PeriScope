make_C <- function() {
  d<-read_pd('C_black_substantia_nigra');stopifnot(nrow(d)==12)
  d$outcome_label<-factor(d$outcome_label,levels=rev(c('QSM L','QSM R','T2* L','T2* R')))
  d$gene<-factor(ifelse(d$gene=='TFRC','TFRC†',d$gene),levels=c('CD22','CD40','TFRC†'))
  ggplot(d,aes(beta,outcome_label,colour=modality))+zero_v()+
    geom_segment(aes(x=ci_lower,xend=ci_upper,yend=outcome_label),linewidth=.4/.pt)+
    geom_point(aes(fill=ifelse(hemisphere=='L',as.character(modality),'open')),shape=21,size=1.1,stroke=.25)+facet_wrap(~gene,nrow=1)+
    scale_fill_manual(values=c(MOD,open='white'),guide='none')+
    scale_colour_manual(values=MOD,guide='none')+
    scale_x_continuous(limits=c(-.15,.21),breaks=c(-.1,.1),labels=num)+
    labs(x='Wald β',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing.x=unit(1,'mm'))
}
