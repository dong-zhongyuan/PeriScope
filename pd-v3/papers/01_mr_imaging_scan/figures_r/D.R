make_D <- function() {
  d<-read_pd('D_tfrc_regions');stopifnot(nrow(d)==32)
  d$region_short<-factor(d$region_short,levels=rev(REG))
  d$yy<-as.numeric(d$region_short)+ifelse(d$hemisphere=='L',.15,-.15)
  ggplot(d,aes(beta,yy,colour=modality))+zero_v()+
    geom_segment(aes(x=ci_lower,xend=ci_upper,yend=yy),linewidth=.4/.pt)+
    geom_point(aes(shape=hemisphere),size=1.15)+facet_wrap(~modality,nrow=1)+
    scale_shape_manual(values=c(L=16,R=1),guide='none')+
    scale_colour_manual(values=MOD,guide='none')+
    scale_y_continuous(breaks=seq_along(REG),labels=rev(REG),expand=c(.04,0))+
    scale_x_continuous(limits=c(-.16,.21),breaks=c(-.1,0,.1),labels=num)+
    labs(x='Wald β',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing.x=unit(1.5,'mm'))
}
