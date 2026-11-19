make_F <- function() {
  d<-read_pd('F_bilateral');stopifnot(nrow(d)==1515*20)
  med<-aggregate(left_minus_right_z~modality,d,median)
  counts<-table(d$modality)
  ggplot(d,aes(left_minus_right_z,modality,fill=modality,colour=modality))+zero_v()+
    geom_violin(orientation='y',scale='width',width=.75,alpha=.2,linewidth=.3/.pt,trim=TRUE)+
    geom_point(position=position_jitter(height=.17,width=0,seed=271),size=.13,alpha=.22,stroke=0)+
    geom_point(data=med,shape=124,size=3.2,colour=INK)+
    scale_fill_manual(values=MOD,guide='none')+scale_colour_manual(values=MOD,guide='none')+
    scale_y_discrete(labels=setNames(paste0(names(counts),'\n',format(as.numeric(counts),big.mark=',')),names(counts)))+
    numeric_x()+labs(x='Left − right z',y=NULL)+
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
}
