make_E <- function() {
  d<-read_pd('E_cross_modality');stopifnot(nrow(d)==1515*8*2)
  d$region_short<-factor(d$region_short,levels=REG)
  hi<-d[d$gene%in%names(shape_gene),]
  ggplot(d,aes(z_qsm,z_t2star))+zero_v()+zero_h()+
    geom_point(size=.12,alpha=.22,colour='#777777',stroke=0)+
    geom_point(data=hi,aes(shape=gene,fill=hemisphere),size=1.1,stroke=.3,colour=INK)+
    facet_wrap(~region_short,nrow=2)+scale_shape_manual(values=shape_gene,guide='none')+
    scale_fill_manual(values=c(L=INK,R='white'),guide='none')+
    scale_x_continuous(limits=c(-10,10),breaks=c(-8,0,8),labels=num,expand=c(.015,0))+
    scale_y_continuous(limits=c(-10,10),breaks=c(-8,0,8),labels=num,expand=c(.015,0))+
    coord_fixed(1)+labs(x='QSM z',y='T2* z')+theme(panel.spacing=unit(1.4,'mm'),strip.text=element_text(size=6.5,face='bold'))
}
