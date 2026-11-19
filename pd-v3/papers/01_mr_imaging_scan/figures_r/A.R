make_A <- function() {
  d<-read_pd('A_landscape'); stopifnot(nrow(d)==68175)
  ch<-read_pd('A_chromosomes')
  labs<-do.call(rbind,lapply(split(d,d$modality),function(x) {
    x<-x[order(-abs(x$z)),];x<-x[!duplicated(x$gene),];head(x,3)
  }))
  labs$label_y<-labs$z-ifelse(labs$z>0,1.6,-1.6)
  ticks<-ch[ch$chr_n%in%c(1,3,5,7,9,11,13,15,17,19,23),]
  ggplot(d,aes(assay_order,z)) +
    geom_rect(data=ch[ch$chr_n%%2==0,],aes(xmin=min-.5,xmax=max+.5,ymin=-10,ymax=10),
              inherit.aes=FALSE,fill='#F4F4F4',colour=NA)+zero_h()+
    geom_point(aes(colour=modality),size=.16,alpha=.35,stroke=0)+
    geom_segment(data=labs,aes(xend=assay_order,yend=label_y),linewidth=.25/.pt,colour='#555555')+
    geom_point(data=labs,aes(fill=modality),shape=21,size=.85,stroke=.15,colour=INK)+
    geom_text(data=labs,aes(y=label_y,label=gene),size=6/.pt,family=FONT,colour=INK)+
    facet_grid(modality~.)+scale_colour_manual(values=MOD,guide='none')+scale_fill_manual(values=MOD,guide='none')+
    scale_x_continuous(breaks=ticks$median,labels=ifelse(ticks$chr_n==23,'X',ticks$chr_n),expand=c(.004,0))+
    scale_y_continuous(limits=c(-10,10),breaks=c(-8,0,8),labels=num,expand=c(0,0))+
    labs(x='Cis-tool chromosome',y='Wald z')+
    theme(panel.spacing.y=unit(1,'mm'),strip.text.y=element_text(angle=0,size=6.5,face='bold'))
}
