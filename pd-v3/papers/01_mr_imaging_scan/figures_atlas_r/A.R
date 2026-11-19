make_A <- function() {
  d<-read_table('A_all_assays');chr<-read_table('A_chromosomes')
  stopifnot(nrow(d)==cfg$assays*2)
  d<-d[d$n_estimated>0,];d$policy_label<-policy_factor(d$policy)
  ticks<-chr[chr$chrom %in% c(seq(1,19,2),23),]
  n<-do.call(rbind,lapply(c('drop','infer'),function(k) {
    data.frame(policy_label=policy_factor(k),x=Inf,y=Inf,
      label=paste0(format(cfg$policies[[k]]$estimated_assays,big.mark=','),' estimated; ',
                   format(cfg$policies[[k]]$unavailable_assays,big.mark=','),' unavailable'))
  }))
  ggplot(d,aes(genome_order,minus_log10_p,colour=policy))+
    geom_rect(data=chr[seq_len(nrow(chr))%%2==0,],inherit.aes=FALSE,
      aes(xmin=min-.5,xmax=max+.5,ymin=-Inf,ymax=Inf),fill='#F3F3F3',colour=NA)+
    geom_point(size=.3,alpha=.65,stroke=0)+
    geom_text(data=n,inherit.aes=FALSE,aes(x=x,y=y,label=label),hjust=1,vjust=1,
      family=FONT,size=6/.pt,colour=INK)+
    facet_wrap(~policy_label,nrow=1)+scale_colour_manual(values=PCOL,guide='none')+
    scale_x_continuous(breaks=ticks$mean,labels=ticks$label,expand=expansion(mult=c(.01,.01)))+
    scale_y_continuous(expand=expansion(mult=c(.02,.2)))+
    labs(x='Cis-tool chromosome (within-chromosome order)',y='Best-pair -log10 P')+
    theme(panel.spacing.x=unit(4,'mm'))
}
