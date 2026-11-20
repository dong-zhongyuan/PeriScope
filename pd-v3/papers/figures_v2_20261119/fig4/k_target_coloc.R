# Figure 4k - colocalisation of the final-target MR pairs: PP.H4 (higher of single-signal and SuSiE) on each protein
# platform that reaches the 0.8 gate (circle = UKB-PPP, diamond = AGES); one row per pair, colour = target.
# Data: target_funnel_final_20261119/final_target_mr_coloc_pairs.csv.
PANEL<-list(id='4k_target_coloc',blocks=c('1x1'),what='Final-target coloc by protein platform')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 f<-read_src(file.path(TF,'final_target_mr_coloc_pairs.csv'),cols=c('UKB_gene','trait','unified_H4_AGES','unified_H4_UKB'))
 f$lab<-unname(SHORT[f$trait]);stopifnot(!anyNA(f$lab))
 f<-f[order(match(f$UKB_gene,c('TFRC','IFNGR1')),-f$unified_H4_UKB),];f$y<-rev(seq_len(nrow(f)))
 L<-rbind(data.frame(y=f$y,g=f$UKB_gene,pf='UKB-PPP',v=f$unified_H4_UKB),data.frame(y=f$y,g=f$UKB_gene,pf='AGES',v=f$unified_H4_AGES))
 L<-L[is.finite(L$v)&L$v>=.8,];stopifnot(all(f$y%in%L$y))
 st<-aggregate(v~y+g,L,max);brk<-min(f$y[f$UKB_gene=='TFRC'])-.5
 P<-ggplot()+
  geom_segment(data=st,aes(x=.8,xend=v,y=y,yend=y),colour=FAINT,linewidth=gglw('emph')*1.6)+
  geom_point(data=L,aes(v,y,shape=pf),fill=FT_COL[L$g],colour='white',stroke=.3,size=ggsz('point'))+
  annotate('segment',x=.8,xend=1,y=brk,yend=brk,colour=FAINT,linewidth=gglw('thin'))+
  scale_shape_manual(values=c(`UKB-PPP`=21,AGES=23),guide='none')+
  scale_x_continuous(limits=c(.8,1),breaks=c(.8,.9,1),labels=c('0.8','0.9','1'),expand=expansion(add=c(0,.012)))+
  scale_y_continuous(breaks=f$y[order(f$y)],labels=f$lab[order(f$y)],expand=expansion(add=.6))+
  labs(x='PP.H4',y=NULL)+theme_fig()+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(colour=unname(FT_COL[f$UKB_gene[order(f$y)]])),
   plot.margin=margin(LETTER_ZONE[['h']]+.8,.6,.6,.6,'mm'))
 W<-b$w;H<-b$h
 kd<-function(pch)list(w=1.6,f=function(x,yc)pointsGrob(unit(x+.8,'mm'),unit(yc,'mm'),pch=pch,size=unit(1.5,'mm'),gp=gpar(col='white',fill=MID,lwd=glwd('thin'))))
 k<-key_row(list(ki(kd(21),'UKB-PPP'),ki(kd(23),'AGES')),6,LETTER_ZONE[['h']]-1.6,h=H)
 gTree(children=gList(gg_grob(P),do.call(gList,k)))}
