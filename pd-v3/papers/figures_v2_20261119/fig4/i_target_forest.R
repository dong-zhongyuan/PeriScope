# Figure 4i - the final targets (CD71/TFRC, CD119/IFNGR1), every MR + colocalisation pair. Left: MR beta with its 95%
# interval per SD of plasma protein; filled = pair-level BH q < 0.05, open = selected by the best-trait/assay rule;
# glyph colour = target, row strip = imaging modality. Right: the unified PP.H4 of the same pair on each protein
# platform that reaches the 0.8 gate (circle = UKB-PPP, diamond = AGES).
# Data: target_funnel_final_20261119/final_target_mr_coloc_pairs.csv.
PANEL<-list(id='4i_final_target_forest',blocks=c('1x2'),what='Final-target MR forest and colocalisation')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 f<-read_src(file.path(TF,'final_target_mr_coloc_pairs.csv'),cols=c('antigen','UKB_gene','trait','phenotype_description',
   'family','beta','se','q_bh_global','unified_H4_AGES','unified_H4_UKB'))
 stopifnot(!anyNA(SHORT[f$trait]));f$lab<-paste0(f$UKB_gene,'  ',unname(SHORT[f$trait]))
 f<-f[order(match(f$UKB_gene,c('TFRC','IFNGR1')),-abs(f$beta/f$se)),];f$y<-rev(seq_len(nrow(f)))
 f$sig<-is.finite(f$q_bh_global)&f$q_bh_global<.05;f$col<-FT_COL[f$UKB_gene]
 YS<-scale_y_continuous(breaks=f$y,labels=f$lab,limits=c(.4,nrow(f)+.6),expand=expansion(0))
 xr<-max(abs(c(f$beta-1.96*f$se,f$beta+1.96*f$se)))*1.05
 P1<-ggplot(f,aes(beta,y))+geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_segment(aes(x=beta-1.96*se,xend=beta+1.96*se,yend=y),colour=f$col,linewidth=gglw('emph'))+
  geom_point(shape=21,fill=ifelse(f$sig,f$col,'white'),colour=f$col,stroke=.45,size=ggsz('big'))+
  YS+scale_x_continuous(labels=num_minus,n.breaks=5)+coord_cartesian(xlim=c(-xr,xr),clip='off')+
  labs(x='MR beta (per SD protein)',y=NULL)+theme_fig()+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(colour=unname(f$col)),
   plot.margin=margin(LETTER_ZONE[['h']],1,.6,.6,'mm'))
 L<-rbind(data.frame(y=f$y,col=f$col,pf='UKB-PPP',v=f$unified_H4_UKB,stringsAsFactors=FALSE),
          data.frame(y=f$y,col=f$col,pf='AGES',v=f$unified_H4_AGES,stringsAsFactors=FALSE))
 L<-L[is.finite(L$v)&L$v>=.8,];stopifnot(all(f$y%in%L$y))
 P2<-ggplot(L,aes(v,y))+geom_vline(xintercept=.8,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_segment(data=aggregate(v~y+col,L,max),aes(x=.8,xend=v,yend=y),colour=FAINT,linewidth=gglw('emph')*1.4)+
  geom_point(aes(shape=pf),fill=L$col,colour='white',stroke=.3,size=ggsz('point'))+
  scale_shape_manual(values=c(`UKB-PPP`=21,AGES=23),guide='none')+
  YS+scale_x_continuous(limits=c(.79,1.01),breaks=c(.8,.9,1),labels=c('0.8','0.9','1.0'))+
  labs(x='PP.H4',y=NULL)+theme_fig()+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),
   plot.margin=margin(LETTER_ZONE[['h']],.6,.6,1,'mm'))
 W<-b$w;H<-b$h;w1<-62
 k<-key_row(list(ki(kg_dot(MID,'white',1.4),'UKB-PPP'),ki(list(w=1.6,f=function(x,yc)ddiam(x+.8,0,1.7,MID,'white',glwd('thin'),h=yc)),'AGES')),
   W-1-key_width(list(ki(kg_dot(MID,'white',1.4),'UKB-PPP'),ki(list(w=1.6,f=function(x,yc)NULL),'AGES'))),2.1,h=H)
 gTree(children=gList(grobTree(gg_grob(P1),vp=viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(w1,'mm'),height=unit(H,'mm'),just=c(0,0))),
  grobTree(gg_grob(P2),vp=viewport(x=unit(w1,'mm'),y=unit(0,'mm'),width=unit(W-w1,'mm'),height=unit(H,'mm'),just=c(0,0))),
  do.call(gList,k)))}
