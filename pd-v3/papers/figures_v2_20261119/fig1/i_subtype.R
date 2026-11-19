# Figure 1i - within-subtype Spearman correlation of predicted vs measured ADT for every non-control antigen, one row per
# held-out donor x subtype (13 rows); points = antigens (10-model means), tick = median, right column = median value.
# Data: figure_redesign_20261004/data/protein_accuracy.csv (read-only).
PANEL<-list(id='1i_subtype_spearman',blocks=c('2x1'),what='Within-subtype protein correlation by donor and subtype')
FRD<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design/manuscript_review/figure_redesign_20261004/data'
make_panel<-function(b){
 d<-read_src(file.path(FRD,'protein_accuracy.csv'),cols=c('target','blood_state','citeseq_subtype','donor','n','is_control','rho'))
 d<-d[!d$is_control&is.finite(d$rho),]
 sub_lab<-c(`CD14 Mono`='CD14 Mono',`CD16 Mono`='CD16 Mono',cDC1='cDC1',cDC2='cDC2',pDC='pDC')
 sub_state<-c(`CD14 Mono`='cMono',`CD16 Mono`='nMono',cDC1='cDC',cDC2='cDC',pDC='pDC')
 ord<-c('CD14 Mono','CD16 Mono','cDC1','cDC2','pDC')
 u<-unique(d[,c('citeseq_subtype','donor')]);u<-u[order(match(u$citeseq_subtype,ord),u$donor),]
 u$row<-rev(seq_len(nrow(u)));u$lab<-paste(u$donor)
 d<-merge(d,u,by=c('citeseq_subtype','donor'))
 set.seed(1);d$y<-d$row+stats::runif(nrow(d),-.32,.32)
 med<-aggregate(rho~row+citeseq_subtype+donor,d,stats::median)
 strip<-aggregate(row~citeseq_subtype,u,range);strip<-data.frame(sub=strip$citeseq_subtype,lo=strip$row[,1],hi=strip$row[,2])
 strip$col<-BLOOD_COL[sub_state[strip$sub]]
 x_med<-1.12;x_strip<- -0.78
 ggplot()+
  geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_point(data=d,aes(rho,y),shape=16,colour=alpha_hex(MID,.35),size=ggsz('unit'))+
  geom_segment(data=med,aes(x=rho,xend=rho,y=row-.42,yend=row+.42),colour=PERISCOPE,linewidth=gglw('mean'))+
  geom_rect(data=strip,aes(xmin=x_strip,xmax=x_strip+.06,ymin=lo-.42,ymax=hi+.42),fill=strip$col,colour=NA)+
  geom_text(data=strip,aes(x=x_strip+.15,y=(lo+hi)/2,label=sub_lab[sub]),family=FONT,size=GT,hjust=.5,vjust=.5,colour=INK,angle=90)+
  scale_y_continuous(breaks=u$row,labels=u$lab,expand=expansion(add=.6))+
  scale_x_continuous(limits=c(-.78,.98),breaks=c(-.5,0,.5),labels=num_minus_acc(.1),expand=expansion(0))+
  labs(x='Within-subtype Spearman',y=NULL)+theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())}
