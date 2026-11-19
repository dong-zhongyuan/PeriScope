# Figure 1j - final targets CD71/TFRC and CD119/IFNGR1: within-group Spearman (fill) and cell count (size) for each
# held-out donor x subtype group, plus the across-group agreement of mean abundance (right column).
# Data: protein_accuracy.csv and diag/protein_across_group_accuracy.csv (read-only).
PANEL<-list(id='1j_final_target_adt',blocks=c('1x2'),what='Final-target protein prediction by donor and subtype')
FRD<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design/manuscript_review/figure_redesign_20261004/data'
make_panel<-function(b){
 d<-read_src(file.path(FRD,'protein_accuracy.csv'),cols=c('target','citeseq_subtype','donor','n','rho'))
 a<-read_src('/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/diag/protein_across_group_accuracy.csv')
 tg<-names(TARGET_COL);d<-d[d$target%in%tg,]
 ord<-c('CD14 Mono','CD16 Mono','cDC1','cDC2','pDC');st<-c(`CD14 Mono`='cMono',`CD16 Mono`='nMono',cDC1='cDC',cDC2='cDC',pDC='pDC')
 u<-unique(d[,c('citeseq_subtype','donor')]);u<-u[order(match(u$citeseq_subtype,ord),u$donor),];u$x<-seq_len(nrow(u));u$state<-st[u$citeseq_subtype]
 d<-merge(d,u);d$y<-match(d$target,rev(tg))
 a<-a[a$target%in%tg,];a$y<-match(a$target,rev(tg));xr<-nrow(u)+1.6
 grp<-aggregate(x~state,u,range);grp<-data.frame(s=grp$state,lo=grp$x[,1],hi=grp$x[,2]);div<-u$x[u$citeseq_subtype=='cDC1']+.5
 ggplot(d,aes(x,y))+
  geom_point(aes(size=n,fill=rho),shape=21,colour='white',stroke=.25)+
  geom_segment(data=grp,aes(x=lo-.35,xend=hi+.35,y=.42,yend=.42),inherit.aes=FALSE,colour=BLOOD_COL[grp$s],linewidth=gglw('emph'))+
  annotate('segment',x=div,xend=div,y=.6,yend=2.4,colour=FAINT,linewidth=gglw('thin'))+
  geom_text(data=grp,aes(x=(lo+hi)/2,y=.2,label=s),inherit.aes=FALSE,family=FONT,size=GT,colour=INK)+
  scale_x_continuous(breaks=u$x,labels=u$donor,limits=c(.4,nrow(u)+.6),expand=expansion(0))+
  scale_y_continuous(breaks=seq_along(tg),labels=target_label(rev(tg)),limits=c(.05,2.95),expand=expansion(0))+
  scale_size_area(max_size=ggsz(4.2),breaks=c(100,1000,5000),labels=num_comma,name='Cells')+
  scale_fill_gradientn(colours=seq_col(seq(0,1,length.out=65),PERISCOPE),limits=c(0,.6),breaks=c(0,.3,.6),labels=num_minus,oob=scales::squish,name='Spearman')+
  coord_cartesian(clip='off')+labs(x=NULL,y=NULL)+theme_fig(legend='right')+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.title=element_text(size=PT),legend.box='vertical',legend.justification=c(0,0),legend.box.margin=margin(0,0,0,1.5,'mm'))+
  guides(fill=guide_bar(FALSE,8),size=guide_legend(override.aes=list(fill=MID)))}
