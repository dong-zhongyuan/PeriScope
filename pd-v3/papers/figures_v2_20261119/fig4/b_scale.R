# Figure 4b - scale of the genetic scan relative to the earlier targeted scan: genes, protein assays, imaging phenotypes and
# MR tests (one palindrome policy), log axis, open = earlier scan, filled = this atlas.
# Data: legacy_comparison.csv; full atlas A_all_assays.csv and phenotypes.csv.
PANEL<-list(id='4b_scale',blocks=c('1x1'),what='Scale of the MR atlas vs the earlier scan')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 L<-read_src(file.path(MR,'legacy_comparison.csv'),cols=c('gene','phenotype'))
 A<-read_src(file.path(OV,'A_all_assays.csv'),cols=c('exposure_id','gene','n_planned','policy'));A<-A[A$policy=='drop',]
 ph<-read_src(file.path(MR,'full_atlas_20261116/full/phenotypes.csv'),cols=c('phenotype'))
 d<-data.frame(item=c('Genes','Assays','Phenotypes','MR tests'),old=c(length(unique(L$gene)),length(unique(L$gene)),length(unique(L$phenotype)),nrow(L)),
  new=c(length(unique(A$gene)),nrow(A),nrow(ph),sum(A$n_planned)));d$y<-rev(seq_len(nrow(d)))
 ggplot(d)+geom_segment(aes(x=old,xend=new,y=y,yend=y),colour=PERISCOPE_LIGHT,linewidth=gglw('emph')*1.6)+
  geom_point(aes(old,y),shape=21,fill='white',colour=MID,stroke=.4,size=ggsz('point'))+
  geom_point(aes(new,y),shape=21,fill=PERISCOPE,colour='white',stroke=.3,size=ggsz('big'))+
  scale_x_log10(breaks=10^c(1,3,5,7),labels=c('10','1,000','100,000','10,000,000'),limits=c(8,4e7))+
  scale_y_continuous(breaks=d$y,labels=d$item,expand=expansion(add=.5))+
  labs(x='Count',y=NULL)+theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(angle=30,hjust=1,vjust=1))}
