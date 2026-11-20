# Figure 4b - scale and yield of the MR atlas at three levels (protein assays, imaging phenotypes, protein-phenotype pairs):
# all (ring), estimable under the primary palindrome policy (light) and carrying a pair at full-family BH q < 0.05 (filled);
# log axis. Data: fig4/data/atlas_scale.csv (counts taken from full_atlas_20261116/full).
PANEL<-list(id='4b_scale',blocks=c('1x1'),what='Scale and yield of the MR atlas')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 d<-read_src(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/data/atlas_scale.csv'),cols=c('level','stage','count'))
 lv<-c('Assays','Phenotypes','Pairs');d$y<-match(d$level,rev(lv));stopifnot(!anyNA(d$y))
 st<-c(all='All',estimable='Estimable',significant='q < 0.05')
 FILL<-c(all=NA,estimable=PERISCOPE_LIGHT,significant=PERISCOPE);COL<-c(all=MID,estimable='white',significant='white')
 SZ<-c(all=ggsz('big'),estimable=ggsz('big'),significant=ggsz('big'))
 d<-d[order(match(d$stage,c('estimable','significant','all'))),]
 rg<-do.call(rbind,lapply(split(d,d$y),function(z)data.frame(y=z$y[1],lo=min(z$count),hi=max(z$count))))
 P<-ggplot()+geom_segment(data=rg,aes(x=lo,xend=hi,y=y,yend=y),colour=FAINT,linewidth=gglw('emph')*1.6)+
  geom_point(data=d,aes(count,y),shape=21,fill=FILL[d$stage],colour=COL[d$stage],stroke=.45,size=SZ[d$stage])+
  scale_x_log10(breaks=10^c(2,4,6),labels=c('100','10k','1M'),limits=c(70,3e7))+
  scale_y_continuous(breaks=seq_along(lv),labels=rev(lv),limits=c(.5,3.5),expand=expansion(0))+
  labs(x='Count',y=NULL)+theme_fig()+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(LETTER_ZONE[['h']]+1.2,.6,.6,.6,'mm'))
 W<-b$w;H<-b$h
 k<-key_row(lapply(names(st),function(s)ki(kg_dot(if(is.na(FILL[[s]]))'white' else FILL[[s]],COL[[s]],1.6),st[[s]])),6,LETTER_ZONE[['h']]-1.6,h=H)
 gTree(children=gList(gg_grob(P),do.call(gList,k)))}
