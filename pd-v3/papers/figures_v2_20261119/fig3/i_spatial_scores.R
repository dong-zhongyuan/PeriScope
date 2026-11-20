# Figure 3i - spatially supported programs in fibrous astrocytes: donor program score (z) per section, control | PD,
# condition bar above. Rows = programs (antigen, blood source). Data: spatial_heat_exact_matrix.csv,
# spatial_heat_annotation.csv (read-only).
PANEL<-list(id='3i_spatial_program_scores',blocks=c('1x1'),what='Spatial program scores per section')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
make_panel<-function(b){
 m<-read_src(file.path(FRD,'spatial_heat_exact_matrix.csv'));rn<-m[[1]];m<-m[,-1];names(m)<-gsub('.','-',names(m),fixed=TRUE)
 sp<-read_src(file.path(FRD,'spatial_fractions.csv'),cols=c('donor','condition'));cond<-tapply(sp$condition,sp$donor,function(z)z[1]);cond<-setNames(as.character(cond),names(cond))
 pp<-prog_parse(rn);lab<-paste0(pp$antigen,' ',pp$blood,' → ',pp$brain,' ',arrow_of(pp$dir))
 L<-do.call(rbind,lapply(seq_along(rn),function(i)data.frame(row=length(rn)-i+1,donor=names(m),v=unlist(m[i,]))))
 L$cond<-cond[L$donor];ord<-c(sort(names(cond)[cond=='Control']),sort(names(cond)[cond=='PD']));L$x<-match(L$donor,ord)+ifelse(L$cond=='PD',.4,0)
 lim<-max(abs(L$v));cb<-unique(L[,c('x','cond')])
 p<-ggplot(L,aes(x,row,fill=v))+geom_tile(width=.92,height=.9)+
  geom_tile(data=cb,aes(x,length(rn)+.75,fill=NULL),fill=DX_COL[cb$cond],width=.92,height=.32,inherit.aes=FALSE)+
  scale_fill_gradientn(colours=c('#3F6E9A','white',PERISCOPE),limits=c(-lim,lim),guide='none')+
  scale_y_continuous(breaks=seq_along(rn),labels=rev(lab),expand=expansion(0))+scale_x_continuous(breaks=NULL,expand=expansion(0))+
  coord_cartesian(ylim=c(.5,length(rn)+1),clip='off')+labs(x=NULL,y=NULL)+theme_fig()+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),plot.margin=margin(LETTER_ZONE[['h']]+.4,.6,4.4,2.4,'mm'))
 W<-b$w;H<-b$h;y1<-H-1.8
 k<-c(list(textGrob('Score',unit(.6,'mm'),unit(1.8,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,fontface='bold',col=INK))),
  key_row(list(ki(kg_grad(c('#3F6E9A','white',PERISCOPE),'−','+',w=10))),.6+sw('Score',TRUE)+1,y1,h=H))
 gTree(children=gList(gg_grob(p),do.call(gList,k)))}
