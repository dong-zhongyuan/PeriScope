# Figure 3h - inferred RNA fractions of fibrous astrocytes (top) and MG-AP (bottom) at the original spot coordinates of
# the ten substantia nigra sections (GSE253975; five control, five PD, condition bar above). Shared 0-1 colour scale on a
# square-root transform. Data: spatial_fractions.csv (read-only).
PANEL<-list(id='3h_spatial_maps',blocks=c('1x4'),what='Spatial RNA fractions of the two glial subtypes')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
make_panel<-function(b){
 d<-read_src(file.path(FRD,'spatial_fractions.csv'),cols=c('donor','condition','pxl_col_in_fullres','pxl_row_in_fullres','astro_fibrous','microglia_mhc2'))
 don<-unique(d[,c('donor','condition')]);don<-don[order(don$condition!='Control',don$donor),];d$donor<-factor(d$donor,don$donor)
 nz<-function(z){x<-z$pxl_col_in_fullres;y<- -z$pxl_row_in_fullres;r<-max(diff(range(x)),diff(range(y)));z$u<-(x-mean(range(x)))/r;z$v<-(y-mean(range(y)))/r;z}
 d<-do.call(rbind,lapply(split(d,d$donor),nz))
 L<-rbind(data.frame(d[,c('donor','condition','u','v')],state='Fibrous astro',f=d$astro_fibrous),data.frame(d[,c('donor','condition','u','v')],state='MG-AP',f=d$microglia_mhc2))
 L$state<-factor(L$state,c('Fibrous astro','MG-AP'))
 cols<-c('#F7F4F0','#E8C9A0','#D9A441','#B4532F','#5E2A1F')
 cb<-don;cb$donor<-factor(cb$donor,levels(d$donor));cb$state<-factor('Fibrous astro',levels(L$state))
 p<-ggplot(L,aes(u,v))+geom_point(aes(colour=sqrt(pmin(pmax(f,0),1))),shape=16,size=ggsz(.42))+
  geom_segment(data=cb,aes(x=-.5,xend=.5,y=.6,yend=.6),inherit.aes=FALSE,colour=unname(DX_COL[cb$condition]),linewidth=gglw('emph')*2.6)+
  facet_grid(state~donor,switch='y')+
  scale_colour_gradientn(colours=cols,limits=c(0,1),guide='none',na.value=MISSING)+
  coord_equal(xlim=c(-.52,.52),ylim=c(-.52,.66),expand=FALSE,clip='off')+
  theme_fig()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),
   strip.text.x=element_blank(),strip.text.y.left=element_text(angle=90,size=PT,face='bold',margin=margin(0,.8,0,0,'mm')),strip.placement='outside',panel.spacing.x=unit(.6,'mm'),panel.spacing.y=unit(1,'mm'),
   plot.margin=margin(LETTER_ZONE[['h']]-.4,14,.4,.6,'mm'))
 W<-b$w;H<-b$h;assign('PANEL_H',H,envir=globalenv())
 gr<-rasterGrob(matrix(rev(grDevices::colorRampPalette(cols)(80)),ncol=1),x=unit(W-11.4,'mm'),y=unit(H/2-3,'mm'),width=unit(1.6,'mm'),height=unit(20,'mm'),interpolate=TRUE)
 tk<-c(0,.25,1);yy<-H/2-3-10+sqrt(tk)*20
 k<-list(gr,rectGrob(unit(W-11.4,'mm'),unit(H/2-3,'mm'),width=unit(1.6,'mm'),height=unit(20,'mm'),gp=gpar(fill=NA,col=INK,lwd=glwd('hair'))))
 for(i in seq_along(tk))k<-c(k,list(textGrob(c('0','0.25','1')[i],unit(W-10.1,'mm'),unit(yy[i],'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,col=INK))))
 k<-c(k,list(textGrob('Fraction',unit(W-11.4,'mm'),unit(H/2-3+11.6,'mm'),just=c('centre','bottom'),gp=gpar(fontfamily=FONT,fontsize=PT,fontface='bold',col=INK))))
 k<-c(k,key_row(list(ki(kg_box(DX_COL[['Control']],w=2.2,ht=1.2),'Ctrl')),W-12.6,H-5.6,h=H),key_row(list(ki(kg_box(DX_COL[['PD']],w=2.2,ht=1.2),'PD')),W-12.6,H-2.6,h=H))
 gTree(children=gList(gg_grob(p),do.call(gList,k)))}
