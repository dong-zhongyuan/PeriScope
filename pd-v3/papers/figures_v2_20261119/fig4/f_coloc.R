# Figure 4f - single-signal colocalisation across all 1,375 protein-imaging regions in genome order. Top: PP.H4 (best of
# the two windows at p12 = 1e-5), H4 >= 0.8 coloured. Bottom: barcode of PP.H4 under the four window x prior settings.
# Data: all_target_delivery_20261119/all_target_region_table.csv.gz (read-only).
PANEL<-list(id='4f_coloc_landscape',blocks=c('1x2'),what='Colocalisation landscape across regions')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 r<-read_src(file.path(MR,'all_target_delivery_20261119/all_target_region_table.csv.gz'),cols=c('region_id','gene','chromosome_label','center_grch37','H4_w1000k_1e5','H4_w500k_1e5','H4_w1000k_1e6','H4_w500k_1e6'))
 r$cn<-suppressWarnings(as.numeric(ifelse(r$chromosome_label%in%c('X','23'),23,r$chromosome_label)));r<-r[order(r$cn,r$center_grch37),];r$x<-seq_len(nrow(r))
 r$h<-pmax(r$H4_w1000k_1e5,r$H4_w500k_1e5,na.rm=TRUE);r$pass<-r$h>=.8
 cb<-aggregate(x~cn,r,function(v)c(min(v),max(v)));cb<-data.frame(cn=cb$cn,lo=cb$x[,1],hi=cb$x[,2]);cb$mid<-(cb$lo+cb$hi)/2
 lb<-cb[cb$cn%in%c(1,5,10,15,20,23),];lb$l<-ifelse(lb$cn==23,'X',lb$cn)
 top<-ggplot(r,aes(x,h))+
  geom_rect(data=cb[cb$cn%%2==0,],aes(xmin=lo-.5,xmax=hi+.5,ymin=-Inf,ymax=Inf),inherit.aes=FALSE,fill=PAPER)+
  geom_hline(yintercept=.8,colour=REF_COL,linewidth=gglw('ref'),linetype=REF_LTY)+
  geom_point(aes(colour=pass),shape=16,size=ggsz('small'))+
  geom_point(data=r[r$gene%in%names(FT_COL)&r$pass,],shape=21,fill=FT_COL[r$gene[r$gene%in%names(FT_COL)&r$pass]],colour='white',stroke=.3,size=ggsz('big'))+
  scale_colour_manual(values=c(`TRUE`=PERISCOPE,`FALSE`=LIGHT),guide='none',na.value=FAINT)+
  scale_x_continuous(breaks=NULL,expand=expansion(0))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','0.5','1'))+
  labs(x=NULL,y='PP.H4')+theme_fig()+theme(axis.line.x=element_blank(),plot.margin=margin(LETTER_ZONE[['h']],.6,0,.6,'mm'))
 top$layers[[1]]$geom_params$na.rm<-TRUE
 st<-c(H4_w1000k_1e5='1 Mb',H4_w500k_1e5='500 kb',H4_w1000k_1e6='1 Mb, strict',H4_w500k_1e6='500 kb, strict')
 B<-do.call(rbind,lapply(names(st),function(k)data.frame(x=r$x,y=match(k,rev(names(st))),v=r[[k]])))
 bot<-ggplot(B,aes(x,y,fill=v))+geom_raster()+
  scale_fill_gradientn(colours=seq_col(seq(0,1,length.out=60),PERISCOPE),limits=c(0,1),guide='none',na.value=MISSING)+
  scale_x_continuous(breaks=lb$mid,labels=lb$l,expand=expansion(0))+scale_y_continuous(breaks=seq_along(st),labels=rev(st),expand=expansion(0))+
  labs(x='Chromosome',y=NULL)+theme_fig()+theme(axis.line=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(.6,.6,.6,.6,'mm'))
 g1<-gg_grob(top);g2<-gg_grob(bot);W<-b$w;H<-b$h;h2<-17
 wmax<-grid::unit.pmax(g1$widths,g2$widths);g1$widths<-wmax;g2$widths<-wmax
 gTree(children=gList(grobTree(g1,vp=viewport(x=unit(0,'mm'),y=unit(h2,'mm'),width=unit(W,'mm'),height=unit(H-h2,'mm'),just=c(0,0))),
  grobTree(g2,vp=viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(W,'mm'),height=unit(h2,'mm'),just=c(0,0)))))}
