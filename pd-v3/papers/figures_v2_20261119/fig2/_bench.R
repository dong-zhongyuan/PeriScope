# shared loader for Figure 2 b-d: human disease-direction benchmark on the locked donors (8 routes). PeriScope = full
# model seeds 42-51 (corrected_benchmark src/eval/per_axis.csv); comparators = manuscript runs, seeds 42-44
# (human_direction_metrics.csv). Read-only sources.
FRD<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design/manuscript_review/figure_redesign_20261004/data'
PAX<-'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/corrected_benchmark_20261119/src/eval/per_axis.csv'
AX_BLOOD<-c(cDC='cDC',classical_mono='cMono',nonclassical_mono='nMono',pDC='pDC');AX_BRAIN<-c(astro='Astro',microglia_mhc2='MG-AP')
M_ORDER<-c('periscope','conditional_mlp','conditional_icnn','shared_label_marginal','cellot','shared_label_mean','ot_ridge')
route_of<-function(axis){p<-do.call(rbind,strsplit(axis,'__'));data.frame(blood=AX_BLOOD[p[,1]],brain=AX_BRAIN[p[,2]],stringsAsFactors=FALSE)}
bench_long<-function(metric='spearman'){
 h<-read_src(file.path(FRD,'human_direction_metrics.csv'),cols=c('axis','seed','method','spearman','pearson','sign_agreement'))
 h<-h[h$method!='periscope',];h$value<-h[[metric]]
 if(metric=='spearman'){p<-read_src(PAX,cols=c('arm','seed','axis','rho'));p<-p[p$arm=='full',];stopifnot(length(unique(p$seed))==10)
  p<-data.frame(axis=p$axis,seed=p$seed,method='periscope',value=p$rho)}else{
  p<-read_src(file.path(FRD,'human_direction_metrics.csv'),cols=c('axis','seed','method',metric));p<-p[p$method=='periscope',];p$value<-p[[metric]];p<-p[,c('axis','seed','method','value')]}
 d<-rbind(h[,c('axis','seed','method','value')],p);d<-cbind(d,route_of(d$axis))
 d$label<-METHOD_LAB[d$method];d$row<-match(d$method,rev(M_ORDER));stopifnot(!anyNA(d$row))
 d$col<-match(d$blood,AX_BLOOD)+ifelse(d$brain=='MG-AP',4.35,0);d}
COLX<-c(1:4,5.35:8.35);MEANX<-9.75
BL<-unname(AX_BLOOD)
div_cols<-function(lim,n=80){f<-scales::div_gradient_pal('#3F6E9A','white',PERISCOPE);v<-seq(lim[1],lim[2],length.out=n);f(scales::rescale_mid(v,to=c(0,1),from=lim,mid=0))}
KEY_BOT<-7.6   # bottom margin (mm) reserved for hand-drawn keys, identical in 2b and 2c so their tile rows align
key_title<-function(t,x,y,h)textGrob(t,unit(x,'mm'),unit(h-y,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,fontface='bold',col=INK))
short2<-function(v)sub('^(−?)0\\.','\\1.',num_minus_acc(.01)(v))
# column header drawn inside the panel (identical in 2b and 2c): rotated blood-state labels, then the brain-state rule and name
Y_TOP<-length(M_ORDER)+.5;Y_EXP<-c(.5,2.25)
col_header<-function(names=TRUE){hx<-data.frame(x=c(2.5,6.85),lab=c('Astro','MG-AP'),col=BRAIN_COL[c('Astro','MG-AP')])
 cl<-data.frame(x=c(COLX,MEANX),lab=c(rep(BL,2),'Mean'))
 list(geom_text(data=cl,aes(x,Y_TOP+.12,label=lab),inherit.aes=FALSE,family=FONT,size=GT,angle=90,hjust=0,vjust=.5,colour=INK),
  geom_segment(data=hx,aes(x=x-1.95,xend=x+1.95,y=Y_TOP+1.22,yend=Y_TOP+1.22),inherit.aes=FALSE,colour=hx$col,linewidth=gglw('emph')*1.4),
  if(names)geom_text(data=hx,aes(x,Y_TOP+1.36,label=lab),inherit.aes=FALSE,family=FONT,size=GT,colour=INK,vjust=0))}
bench_theme<-function()theme_fig()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),
 plot.margin=margin(LETTER_ZONE[['h']]+.4,.6,KEY_BOT,.6,'mm'))
# ggplot2 4 gives the (blank) bottom-axis row a share of the free space; give it all to the panel so 2b and 2c fill and align
fix_null_rows<-function(G){pr<-unique(G$layout$t[grepl('^panel',G$layout$name)])
 for(i in seq_along(G$heights))if(!(i%in%pr)&&grepl('null',as.character(G$heights[i])))G$heights[i]<-unit(0,'mm');G}
# fixed plotting window: the in-panel header layers must not widen the y range (2b and 2c share it)
bench_coord<-function()coord_cartesian(xlim=c(.5,MEANX+.5),ylim=c(.5,length(M_ORDER)+Y_EXP[2]),expand=FALSE,clip='off')
