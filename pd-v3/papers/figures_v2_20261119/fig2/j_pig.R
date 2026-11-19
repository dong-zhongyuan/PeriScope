# Figure 2j - pig LPS prefrontal-cortex response. Left: measured vs PeriScope-predicted LPS-minus-saline effect for the
# 3,730 one-to-one orthologs (points = genes; line = linear fit). Right: Spearman by method; ridge = 1,000
# treatment-stratified paired animal bootstraps, glyph = full-data estimate. Constant predictors (label mean and
# marginal) have undefined correlation and are not shown. Data: pig_external_validation_20261003 (read-only).
PANEL<-list(id='2j_pig_lps',blocks=c('1x2'),what='Pig LPS response prediction')
P0<-'/public/home/mengxl/dzy/pd_product_assets/results/pig_external_validation_20261003'
make_panel<-function(b){
 g<-read_src(file.path(P0,'gene_level_group_effects.csv'),cols=c('gene','observed_brain_effect','periscope'))
 r<-read_src(file.path(P0,'response_metrics.csv'),cols=c('method','condition','subset','spearman'));r<-r[r$condition==0&r$subset=='all_one2one'&is.finite(r$spearman),]
 bt<-read_src(file.path(P0,'paired_animal_bootstrap.csv'),cols=c('replicate','method','spearman'));bt<-bt[is.finite(bt$spearman),]
 rho<-r$spearman[r$method=='periscope'];stopifnot(nrow(g)==3730)
 sc<-ggplot(g,aes(periscope*1e3,observed_brain_effect))+
  geom_hline(yintercept=0,colour=ZERO_COL,linewidth=gglw('thin'))+geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('thin'))+
  geom_point(shape=16,size=ggsz('unit'),colour=alpha_hex(PERISCOPE,.28))+
  geom_smooth(method='lm',formula=y~x,se=FALSE,colour=INK,linewidth=gglw('data'))+
  scale_x_continuous(labels=num_minus)+scale_y_continuous(labels=num_minus_acc(.1))+
  labs(x='1,000 \u00d7 predicted',y='Measured')+theme_fig()
 ord<-c('periscope','conditional_mlp','cellot','conditional_icnn','ot_ridge');bt<-bt[bt$method%in%ord,];r<-r[r$method%in%ord,]
 r$row<-match(r$method,rev(ord));r$label<-METHOD_LAB[r$method]
 dens<-do.call(rbind,lapply(ord,function(m){v<-bt$spearman[bt$method==m];d<-stats::density(v,n=256,adjust=1.1);data.frame(method=m,x=d$x,h=d$y/max(d$y)*.78,row=match(m,rev(ord)))}))
 dens$fill<-ifelse(dens$method=='periscope',PERISCOPE_LIGHT,FAINT)
 fo<-ggplot()+
  geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_ribbon(data=dens,aes(x,ymin=row-.3,ymax=row-.3+h,group=method,fill=fill),colour=NA)+
  geom_point(data=r,aes(spearman,row-.3),shape=unname(METHOD_SHAPE[r$label]),fill=unname(METHOD_FILL[r$label]),colour=unname(METHOD_COL[r$label]),size=ggsz('point')*1.15,stroke=.3)+
  scale_fill_identity()+
  scale_y_continuous(breaks=seq_along(ord)-.3,labels=rev(METHOD_LAB[ord]),expand=expansion(add=c(.35,.6)))+
  scale_x_continuous(breaks=c(-.2,0,.2,.4),labels=num_minus_acc(.1))+
  labs(x='Spearman',y=NULL)+theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
 W<-b$w;H<-b$h;w1<-37
 gTree(children=gList(grobTree(gg_grob(sc),vp=viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(w1,'mm'),height=unit(H,'mm'),just=c(0,0))),
  grobTree(gg_grob(fo),vp=viewport(x=unit(w1+1,'mm'),y=unit(0,'mm'),width=unit(W-w1-1,'mm'),height=unit(H,'mm'),just=c(0,0)))))}
