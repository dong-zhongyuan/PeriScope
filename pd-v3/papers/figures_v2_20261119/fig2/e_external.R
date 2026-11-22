# Figure 2e - independent human cohorts (GSE157783, GSE184950): one card per cohort, one row per method.
# Small points = Astro and MG-AP disease-direction Spearman, glyph = their mean.
# of the mean. All methods = mean of seeds 42-44. Data: TABLE_benchmark_key_summary.csv.
PANEL<-list(id='2e_external_cohorts',blocks=c('1x3'),what='External human cohorts by method')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/_bench.R'))
make_panel<-function(b){
 k<-read_src('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/corrected_benchmark_20261119/TABLE_benchmark_key_summary.csv')
 pick<-c(periscope='PeriScope seeds 42-44, mean of seeds',conditional_mlp='conditional_mlp, mean of 3 seeds',conditional_icnn='conditional_icnn, mean of 3 seeds',
  shared_label_marginal='shared_label_marginal, mean of 3 seeds',cellot='cellot, mean of 3 seeds',shared_label_mean='shared_label_mean, mean of 3 seeds',ot_ridge='ot_ridge, mean of 3 seeds')
 k<-k[match(pick,k$row),];stopifnot(!anyNA(k$row));k$method<-names(pick)
 co<-c(smajic='Smaji\u0107',gse184950='Wang')
 L<-do.call(rbind,lapply(names(co),function(c)data.frame(method=k$method,cohort=co[[c]],astro=k[[paste0(c,'_astro')]],mg=k[[paste0(c,'_mhc2')]],mean=k[[paste0(c,'_mean')]])))
 L$row<-match(L$method,rev(M_ORDER));L$label<-METHOD_LAB[L$method];L$cohort<-factor(L$cohort,co)
 xr<-range(c(L$astro,L$mg,0));xr<-c(floor(xr[1]*10)/10,ceiling(xr[2]*10)/10)
 P<-ggplot(L,aes(y=row))+
  geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_segment(aes(x=pmin(astro,mg),xend=pmax(astro,mg),yend=row),colour=LIGHT,linewidth=gglw('data'))+
  geom_point(aes(x=astro),shape=16,colour=BRAIN_COL[['Astro']],size=ggsz('small'))+
  geom_point(aes(x=mg),shape=16,colour=BRAIN_COL[['MG-AP']],size=ggsz('small'))+
  geom_point(aes(x=mean,shape=label,fill=label,colour=label),size=ggsz('big'),stroke=gglw('thin')*1.6)+
  facet_wrap(~cohort,nrow=1)+
  scale_shape_manual(values=METHOD_SHAPE,guide='none')+scale_fill_manual(values=METHOD_FILL,guide='none')+scale_colour_manual(values=METHOD_COL,guide='none')+
  scale_y_continuous(breaks=seq_along(M_ORDER),labels=rev(METHOD_LAB[M_ORDER]),limits=c(.5,length(M_ORDER)+.5),expand=expansion(0))+
  scale_x_continuous(limits=xr,breaks=seq(xr[1],xr[2],.1)[abs(seq(xr[1],xr[2],.1)*10)%%2<1e-9],labels=num_minus_acc(.1),expand=expansion(mult=.03))+
  labs(x='Spearman',y=NULL)+theme_fig()+
  theme(panel.spacing=unit(3,'mm'),strip.text=element_text(size=PT,face='bold',hjust=0,margin=margin(0,0,.4,0,'mm')),axis.line.y=element_blank(),axis.ticks.y=element_blank(),
   plot.margin=margin(LETTER_ZONE[['h']]-.6,.6,.6,.6,'mm'))
 W<-b$w;H<-b$h;y1<-H-1.6
 kk<-key_row(list(ki(kg_dot(BRAIN_COL[['Astro']],NA,D[['small']]*1.2),'Astro'),ki(kg_dot(BRAIN_COL[['MG-AP']],NA,D[['small']]*1.2),'MG-AP')),W-24,LETTER_ZONE[['h']]-1.6,h=H)
 gTree(children=gList(gg_grob(P),do.call(gList,kk)))}
