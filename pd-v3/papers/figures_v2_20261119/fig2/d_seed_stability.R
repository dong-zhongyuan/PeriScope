# Figure 2d - seed stability on each route: the ten PeriScope seeds (42-51) as points, the best comparator's
# seed-mean on that route as a grey tick, right column = seeds above that tick. Data as 2b.
PANEL<-list(id='2d_seed_stability',blocks=c('1x1'),what='Seed stability per route against the best comparator')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/_bench.R'))
make_panel<-function(b){d<-bench_long('spearman');d$route<-paste(d$blood,d$brain)
 ord<-as.vector(outer(BL,c('Astro','MG-AP'),paste));d$y<-match(d$route,rev(ord))
 p<-d[d$method=='periscope',];c_<-aggregate(value~method+route+y,d[d$method!='periscope',],mean)
 best<-do.call(rbind,lapply(split(c_,c_$route),function(z)z[which.max(z$value),]))
 kk<-do.call(rbind,lapply(split(p,p$route),function(z){bv<-best$value[best$route==z$route[1]];data.frame(route=z$route[1],y=z$y[1],k=sum(z$value>bv))}))
 set.seed(2);p$yj<-p$y+stats::runif(nrow(p),-.22,.22)
 bc<-BRAIN_COL[sub('^\\S+ ','',ord)];lab<-sub(' (Astro|MG-AP)$','',rev(ord))
 ggplot()+
  geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_segment(data=best,aes(x=value,xend=value,y=y-.38,yend=y+.38),colour=MID,linewidth=gglw('emph'))+
  geom_point(data=p,aes(value,yj),shape=21,fill=PERISCOPE,colour='white',stroke=.2,size=ggsz('point'))+
  annotate('segment',x=-.4,xend=-.4,y=c(.6,4.6),yend=c(4.4,8.4),colour=BRAIN_COL[c('MG-AP','Astro')],linewidth=gglw('emph')*1.6)+
  scale_y_continuous(breaks=1:8,labels=lab,expand=expansion(add=.5))+
  scale_x_continuous(limits=c(-.4,.55),breaks=c(0,.5),labels=c('0','0.5'),expand=expansion(0))+
  coord_cartesian(clip='off')+labs(x='Spearman',y=NULL)+theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())}
