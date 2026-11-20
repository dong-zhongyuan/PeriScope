# Figure 2d - seed stability on each route: the ten PeriScope seeds (42-51) as points, their mean as a bar and their
# range as a line. Data as 2b.
PANEL<-list(id='2d_seed_stability',blocks=c('1x1'),what='Seed stability of PeriScope on each route')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/_bench.R'))
make_panel<-function(b){d<-bench_long('spearman');d$route<-paste(d$blood,d$brain)
 ord<-as.vector(outer(BL,c('Astro','MG-AP'),paste))
 p<-d[d$method=='periscope',];stopifnot(all(table(p$route)==10),length(unique(p$route))==8)
 p$y<-match(p$route,rev(ord))
 s<-do.call(rbind,lapply(split(p,p$route),function(z)data.frame(y=z$y[1],m=mean(z$value),lo=min(z$value),hi=max(z$value))))
 set.seed(2);p$yj<-p$y+stats::runif(nrow(p),-.22,.22)
 lab<-sub(' (Astro|MG-AP)$','',rev(ord))
 ggplot()+
  geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_segment(data=s,aes(x=lo,xend=hi,y=y,yend=y),colour=PERISCOPE_LIGHT,linewidth=gglw('emph'))+
  geom_point(data=p,aes(value,yj),shape=21,fill=PERISCOPE,colour='white',stroke=.2,size=ggsz('point'))+
  geom_segment(data=s,aes(x=m,xend=m,y=y-.38,yend=y+.38),colour=INK,linewidth=gglw('mean'))+
  annotate('segment',x=-.06,xend=-.06,y=c(.6,4.6),yend=c(4.4,8.4),colour=BRAIN_COL[c('MG-AP','Astro')],linewidth=gglw('emph')*1.6)+
  scale_y_continuous(breaks=1:8,labels=lab,expand=expansion(add=.5))+
  scale_x_continuous(limits=c(-.06,.45),breaks=c(0,.2,.4),labels=c('0','0.2','0.4'),expand=expansion(0))+
  coord_cartesian(clip='off')+labs(x='Spearman',y=NULL)+theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())}
