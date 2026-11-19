# Figure 2g - structural ablation on the six cDC/monocyte-to-glia routes (manuscript Fig. 3D design, published
# checkpoints, seeds 42-44): small points = seeds, glyph = mean, rail = SD; last group = six-route mean.
# Data: manuscript_readout_ablation_20261119/routes_by_arm_seed.csv (checkpoint_set == published).
PANEL<-list(id='2g_route_ablation',blocks=c('1x2'),what='Structural ablation by route')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/_arms.R'))
make_panel<-function(b){
 d<-read_src(file.path(RES,'manuscript_readout_ablation_20261119/routes_by_arm_seed.csv'),cols=c('checkpoint_set','arm','seed','route','rho'))
 d<-d[d$checkpoint_set=='published'&d$arm%in%names(ARM_LAB)&d$seed%in%42:44,];stopifnot(nrow(d)==3*3*6)
 rt<-c('cDC_to_Astro','cMono_to_Astro','nMono_to_Astro','cDC_to_MG-AP','cMono_to_MG-AP','nMono_to_MG-AP')
 mn<-aggregate(rho~arm+seed,d,mean);mn$route<-'mean';d<-rbind(d[,c('arm','seed','route','rho')],mn[,c('arm','seed','route','rho')])
 d$g<-match(d$route,c(rt,'mean'));d$g<-d$g+ifelse(d$g>3,.35,0)+ifelse(d$g>6,.35,0)
 d$arm<-ARM_LAB[d$arm];d$off<-c(Full=-.27,`Random pairing`=0,`Coupling off`=.27)[d$arm];d$x<-d$g+d$off
 s<-aggregate(rho~arm+g+off+x,d,function(v)c(m=mean(v),s=stats::sd(v)));s<-data.frame(s[,1:4],m=s$rho[,'m'],sd=s$rho[,'s'])
 gx<-sort(unique(d$g));lab<-c('cDC','cMono','nMono','cDC','cMono','nMono','Mean')
 br<-data.frame(x0=c(gx[1],gx[4])-.4,x1=c(gx[3],gx[6])+.4,col=BRAIN_COL[c('Astro','MG-AP')],lab=c('Astro','MG-AP'))
 yl<-c(-.3,.4)
 P<-ggplot()+
  geom_hline(yintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_linerange(data=s,aes(x=x,ymin=m-sd,ymax=m+sd,colour=arm),linewidth=gglw('emph'))+
  geom_point(data=d,aes(x+.09,rho,colour=arm),shape=16,size=ggsz('small'),alpha=.7)+
  geom_point(data=s,aes(x,m,shape=arm,fill=arm,colour=arm),size=ggsz('point')*1.15,stroke=.3)+
  geom_segment(data=br,aes(x=x0,xend=x1,y=yl[1]-.075,yend=yl[1]-.075),colour=br$col,linewidth=gglw('emph')*1.4)+
  geom_text(data=br,aes((x0+x1)/2,yl[1]-.1,label=lab),family=FONT,size=GT,vjust=1,colour=INK)+
  scale_colour_manual(values=ARM_COL,guide='none')+scale_fill_manual(values=ARM_COL,guide='none')+scale_shape_manual(values=ARM_SHAPE,guide='none')+
  scale_x_continuous(breaks=gx,labels=lab,expand=expansion(add=.45))+
  scale_y_continuous(breaks=seq(-.2,.4,.2),labels=num_minus_acc(.1),expand=expansion(0))+
  coord_cartesian(ylim=yl,clip='off')+labs(x=NULL,y='Spearman')+theme_fig()+
  theme(axis.ticks.x=element_blank(),plot.margin=margin(LETTER_ZONE[['h']]+.4,.6,4.4,.6,'mm'))
 gg_grob(P)}
