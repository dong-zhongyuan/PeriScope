# Figure 2h - program-direction ablation: median cosine between each seed's program responses and the confirmation-seed
# (47-51) reference, for all programs and per recipient state; full model seeds 42-51 vs random pairing seeds 52-56.
# Numbers = change of the random-pairing mean relative to the full-model mean.
# Data: random_pairing_20261119/evaluation/tables/per_seed_metrics.csv and arm_vs_full_pool.csv.
PANEL<-list(id='2h_program_ablation',blocks=c('1x1'),what='Program-direction ablation, random pairing vs full')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/_arms.R'))
make_panel<-function(b){
 s<-read_src(file.path(RES,'random_pairing_20261119/evaluation/tables/per_seed_metrics.csv'))
 a<-read_src(file.path(RES,'random_pairing_20261119/evaluation/tables/arm_vs_full_pool.csv'),cols=c('arm_set','metric','pct_diff_vs_full'))
 s<-s[(s$arm=='full'&s$seed%in%42:51)|(s$arm=='random_pairing'&s$seed%in%52:56),];stopifnot(sum(s$arm=='full')==10,sum(s$arm=='random_pairing')==5)
 mets<-c(p_cos_median='All',p_cos_median_astro='Astro',p_cos_median_microglia_homeostatic='MG-H',p_cos_median_microglia_mhc2='MG-AP')
 L<-do.call(rbind,lapply(names(mets),function(m)data.frame(state=mets[[m]],metric=m,arm=ARM_LAB[s$arm],seed=s$seed,v=s[[m]])))
 L$g<-match(L$state,mets);L$x<-L$g+ifelse(L$arm=='Full',-.2,.2)
 set.seed(3);L$xj<-L$x+stats::runif(nrow(L),-.07,.07)
 m<-aggregate(v~arm+g+x,L,mean)
 pc<-a[a$arm_set=='random_pairing_new'&a$metric%in%names(mets),];pc$g<-match(pc$metric,names(mets))
 pc$lab<-paste0(sub('-','−',formatC(pc$pct_diff_vs_full,format='f',digits=1)),'%')
 P<-ggplot()+
  geom_hline(yintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_point(data=L,aes(xj,v,colour=arm),shape=16,size=ggsz('small'),alpha=.75)+
  geom_segment(data=m,aes(x=x-.15,xend=x+.15,y=v,yend=v,colour=arm),linewidth=gglw('mean')*1.2)+
  scale_colour_manual(values=ARM_COL,guide='none')+
  scale_x_continuous(breaks=seq_along(mets),labels=unname(mets),expand=expansion(add=.55))+
  scale_y_continuous(limits=c(-.2,.65),breaks=c(-.2,0,.2,.4),labels=num_minus_acc(.1),expand=expansion(0))+
  coord_cartesian(clip='off')+labs(x=NULL,y='Program cosine')+theme_fig()+theme(axis.ticks.x=element_blank(),axis.text.x=element_text(angle=40,hjust=1,vjust=1))
 gg_grob(P)}
