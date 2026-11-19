# Figure 2i - random pairing (seeds 52-56) relative to the full-model mean (seeds 42-51) on six readouts: points = seeds,
# bar = mean, grey band = the full model's own seed range, number = mean change. Higher is better for the two Spearman
# readouts, lower for the four error readouts. Data: arm_vs_full_pool.csv (arm_set random_pairing_new).
PANEL<-list(id='2i_metric_ablation',blocks=c('1x1'),what='Readout changes under random pairing')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/_arms.R'))
make_panel<-function(b){
 a<-read_src(file.path(RES,'random_pairing_20261119/evaluation/tables/arm_vs_full_pool.csv'))
 mets<-c(f_all_mean='Route ρ',h_rho_all='Axis ρ',h_contrast_mse='Contrast MSE',h_mean_mse='Mean MSE',h_mmd='MMD',h_energy='Energy')
 a<-a[a$arm_set=='random_pairing_new'&a$metric%in%names(mets),];a<-a[match(names(mets),a$metric),];stopifnot(nrow(a)==6)
 a$y<-rev(seq_len(nrow(a)));a$lab<-mets[a$metric]
 pts<-do.call(rbind,lapply(seq_len(nrow(a)),function(i){v<-as.numeric(strsplit(a$arm_values[i],';')[[1]]);data.frame(y=a$y[i],p=(v/a$full_mean[i]-1)*100)}))
 a$lo<-(a$full_min/a$full_mean-1)*100;a$hi<-(a$full_max/a$full_mean-1)*100
 a$worse<-ifelse(a$better=='high',a$pct_diff_vs_full<0,a$pct_diff_vs_full>0)
 a$txt<-paste0(ifelse(a$pct_diff_vs_full>0,'+','−'),formatC(abs(a$pct_diff_vs_full),format='f',digits=1),'%')
 P<-ggplot()+
  geom_rect(data=a,aes(xmin=lo,xmax=hi,ymin=y-.32,ymax=y+.32),fill=PERISCOPE_FAINT,colour=NA)+
  geom_vline(xintercept=0,colour=PERISCOPE,linewidth=gglw('ref'))+
  geom_point(data=pts,aes(p,y),shape=24,fill=ARM_COL[['Random pairing']],colour='white',stroke=.2,size=ggsz('point'))+
  geom_segment(data=a,aes(x=pct_diff_vs_full,xend=pct_diff_vs_full,y=y-.36,yend=y+.36),colour=INK,linewidth=gglw('mean'))+
  scale_y_continuous(breaks=a$y,labels=a$lab,expand=expansion(add=.55))+
  scale_x_continuous(breaks=c(-30,0,30),labels=num_minus,expand=expansion(0))+
  coord_cartesian(xlim=c(-42,32),clip='off')+labs(x='Change vs full (%)',y=NULL)+theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(LETTER_ZONE[['h']],1.2,.6,.6,'mm'))
 gg_grob(P)}
