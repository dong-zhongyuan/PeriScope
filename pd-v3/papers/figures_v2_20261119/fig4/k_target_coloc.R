# Figure 4k - colocalisation of each final-target MR pair on the two protein platforms: unified PP.H4 (higher of single-
# signal and SuSiE) on AGES (x) and UKB-PPP (y); dashed = 0.8 gate (pass on either platform); colour = target.
# Data: target_funnel_final_20261119/final_target_mr_coloc_pairs.csv.
PANEL<-list(id='4k_target_coloc',blocks=c('1x1'),what='Final-target coloc on AGES vs UKB-PPP')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 f<-read_src(file.path(TF,'final_target_mr_coloc_pairs.csv'),cols=c('UKB_gene','unified_H4_AGES','unified_H4_UKB'))
 f$a<-ifelse(is.finite(f$unified_H4_AGES),f$unified_H4_AGES,0);f$u<-ifelse(is.finite(f$unified_H4_UKB),f$unified_H4_UKB,0)
 ggplot(f,aes(a,u))+
  geom_hline(yintercept=.8,colour=REF_COL,linewidth=gglw('ref'),linetype=REF_LTY)+geom_vline(xintercept=.8,colour=REF_COL,linewidth=gglw('ref'),linetype=REF_LTY)+
  geom_point(shape=21,fill=FT_COL[f$UKB_gene],colour='white',stroke=.3,size=ggsz('big'))+
  scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','0.5','1'))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','0.5','1'))+
  coord_equal()+labs(x='AGES PP.H4',y='UKB-PPP PP.H4')+theme_fig()}
