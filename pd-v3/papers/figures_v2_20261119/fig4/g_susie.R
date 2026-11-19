# Figure 4g - SuSiE-coloc vs single-signal coloc for every pair where SuSiE ran, per protein platform: x = single-signal
# PP.H4, y = SuSiE PP.H4 for the credible-set pair carrying the MR instrument; colour = instrument inside the protein
# credible set. Dashed lines = 0.8. Data: all_target_susie_coloc_20261119/susie_coloc_pairs.csv.gz (read-only).
PANEL<-list(id='4g_susie',blocks=c('1x2'),what='SuSiE-coloc vs single-signal coloc')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 s<-read_src(file.path(MR,'all_target_susie_coloc_20261119/susie_coloc_pairs.csv.gz'),cols=c('platform','status','H4_single','H4_susie_mr_signal','H4_susie_any','instrument_in_signal'))
 s<-s[s$status=='ran'&is.finite(s$H4_single),];s$y<-ifelse(is.finite(s$H4_susie_mr_signal),s$H4_susie_mr_signal,s$H4_susie_any);s<-s[is.finite(s$y),]
 s$ins<-s$instrument_in_signal%in%c(TRUE,'True');s$platform<-factor(ifelse(s$platform=='UKB','UKB-PPP',s$platform),c('UKB-PPP','AGES'))
 p<-ggplot(s,aes(H4_single,y))+
  geom_abline(slope=1,intercept=0,colour=GUIDE,linewidth=gglw('thin'))+
  geom_hline(yintercept=.8,colour=REF_COL,linewidth=gglw('ref'),linetype=REF_LTY)+geom_vline(xintercept=.8,colour=REF_COL,linewidth=gglw('ref'),linetype=REF_LTY)+
  geom_point(aes(colour=ins),shape=16,size=ggsz('small'),alpha=.75)+
  facet_wrap(~platform,nrow=1)+
  scale_colour_manual(values=c(`TRUE`=PERISCOPE,`FALSE`=LIGHT),guide='none')+
  scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','0.5','1'))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','0.5','1'))+
  coord_equal()+labs(x='Single-signal PP.H4',y='SuSiE PP.H4')+theme_fig()+
  theme(strip.text=element_text(size=PT,face='bold',margin=margin(0,0,.5,0,'mm')),panel.spacing=unit(3,'mm'),plot.margin=margin(LETTER_ZONE[['h']]-.6,22,.6,.6,'mm'))
 W<-b$w;H<-b$h;k<-c(key_row(list(ki(kg_dot(PERISCOPE,NA,1.4),'In set')),W-20.5,H/2-2,h=H),key_row(list(ki(kg_dot(LIGHT,NA,1.4),'Outside')),W-20.5,H/2+1.6,h=H))
 gTree(children=gList(gg_grob(p),do.call(gList,k)))}
