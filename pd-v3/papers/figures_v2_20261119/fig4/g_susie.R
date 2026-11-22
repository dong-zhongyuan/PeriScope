# Figure 4g - SuSiE-coloc vs single-signal coloc for every pair where SuSiE ran, on both protein platforms: x =
# single-signal PP.H4, y = SuSiE PP.H4 for the credible-set pair carrying the MR instrument; colour = instrument
# inside the protein credible set, shape = platform. Dashed lines = 0.8.
# Data: all_target_susie_coloc_20261119/susie_coloc_pairs.csv.gz (read-only).
PANEL<-list(id='4g_susie',blocks=c('1x1','1x2'),what='SuSiE-coloc vs single-signal coloc')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 s<-read_src(file.path(MR,'all_target_susie_coloc_20261119/susie_coloc_pairs.csv.gz'),cols=c('platform','status','H4_single','H4_susie_mr_signal','H4_susie_any','instrument_in_signal'))
 s<-s[s$status=='ran'&is.finite(s$H4_single),];s$y<-ifelse(is.finite(s$H4_susie_mr_signal),s$H4_susie_mr_signal,s$H4_susie_any);s<-s[is.finite(s$y),]
 s$ins<-s$instrument_in_signal%in%c(TRUE,'True');s$platform<-factor(ifelse(s$platform=='UKB','UKB-PPP',s$platform),c('UKB-PPP','AGES'))
 base<-function(p)p+geom_abline(slope=1,intercept=0,colour=GUIDE,linewidth=gglw('thin'))+
  geom_hline(yintercept=.8,colour=REF_COL,linewidth=gglw('ref'),linetype=REF_LTY)+geom_vline(xintercept=.8,colour=REF_COL,linewidth=gglw('ref'),linetype=REF_LTY)+
  scale_colour_manual(values=c(`TRUE`=PERISCOPE,`FALSE`=LIGHT),guide='none')+
  scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','0.5','1'))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','0.5','1'))+
  labs(x='Single-signal PP.H4',y='SuSiE PP.H4')+theme_fig()
 W<-b$w;H<-b$h
 if(b$name=='1x1'){
  p<-base(ggplot(s,aes(H4_single,y)))+geom_point(aes(colour=ins,shape=platform),size=ggsz('small'),stroke=.3,alpha=.8)+
   scale_shape_manual(values=c(`UKB-PPP`=16,AGES=18),guide='none')+coord_equal()+
   theme(plot.margin=margin(LETTER_ZONE[['h']]+5.4,.6,.6,.6,'mm'))
  k<-c(key_row(list(ki(kg_dot(PERISCOPE,NA,1.4),'In set'),ki(kg_dot(LIGHT,NA,1.4),'Outside')),10.5,LETTER_ZONE[['h']]-1.3,h=H),
       key_row(list(ki(list(w=1.5,f=function(x,yc)circleGrob(unit(x+.75,'mm'),unit(yc,'mm'),r=unit(.7,'mm'),gp=gpar(fill=MID,col=NA))),'UKB-PPP'),
                    ki(list(w=1.7,f=function(x,yc)ddiam(x+.85,0,1.9,MID,NA,glwd('thin'),h=yc)),'AGES')),10.5,LETTER_ZONE[['h']]+2.1,h=H))
  return(gTree(children=gList(gg_grob(p),do.call(gList,k))))}
 p<-base(ggplot(s,aes(H4_single,y)))+geom_point(aes(colour=ins),shape=16,size=ggsz('small'),alpha=.75)+facet_wrap(~platform,nrow=1)+coord_equal()+
  theme(strip.text=element_text(size=PT,face='bold',margin=margin(0,0,.5,0,'mm')),panel.spacing=unit(3,'mm'),plot.margin=margin(LETTER_ZONE[['h']]-.6,22,.6,.6,'mm'))
 k<-c(key_row(list(ki(kg_dot(PERISCOPE,NA,1.4),'In set')),W-20.5,H/2-2,h=H),key_row(list(ki(kg_dot(LIGHT,NA,1.4),'Outside')),W-20.5,H/2+1.6,h=H))
 gTree(children=gList(gg_grob(p),do.call(gList,k)))}
