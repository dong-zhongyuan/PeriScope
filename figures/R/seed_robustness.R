# cp panel_grid/source/round2_panels.R, then retain/edit only the original robustness block.
# Load the cp-edited base_main.R before calling; tpl_7_G is the actual original factory.
make_seed_robustness <- function(input, panel_letter='F', font_pt=9) {
short<-function(x){x<-gsub('raw::','',x,fixed=TRUE);x<-gsub('relative::','Rel. ',x,fixed=TRUE);x<-gsub('nonclassical_mono','nMono',x);x<-gsub('classical_mono','cMono',x);x<-gsub('microglia_homeostatic','MG-H',x);x<-gsub('microglia_mhc2','MG-AP',x);x<-gsub('astro','Astro',x);x<-gsub('__x__',' > ',x);gsub('__',' | ',x)}
# 4E: leave-one-extension-seed-out significance retention, a sensitivity question.
a<-read.csv(input,check.names=FALSE,stringsAsFactors=FALSE);a$state<-ifelse(a$q_within_omission<=.05,'q <= 0.05','q > 0.05');a$state<-factor(a$state,c('q <= 0.05','q > 0.05'));a$count<-1;d<-aggregate(count~combo+state,a,sum);grid<-expand.grid(combo=unique(a$combo),state=levels(a$state));d<-merge(grid,d,all.x=TRUE);d$count[is.na(d$count)]<-0;d$percent<-100*d$count/5;d$state<-factor(d$state,levels(a$state));rate<-aggregate(count~combo,d[d$state=='q <= 0.05',],sum);ord<-sort(unique(d$combo));d$condition<-factor(short(d$combo),rev(short(ord)))
p<-tpl_7_G(d)+coord_flip()+scale_fill_manual(values=c('q <= 0.05'='#42988a','q > 0.05'='#e6e9e9'),name='Omission result')+scale_y_continuous(limits=c(0,100),oob=scales::squish,breaks=seq(0,100,25),expand=expansion(mult=c(0,.02)))+labs(x=NULL,y='Leave-one-seed-out runs (%)')+theme(axis.text.y=element_text(size=font_pt),legend.position='bottom',legend.text=element_text(size=font_pt))
p$layers[[2]]$mapping$label<-aes(label=ifelse(count>0,paste0(count,'/5'),''))$label;p$layers[[2]]$aes_params$size<-font_pt/ggplot2::.pt
p+labs(tag=panel_letter)+theme_classic(base_size=font_pt,base_family='Helvetica')+theme(legend.position='bottom',axis.title=element_text(size=font_pt+1,face='bold'),plot.tag=element_text(size=font_pt+4,face='bold'),plot.title=element_blank(),plot.subtitle=element_blank(),plot.caption=element_blank())

}
