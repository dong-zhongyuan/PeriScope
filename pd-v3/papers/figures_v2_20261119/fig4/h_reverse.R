# Figure 4h - reverse MR (imaging -> protein), sensitivity to protein-cis-region instruments: the share of reverse
# estimates with P < 0.05 among those carrying no cis instrument, among those carrying one, and after removing it,
# for discovery and replication imaging effects under fixed- and random-effects models. At 1x2 a coverage tile of
# region status by scope is drawn beside it. Data: all_target_delivery_20261119/reverse_cis_summary.csv and
# reverse_status_table.csv (read-only).
PANEL<-list(id='4h_reverse_mr',blocks=c('1x1','1x2'),what='Reverse MR cis sensitivity')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
cis_rows<-function(){
 c_<-read_src(file.path(MR,'all_target_delivery_20261119/reverse_cis_summary.csv'));c_<-c_[c_$scope=='all'&c_$palindrome_policy=='drop',]
 R<-do.call(rbind,lapply(seq_len(nrow(c_)),function(i)do.call(rbind,lapply(c('fixed','random'),function(m)data.frame(set=c_$imaging_effect_set[i],model=m,
  without=c_[[paste0('frac_p_lt_0p05_without_cis_',m)]][i],with=c_[[paste0('frac_p_lt_0p05_with_cis_',m)]][i],after=c_[[paste0('frac_p_lt_0p05_after_removal_',m)]][i])))))
 R$lab<-paste(tools::toTitleCase(R$set),R$model);R$y<-rev(seq_len(nrow(R)));R}
dumbbell<-function(R,lab_inside=FALSE,xmax=32){
 p<-ggplot(R)+geom_vline(xintercept=5,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_segment(aes(x=without*100,xend=with*100,y=y,yend=y),colour=FAINT,linewidth=gglw('emph')*1.5)+
  geom_point(aes(without*100,y),shape=21,fill=MID,colour='white',stroke=.3,size=ggsz('point'))+
  geom_point(aes(with*100,y),shape=21,fill=PERISCOPE,colour='white',stroke=.3,size=ggsz('point'))+
  geom_point(aes(after*100,y),shape=21,fill='white',colour=PERISCOPE,stroke=.4,size=ggsz('point'))+
  scale_x_continuous(limits=c(0,xmax),breaks=c(0,10,20,30))+labs(x='P < 0.05 (%)',y=NULL)+theme_fig()+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
 if(lab_inside)p+geom_text(aes(xmax,y+.34,label=lab),hjust=1,vjust=0,size=GT,colour=INK,family=FONT)+
   scale_y_continuous(breaks=NULL,limits=c(.6,nrow(R)+.9),expand=expansion(0))+theme(axis.text.y=element_blank())
 else p+scale_y_continuous(breaks=R$y,labels=R$lab,expand=expansion(add=.6))}
key_items<-function()list(ki(kg_dot(MID,'white',1.5),'No cis'),ki(kg_dot(PERISCOPE,'white',1.5),'Cis'),
 ki(list(w=1.5,f=function(x,yc)circleGrob(unit(x+.75,'mm'),unit(yc,'mm'),r=unit(.7,'mm'),gp=gpar(fill='white',col=PERISCOPE,lwd=glwd('data')))),'Removed'))
make_panel<-function(b){
 W<-b$w;H<-b$h;R<-cis_rows()
 if(b$name=='1x1'){
  p<-dumbbell(R,lab_inside=TRUE)+theme(plot.margin=margin(LETTER_ZONE[['h']]+2.2,1.2,.6,.6,'mm'))
  k<-key_row(key_items(),4.6,LETTER_ZONE[['h']]-1.4,gap=1.5,h=H)
  return(gTree(children=gList(gg_grob(p),do.call(gList,k))))}
 t<-read_src(file.path(MR,'all_target_delivery_20261119/reverse_status_table.csv'));t<-t[t$palindrome_policy=='drop'&t$imaging_effect_set=='discovery',]
 sc<-c(all='All',model_linked='Model-linked',chrX_protein_region='chrX');st<-c(estimated='Estimated',not_estimated_no_usable_instrument='No instrument',not_estimated_no_published_peak='No peak')
 t<-t[t$scope%in%names(sc)&t$status%in%names(st),];t$f<-ave(t$n_regions,t$scope,FUN=function(v)v/sum(v))
 t$x<-match(t$status,names(st));t$y<-match(t$scope,rev(names(sc)))
 tl<-ggplot(t,aes(x,y,fill=f))+geom_tile(width=.94,height=.9)+
  scale_fill_gradientn(colours=seq_col(seq(0,1,length.out=60),MID),limits=c(0,1),guide='none')+
  scale_x_continuous(breaks=seq_along(st),labels=st,position='top',expand=expansion(0))+scale_y_continuous(breaks=seq_along(sc),labels=rev(sc),expand=expansion(0))+
  labs(x=NULL,y=NULL)+theme_fig()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x.top=element_text(angle=35,hjust=0,vjust=0),plot.margin=margin(LETTER_ZONE[['h']],1,6,.6,'mm'))
 dm<-dumbbell(R)+theme(plot.margin=margin(LETTER_ZONE[['h']],.6,6,.6,'mm'))
 w1<-40
 k<-key_row(key_items(),w1+3,H-1.8,h=H)
 k<-c(k,list(key_title('Share',1,H-1.8,H)),key_row(list(ki(kg_grad(seq_col(seq(0,1,length.out=60),MID),'0','1',w=10))),1+sw('Share',TRUE)+1.2,H-1.8,h=H))
 gTree(children=gList(grobTree(gg_grob(tl),vp=viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(w1,'mm'),height=unit(H,'mm'),just=c(0,0))),
  grobTree(gg_grob(dm),vp=viewport(x=unit(w1+1,'mm'),y=unit(0,'mm'),width=unit(W-w1-1,'mm'),height=unit(H,'mm'),just=c(0,0))),do.call(gList,k)))}
