# Figure 2k - mouse microglial LPS response in two held-out studies (Boles, Shemer). Left: measured vs predicted group
# mean response for 3,381 genes (points = genes). Right: per-mouse MSE of the zero-response baseline and of PeriScope,
# as % of the study's mean baseline MSE; thin lines = mice, thick = study means; number = mean reduction.
# Data: mouse_followup/reference_animal_errors.csv and reference_heldout_predictions.npz (group means exported to
# figures_v2_20261119/fig2/data/mouse_gene_group_means.csv).
PANEL<-list(id='2k_mouse_lps',blocks=c('1x2'),what='Mouse microglial LPS response prediction')
M0<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/diag/periscope_panel_correlation_20261110/payload/pd_product_figure_design/manuscript_review/periscope_panel_correlation_20261110/mouse_followup'
make_panel<-function(b){
 g<-read_src(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/data/mouse_gene_group_means.csv'),cols=c('study','gene','predicted','observed'))
 e<-read_src(file.path(M0,'reference_animal_errors.csv'),cols=c('test_study','animal','zero_mse','mse'))
 r<-do.call(rbind,lapply(split(g,g$study),function(z)data.frame(study=z$study[1],r=stats::cor(z$predicted,z$observed))))
 sc<-ggplot(g,aes(predicted,observed))+
  geom_abline(slope=1,intercept=0,colour=GUIDE,linewidth=gglw('thin'),linetype=REF_LTY)+
  geom_point(shape=16,size=ggsz('unit'),colour=alpha_hex(PERISCOPE,.3))+
  facet_wrap(~study,nrow=1)+scale_x_continuous(labels=num_minus)+scale_y_continuous(labels=num_minus)+
  labs(x='Predicted',y='Measured')+theme_fig()+theme(strip.text=element_text(size=PT,face='bold',margin=margin(0,0,.5,0,'mm')),panel.spacing=unit(2,'mm'))
 e$base<-ave(e$zero_mse,e$test_study,FUN=mean)
 L<-rbind(data.frame(study=e$test_study,animal=e$animal,x=1,v=e$zero_mse/e$base*100),data.frame(study=e$test_study,animal=e$animal,x=2,v=e$mse/e$base*100))
 L$x<-L$x+ifelse(L$study=='Shemer',2.4,0);L$id<-paste(L$study,L$animal)
 m<-aggregate(v~study+x,L,mean)
 red<-do.call(rbind,lapply(split(e,e$test_study),function(z){s<-z$test_study[1];x2<-if(s=='Shemer')4.4 else 2
  data.frame(study=s,x=x2,y=m$v[m$study==s&m$x==x2],pct=(1-mean(z$mse)/mean(z$zero_mse))*100)}))
 red$lab<-paste0('−',formatC(red$pct,format='f',digits=2),'%')
 pl<-ggplot()+
  geom_line(data=L,aes(x,v,group=id),colour=LIGHT,linewidth=gglw('thin'))+
  geom_point(data=L,aes(x,v),shape=16,colour=LIGHT,size=ggsz('small'))+
  geom_line(data=m,aes(x,v,group=study),colour=PERISCOPE,linewidth=gglw('mean'))+
  geom_point(data=m,aes(x,v),shape=21,fill=PERISCOPE,colour='white',stroke=.3,size=ggsz('point')*1.2)+
  annotate('text',x=c(1.5,3.9),y=130,label=c('Boles','Shemer'),family=FONT,size=GT,fontface='bold')+
  scale_x_continuous(breaks=c(1,2,3.4,4.4),labels=c('Zero','PeriScope','Zero','PeriScope'),expand=expansion(add=.45))+
  scale_y_continuous(limits=c(0,136),breaks=c(0,50,100),expand=expansion(0))+
  labs(x=NULL,y='MSE (% of zero)')+theme_fig()+theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))
 W<-b$w;H<-b$h;w1<-48
 gTree(children=gList(grobTree(gg_grob(sc),vp=viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(w1,'mm'),height=unit(H,'mm'),just=c(0,0))),
  grobTree(gg_grob(pl),vp=viewport(x=unit(w1+1,'mm'),y=unit(0,'mm'),width=unit(W-w1-1,'mm'),height=unit(H,'mm'),just=c(0,0)))))}
