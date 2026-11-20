# Figure 1h - protein prediction error. Left: overall and donor x subtype-balanced MSE for the four predictors (slope
# lines; PeriScope coloured; end labels spread to one line apart). Right: PeriScope MSE per antigen for every held-out
# donor x subtype group (median point, interquartile rail), grouped by subtype. Data: protein_bridge.csv,
# protein_accuracy.csv (read-only).
PANEL<-list(id='1h_protein_mse',blocks=c('1x2'),what='Protein MSE by predictor and by donor x subtype')
FRD<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design/manuscript_review/figure_redesign_20261004/data'
spread<-function(y,gap){o<-order(y);z<-y[o];for(i in seq_along(z)[-1])z[i]<-max(z[i],z[i-1]+gap);z<-z-(mean(z)-mean(y[o]));out<-y;out[o]<-z;out}
make_panel<-function(b){
 m<-read_src(file.path(FRD,'protein_bridge.csv'),cols=c('method','overall_mse','donor_subtype_macro_mse'))
 key<-c(ridge='ridge',coarse_state_mean='state_mean',subtype_mean='subtype_mean',periscope_42='periscope');m$key<-key[m$method]
 L<-rbind(data.frame(key=m$key,x=1,mse=m$overall_mse),data.frame(key=m$key,x=2,mse=m$donor_subtype_macro_mse))
 L$lab<-METHOD_LAB[L$key]
 yr<-range(L$mse);W<-b$w;H<-b$h;w1<-37
 ph<-H-LETTER_ZONE[['h']]-.6-LINE_MM-1.6-2*LINE_MM-2   # approx. panel height (mm) above the axis and the 2-row legend
 gap<-diff(yr)*1.08/(1-.12)*LINE_MM*1.05/ph
 E<-L[L$x==2,];E$ly<-spread(E$mse,gap)
 p1<-ggplot(L,aes(x,mse,group=lab))+
  geom_line(aes(colour=lab),linewidth=gglw('data'))+geom_point(aes(shape=lab,fill=lab,colour=lab),size=ggsz('point'),stroke=gglw('thin')*1.5)+
  scale_colour_manual(values=METHOD_COL)+scale_fill_manual(values=METHOD_FILL)+scale_shape_manual(values=METHOD_SHAPE)+
  scale_x_continuous(breaks=1:2,labels=c('Overall','Balanced'),limits=c(.8,2.2),expand=expansion(0))+
  scale_y_continuous(labels=num_minus_acc(.01),expand=expansion(mult=.06))+coord_cartesian(clip='off')+
  labs(x=NULL,y='Protein MSE')+theme_fig(legend='bottom')+theme(legend.key.spacing.x=unit(.8,'mm'),legend.box.margin=margin(0,0,0,-6,'mm'))+
  guides(colour=guide_legend(nrow=2),shape=guide_legend(nrow=2),fill=guide_legend(nrow=2))
 a<-read_src(file.path(FRD,'protein_accuracy.csv'),cols=c('target','citeseq_subtype','donor','n','is_control','mse'))
 a<-a[!a$is_control&is.finite(a$mse),]
 st<-c(`CD14 Mono`='cMono',`CD16 Mono`='nMono',cDC1='cDC',cDC2='cDC',pDC='pDC');ord<-c('CD14 Mono','CD16 Mono','cDC1','cDC2','pDC')
 slab<-c(`CD14 Mono`='CD14 Mono',`CD16 Mono`='CD16 Mono',cDC1='cDC1',cDC2='cDC2',pDC='pDC')
 s<-do.call(rbind,lapply(split(a,list(a$citeseq_subtype,a$donor),drop=TRUE),function(z)data.frame(row.names=NULL,sub=z$citeseq_subtype[1],donor=z$donor[1],
  med=stats::median(z$mse),lo=unname(stats::quantile(z$mse,.25)),hi=unname(stats::quantile(z$mse,.75)))))
 s$state<-st[s$sub];s$gx<-match(s$sub,ord);s$x<-6-s$gx-(match(s$donor,c('P6','P7','P8'))-2)*.27
 p2<-ggplot(s,aes(med,x))+
  geom_linerange(aes(xmin=lo,xmax=hi,colour=state),linewidth=gglw('emph')*1.6)+
  geom_point(aes(fill=state),shape=21,colour='white',stroke=.3,size=ggsz('big'))+
  scale_colour_manual(values=BLOOD_COL)+scale_fill_manual(values=BLOOD_COL)+
  scale_y_continuous(breaks=6-seq_along(ord),labels=slab[ord],expand=expansion(add=.45))+scale_x_continuous(labels=num_minus_acc(.1),expand=expansion(mult=.04))+
  labs(y=NULL,x='MSE per antigen')+theme_fig()+
  theme(axis.ticks.y=element_blank(),axis.line.y=element_blank())
 g1<-gg_grob(p1);g2<-gg_grob(p2)
 gTree(children=gList(grobTree(g1,vp=viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(w1,'mm'),height=unit(H,'mm'),just=c(0,0))),
  grobTree(g2,vp=viewport(x=unit(w1+1,'mm'),y=unit(0,'mm'),width=unit(W-w1-1,'mm'),height=unit(H,'mm'),just=c(0,0)))))}
