rings_grid<-function(){
 d<-read.csv(file.path(ROOT,'pd_product_figure_design/template_semantic_review/data/input_1A.csv'));d$ring<-factor(d$ring,c('tissue','condition','cohort'));d$category<-factor(d$category)
 cols<-c(setNames(c('#b46d44','#dfac88','#f4d6bd'),as.character(d$category[d$ring=='tissue'])),setNames(c('#326b62','#498879','#64a08b','#7eb39b','#96c4ad','#b0d2c0','#c7e1d3','#e0eee7'),as.character(d$category[d$ring=='condition'])),setNames(c('#355f8d','#4c7ba3','#6998b8','#89b3c9','#afd0dd','#d3e6ed'),as.character(d$category[d$ring=='cohort'])))
 
 q<-tpl_1_A(d,rings=d)+aes(x=radius)+scale_x_continuous(limits=c(0,4.5),expand=c(0,0))+scale_y_continuous(limits=c(0,1),expand=c(0,0),oob=scales::squish)+scale_fill_manual(values=cols)+theme_void()+theme(legend.position='none',plot.margin=margin(0,0,0,0))+labs(subtitle=NULL,title=NULL,caption=NULL)
 # Ring geometry is reflowed within a square device; leaders terminate at the relevant ring.
 out<-ggdraw()+draw_plot(q,-.005,.135,.57,.713)+draw_label(format(sum(d$count[d$ring=='cohort']),big.mark=','),.280,.514,size=13,fontface='bold')+draw_label('observations',.280,.472,size=12)
 keys<-c('cohort','condition','tissue');titles<-c('Cohort (outer)' ,'Diagnosis (middle)' ,'Observation (inner)' )
 heads<-c(.925,.615,.225);starts<-c(.878,.568,.178)
 anchors<-list(c(.424,.71),c(.446,.496),c(.362,.385))
 for(k in 1:3){
  # The three leaders indicate ring semantics, not individual segments.
  a<-anchors[[k]];elbow<-c(.53,heads[k]);out<-out+draw_line(c(a[1],elbow[1],.55),c(a[2],elbow[2],heads[k]),color='#555555',size=.5)+draw_grob(pointsGrob(x=.5,y=.5,pch=21,gp=gpar(fill='white',col='#555555',fontsize=3)),x=a[1]-.004,y=a[2]-.004,width=.008,height=.008)
  z<-d[d$ring==keys[k],];x<-.57;out<-out+draw_label(titles[k],x,heads[k],hjust=0,size=12,fontface='bold');lab<-sub('^[^:]+: ','',as.character(z$category));lab<-sub('^gse','GSE',lab);lab[lab=='Blood']<-'Blood cells'
  for(j in seq_len(nrow(z))){y<-starts[k]-(j-1)*.035;out<-out+draw_grob(rectGrob(gp=gpar(fill=cols[as.character(z$category[j])],col=NA)),x=x,y=y-.006,width=.012,height=.012)+draw_label(lab[j],x+.023,y,hjust=0,size=12)+draw_label(sprintf('%.1f%%',100*z$weight[j]),.981,y,hjust=1,size=12)}
 }
 out
}
