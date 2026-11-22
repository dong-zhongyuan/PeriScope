# Figure 4d - effect profiles across imaging categories for the top-ranked assays (atlas display rows): best signed MR z
# per assay x imaging category, diverging, saturating at the 98th percentile of |z|; empty = not estimated. Data: B_display_rows.csv, B_profiles.csv, categories.csv.
PANEL<-list(id='4d_category_heatmap',blocks=c('2x2'),what='Assay x imaging-category effect profiles')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 P<-read_src(file.path(OV,'B_profiles.csv'),cols=c('gene','imaging_family','is_qc','best_z','display_row','policy'));P<-P[P$policy=='drop'&!P$is_qc,]
 ca<-read_src(file.path(OV,'categories.csv'),cols=c('family','label','is_qc','family_order'));ca<-ca[!ca$is_qc,];ca<-ca[order(ca$family_order),];ab<-c('Reg. volume'='Volume','Cort. area'='Cortical area','Reg. intensity'='Intensity','WM diff.'='WM diffusivity','Rest ampl.'='Rest amplitude','Rest connect.'='Connectivity');ca$label<-ifelse(ca$label%in%names(ab),ab[ca$label],ca$label);ca$label<-gsub('.','',ca$label,fixed=TRUE)
 P<-P[P$imaging_family%in%ca$family,];P$x<-match(P$imaging_family,ca$family)
 rows<-unique(P[,c('gene','display_row')]);rows<-rows[order(rows$display_row),];P$y<-match(P$display_row,rev(rows$display_row))
 lim<-stats::quantile(abs(P$best_z),.98,na.rm=TRUE)
 fr<-data.frame(y=match(rows$gene[rows$gene%in%names(FT_COL)],rev(rows$gene)),g=rows$gene[rows$gene%in%names(FT_COL)])
 p<-ggplot(P,aes(x,y,fill=best_z))+geom_tile(width=.92,height=.9)+
  geom_polygon(data=do.call(rbind,lapply(seq_len(nrow(fr)),function(i)data.frame(g=fr$g[i],x=c(.45,nrow(ca)+.55,nrow(ca)+.55,.45),y=fr$y[i]+c(-.5,-.5,.5,.5)))),aes(x,y,group=g),inherit.aes=FALSE,fill=NA,colour=unname(FT_COL[fr$g[1]]),linewidth=gglw('emph'))+
  scale_fill_gradientn(colours=c('#3F6E9A','white',PERISCOPE),limits=c(-lim,lim),oob=scales::squish,guide='none',na.value=MISSING)+
  scale_x_continuous(breaks=seq_len(nrow(ca)),labels=ca$label,position='top',expand=expansion(add=.5))+
  scale_y_continuous(breaks=seq_len(nrow(rows)),labels=rev(rows$gene),expand=expansion(add=.5))+
  labs(x=NULL,y=NULL)+theme_fig()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x.top=element_text(angle=60,hjust=0,vjust=0),axis.text.y=element_text(face='italic'),
   plot.margin=margin(LETTER_ZONE[['h']],1.5,5,.6,'mm'))
 W<-b$w;H<-b$h;y1<-H-2.2;x<-20
 k<-list(key_title('MR z',x,y1,H));k<-c(k,key_row(list(ki(kg_grad(c('#3F6E9A','white',PERISCOPE),paste0('\u2264 ',fnum(-lim,1)),paste0('\u2265 ',fnum(lim,1)),w=16))),x+sw('MR z',TRUE)+1.2,y1,h=H))
 gTree(children=gList(gg_grob(p),do.call(gList,k)))}
