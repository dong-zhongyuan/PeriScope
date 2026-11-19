make_B <- function() {
  d<-read_pd('B_matrix');cols<-read_pd('B_columns')
  stopifnot(nrow(d)==length(unique(d$exposure_id))*45)
  rows<-d[!duplicated(d$assay_row),];rows<-rows[order(rows$assay_row),]
  coding<-tolower(as.character(rows$pav_coding_consequence))=='true'
  rows$label<-paste0(rows$gene,ifelse(coding,'†',''))
  d$pid<-factor(d$phenotype,levels=cols$phenotype)
  ggplot(d,aes(pid,assay_row,fill=z))+geom_tile(width=.96,height=.94)+
    facet_grid(.~modality,scales='free_x',space='free_x')+
    scale_y_reverse(breaks=rows$assay_row,labels=rows$label,expand=c(0,0))+
    scale_x_discrete(labels=setNames(cols$idp_label,cols$phenotype),expand=c(0,0))+
    scale_fill_gradient2(low='#0072B2',mid='white',high='#D55E00',midpoint=0,limits=c(-9,9),
                         breaks=c(-8,0,8),labels=num,name='Wald z',
                         guide=guide_colourbar(display='rectangles',barwidth=unit(1.5,'mm'),barheight=unit(18,'mm')))+
    labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),
      axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=6,margin=margin(t=.4,unit='mm')),
      axis.text.y=element_text(size=6),panel.spacing.x=unit(1.2,'mm'),
      strip.text=element_text(size=6.5,face='bold'),legend.position='right')
}
