make_B <- function() {
  d<-read_table('B_profiles');rows<-read_table('B_display_rows')
  families<-cats[!cats$is_qc,];L<-cfg$matrix_z_limit
  stopifnot(nrow(d)==nrow(rows)*nrow(families)*2,!any(d$is_qc))
  d$column<-match(d$imaging_family,families$family)
  d$policy_label<-policy_factor(d$policy)
  labels<-paste0(rows$gene,ifelse(duplicated(rows$gene)|duplicated(rows$gene,fromLast=TRUE),
                                 paste0(' [',rows$assay,']'),''))
  valid<-d[d$n_estimated>0,];missing<-d[d$n_estimated==0,]
  p<-ggplot(d,aes(column,display_row,fill=best_z))+
    geom_tile(width=.96,height=.94)+
    geom_segment(data=valid,aes(x=column+z_min/L*.4,xend=column+z_max/L*.4,
      y=display_row+.30,yend=display_row+.30),inherit.aes=FALSE,linewidth=.35/.pt,colour='#252525')+
    geom_segment(data=valid,aes(x=column,xend=column,y=display_row+.23,yend=display_row+.37),
      inherit.aes=FALSE,linewidth=.25/.pt,colour='#252525')+
    geom_point(data=missing,aes(column,display_row),inherit.aes=FALSE,shape=4,size=.9,stroke=.2)+
    facet_wrap(~policy_label,nrow=1)+
    scale_y_reverse(breaks=rows$display_row,labels=labels,expand=c(0,0))+
    scale_x_continuous(breaks=seq_len(nrow(families)),
      labels=paste0(families$label,' (',families$n_phenotypes,')'),expand=c(0,0))+
    scale_fill_gradient2(low='#0072B2',mid='white',high='#D55E00',midpoint=0,
      limits=c(-L,L),na.value='#E6E6E6',breaks=c(-L,0,L),name='Best-pair z',
      guide=guide_colourbar(display='rectangles',barwidth=unit(1.5,'mm'),barheight=unit(18,'mm')))+
    labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),
      axis.text.x=element_text(angle=if(nrow(families)>4)90 else 0,hjust=if(nrow(families)>4)1 else .5,
        vjust=if(nrow(families)>4).5 else 1,size=6),
      axis.text.y=element_text(size=6),panel.spacing.x=unit(3,'mm'),legend.position='right')
  if(nrow(families)<=4) p<-p+geom_text(data=valid,
    aes(label=formatC(best_z,format='f',digits=1)),family=FONT,size=6/.pt,colour=INK)
  p
}
