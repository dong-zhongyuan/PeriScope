args<-commandArgs(TRUE);stopifnot(length(args)==1)
arg0<-commandArgs(FALSE);script<-sub('^--file=','',arg0[grep('^--file=',arg0)][1])
HERE<-dirname(normalizePath(script));OUT<-normalizePath(args[1]);DATA<-file.path(OUT,'plot_data')
for(d in c('panels','qa'))dir.create(file.path(OUT,d),showWarnings=FALSE,recursive=TRUE)
source(file.path(HERE,'common.R'))
layout<-read.csv(file.path(HERE,'layout.csv'));W<-174;H<-225
stopifnot(all(layout$x_mm>=1),all(layout$y_mm>=1),all(layout$x_mm+layout$w_mm<=W-1),all(layout$y_mm+layout$h_mm<=H-1))
for(i in 1:(nrow(layout)-1))for(j in (i+1):nrow(layout)) {
  a<-layout[i,];b<-layout[j,]
  stopifnot(!(a$x_mm<b$x_mm+b$w_mm && b$x_mm<a$x_mm+a$w_mm && a$y_mm<b$y_mm+b$h_mm && b$y_mm<a$y_mm+a$h_mm))
}
grobs<-list()
for(i in seq_len(nrow(layout))) {
  l<-layout[i,];source(file.path(HERE,paste0(l$panel,'.R')))
  grobs[[l$panel]]<-wrap_panel(get(paste0('make_',l$panel))(),l$panel,l$title,l$w_mm,l$h_mm)
  cat('Rendered',l$panel,'\n')
}
draw_page<-function() {
  grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
  for(i in seq_len(nrow(layout))) {
    l<-layout[i,];pushViewport(viewport(x=unit(l$x_mm,'mm'),y=unit(H-l$y_mm,'mm'),width=unit(l$w_mm,'mm'),height=unit(l$h_mm,'mm'),just=c('left','top'),clip='on'))
    grid.draw(grobs[[l$panel]]);popViewport()
  }
}
open_pdf(file.path(OUT,'Figure1_expanded_MR.pdf'),W,H);draw_page();dev.off()
ragg::agg_png(file.path(OUT,'Figure1_expanded_MR.png'),width=W,height=H,units='mm',res=600,background='white');draw_page();dev.off()
ragg::agg_png(file.path(OUT,'Figure1_expanded_MR_preview.png'),width=W,height=H,units='mm',res=150,background='white');draw_page();dev.off()
# Supporting figure: all 135 estimates, fixed declared phenotype order, no ranking.
d<-read_pd('S1_historical');cols<-read_pd('B_columns');stopifnot(nrow(d)==135)
d$gene<-factor(ifelse(d$gene=='TFRC','TFRC†',d$gene),levels=c('CD22','CD40','TFRC†'))
d$label<-paste(as.character(d$modality),d$idp_label)
labs<-d[!duplicated(d$idp_order),];labs<-labs[order(labs$idp_order),]
p<-ggplot(d,aes(beta,idp_order,colour=modality))+zero_v()+
  geom_segment(aes(x=ci_lower,xend=ci_upper,yend=idp_order),linewidth=.4/.pt)+geom_point(size=1.3)+
  facet_wrap(~gene,nrow=1,scales='free_x')+scale_colour_manual(values=MOD,guide='none')+
  scale_y_reverse(breaks=labs$idp_order,labels=labs$label,expand=c(.01,0))+
  numeric_x()+labs(x='Imaging SD / protein analysis unit',y=NULL)+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing.x=unit(4,'mm'),plot.margin=margin(2,2,2,1,'mm'))
ragg::agg_capture(width=174,height=200,units='mm',res=300);sg<-ggplotGrob(p);dev.off()
stopifnot(all(text_audit(sg)$pt>=6-1e-6))
open_pdf(file.path(OUT,'historical_proteins_complete_panel.pdf'),174,200);grid.newpage();grid.draw(sg);dev.off()
ragg::agg_png(file.path(OUT,'historical_proteins_complete_panel.png'),width=174,height=200,units='mm',res=600,background='white');grid.newpage();grid.draw(sg);dev.off()
ragg::agg_png(file.path(OUT,'historical_proteins_complete_panel_preview.png'),width=174,height=200,units='mm',res=150,background='white');grid.newpage();grid.draw(sg);dev.off()
write.csv(layout,file.path(OUT,'qa','layout_mm.csv'),row.names=FALSE)
write_json(list(page_mm=c(W,H),panels=nrow(layout),plotting_areas=sum(layout$areas),minimum_font_pt=6,
                slot_coverage=sum(layout$w_mm*layout$h_mm)/(172*223),simple_bar_charts=0,
                layout_bounds='PASS',slot_overlap='PASS',text_sizes='PASS'),file.path(OUT,'qa','layout_checks.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'qa','R_sessionInfo.txt'))
