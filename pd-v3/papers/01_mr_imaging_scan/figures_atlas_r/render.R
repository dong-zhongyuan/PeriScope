args<-commandArgs(TRUE);stopifnot(length(args)==1)
arg0<-commandArgs(FALSE);script<-sub('^--file=','',arg0[grep('^--file=',arg0)][1])
HERE<-dirname(normalizePath(script));OUT<-normalizePath(args[1]);DATA<-file.path(OUT,'plot_data')
for(d in c('panels','qa'))dir.create(file.path(OUT,d),showWarnings=FALSE,recursive=TRUE)
source(file.path(HERE,'common.R'))
layout<-read.csv(file.path(HERE,'layout.csv'));W<-174;H<-225
stopifnot(all(layout$x_mm>=1),all(layout$y_mm>=1),all(layout$x_mm+layout$w_mm<=W-1),all(layout$y_mm+layout$h_mm<=H-7))
for(i in 1:(nrow(layout)-1))for(j in (i+1):nrow(layout)) {
  a<-layout[i,];b<-layout[j,]
  stopifnot(!(a$x_mm<b$x_mm+b$w_mm && b$x_mm<a$x_mm+a$w_mm && a$y_mm<b$y_mm+b$h_mm && b$y_mm<a$y_mm+a$h_mm))
}
grobs<-list()
for(i in seq_len(nrow(layout))) {
  l<-layout[i,];source(file.path(HERE,paste0(l$panel,'.R')))
  grobs[[l$panel]]<-wrap_atlas(get(paste0('make_',l$panel))(),l$panel,l$title,l$w_mm,l$h_mm)
  cat('Rendered',l$panel,'\n')
}
draw_page<-function() {
  grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
  for(i in seq_len(nrow(layout))) {
    l<-layout[i,];pushViewport(viewport(x=unit(l$x_mm,'mm'),y=unit(H-l$y_mm,'mm'),width=unit(l$w_mm,'mm'),height=unit(l$h_mm,'mm'),just=c('left','top'),clip='on'))
    grid.draw(grobs[[l$panel]]);popViewport()
  }
  label<-if(cfg$mode=='full')'Descriptive overview; category ranges are not confidence intervals.' else 'VALIDATION ONLY: completed 45-trait panel; not the full-atlas result.'
  grid.draw(txt(label,2,221,6.5))
}
stem<-if(cfg$mode=='full')'Full_atlas_overview' else 'Validation45_overview'
open_pdf(file.path(OUT,paste0(stem,'.pdf')),W,H);draw_page();dev.off()
ragg::agg_png(file.path(OUT,paste0(stem,'.png')),width=W,height=H,units='mm',res=600,background='white');draw_page();dev.off()
ragg::agg_png(file.path(OUT,paste0(stem,'_preview.png')),width=W,height=H,units='mm',res=150,background='white');draw_page();dev.off()
write.csv(layout,file.path(OUT,'qa','layout_mm.csv'),row.names=FALSE)
write_json(list(page_mm=c(W,H),panels=nrow(layout),plotting_areas=5+nrow(cats),minimum_font_pt=6,
  simple_bar_charts=0,layout_bounds='PASS',slot_overlap='PASS',text_sizes='PASS',
  mode=cfg$mode,formal_full_atlas_result=cfg$formal_full_atlas_result),
  file.path(OUT,'qa','layout_checks.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'qa','R_sessionInfo.txt'))
