suppressPackageStartupMessages({library(ggplot2); library(grid); library(gtable); library(jsonlite)})
FONT <- 'Arial'
MOD <- c(Volume='#758AA2', QSM='#3A8B84', `T2*`='#C48A53')
INK <- '#252525'; ZERO <- '#BDBDBD'
num <- function(x) sub('^-','−', format(x, trim=TRUE, scientific=FALSE))
BASE <- theme_classic(base_family=FONT, base_size=6) + theme(
  axis.title=element_text(size=7), axis.text=element_text(size=6, colour=INK),
  axis.line=element_line(linewidth=.4/.pt, colour=INK),
  axis.ticks=element_line(linewidth=.3/.pt), axis.ticks.length=unit(1,'mm'),
  axis.title.x=element_text(margin=margin(t=1,unit='mm')),
  axis.title.y=element_text(margin=margin(r=1,unit='mm')),
  strip.background=element_blank(), strip.text=element_text(size=6.5,face='bold'),
  panel.spacing=unit(2,'mm'), plot.margin=margin(5.5,1,1,1,'mm'),
  legend.text=element_text(size=6),legend.title=element_text(size=6.5),
  legend.key.height=unit(2.5,'mm'),legend.key.width=unit(4,'mm'),
  legend.margin=margin(0,0,0,0,'mm'),legend.box.spacing=unit(1,'mm'))
theme_set(BASE)
read_pd <- function(name) {
  d <- read.csv(file.path(DATA,paste0(name,'.csv')),check.names=FALSE)
  if('modality'%in%names(d)) d$modality <- factor(d$modality,levels=names(MOD))
  d
}
REG <- c('THA','CAU','PUT','PAL','HIP','AMY','ACC','SN')
shape_gene <- c(TFRC=21,MRC1=24,LRRC37A2=23)
zero_v <- function() geom_vline(xintercept=0,linewidth=.3/.pt,colour=ZERO,linetype='dashed')
zero_h <- function() geom_hline(yintercept=0,linewidth=.3/.pt,colour=ZERO,linetype='dashed')
numeric_x <- function(...) scale_x_continuous(labels=num,...)
numeric_y <- function(...) scale_y_continuous(labels=num,...)
open_pdf <- function(file,w,h) {
  if(capabilities('aqua')) grDevices::quartz(file=file,type='pdf',width=w/25.4,height=h/25.4,family=FONT,bg='white')
  else grDevices::cairo_pdf(file,width=w/25.4,height=h/25.4,family=FONT,bg='white')
}
txt <- function(label,x,y,size=6,face='plain',just=c('left','top')) {
  stopifnot(size>=6)
  textGrob(label,x=unit(x,'mm'),y=unit(1,'npc')-unit(y,'mm'),just=just,
           gp=gpar(fontfamily=FONT,fontsize=size,fontface=face,col=INK))
}
text_audit <- function(g) {
  found<-list()
  walk<-function(x,fs=12,cex=1) {
    if(!is.null(x$gp$fontsize)) fs<-x$gp$fontsize
    if(!is.null(x$gp$cex)) cex<-cex*x$gp$cex
    if(inherits(x,'text')) found[[length(found)+1]]<<-data.frame(label=as.character(x$label),pt=rep_len(fs*cex,length(x$label)))
    if(inherits(x,'gtable')) for(ch in x$grobs) walk(ch,fs,cex)
    else if(!is.null(x$children)) for(ch in x$children) walk(ch,fs,cex)
  }
  walk(g); if(length(found)) do.call(rbind,found) else data.frame(label=character(),pt=numeric())
}
wrap_panel <- function(p,id,title,w,h) {
  stopifnot(is.null(p$labels$title),is.null(p$labels$subtitle),is.null(p$labels$caption))
  ragg::agg_capture(width=w,height=h,units='mm',res=300)
  g<-ggplotGrob(p); dev.off()
  keys<-list()
  mark<-function(x,y,pch,fill='white') pointsGrob(x=unit(x,'mm'),y=unit(1,'npc')-unit(y,'mm'),
      pch=pch,size=unit(1.1,'mm'),gp=gpar(col=INK,fill=fill,lwd=.3))
  box<-function(x,y,ww,hh,fill,col=NA) rectGrob(x=unit(x,'mm'),y=unit(1,'npc')-unit(y,'mm'),
      width=unit(ww,'mm'),height=unit(hh,'mm'),just=c('left','top'),gp=gpar(fill=fill,col=col,lwd=.3))
  if(id=='D')keys<-list(mark(43,2.5,21,INK),txt('L',44.5,1.5),mark(51,2.5,21),txt('R',52.5,1.5))
  if(id=='E')keys<-list(mark(42,2.5,21),txt('TFRC',43.5,1.5),mark(54,2.5,24),txt('MRC1',55.5,1.5),
                        mark(66,2.5,23),txt('LRRC37A2',67.5,1.5))
  if(id=='H')keys<-list(box(43,1,3,1.5,'#B59DCB'),box(43,2.8,3,1.5,'#B59DCB'),
                        txt('Protein',47,.5),txt('Imaging',47,2.8),
                        box(45,h-7.8,2,2,'white',INK),txt('Unresolved 42',48,h-7.8))
  result<-gTree(children=gList(rectGrob(gp=gpar(fill='white',col=NA)),g,
                              txt(tolower(id),.6,.8,8,'bold'),txt(title,4.6,1.2,6.5,'bold'),
                              gTree(children=do.call(gList,keys))))
  audit<-text_audit(result);stopifnot(all(audit$pt>=6-1e-6))
  write.csv(audit,file.path(OUT,'qa',paste0('text_',id,'.csv')),row.names=FALSE)
  for(suffix in c('pdf','png')) {
    f<-file.path(OUT,'panels',paste0(id,'.',suffix))
    if(suffix=='pdf')open_pdf(f,w,h) else ragg::agg_png(f,width=w,height=h,units='mm',res=600,background='white')
    grid.newpage();grid.draw(result);dev.off()
  }
  saveRDS(result,file.path(OUT,'panels',paste0(id,'.rds')))
  result
}
