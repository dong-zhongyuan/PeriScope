source(file.path(HERE,'..','figures_r','common.R'))
POLICY <- c(drop='Primary', infer='EAF-assisted')
PCOL <- c(drop='#758AA2', infer='#80609B')
PLINE <- c(drop='solid', infer='dashed')
cfg <- jsonlite::fromJSON(file.path(DATA,'CONFIG.json'))
cats <- read.csv(file.path(DATA,'categories.csv'),check.names=FALSE)
as_csv_flag <- function(x) {
  stopifnot(all(tolower(as.character(x)) %in% c('true','false')))
  tolower(as.character(x))=='true'
}
cats$is_qc <- as_csv_flag(cats$is_qc)
read_table <- function(name) {
  d<-read.csv(file.path(DATA,paste0(name,'.csv')),check.names=FALSE)
  if('is_qc' %in% names(d))d$is_qc<-as_csv_flag(d$is_qc)
  d
}
policy_factor <- function(x) factor(x,levels=c('drop','infer'),labels=unname(POLICY))
wrap_atlas <- function(p,id,title,w,h) {
  stopifnot(is.null(p$labels$title),is.null(p$labels$subtitle),is.null(p$labels$caption))
  ragg::agg_capture(width=w,height=h,units='mm',res=300)
  g<-ggplotGrob(p);dev.off()
  key<-if(id=='B')gTree(children=gList(
    segmentsGrob(x0=unit(135,'mm'),x1=unit(141,'mm'),y0=unit(1,'npc')-unit(2.5,'mm'),
      y1=unit(1,'npc')-unit(2.5,'mm'),gp=gpar(col=INK,lwd=.4)),
    segmentsGrob(x0=unit(138,'mm'),x1=unit(138,'mm'),y0=unit(1,'npc')-unit(1.8,'mm'),
      y1=unit(1,'npc')-unit(3.2,'mm'),gp=gpar(col=INK,lwd=.3)),
    txt('z range',142,1.2,6),txt('x missing',158,1.2,6))) else nullGrob()
  result<-gTree(children=gList(rectGrob(gp=gpar(fill='white',col=NA)),g,
    txt(tolower(id),.6,.8,8,'bold'),txt(title,4.6,1.2,6.5,'bold'),key))
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
