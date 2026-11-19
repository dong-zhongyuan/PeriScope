# render_panel.R - render one panel script at one block size.
# Usage (after `. _system/env.sh`):  Rscript _system/render_panel.R <panel.R> <block RxC> [<out_dir>]
# Writes <out_dir>/<id>_<block>.pdf (vector, page = exact block), .png (600 dpi), _preview.png (150 dpi),
# .manifest.json (sources read, grob text audit, bar audit) and .check.json (check_panel.py), and prints the check.
# Panel script contract:  PANEL <- list(id = '1b', blocks = c('2x1', '1x2'), what = '...')  and
#   make_panel(b) -> a grob or ggplot laid out for block b (b$w, b$h in mm; PANEL_W / PANEL_H are also set). The code
#   chooses positions and sizes for each block; nothing is scaled and every text stays 6 pt.
args<-commandArgs(TRUE)
if(length(args)<2)stop('usage: Rscript render_panel.R <panel.R> <block> [<out_dir>]',call.=FALSE)
SYS_DIR<-Sys.getenv('FIGV2_SYS');if(!nzchar(SYS_DIR))stop('source _system/env.sh first',call.=FALSE)
source(file.path(SYS_DIR,'common.R'))
script<-normalizePath(args[1],mustWork=TRUE);b<-block(args[2])
out<-if(length(args)>=3)args[3]else file.path(dirname(script),'panels')
dir.create(out,showWarnings=FALSE,recursive=TRUE);out<-normalizePath(out)
PANEL_W<-b$w;PANEL_H<-b$h;BLOCK<-b$name;SLOT_W<-PANEL_W;SLOT_H<-PANEL_H
src_reset();bars_reset();PANEL<-NULL;make_panel<-NULL
source(script,local=globalenv())
if(!is.list(PANEL)||!is.function(make_panel))stop(script,' must define PANEL (list) and make_panel(b)',call.=FALSE)
id<-PANEL$id
if(!is.null(PANEL$blocks)&&!b$name%in%PANEL$blocks)stop(sprintf('panel %s has no layout for block %s (layouts: %s)',id,b$name,paste(PANEL$blocks,collapse=' ')),call.=FALSE)
op<-meas_open();res<-make_panel(b);if(inherits(res,'ggplot'))res<-gg_grob(res,id);if(op)dev.off()
if(!inherits(res,'grob'))stop('make_panel() must return a grob or a ggplot',call.=FALSE)
final<-gTree(children=gList(res),name=paste0('panel_',id))
audit<-text_audit(final)
base<-file.path(out,paste0(id,'_',b$name))
raw<-tempfile(fileext='.pdf')
open_block_pdf(raw,b);grid.draw(final);invisible(dev.off());pdf_crop(raw,paste0(base,'.pdf'),b$w,b$h);unlink(raw)
open_block_png(paste0(base,'.png'),b,600);grid.draw(final);invisible(dev.off())
open_block_png(paste0(base,'_preview.png'),b,150);grid.draw(final);invisible(dev.off())
man<-list(panel=id,block=b$name,size_mm=c(b$w,b$h),what=PANEL$what%||%'',script=script,rendered=format(Sys.time(),'%Y-%m-%dT%H:%M:%S'),
 r=R.version.string,ggplot2=as.character(packageVersion('ggplot2')),sources=src_list(),bars=bars_list(),text=audit,
 files=list(pdf=paste0(base,'.pdf'),png=paste0(base,'.png'),preview=paste0(base,'_preview.png')))
jsonlite::write_json(man,paste0(base,'.manifest.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA,null='null')
st<-system2(Sys.getenv('FIGV2_PY'),c(shQuote(file.path(SYS_DIR,'check_panel.py')),shQuote(paste0(base,'.pdf')),'--png',shQuote(paste0(base,'.png')),
 '--manifest',shQuote(paste0(base,'.manifest.json'))),stdout=TRUE,stderr=TRUE)
cat(st,sep='\n')
if(!is.null(attr(st,'status'))&&attr(st,'status')!=0)quit(status=1)
