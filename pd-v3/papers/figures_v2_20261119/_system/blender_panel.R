# R side of the Blender schematics: place the Blender render (raw RGBA written by to_rgba.py) at the exact block size and
# add the layout's labels as 7-pt vector text. Label fields: t, x, y (mm from block top), hj (left/centre/right), bold, col.
read_rgba<-function(f){dm<-jsonlite::fromJSON(src_log(paste0(f,'.json')));v<-as.integer(readBin(src_log(f),'raw',n=dm$width*dm$height*4))
 matrix(grDevices::rgb(v[c(TRUE,FALSE,FALSE,FALSE)],v[c(FALSE,TRUE,FALSE,FALSE)],v[c(FALSE,FALSE,TRUE,FALSE)],v[c(FALSE,FALSE,FALSE,TRUE)],maxColorValue=255),nrow=dm$height,byrow=TRUE)}
blender_panel<-function(b,dir,stem){W<-b$w;H<-b$h;assign('PANEL_H',H,envir=globalenv())
 L<-jsonlite::fromJSON(src_log(file.path(dir,paste0('layout_',stem,'.json'))),simplifyVector=FALSE);stopifnot(abs(L$W-W)<1e-6,abs(L$H-H)<1e-6)
 g<-list(rasterGrob(read_rgba(file.path(dir,paste0(stem,'_art.rgba'))),x=unit(0,'mm'),y=unit(0,'mm'),width=unit(W,'mm'),height=unit(H,'mm'),just=c(0,0),interpolate=TRUE))
 for(l in L$labels){hj<-if(is.null(l$hj))'centre' else l$hj;g<-c(g,list(dtxt(l$t,l$x,l$y,just=c(hj,'centre'),bold=isTRUE(l$bold),col=if(is.null(l$col))INK else l$col)))}
 dgroup(g)}
