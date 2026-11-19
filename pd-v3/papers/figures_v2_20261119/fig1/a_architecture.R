# Figure 1a - PeriScope architecture. Art rendered in Blender (fig1/blender/render_1a.py, Cycles, assets from
# pd_product_figure_design/workflow_assets) at the exact block size; labels added here as 7-pt vector text from the same
# layout file (fig1/blender/layout_1a.json). Main line: blood RNA -> encoder -> blood state (locked) -> coupled potential
# -> brain state -> decoder -> brain RNA -> glial programs -> dose response. Branches: the blood state feeds the protein
# head, whose ADT output (supervised by same-cell CITE-seq ADT) enters the potential as the dose query; observed brain
# snRNA supervises the predicted brain RNA through the Sinkhorn loss.
PANEL<-list(id='1a_architecture',blocks=c('1x3'),what='Model architecture (Blender art, vector labels)')
make_panel<-function(b){W<-b$w;H<-b$h;assign('PANEL_H',H,envir=globalenv())
 D<-file.path(Sys.getenv('FIGV2_ROOT'),'fig1/blender');L<-jsonlite::fromJSON(src_log(file.path(D,'layout_1a.json')),simplifyVector=FALSE)
 stopifnot(abs(L$W-W)<1e-6,abs(L$H-H)<1e-6)
 dm<-jsonlite::fromJSON(src_log(file.path(D,'1a_art.rgba.json')));rw<-readBin(src_log(file.path(D,'1a_art.rgba')),'raw',n=dm$width*dm$height*4)
 v<-as.integer(rw);img<-matrix(grDevices::rgb(v[c(TRUE,FALSE,FALSE,FALSE)],v[c(FALSE,TRUE,FALSE,FALSE)],v[c(FALSE,FALSE,TRUE,FALSE)],v[c(FALSE,FALSE,FALSE,TRUE)],maxColorValue=255),nrow=dm$height,byrow=TRUE)
 g<-list(rasterGrob(img,x=unit(0,'mm'),y=unit(0,'mm'),width=unit(W,'mm'),height=unit(H,'mm'),just=c(0,0),interpolate=TRUE))
 for(l in L$labels){hj<-if(is.null(l$hj))'centre' else l$hj;g<-c(g,list(dtxt(l$t,l$x,l$y,just=c(hj,'centre'))))}
 dgroup(g)}
