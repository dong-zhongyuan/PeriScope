# Figure 5a - antigen funnel (Blender art from fig5/blender/layout_5a.json, built by prep_5a.py from
# target_funnel_final_20261119/antigen_gate_matrix.csv): one disc per gate holding every antigen that reaches it, one
# sphere per antigen (grey passes, pink removed at this gate, CD71 and CD119 in their colours); disc area follows the
# count; last disc = the two final targets. Labels and key are 7-pt vector text.
PANEL<-list(id='5a_antigen_funnel',blocks=c('1x4'),what='Antigen funnel through all gates (Blender)')
source(file.path(Sys.getenv('FIGV2_ROOT'),'_system/blender_panel.R'))
make_panel<-function(b){D<-file.path(Sys.getenv('FIGV2_ROOT'),'fig5/blender');g<-blender_panel(b,D,'5a')
 K<-jsonlite::fromJSON(file.path(D,'layout_5a.json'))$keys;H<-b$h
 k<-key_row(list(ki(kg_dot(K$keep,NA,1.6),'Retained'),ki(kg_dot(K$gone,NA,1.6),'Removed'),ki(kg_dot(K$CD71,NA,1.6),'CD71'),ki(kg_dot(K$CD119,NA,1.6),'CD119')),1,H-1.9,gap=2.4,h=H)
 # real-world outcome marker: schematic cumulative-incidence curves (treated below control),
 # set apart from the funnel by a dashed connector because this check follows the funnel and is not a gate
 x0<-151.6;x1<-170.2;yb<-11.6;xs<-c(x0,154,154,157,157,160,160,163,163,166.5,166.5,x1)
 km<-list(dline(x0,3.7,x0,yb,col=GUIDE,lwd=glwd('hair'),h=H),dline(x0,yb,x1,yb,col=GUIDE,lwd=glwd('hair'),h=H),
  dpath(xs,c(11.4,11.4,10.3,10.3,9.1,9.1,7.9,7.9,6.7,6.7,5.5,5.5),col=MID,lwd=glwd('data'),h=H),
  dpath(xs,c(11.4,11.4,10.9,10.9,10.3,10.3,9.8,9.8,9.3,9.3,8.8,8.8),col=LIGHT,lwd=glwd('data'),h=H),
  dcurve(166.4,16.4,164.6,12.6,curv=.35,col=GUIDE,lwd=glwd('hair'),lty=2,arrow=TRUE,h=H),
  dtxt('Outcomes',160.9,2.6,col=INK,h=H))
 gTree(children=gList(g,do.call(gList,k),do.call(gList,km)))}
