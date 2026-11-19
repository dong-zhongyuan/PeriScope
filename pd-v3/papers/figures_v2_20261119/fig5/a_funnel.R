# Figure 5a - antigen funnel (Blender art from fig5/blender/layout_5a.json, built by prep_5a.py from
# target_funnel_final_20261119/antigen_gate_matrix.csv): one disc per gate holding every antigen that reaches it, one
# sphere per antigen (grey passes, pink removed at this gate, CD71 and CD119 in their colours); disc area follows the
# count; last disc = the two final targets. Labels and key are 7-pt vector text.
PANEL<-list(id='5a_antigen_funnel',blocks=c('1x4'),what='Antigen funnel through all gates (Blender)')
source(file.path(Sys.getenv('FIGV2_ROOT'),'_system/blender_panel.R'))
make_panel<-function(b){D<-file.path(Sys.getenv('FIGV2_ROOT'),'fig5/blender');g<-blender_panel(b,D,'5a')
 K<-jsonlite::fromJSON(file.path(D,'layout_5a.json'))$keys;H<-b$h
 k<-key_row(list(ki(kg_dot(K$keep,NA,1.6),'Passes'),ki(kg_dot(K$gone,NA,1.6),'Removed'),ki(kg_dot(K$CD71,NA,1.6),'CD71'),ki(kg_dot(K$CD119,NA,1.6),'CD119')),1,H-1.9,gap=2.4,h=H)
 gTree(children=gList(g,do.call(gList,k)))}
