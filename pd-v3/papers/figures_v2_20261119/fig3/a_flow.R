# Figure 3a - program readout flow (Blender art: fig3/blender/layout_3a.json): measured-dose protein queries -> response
# slope -> discovery seeds (42-46) and confirmation seeds (47-51) -> frozen program (lock) -> matched brain states and
# spatial sections. Labels are 7-pt vector text.
PANEL<-list(id='3a_readout_flow',blocks=c('1x1'),what='Program readout flow (Blender)')
source(file.path(Sys.getenv('FIGV2_ROOT'),'_system/blender_panel.R'))
make_panel<-function(b)blender_panel(b,file.path(Sys.getenv('FIGV2_ROOT'),'fig3/blender'),'3a')
