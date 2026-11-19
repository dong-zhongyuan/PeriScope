# Figure 4a - genetic design (Blender art: fig4/blender/layout_4a.json): cis-pQTL instrument -> plasma protein -> brain
# imaging (MR, coloured), reverse MR (dashed), and colocalisation of the protein and imaging association signals at the
# same locus. Labels are 7-pt vector text.
PANEL<-list(id='4a_design',blocks=c('1x1'),what='Genetic analysis design (Blender)')
source(file.path(Sys.getenv('FIGV2_ROOT'),'_system/blender_panel.R'))
make_panel<-function(b)blender_panel(b,file.path(Sys.getenv('FIGV2_ROOT'),'fig4/blender'),'4a')
