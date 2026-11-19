# Figure 2a - evaluation design (Blender art: fig2/blender/layout_2a.json, rendered by _system/blender_render.py).
# Top: eight routes from four blood source states to two brain recipient states. Bottom: evaluation settings with one
# cube per held-out unit - locked human brain donors (PD, control), independent cohorts, pig animals (saline, LPS) and
# LPS-treated mice (two studies). Labels are 7-pt vector text.
PANEL<-list(id='2a_design',blocks=c('1x1'),what='Evaluation design (Blender)')
source(file.path(Sys.getenv('FIGV2_ROOT'),'_system/blender_panel.R'))
make_panel<-function(b)blender_panel(b,file.path(Sys.getenv('FIGV2_ROOT'),'fig2/blender'),'2a')
