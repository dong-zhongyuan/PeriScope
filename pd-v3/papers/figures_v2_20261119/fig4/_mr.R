# shared paths for Figure 4 (read-only result tables)
MR<-'/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results'
OV<-file.path(MR,'full_atlas_20261116/full/overview_figure/plot_data')
TF<-'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/target_funnel_final_20261119'
FT_GENE<-c(TFRC='CD71',IFNGR1='CD119');FT_COL<-c(TFRC=unname(TARGET_COL['CD71']),IFNGR1=unname(TARGET_COL['CD119']))
key_title<-function(t,x,y,h)textGrob(t,unit(x,'mm'),unit(h-y,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,fontface='bold',col=INK))
