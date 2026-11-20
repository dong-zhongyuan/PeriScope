# shared paths for Figure 4 (read-only result tables)
MR<-'/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results'
OV<-file.path(MR,'full_atlas_20261116/full/overview_figure/plot_data')
TF<-'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/target_funnel_final_20261119'
FT_GENE<-c(TFRC='CD71',IFNGR1='CD119');FT_COL<-c(TFRC=unname(TARGET_COL['CD71']),IFNGR1=unname(TARGET_COL['CD119']))
key_title<-function(t,x,y,h)textGrob(t,unit(x,'mm'),unit(h-y,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,fontface='bold',col=INK))
# short phenotype names shared by 4i and 4k (SN = substantia nigra, GP = globus pallidus; DMN-temporal = resting-state
# connectivity between the default mode network and bilateral anterior/medial temporal lobe, UKB ICA100 edge 868)
SHORT<-c(QSM_Left_SN='QSM SN L',QSM_Right_pallidum='QSM GP R',QSM_Left_pallidum='QSM GP L',`T2*_Left_SN`='T2* SN L',
 IDP_SWI_T2star_left_pallidum='T2* GP L',IDP_T1_FAST_ROIs_R_thalamus='GM thalamus R',`rfMRI connectivity (ICA100 edge 868)`='DMN–temporal')
