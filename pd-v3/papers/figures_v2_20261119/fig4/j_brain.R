# Figure 4j - where in the brain the two final targets act. Four sections of the 1 mm MNI152 template.
# CD71/TFRC: the substantia nigra and the pallidum, the structures carrying its MR estimates with
# colocalisation support (Talairach Brodmann and Harvard-Oxford masks). CD119/IFNGR1: the two nodes of
# its only imaging hit, one edge of the UK Biobank resting-state network matrix (ICA100 edge 868);
# UK Biobank's edge list gives that edge as good nodes 7 and 43, which its good-component list maps to
# group-ICA components 8 and 45 of rfMRI_ICA_d100.nii.gz. Each component is drawn above z = 8 on the
# section where it is largest. Sections are rendered by fig4/brain/render_brain.py and placed as raw
# RGBA at the exact size, as the Blender art is.
PANEL<-list(id='4j_target_brain',blocks=c('1x4'),what='Brain sections of the two final targets')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
source(file.path(Sys.getenv('FIGV2_ROOT'),'_system/blender_panel.R'))
BR<-file.path(Sys.getenv('FIGV2_ROOT'),'fig4/brain')
SEC<-list(list(s='sn', lab='Substantia nigra',       m=c(L=2.2,R=-2.2)),
          list(s='gp', lab='Globus pallidus',        m=c(L=2.2,R=-2.2)),
          list(s='dmn',lab='Default mode network',   m=c(P=2.2,A=-2.2)),
          list(s='amt',lab='Anterior–medial temporal',m=c(L=2.2,R=-2.2)))
make_panel<-function(b){
 W<-b$w;H<-b$h;M<-read_src(file.path(BR,'brain.render.json'))
 IW<-M$display_mm[1];IH<-M$display_mm[2];Y0<-4.6;GAP<-(W-2*1.05-4*IW)/3
 x0<-1.05+(seq_along(SEC)-1)*(IW+GAP)
 g<-list()
 for(i in seq_along(SEC)){S<-SEC[[i]];u<-M$sections[[S$s]]
  g<-c(g,list(rasterGrob(read_rgba(file.path(BR,paste0(S$s,'.rgba'))),x=unit(x0[i],'mm'),y=unit(H-Y0-IH,'mm'),
    width=unit(IW,'mm'),height=unit(IH,'mm'),just=c(0,0),interpolate=TRUE)))
  for(k in seq_along(S$m))g<-c(g,list(dtxt(names(S$m)[k],x0[i]+if(S$m[[k]]>0)S$m[[k]] else IW+S$m[[k]],Y0+2.1,h=H)))
  g<-c(g,list(dtxt(sprintf('%s  %s %s',S$lab,u$plane,fnum(u$coord,0)),x0[i]+IW/2,Y0+IH+2.4,h=H)))}
 k<-list(ki(kg_box(unname(TARGET_COL['CD71']),NA,2.2,2.2),'CD71/TFRC'),
         ki(kg_box(unname(TARGET_COL['CD119']),NA,2.2,2.2),'CD119/IFNGR1'))
 dgroup(c(g,key_row(k,W-1-key_width(k),2.4,h=H)))}
