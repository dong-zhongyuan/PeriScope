work <- Sys.getenv('SPATIAL_WORK_DIR','/public/home/mengxl/dzy/pd_product/_work/spatial_subtype_20261003')
out <- Sys.getenv('SPATIAL_OUTPUT_DIR','/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/spatial_subtype_20261003')
.libPaths(c(file.path(work,'R_library'),.libPaths()))
library(ggplot2)
dir.create(file.path(out,'preview'),showWarnings=FALSE)
theme_set(theme_classic(base_size=14)+theme(plot.title=element_text(face='bold'),strip.background=element_rect(fill='grey95',colour='black'),strip.text=element_text(face='bold')))
states<-c('astro_fibrous','astro_protoplasmic','microglia_homeostatic','microglia_mhc2')
labels<-c(astro_fibrous='Fibrous astrocytes',astro_protoplasmic='Protoplasmic astrocytes',microglia_homeostatic='Homeostatic microglia',microglia_mhc2='MHC-II microglia')
truth<-read.csv(file.path(out,'inputs/qc_truth.csv'),row.names=1,check.names=FALSE)
est<-read.csv(file.path(out,'samples/QC/cell_state_RNA_fractions.csv'),row.names=1,check.names=FALSE)
q<-do.call(rbind,lapply(states,function(st)data.frame(state=labels[st],truth=truth[rownames(est),st],estimated=est[,st])))
p<-ggplot(q,aes(truth,estimated))+geom_abline(slope=1,intercept=0,colour='grey60',linewidth=.6)+geom_point(size=1.8,alpha=.55,colour='#2166AC')+
  facet_wrap(~state,nrow=2)+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(x='Known RNA fraction',y='Estimated RNA fraction',title='Cell-state recovery in held-out reference mixtures')
ggsave(file.path(out,'preview/state_recovery.png'),p,width=10,height=9,dpi=220,bg="white")
ggsave(file.path(out,'preview/state_recovery.pdf'),p,width=10,height=9,bg="white")
if('QC' %in% commandArgs(trailingOnly=TRUE)) quit(save='no')
meta<-read.csv(file.path(out,'inputs/spatial_metadata.csv'),row.names=1)
all<-do.call(rbind,lapply(unique(meta$sample_gsm),function(id){
  w<-read.csv(file.path(out,'samples',id,'cell_state_RNA_fractions.csv'),row.names=1,check.names=FALSE)
  m<-meta[rownames(w),];m$section<-paste(m$donor,m$condition)
  m$x<-m$x-mean(m$x);m$y<-m$y-mean(m$y)
  cbind(m,w)
}))
for(st in states){
  d<-all;d$fraction<-d[[st]]
  p<-ggplot(d,aes(x,-y,colour=fraction))+geom_point(size=1.1,shape=16)+
    scale_colour_viridis_c(option='magma',limits=c(0,max(d$fraction)),name='RNA fraction')+
    facet_wrap(~section,ncol=5)+theme_void(base_size=14)+
    theme(strip.text=element_text(face='bold',size=12),legend.title=element_text(size=13),plot.title=element_text(face='bold',size=18))+
    labs(title=paste(labels[st],': spatial RNA contribution'))
  # Centering changes position only; common scales and equal coordinates retain
  # image-coordinate geometry without stretching individual point clouds.
  p<-p+coord_equal()
  ggsave(file.path(out,'preview',paste0(st,'_maps.png')),p,width=16,height=7,dpi=220,bg="white")
  ggsave(file.path(out,'preview',paste0(st,'_maps.pdf')),p,width=16,height=7,bg="white")
}
