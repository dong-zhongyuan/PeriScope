args <- commandArgs(trailingOnly=TRUE)
ROOT <- normalizePath(args[1])
DATA <- file.path(ROOT,'papers/01_mr_imaging_scan/results/public_cd40_20261116')
OUT <- file.path(DATA,'figures')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(OUT,'qa'),showWarnings=FALSE)
source(file.path(ROOT,'papers/01_mr_imaging_scan/figures_r/common.R'))
readz <- function(name) read.csv(gzfile(file.path(DATA,name)),check.names=FALSE)
protein <- readz('protein_region_all_rows.csv.gz')
regional <- data.frame(pos=protein$pos37,z=abs(protein$beta_x/protein$se_x),trait='SCALLOP CD40')
compare <- list()
for(idp in c('0017','0018')) {
  hemi <- if(idp=='0017') 'Left pallidum' else 'Right pallidum'
  d <- readz(paste0(idp,'_imaging_region_all_rows.csv.gz'))
  regional <- rbind(regional,data.frame(pos=d$pos37,z=abs(d$beta_y/d$se_y),trait=hemi))
  s <- readz(paste0(idp,'_shared_snps.csv.gz'))
  compare[[idp]] <- data.frame(z_protein=abs(s$beta_x/s$se_x),z_imaging=abs(s$beta_y/s$se_y),
                       position_Mb=s$pos37/1e6,hemi=hemi,lead=s$variant_y=='rs4810485')
}
regional$trait <- factor(regional$trait,levels=c('SCALLOP CD40','Left pallidum','Right pallidum'))
regional$view <- '±1 Mb'
zoom <- regional[abs(regional$pos-44747947)<=100000,];zoom$view <- '±100 kb'
regional <- rbind(regional,zoom);regional$view <- factor(regional$view,levels=c('±1 Mb','±100 kb'))
lead <- regional[regional$pos==44747947,]
a <- ggplot(regional,aes(pos/1e6,z)) +
  geom_point(colour='#758AA2',size=.22,alpha=.45,stroke=0) +
  geom_vline(xintercept=44.747947,linewidth=.3/.pt,linetype='dashed',colour='#C7C7C7')+
  geom_point(data=lead,shape=23,fill=INK,colour=INK,size=1.15,stroke=.18)+
  facet_grid(trait~view,scales='free')+
  labs(x='Chromosome 20 position (GRCh37, Mb)',y='Absolute association z')+
  theme(strip.text.y=element_text(angle=0,hjust=0),panel.spacing.y=unit(2,'mm'),panel.spacing.x=unit(5,'mm'))
bdata <- do.call(rbind,compare)
b <- ggplot(bdata,aes(z_protein,z_imaging))+
  geom_point(aes(colour=position_Mb),size=.25,alpha=.6,stroke=0)+
  geom_point(data=bdata[bdata$lead,],shape=23,fill=INK,colour=INK,size=1.3,stroke=.2)+
  scale_colour_gradientn(colours=c('#719DB9','#E5E5E5','#A66868'),name='Position (Mb)')+
  guides(colour=guide_colourbar(display='rectangles'))+
  facet_wrap(~hemi,nrow=1)+labs(x='CD40 absolute z',y='Imaging absolute z')+
  theme(legend.position='bottom',legend.key.width=unit(8,'mm'),panel.spacing=unit(3,'mm'))
s <- read.csv(file.path(DATA,'coloc_prior_sensitivity.csv'),colClasses=c(idp='character'),check.names=FALSE)
s$hemi <- ifelse(s$idp=='0017','Left pallidum','Right pallidum')
s$window <- factor(ifelse(s$window_bp==500000,'±500 kb','±1 Mb'),levels=c('±500 kb','±1 Mb'))
c <- ggplot(s,aes(p12_requested,PP.H4.abf,linetype=window))+
  geom_vline(xintercept=1e-5,linewidth=.3/.pt,colour='#C7C7C7')+
  geom_line(linewidth=.6/.pt,colour='#3A8B84')+
  geom_point(data=s[abs(s$p12_requested-1e-5)<1e-12,],size=1,colour='#3A8B84')+
  scale_x_log10(breaks=c(1e-6,1e-5,3e-5),labels=c('1e−6','1e−5','3e−5'))+
  scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=c(.01,.02)))+
  scale_linetype_manual(values=c('solid','dashed'))+
  facet_wrap(~hemi,nrow=1)+labs(x='Prior probability of a shared variant',y='Posterior probability H4',linetype=NULL)+
  theme(legend.position='bottom',panel.spacing=unit(3,'mm'))
panels <- list(a,b,c)
layout <- data.frame(id=c('a','b','c'),x=c(1,1,89),y=c(1,102,102),w=c(172,84,84),h=c(97,67,67),
                     title=c('CD40 regional association profiles','Shared SNP association pattern','Sensitivity to shared-variant prior'))
grobs <- list()
for(i in seq_len(nrow(layout))) {
  r <- layout[i,]
  ragg::agg_capture(width=r$w,height=r$h,units='mm',res=300)
  g <- ggplotGrob(panels[[i]]);dev.off()
  g <- gTree(children=gList(g,txt(r$id,.6,.8,8,'bold'),txt(r$title,4.6,1.2,6.5,'bold')))
  audit <- text_audit(g);stopifnot(all(audit$pt>=6-1e-6))
  write.csv(audit,file.path(OUT,'qa',paste0('text_',r$id,'.csv')),row.names=FALSE)
  grobs[[i]] <- g
}
draw <- function() {
 grid.newpage()
 for(i in seq_len(nrow(layout))) {
  r <- layout[i,]
  pushViewport(viewport(x=unit(r$x,'mm'),y=unit(170-r$y,'mm'),width=unit(r$w,'mm'),height=unit(r$h,'mm'),just=c('left','top')))
  grid.draw(grobs[[i]]);popViewport()
 }
}
open_pdf(file.path(OUT,'CD40_public_coloc.pdf'),174,170);draw();dev.off()
for(dpi in c(600,150)) {
 name <- if(dpi==600) 'CD40_public_coloc.png' else 'CD40_public_coloc_preview.png'
 ragg::agg_png(file.path(OUT,name),width=174,height=170,units='mm',res=dpi,background='white');draw();dev.off()
}
write.csv(layout,file.path(OUT,'qa/layout_mm.csv'),row.names=FALSE)
write_json(list(page_mm=c(174,170),panels=3,plotting_areas=10,minimum_font_pt=6,
                simple_bar_charts=0,regional_rows=nrow(regional),shared_snp_points=nrow(bdata),
                colour='Position, not LD',lead_marker='Historical rs4810485',
                assumptions='H4 uses single-causal ABF; region and prior sensitivities shown'),
           file.path(OUT,'qa/FIGURE_AUDIT.json'),pretty=TRUE,auto_unbox=TRUE)
