source(file.path(Sys.getenv('PD_PROJECT_ROOT',unset=getwd()),'pd_product_figure_design/panel_grid/source/paths.R'))
suppressPackageStartupMessages({library(ggplot2);library(grid);library(cowplot);library(jsonlite);library(ggrepel);library(ggalluvial);library(ComplexHeatmap);library(circlize)})
ROOT<-normalizePath('pd_product_figure_design/panel_grid');setwd('pd_product_figure_design/template_semantic_review')
map<-fromJSON(file.path(ROOT,'grid.json'));theme_set(theme_classic(base_size=14,base_family='Helvetica'))
# Text remains measured in points, while the plot layout is recomputed at the new device size.
style<-function(p){
 fields<-c('text','axis.text','axis.text.x','axis.text.y','axis.title','axis.title.x','axis.title.y','legend.text','legend.title','strip.text','strip.text.x','strip.text.y')
 th<-p$theme
 for(k in fields){
  old<-th[[k]]; if(inherits(old,'element_blank') || (grepl('^axis.title[.]',k)&&inherits(th[['axis.title']],'element_blank')) || (grepl('^axis.text[.]',k)&&inherits(th[['axis.text']],'element_blank')))next
  size<-if(grepl('axis.title',k))16 else 14
  if(is.null(old))old<-element_text()
  old$size<-size;old$family<-'Helvetica';if(grepl('title|strip',k))old$face<-'bold'
  th[[k]]<-old
 }
 p$theme<-th
 p<-p+theme(plot.margin=margin(30,18,16,18),panel.spacing=unit(7,'mm'),legend.key.size=unit(5,'mm'),legend.box.spacing=unit(3,'mm'))
 for(i in seq_along(p$layers)){
  if(inherits(p$layers[[i]]$geom,c('GeomText','GeomTextRepel','GeomLabelRepel'))){p$layers[[i]]$aes_params$size<-14/ggplot2::.pt;p$layers[[i]]$aes_params$family<-'Helvetica'}
 }
 p
}
# Composite grobs are relaid out by grid; semantic text sizes are unified inside nested legends/axes.
fixg<-function(g,context=""){
 context<-paste(context,g$name)
 if(inherits(g,'text')&&!is.null(g$gp$fontsize)){
  z<-g$gp$fontsize;g$gp$fontsize<-if(grepl('axis.title',context))16 else ifelse(z>=18,16,14);g$gp$fontfamily<-'Helvetica'
 }
 if(!is.null(g$grobs))g$grobs<-lapply(g$grobs,function(a)fixg(a,context))
 if(!is.null(g$children))g$children<-do.call(gList,lapply(g$children,function(a)fixg(a,context)))
 g
}
letter<-function(id)grid.text(toupper(substr(id,2,2)),unit(5,'pt'),unit(1,'npc')-unit(10,'pt'),just=c('left','top'),gp=gpar(fontsize=20,fontface='bold',fontfamily='Helvetica'))
selected<-Sys.getenv('PANEL_IDS','');ids<-if(nchar(selected))strsplit(selected,',')[[1]]else map$id
for(i in seq_len(nrow(map))){
 r<-map[i,];if(is.na(r$old)||r$id%in%c('2d','3e','5c')||!r$id%in%ids)next
 p<-readRDS(file.path(ROOT,'objects',paste0(r$old,'.rds')))
 composite<-r$id%in%c('1b','1c','1g')
 if(!composite)p<-style(p)
 if(r$id=='1b'){source(file.path(ROOT,'source/rings_grid.R'));p<-rings_grid();composite<-FALSE;theme_set(theme_classic(base_size=14,base_family='Helvetica'))}
 if(r$id=='1c'){source(file.path(ROOT,'source/umap_grid.R'));p<-umap_grid();composite<-FALSE;theme_set(theme_classic(base_size=14,base_family='Helvetica'))}
 if(r$id=='3b'){sc<-p$scales$get_scales('fill');sc$breaks<-c(sc$limits[1],0,sc$limits[2]);sc$labels<-function(x)sprintf('%.2f',x);p<-p+guides(fill=guide_colorbar(barheight=unit(38,'mm')))}
 if(r$id=='2b'){p<-p+theme(legend.position='bottom');cs<-p$scales$get_scales('colour');cs$labels<-c('Held-out ID 5'='P6','Held-out ID 6'='P7','Held-out ID 7'='P8')}
 if(r$id=='1d')p<-p+facet_grid(~stage,scales='free_x',space='free_x',labeller=labeller(stage=c('Normal'='Control','Early PD'='Early PD','Late PD'='Late PD')))
 if(r$id=='3f')p<-p+scale_y_discrete(labels=function(x)vapply(x,function(t)paste(strwrap(t,25),collapse='\n'),''))+theme(legend.position='bottom')+guides(color=guide_colorbar(title.position='top',barwidth=unit(45,'mm')))
 if(r$id%in%c('1h','5h','5i'))p<-p+labs(x=NULL,y=NULL)+facet_wrap(~region,ncol=5,labeller=labeller(region=function(x)gsub(' | ','\n',x,fixed=TRUE)))+theme(strip.text=element_text(size=14,face='bold'),panel.spacing=unit(4,'mm'))
 if(r$id=='5i'){sc<-p$scales$get_scales('colour');sc$breaks<-c(-.9,1,10,100);p<-p+guides(colour=guide_colorbar(title.position='top',barwidth=unit(70,'mm'),barheight=unit(4,'mm')))}
 if(r$id=='2f')p<-p+scale_y_continuous(limits=c(.54,.612),breaks=seq(.54,.61,.01),expand=c(0,0))
 if(r$id=='3d'){p<-p+labs(x='Discovery mean gene slope',y='Extension mean\ngene slope',color=NULL)+guides(color=guide_legend(nrow=1))+theme(legend.position='bottom',legend.box='horizontal',legend.margin=margin(0,0,0,0));p$layers[[3]]$geom_params$seed<-6674;p$layers[[3]]$geom_params$max.overlaps<-Inf;p$layers[[3]]$geom_params$max.time<-2}
 if(r$id=='3g')p<-p+scale_y_continuous(expand=expansion(mult=c(.12,.28)))
 if(r$id=='5f')p<-p+labs(y='Association r and\ndonor-bootstrap interval')
 if(r$id=='1f'){
  p$layers[[3]]$aes_params$size<-14/ggplot2::.pt
  p$layers[[3]]$mapping$label<-aes(label=ifelse(after_stat(count)>=3000,gsub('Monocyte_','Monocyte\n',as.character(after_stat(stratum))),''))$label
  p$layers[[4]]$aes_params$size<-14/ggplot2::.pt
  p<-p+theme(plot.margin=margin(30,24,20,60))
 }
 if(r$id%in%c('4b','6c')){p$layers[[1]]$aes_params$size<-8.7;p$layers[[2]]$aes_params$size<-12/ggplot2::.pt}
 if(r$id%in%c('2c','3g','4f','1d','6b','6e'))for(k in seq_along(p$layers))if(inherits(p$layers[[k]]$geom,'GeomText'))p$layers[[k]]$aes_params$size<-12/ggplot2::.pt
 pdf(file.path(PANEL_BUILD,paste0(r$id,'.pdf')),width=r$width_in,height=r$height_in,family='Helvetica',useDingbats=FALSE)
 if(composite){grid.newpage();grid.draw(fixg(ggplotGrob(p)))}else print(p)
 letter(r$id);dev.off();cat('DRAW',r$id,'\n')
}
# Direct calls to copied native matrix/pathway renderers; regenerate geometry at device dimensions.
if(any(c('2d','3e')%in%ids)){
 source('R/copied/base_helpers.R');source('R/copied/base_main.R');source('R/copied/base_supplement.R');source('R/copied/calibrated_panels.R')
 used<-list();ne<-new.env(parent=globalenv());sys.source('R/copied/render_native.R',ne);sys.source('R/copied/fig4_corrections.R',ne);ne$DATA_CACHE<-list();ne$rd<-function(key)ne$DATA_CACHE[[key]]
 rd<-function(f)read.csv(f,check.names=FALSE,stringsAsFactors=FALSE)
 short<-function(x){x<-gsub('nonclassical_mono','nMono',x);x<-gsub('classical_mono','cMono',x);x<-gsub('microglia_homeostatic','MG-H',x);x<-gsub('microglia_mhc2','MG-AP',x);x<-gsub('astro','Astro',x);x<-gsub('__x__',' > ',x);gsub('__',' | ',x)}
 src<-readLines(file.path(ROOT,'source/captured_render.R'))
 if('2d'%in%ids){
  a<-grep('^d<-rd.*direction_seeds',src);b<-grep('^saveRDS.*native2C',src)
  eval(parse(text=src[a:(b-1)]))
  oldtxt<-ne$txt;ne$txt<-function(t,x,y,size=14,col='#151515',rot=0,just='centre',face='plain'){if(t=='C'&&x==12)return(invisible(NULL));oldtxt(t,x,y,if(face=='bold')16 else 14,col,rot,just,face)}
  pdf(file.path(PANEL_BUILD,'2d.pdf'),width=10,height=10,family='Helvetica',useDingbats=FALSE);nativeC();letter('2d');dev.off();ne$txt<-oldtxt
 }
 if('3e'%in%ids){
  a<-grep('^terms3<-',src);b<-grep('^saveRDS.*native3D',src);eval(parse(text=src[a:(b-1)]))
  oldtxt<-ne$txt;ne$txt<-function(t,x,y,size=14,col='#151515',rot=0,just='centre',face='plain'){if(t==ne$TERM_LABELS[2]){y<-800;x<-965};if(t==ne$TERM_LABELS[1])y<-800;oldtxt(t,x,y,if(face=='bold')16 else 14,col,rot,just,face)}
  pdf(file.path(PANEL_BUILD,'3e.pdf'),width=10,height=10,family='Helvetica',useDingbats=FALSE)
  grid.newpage();pushViewport(viewport(x=.5,y=.49,width=.94,height=.94,xscale=c(520,1060),yscale=c(1300,760)));ne$draw_pathway6();popViewport();letter('3e');dev.off()
 }
}
if('5c'%in%ids){
 setwd('../grouped_figures')
 src<-readLines(file.path(ROOT,'source/copied_grouped.R'))
 short<-function(x){x<-gsub('nonclassical_mono','nMono',x);x<-gsub('classical_mono','cMono',x);x<-gsub('microglia_homeostatic','MG-H',x);x<-gsub('microglia_mhc2','MG-AP',x);x<-gsub('astro','Astro',x);x<-gsub('__x__',' > ',x);gsub('__',' | ',x)}
 ev<-read.csv('../data/program_evidence.csv');oo<-ev$combo[!is.na(ev$spatial_meta_q)][order(ev$spatial_meta_q[!is.na(ev$spatial_meta_q)])]
 status<-function(x)ifelse(is.na(x),'Not tested',ifelse(x<=.05,'Pass','No pass'));statuscols<-c(Pass='#398b82','No pass'='#d4d9dc','Not tested'='#e4b75e');lg<-list(labels_gp=gpar(fontsize=14),title_gp=gpar(fontsize=14,fontface='bold'))
 code<-src[44:47];code<-gsub('fontsize=13','fontsize=14',code,fixed=TRUE);code<-gsub('fontsize=11','fontsize=12',code,fixed=TRUE);eval(parse(text=code))
 pdf(file.path(PANEL_BUILD,'5c.pdf'),width=10,height=10,family='Helvetica',useDingbats=FALSE);draw(hm,padding=unit(c(9,7,14,7),'mm'),merge_legends=TRUE);letter('5c');dev.off()
}
