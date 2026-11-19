# Reuse existing manuscript grid, typography and export conventions.
args<-commandArgs(trailingOnly=TRUE);ROOT<-normalizePath(args[1])
DATA<-file.path(ROOT,'papers/03_world_model_integration/results/perturbation_readout_20261117')
OUT<-file.path(DATA,'figures')
source(file.path(ROOT,'papers/01_mr_imaging_scan/figures_r/common.R'))
stopifnot(isTRUE(fromJSON(file.path(DATA,'INDEPENDENT_READOUT_VERIFICATION.json'))$complete),isTRUE(fromJSON(file.path(DATA,'INDEPENDENT_COMPARISON_VERIFICATION.json'))$complete))
dir.create(file.path(OUT,'qa'),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(OUT,'plot_data'),showWarnings=FALSE)
effects<-read.csv(file.path(DATA,'target_gene_readouts.csv'))
effects$target_gene<-factor(effects$target_gene,levels=c('CD86','IFNGR1','MERTK'))
effects$guide<-ifelse(effects$kind=='gene','Pooled TSS',sub('.*_tss_','',effects$label))
effects$guide<-factor(effects$guide,levels=c('g1','g2','g3','Pooled TSS'))
effects$well<-match(effects$sample,sort(unique(effects$sample)))
effects$position<-as.numeric(effects$guide)+(effects$well-4)*.04+ifelse(effects$assignment=='author_inclusive',-.015,.015)
effects$assignment<-factor(effects$assignment,levels=c('author_inclusive','single_positive_guide'),labels=c('Author inclusion','Single guide'))
a<-ggplot(effects,aes(position,target_gene_delta,colour=assignment,shape=assignment))+zero_h()+geom_point(size=.85,stroke=.3)+facet_wrap(~target_gene,nrow=1)+scale_x_continuous(breaks=1:4,labels=levels(effects$guide))+scale_colour_manual(values=c('#758AA2','#A66868'))+scale_shape_manual(values=c(16,1))+labs(x='TSS guide (seven capture wells per guide)',y='Target RNA change\n(mean log1p counts per 10,000)',colour=NULL,shape=NULL)+theme(legend.position='bottom',axis.text.x=element_text(angle=30,hjust=1))
d<-read.csv(gzfile(file.path(DATA,'all_direction_comparisons.csv.gz')))
d<-d[d$kind=='gene'&d$gene_scope=='exclude_perturbed_gene',]
d$target_gene<-factor(d$target_gene,levels=c('CD86','IFNGR1','MERTK'))
corlim<-max(.05,ceiling(max(abs(d$spearman),na.rm=TRUE)*20)/20)
axes<-c('cDC__x__microglia_homeostatic','cDC__x__microglia_mhc2','classical_mono__x__microglia_homeostatic','classical_mono__x__microglia_mhc2','nonclassical_mono__x__microglia_homeostatic','nonclassical_mono__x__microglia_mhc2','pDC__x__microglia_homeostatic','pDC__x__microglia_mhc2')
axlabs<-c('cDC → homeo','cDC → MHC-II','cMono → homeo','cMono → MHC-II','ncMono → homeo','ncMono → MHC-II','pDC → homeo','pDC → MHC-II')
m<-aggregate(spearman~target_gene+axis+sample,d[d$assignment=='author_inclusive',],median)
m<-merge(expand.grid(target_gene=levels(d$target_gene),axis=axes,sample=sort(unique(d$sample))),m,by=c('target_gene','axis','sample'),all.x=TRUE);m$target_gene<-factor(m$target_gene,levels=c('CD86','IFNGR1','MERTK'))
m$axis<-factor(m$axis,levels=rev(axes));m$well<-paste0('S',match(m$sample,sort(unique(m$sample))))
b<-ggplot(m,aes(well,axis,fill=spearman))+geom_tile(width=.93,height=.9)+facet_wrap(~target_gene,nrow=1)+scale_y_discrete(labels=setNames(axlabs,axes))+scale_fill_gradient2(low='#2166AC',mid='#F7F7F7',high='#B35806',limits=c(-corlim,corlim),breaks=c(-corlim,0,corlim),name='Spearman ρ',na.value='#BDBDBD')+guides(fill=guide_colourbar(display='rectangles',title.position='top',title.hjust=.5))+labs(x='Capture well; median across ten frozen model seeds',y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='bottom',legend.key.width=unit(10,'mm'))
key<-c('target_gene','axis','sample','seed')
pair<-merge(d[d$assignment=='author_inclusive',c(key,'spearman')],d[d$assignment=='single_positive_guide',c(key,'spearman')],by=key,suffixes=c('_author','_single'),all=TRUE)
pair$blood<-sub('__x__.*','',pair$axis);pair$state<-ifelse(grepl('homeostatic',pair$axis),'Homeostatic','MHC-II')
pair$blood<-factor(pair$blood,levels=c('cDC','classical_mono','nonclassical_mono','pDC'),labels=c('cDC','cMono','ncMono','pDC'))
c<-ggplot(pair,aes(spearman_author,spearman_single,colour=blood,shape=state))+geom_abline(intercept=0,slope=1,colour='#BBBBBB',linewidth=.3/.pt)+geom_vline(xintercept=0,colour='#E0E0E0',linewidth=.3/.pt)+geom_hline(yintercept=0,colour='#E0E0E0',linewidth=.3/.pt)+geom_point(size=.65,alpha=.65,stroke=.25)+facet_wrap(~target_gene,nrow=1)+scale_colour_manual(values=c('#758AA2','#3A8B84','#C48A53','#A66868'))+scale_shape_manual(values=c(16,1))+scale_x_continuous(limits=c(-corlim,corlim),breaks=c(-corlim,0,corlim))+scale_y_continuous(limits=c(-corlim,corlim),breaks=c(-corlim,0,corlim))+coord_equal()+labs(x='Direction correlation: author guide inclusion',y='Direction correlation: single guide',colour='Blood model input',shape='Brain model state')+theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(0,0,0,0),legend.box.spacing=unit(.5,'mm'))
panels<-list(a,b,c);layout<-data.frame(id=c('a','b','c'),x=1,y=c(1,69,143),w=172,h=c(65,71,70),title=c('Measured microglial TSS knockdown responses','Model / experiment directions across all original microglial axes','Sensitivity to co-occurring perturbation guides'))
grobs<-list()
for(i in 1:3){rec<-layout[i,];ragg::agg_capture(width=rec$w,height=rec$h,units='mm',res=300);g<-ggplotGrob(panels[[i]]);dev.off();g<-gTree(children=gList(g,txt(rec$id,.6,.8,8,'bold'),txt(rec$title,4.6,1.2,6.5,'bold')));audit<-text_audit(g);stopifnot(all(audit$pt>=6-1e-6));write.csv(audit,file.path(OUT,'qa',paste0('text_',rec$id,'.csv')),row.names=FALSE);grobs[[i]]<-g}
draw<-function(){
 grid.newpage()
 for(i in 1:3){rec<-layout[i,];pushViewport(viewport(x=unit(rec$x,'mm'),y=unit(225-rec$y,'mm'),width=unit(rec$w,'mm'),height=unit(rec$h,'mm'),just=c('left','top')));grid.draw(grobs[[i]]);popViewport()}
 for(j in 1:3){label<-c('Experimental perturbation: microglial CRISPRi. Model control: blood-cell ADT. Interventions are unmatched.','Seven capture wells share one iPSC clone. Correlations are descriptive; no biological confidence interval.',paste0('Empty target groups remain blank: ',sum(effects$n_cells==0|effects$n_control_cells==0),' / ',nrow(effects),'. No effect-based guide selection.'))[j];grid.text(label,x=unit(4,'mm'),y=unit(9-(j-1)*3.5,'mm'),just='left',gp=gpar(fontfamily='Arial',fontsize=6))}
}
open_pdf(file.path(OUT,'Public_perturbation_direction_comparison.pdf'),174,225);draw();dev.off()
for(dpi in c(600,160)){name<-if(dpi==600)'Public_perturbation_direction_comparison.png' else 'Public_perturbation_direction_comparison_preview.png';ragg::agg_png(file.path(OUT,name),width=174,height=225,units='mm',res=dpi,background='white');draw();dev.off()}
write.csv(effects,file.path(OUT,'plot_data/target_effects.csv'),row.names=FALSE);write.csv(m,file.path(OUT,'plot_data/well_correlations.csv'),row.names=FALSE);write.csv(pair,file.path(OUT,'plot_data/assignment_sensitivity.csv'),row.names=FALSE);write.csv(layout,file.path(OUT,'qa/layout_mm.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'qa/R_sessionInfo.txt'))
write_json(list(page_mm=c(174,225),minimum_font_pt=6,correlation_axis_limit=corlim,panels=3,plot_areas=9,simple_bar_charts=0,target_readout_rows=nrow(effects),target_readout_points=sum(is.finite(effects$target_gene_delta)),undefined_target_groups=sum(!is.finite(effects$target_gene_delta)),correlation_tiles=nrow(m),sensitivity_rows=nrow(pair),sensitivity_points=sum(is.finite(pair$spearman_author)&is.finite(pair$spearman_single)),direction_scope='Pooled TSS guide union; excludes the perturbed gene from correlation; all eight model axes; intervention mismatch explicit',visual_inspection_required=TRUE),file.path(OUT,'qa/FIGURE_AUDIT.json'),pretty=TRUE,auto_unbox=TRUE)
