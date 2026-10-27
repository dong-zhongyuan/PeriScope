umap_grid<-function(){
theme_set(theme_classic(base_size=12,base_family='Helvetica'));lookup<-fromJSON(file.path(ROOT,'pd_product_figure_design/template_semantic_review/inputs/brain_diagnosis_lookup.json'));rd<-function(f)read.csv(file.path(ROOT,'pd_product_figure_design/template_semantic_review',f),stringsAsFactors=FALSE);typo<-theme(text=element_text(family='Helvetica',size=12),legend.text=element_text(size=12),legend.title=element_text(size=12,face='bold'))
# 1B: shared expression coordinates, two annotation rows with actual separators/axes.
d<-rd('data/real_cell_embedding.csv');d$condition<-vapply(d$donor,function(w)lookup[[w]][1],'');d$condition[d$condition=='Ctrl']<-'Control';d$cell_type<-factor(d$cell_type);levels(d$cell_type)<-c(Astrocyte='Astro',B_cell='B cell',DA_neuron='DA neuron',Endothelial='Endothelial',Microglia='Microglia',`Mono-macro`='Mono/macro',Neuron_other='Other neuron',Oligodendrocyte='Oligo',OPC='OPC')[levels(d$cell_type)];d$split<-factor(d$split,c('Train','Held out'))
limits1<-range(d$x);limits2<-range(d$y)
cellp<-function(z,color,legtitle){q<-emb_plot(z,color=color,legend=TRUE)+typo+coord_equal(xlim=limits1,ylim=limits2,clip='off')+scale_color_manual(values=setNames(pal[seq_len(length(unique(d[[color]])))],sort(unique(d[[color]]))),name=legtitle)+theme(plot.margin=margin(1,1,1,1),legend.text=element_text(size=12),legend.position='none',axis.text=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title=element_blank());q$layers[[1]]$aes_params$size<-.7;q}
ps<-list();for(v in c('condition','cell_type'))for(s in levels(d$split))ps[[length(ps)+1]]<-cellp(d[d$split==s,],v,if(v=='condition')'Condition'else'Cell type')
p<-ggdraw()+draw_plot(plot_grid(plotlist=ps,ncol=2),.08,.38,.9,.52)+draw_line(c(.08,.98),c(.64,.64),color='#686868',size=.65)+draw_line(c(.53,.53),c(.38,.90),color='#686868',size=.65)+draw_label('Train',.305,.96,fontface='bold',size=13)+draw_label('Held out',.755,.96,fontface='bold',size=13)+draw_label('Condition',.04,.77,angle=90,fontface='bold',size=13)+draw_label('Cell type',.025,.51,angle=90,fontface='bold',size=13)+draw_line(c(.72,.79),c(.30,.30),arrow=arrow(length=unit(2,'mm')),size=.6)+draw_line(c(.72,.72),c(.30,.36),arrow=arrow(length=unit(2,'mm')),size=.6)+draw_label('UMAP 1',.875,.30,size=12)+draw_label('UMAP 2',.83,.365,size=12)
for(v in c('condition','cell_type')){
 vals<-sort(unique(d[[v]]));cols<-setNames(pal[seq_along(vals)],vals)
 yy<-if(v=='condition').25 else .19
 for(i in seq_along(vals)){x<-.065+((i-1)%%3)*.32;y<-yy-floor((i-1)/3)*.063;p<-p+draw_grob(grid::pointsGrob(x=.5,y=.5,pch=16,size=unit(1.7,'mm'),gp=gpar(col=cols[[as.character(vals[i])]])),x-.01,y-.01,.02,.02)+draw_label(as.character(vals[i]),x+.025,y,hjust=0,size=12)}
}
p
}
