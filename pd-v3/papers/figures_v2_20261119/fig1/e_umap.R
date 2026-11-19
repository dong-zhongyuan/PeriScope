# Figure 1e - brain RNA embedding: training cells as a light grey base, held-out cells coloured on top; left by
# diagnosis, right by cell type with its key. Data: template_semantic_review/data/real_cell_embedding.csv (read-only).
PANEL<-list(id='1e_brain_umap',blocks=c('1x2'),what='Brain RNA embedding, training and held-out cells')
T0<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design/template_semantic_review'
CT_COL<-c(Astro=unname(BRAIN_COL['Astro']),Microglia=unname(BRAIN_COL['MG-AP']),`Mono/macro`='#D98C6B',Oligo='#7E8C9E',OPC='#A9B6C4',
 `DA neuron`='#5B7F5A',`Other neuron`='#9DB59A',Endothelial='#8A6E9E',Pericyte='#C2AECF',`B cell`='#BFA77A',`T cell`='#D9C7A0')
make_panel<-function(b){
 d<-read_src(file.path(T0,'data/real_cell_embedding.csv'),cols=c('donor','condition','cell_type','split','x','y'))
 lk<-jsonlite::fromJSON(src_log(file.path(T0,'inputs/brain_diagnosis_lookup.json')))
 d$dx<-vapply(d$donor,function(w)lk[[w]][1],'');d$dx[d$dx=='Ctrl']<-'Control'
 map<-c(Astrocyte='Astro',B_cell='B cell',DA_neuron='DA neuron',Endothelial='Endothelial',Microglia='Microglia',Mono_macro='Mono/macro',T_cell='T cell',
  Neuron_other='Other neuron',Oligodendrocyte='Oligo',OPC='OPC',Pericyte='Pericyte')
 d$ct<-map[d$cell_type];stopifnot(!anyNA(d$ct))
 tr<-d[d$split=='Train',];ho<-d[d$split!='Train',];set.seed(1);ho<-ho[sample(nrow(ho)),]
 xl<-range(d$x);yl<-range(d$y)
 base<-function()ggplot()+geom_point(data=tr,aes(x,y),shape=16,colour=FAINT,size=ggsz('point'))+
  scale_x_continuous(limits=xl,expand=expansion(mult=.03))+scale_y_continuous(limits=yl,expand=expansion(mult=.03))+theme_fig(legend='none')+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),
   plot.margin=margin(0,0,0,0))
 dxs<-intersect(names(DX_COL),unique(ho$dx))
 p1<-base()+geom_point(data=ho,aes(x,y,colour=dx),shape=16,size=ggsz('small'))+scale_colour_manual(values=DX_COL[dxs],breaks=dxs)+
  theme(legend.position='bottom',legend.key.size=unit(2.2,'mm'),legend.key.spacing.x=unit(1,'mm'),legend.box.spacing=unit(.8,'mm'))+
  guides(colour=guide_legend(nrow=1,override.aes=list(size=ggsz('big'))))
 cts<-names(CT_COL)[names(CT_COL)%in%ho$ct|names(CT_COL)%in%tr$ct]
 p2<-base()+geom_point(data=ho,aes(x,y,colour=ct),shape=16,size=ggsz('small'))+scale_colour_manual(values=CT_COL[cts],breaks=cts)+
  theme(legend.position='right',legend.key.size=unit(2.2,'mm'),legend.key.spacing.y=unit(.35,'mm'),legend.box.spacing=unit(1.2,'mm'))+
  guides(colour=guide_legend(ncol=1,override.aes=list(size=ggsz('big'))))
 g1<-gg_grob(p1);g2<-gg_grob(p2);W<-b$w;H<-b$h;top<-LETTER_ZONE[['h']]+.4;w1<-(W-3)*.42
 gTree(children=gList(
  grobTree(g1,vp=viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(w1,'mm'),height=unit(H-top,'mm'),just=c(0,0))),
  grobTree(g2,vp=viewport(x=unit(w1+3,'mm'),y=unit(0,'mm'),width=unit(W-w1-3,'mm'),height=unit(H,'mm'),just=c(0,0)))))}
