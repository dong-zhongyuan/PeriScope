# Figure 5g - gene membership across the ten final-target programs: one column per gene of their union (sorted by the
# number of programs containing it), one row per program; top strip = programs per gene.
# Data: program members from optimization_20261003/program_gene_rankings.json (fig5/data/final_program_genes.csv).
PANEL<-list(id='5g_gene_upset',blocks=c('2x2'),what='Gene membership barcode of final-target programs')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig5/_fun.R'))
make_panel<-function(b){fp<-final_prog_genes();G<-strsplit(fp$genes,';');names(G)<-fp$combo;U<-unique(unlist(G))
 M<-sapply(G,function(z)U%in%z);cnt<-rowSums(M);first<-apply(M,1,function(r)which(r)[1]);o<-order(-cnt,first,U);U<-U[o];M<-M[o,,drop=FALSE];cnt<-cnt[o]
 L<-do.call(rbind,lapply(seq_len(ncol(M)),function(j)data.frame(x=which(M[,j]),y=ncol(M)-j+1,col=TARGET_COL[as.character(fp$antigen[j])])))
 top<-data.frame(x=seq_along(U),n=cnt);ny<-ncol(M)
 p<-ggplot()+geom_tile(data=L,aes(x,y),fill=L$col,width=1,height=.8)+
  geom_segment(data=top,aes(x=x,xend=x,y=ny+.8,yend=ny+.8+n*.45),colour=MID,linewidth=gglw('hair')*1.6)+
  scale_y_continuous(breaks=ny:1,labels=paste(fp$antigen,fp$lab),expand=expansion(add=c(.5,.3)))+scale_x_continuous(expand=expansion(add=.5),breaks=NULL)+
  coord_cartesian(clip='off')+labs(x='Genes (union)',y=NULL)+theme_fig()+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(colour=unname(TARGET_COL[as.character(fp$antigen)])))
 gg_grob(p)}
