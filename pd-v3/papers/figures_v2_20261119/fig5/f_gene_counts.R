# Figure 5f - size of the ten final-target programs: member genes per program (filled) and genes found in no other of the
# ten programs (open); rows grouped by target. Data: optimization_20261003/program_gene_rankings.json (first n_genes ranked genes = members), exported to fig5/data/final_program_genes.csv.
PANEL<-list(id='5f_gene_counts',blocks=c('1x2'),what='Gene counts of the final-target programs')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig5/_fun.R'))
make_panel<-function(b){fp<-final_prog_genes();G<-strsplit(fp$genes,';');tab<-table(unlist(G))
 fp$uniq<-vapply(G,function(z)sum(tab[z]==1),0);fp$y<-rev(seq_len(nrow(fp)));fp$col<-TARGET_COL[as.character(fp$antigen)];fp$lab<-paste(fp$antigen,fp$lab)
 ggplot(fp,aes(y=y))+geom_segment(aes(x=0,xend=n_genes,yend=y),colour=FAINT,linewidth=gglw('emph'))+
  geom_point(aes(x=n_genes),shape=21,fill=fp$col,colour='white',stroke=.3,size=ggsz('big'))+
  geom_point(aes(x=uniq),shape=21,fill='white',colour=fp$col,stroke=.5,size=ggsz('point'))+
  scale_y_continuous(breaks=fp$y,labels=fp$lab,expand=expansion(add=.6))+scale_x_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.05)))+
  labs(x='Genes',y=NULL)+theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(colour=unname(fp$col)))}
