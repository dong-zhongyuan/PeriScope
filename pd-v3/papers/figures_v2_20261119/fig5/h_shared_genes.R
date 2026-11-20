# Figure 5h - genes shared most widely across the ten final-target programs (in at least four): gene x program membership;
# rows = programs grouped by target. Data: program members from optimization_20261003/program_gene_rankings.json (fig5/data/final_program_genes.csv).
PANEL<-list(id='5h_shared_genes',blocks=c('1x2'),what='Most widely shared genes of final-target programs')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig5/_fun.R'))
make_panel<-function(b){fp<-final_prog_genes();G<-strsplit(fp$genes,';');U<-unique(unlist(G))
 M<-sapply(G,function(z)U%in%z);cnt<-rowSums(M);k<-4;while(sum(cnt>=k)>26)k<-k+1;keep<-cnt>=k;U<-U[keep];M<-M[keep,,drop=FALSE];cnt<-cnt[keep]
 o<-order(-cnt,U);U<-U[o];M<-M[o,,drop=FALSE];ny<-ncol(M)
 L<-do.call(rbind,lapply(seq_len(ny),function(j)data.frame(x=seq_along(U),y=ny-j+1,m=M[,j],col=TARGET_COL[as.character(fp$antigen[j])])))
 ggplot(L,aes(x,y))+geom_tile(fill=ifelse(L$m,L$col,FAINT),width=.9,height=.85)+
  scale_x_continuous(breaks=seq_along(U),labels=U,position='top',expand=expansion(add=.5))+
  scale_y_continuous(breaks=ny:1,labels=paste(fp$antigen,fp$lab),expand=expansion(add=.5))+
  labs(x=NULL,y=NULL)+theme_fig()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x.top=element_text(angle=90,hjust=0,vjust=.5,face='italic'),
   axis.text.y=element_text(colour=unname(TARGET_COL[as.character(fp$antigen)])))}
