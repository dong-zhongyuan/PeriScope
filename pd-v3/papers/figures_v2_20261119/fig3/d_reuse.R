# Figure 3d - how many programs contain each gene (2,503 distinct genes across the program library): genes binned on a
# log axis, stems + dots = gene counts per bin. Data: gene_reuse.csv.
PANEL<-list(id='3d_gene_multiplicity',blocks=c('1x1'),what='Gene multiplicity across programs')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
make_panel<-function(b){
 r<-read_src(file.path(FRD,'gene_reuse.csv'),cols=c('gene','n_programs'))
 br<-10^seq(0,log10(max(r$n_programs))+.05,length.out=19);h<-graphics::hist(r$n_programs,breaks=c(.999,br[-1]),plot=FALSE)
 d<-data.frame(x=sqrt(head(c(.999,br[-1]),-1)*br[-1]),n=h$counts);d<-d[d$n>0,]
 ggplot(d,aes(x,n))+geom_segment(aes(xend=x,y=0,yend=n),colour=LIGHT,linewidth=gglw('emph'))+
  geom_point(shape=21,fill=PERISCOPE,colour='white',stroke=.2,size=ggsz('point'))+
  scale_x_log10(breaks=c(1,10,100,1000),labels=num_comma)+scale_y_continuous(labels=num_comma,expand=expansion(mult=c(0,.06)))+
  labs(x='Programs per gene',y='Genes')+theme_fig()}
