# Figure 4c - MR atlas across the genome: best -log10 P over all imaging phenotypes for each of the 2,941 UKB-PPP assays,
# at the cis-instrument position; alternating chromosome shades; final targets coloured. Data: A_all_assays.csv (drop policy).
PANEL<-list(id='4c_atlas',blocks=c('1x2'),what='Genome-wide MR atlas, best P per assay')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 A<-read_src(file.path(OV,'A_all_assays.csv'),cols=c('gene','genome_order','chrom','minus_log10_p','policy'));A<-A[A$policy=='drop'&is.finite(A$minus_log10_p),]
 ch<-read_src(file.path(OV,'A_chromosomes.csv'),cols=c('chrom','min','max','mean','label'))
 A$shade<-ifelse(as.numeric(A$chrom)%%2==0,'#9A9A9A','#C9C9C9');ft<-A[A$gene%in%names(FT_COL),]
 lb<-ch[ch$label%in%c('1','5','10','15','20','X','23'),]
 ggplot(A,aes(genome_order,minus_log10_p))+
  geom_point(colour=A$shade,shape=16,size=ggsz('small'))+
  geom_point(data=ft,shape=21,fill=FT_COL[ft$gene],colour='white',stroke=.3,size=ggsz('big'))+
  geom_text(data=ft,aes(label=gene),family=FONT,size=GT,fontface='italic',hjust=-.25,colour=FT_COL[ft$gene])+
  scale_x_continuous(breaks=lb$mean,labels=lb$label,expand=expansion(mult=.01))+scale_y_continuous(expand=expansion(mult=c(0,.05)))+
  labs(x='Chromosome',y='−log10 P')+theme_fig()}
