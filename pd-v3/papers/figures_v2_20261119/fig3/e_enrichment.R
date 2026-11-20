# Figure 3e - enrichment of two discovery queries. Left: the 21 overlap genes of the CD22 cDC -> MG-AP down query; dots
# mark membership in protein modification (grey) and phosphorylation (colour) - all 37 term-gene links. Right: six
# leading GO biological-process terms of the CD40 cDC -> Astro up query; x = overlap genes, colour = -log10 q.
# Data: discovery_go.csv, cd40_functional_go.csv (read-only).
PANEL<-list(id='3e_enrichment',blocks=c('1x2'),what='Enrichment of two discovery queries')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
make_panel<-function(b){W<-b$w;H<-b$h;assign('PANEL_H',H,envir=globalenv())
 go<-read_src(file.path(FRD,'discovery_go.csv'),cols=c('term','genes'));stopifnot(nrow(go)==2)
 s<-lapply(strsplit(go$genes,';'),sort);gs<-sort(unique(unlist(s)));stopifnot(sum(lengths(s))==37)
 tc<-c(MID,PERISCOPE);g<-list();nc<-3;nr<-ceiling(length(gs)/nc);top<-LETTER_ZONE[['h']]+2.2;rh<-2.95;cw<-14.6
 for(i in seq_along(gs)){cc<-(i-1)%/%nr;rr<-(i-1)%%nr;x<-.6+cc*cw;y<-top+rr*rh
  for(k in 1:2)if(gs[i]%in%s[[k]])g<-c(g,list(circleGrob(unit(x+.8+(k-1)*1.9,'mm'),unit(H-y,'mm'),r=unit(.75,'mm'),gp=gpar(fill=tc[k],col=NA))))
  g<-c(g,list(dtxt(gs[i],x+4.6,y,just=c('left','centre'),italic=TRUE)))}
 ky<-top+nr*rh+1.4
 g<-c(g,key_row(list(ki(kg_dot(tc[1],NA,1.5),'Modification'),ki(kg_dot(tc[2],NA,1.5),'Phosphorylation')),.6,ky,h=H))
 cd<-read_src(file.path(FRD,'cd40_functional_go.csv'),cols=c('term','q','overlap'))
 cd$t<-sub(' \\(GO:[0-9]+\\)$','',cd$term);cd$t<-paste0(toupper(substr(cd$t,1,1)),tolower(substring(cd$t,2)))
 sh<-c('Modulation of chemical synaptic transmission'='Synaptic modulation','Chemical synaptic transmission'='Synaptic transmission','Anterograde trans-synaptic signaling'='Anterograde signaling','Positive regulation of synaptic transmission'='Synaptic upregulation','Synapse organization'='Synapse organization','Regulation of trans-synaptic signaling'='Trans-synaptic regulation');stopifnot(all(cd$t%in%names(sh)));cd$t<-unname(sh[cd$t]);cd<-cd[order(cd$q,decreasing=TRUE),];cd$y<-seq_len(nrow(cd));cd$lq<- -log10(cd$q)
 p<-ggplot(cd,aes(overlap,y))+geom_segment(aes(x=0,xend=overlap,yend=y),colour=FAINT,linewidth=gglw('emph'))+
  geom_point(aes(fill=lq),shape=21,colour='white',stroke=.2,size=ggsz('big'))+
  scale_fill_gradientn(colours=seq_col(seq(0,1,length.out=60),PERISCOPE),limits=c(0,max(cd$lq)),guide='none')+
  scale_y_continuous(breaks=cd$y,labels=cd$t,expand=expansion(add=.6))+scale_x_continuous(limits=c(0,max(cd$overlap)*1.08),breaks=c(0,10,20),expand=expansion(0))+
  labs(x='Overlap',y=NULL)+theme_fig()+theme(axis.text.y=element_text(lineheight=.9),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(LETTER_ZONE[['h']],.6,.6,.6,'mm'))
 x0<-45.5
 kk<-list(key_title('−log10 q',x0+1,H-1.6,H));kk<-c(kk,key_row(list(ki(kg_grad(seq_col(seq(0,1,length.out=60),PERISCOPE),'0',formatC(floor(max(cd$lq)),format='d'),w=10))),x0+1+sw('−log10 q',TRUE)+1.2,H-1.6,h=H))
 gTree(children=gList(dgroup(g),grobTree(gg_grob(p),vp=viewport(x=unit(x0,'mm'),y=unit(0,'mm'),width=unit(W-x0,'mm'),height=unit(H,'mm'),just=c(0,0))),do.call(gList,kk)))}
key_title<-function(t,x,y,h)textGrob(t,unit(x,'mm'),unit(h-y,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,fontface='bold',col=INK))
