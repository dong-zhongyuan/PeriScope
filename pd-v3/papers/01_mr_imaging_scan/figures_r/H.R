make_H <- function() {
  d<-read_pd('H_palindromes');d$retained<-tolower(as.character(d$retained))=='true'
  stopifnot(nrow(d)==269,sum(d$retained)==227)
  d$x<-(d$tool_order-1)%%27+1;d$y<-10-floor((d$tool_order-1)/27)
  tiles<-rbind(data.frame(x=d$x,y0=d$y+.02,y1=d$y+.43,maf=d$maf_protein),
               data.frame(x=d$x,y0=d$y-.43,y1=d$y-.02,maf=d$maf_imaging))
  ggplot(tiles)+geom_rect(aes(xmin=x-.43,xmax=x+.43,ymin=y0,ymax=y1,fill=maf),colour=NA)+
    geom_rect(data=d[!d$retained,],aes(xmin=x-.47,xmax=x+.47,ymin=y-.47,ymax=y+.47),
              inherit.aes=FALSE,fill=NA,colour=INK,linewidth=.3/.pt)+
    scale_fill_gradient(low='white',high='#78549D',limits=c(0,.5),breaks=c(0,.25,.5),name='Minor-allele frequency',
                        guide=guide_colourbar(display='rectangles',direction='horizontal',barwidth=unit(28,'mm'),barheight=unit(2,'mm')))+
    scale_x_continuous(limits=c(.5,27.5),expand=c(0,0))+
    scale_y_continuous(limits=c(.5,10.5),expand=c(0,0))+
    labs(x=NULL,y=NULL)+
    theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),legend.position='bottom',
          legend.title.position='top',legend.justification='left',legend.margin=margin(0,0,0,2.5,'mm'))
}
