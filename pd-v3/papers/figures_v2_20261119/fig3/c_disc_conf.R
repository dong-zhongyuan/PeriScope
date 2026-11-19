# Figure 3c - discovery vs confirmation gene slopes (x 10^3) for the CD22 cDC -> MG-AP down program, all 4,741 genes;
# program members coloured, the four most extreme members named. Data: gene_responses.csv.
PANEL<-list(id='3c_discovery_confirmation',blocks=c('1x1'),what='Discovery vs confirmation gene slopes')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
make_panel<-function(b){
 gr<-read_src(file.path(FRD,'gene_responses.csv'),cols=c('combo','gene','discovery','confirmation','member'))
 d<-gr[gr$combo=='raw::CD22__cDC__x__microglia_mhc2__down',];d$x<-d$discovery*1e3;d$y<-d$confirmation*1e3
 m<-d[d$member,];lab<-m[order(m$x),][1:3,];yr<-diff(range(d$y));o<-order(lab$y);ly<-lab$y[o];for(i in 2:3)ly[i]<-max(ly[i],ly[i-1]+.12*yr);lab$ly[o]<-ly;lab$lx<-min(d$x)+.02*diff(range(d$x))
 ggplot(d,aes(x,y))+
  geom_hline(yintercept=0,colour=ZERO_COL,linewidth=gglw('thin'))+geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('thin'))+
  geom_point(data=d[!d$member,],shape=16,colour=alpha_hex(MID,.25),size=ggsz('unit'))+
  geom_point(data=m,shape=21,fill=PERISCOPE,colour='white',stroke=.15,size=ggsz('small')*1.2)+
  geom_segment(data=lab,aes(x=x,y=y,xend=lx+.12*diff(range(d$x)),yend=ly-.07*yr),colour=GUIDE,linewidth=gglw('hair'))+
  geom_text(data=lab,aes(lx,ly,label=gene),family=FONT,size=GT,hjust=0,vjust=0,colour=INK,fontface='italic')+
  scale_x_continuous(labels=num_minus,n.breaks=4)+scale_y_continuous(labels=num_minus,n.breaks=4)+
  labs(x='1,000 \u00d7 discovery slope',y='1,000 \u00d7 confirmation')+theme_fig()}
