# Figure 3j - CD71 -> MG-AP program scores in GSE157783 donors (six control, five PD), one group per blood source;
# points = donors, bar = mean, rail = 95% interval of the mean. Data: matched_independent_donor_scores.csv.
PANEL<-list(id='3j_donor_scores',blocks=c('1x1'),what='CD71 program donor scores in GSE157783')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
make_panel<-function(b){
 d<-read_src(file.path(FRD,'matched_independent_donor_scores.csv'),cols=c('combo','donor','condition','score'));d<-d[grepl('^raw::CD71__',d$combo),]
 pp<-prog_parse(unique(d$combo));d$blood<-pp$blood[match(d$combo,pp$combo)];d$g<-match(d$blood,B_LAB)
 d$condition[d$condition=='Ctrl']<-'Control';d$x<-d$g+ifelse(d$condition=='PD',.2,-.2)
 set.seed(7);d$xj<-d$x+stats::runif(nrow(d),-.06,.06)
 s<-do.call(rbind,lapply(split(d,list(d$g,d$condition),drop=TRUE),function(z){m<-mean(z$score);se<-stats::sd(z$score)/sqrt(nrow(z));tq<-stats::qt(.975,nrow(z)-1)
  data.frame(x=z$x[1],condition=z$condition[1],m=m,lo=m-tq*se,hi=m+tq*se)}))
 ggplot()+
  geom_hline(yintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_linerange(data=s,aes(x=x+.1,ymin=lo,ymax=hi,colour=condition),linewidth=gglw('emph'))+
  geom_segment(data=s,aes(x=x-.12,xend=x+.12,y=m,yend=m,colour=condition),linewidth=gglw('mean'))+
  geom_point(data=d,aes(xj,score,fill=condition),shape=21,colour='white',stroke=.2,size=ggsz('point'))+
  scale_colour_manual(values=DX_COL[c('Control','PD')],guide='none')+scale_fill_manual(values=DX_COL[c('Control','PD')],guide='none')+
  scale_x_continuous(breaks=sort(unique(d$g)),labels=B_LAB[sort(unique(d$g))],expand=expansion(add=.4))+scale_y_continuous(labels=num_minus,n.breaks=4)+
  labs(x=NULL,y='Program score')+theme_fig()+theme(axis.ticks.x=element_blank(),axis.text.x=element_text(angle=35,hjust=1,vjust=1))}
