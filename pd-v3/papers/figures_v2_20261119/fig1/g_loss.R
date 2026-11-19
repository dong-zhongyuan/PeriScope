# Figure 1g - protein-head training loss for full-model seeds 42-51: one thin line per seed, seed mean thick,
# open circle at each seed's stopping step. Data: figure_redesign_20261004/data/training.csv (read-only).
PANEL<-list(id='1g_training_loss',blocks=c('1x1'),what='Protein-head training loss, 10 seeds')
FRD<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design/manuscript_review/figure_redesign_20261004/data'
make_panel<-function(b){
 d<-read_src(file.path(FRD,'training.csv'),cols=c('step','protein','seed'))
 d$k<-d$step/1000
 end<-do.call(rbind,lapply(split(d,d$seed),function(x)x[which.max(x$step),]))
 m<-aggregate(protein~step,d,mean);m<-m[m$step<=min(end$step),];m$k<-m$step/1000
 ggplot(d,aes(k,protein,group=seed))+
  geom_line(colour=PERISCOPE_LIGHT,linewidth=gglw('thin'))+
  geom_line(data=m,aes(k,protein),inherit.aes=FALSE,colour=PERISCOPE,linewidth=gglw('mean'))+
  geom_point(data=end,aes(k,protein),shape=21,fill='white',colour=PERISCOPE,size=ggsz('small'),stroke=gglw('thin')/.25*.25)+
  scale_y_log10(breaks=c(.3,.5,1,2),labels=num_minus,expand=expansion(mult=c(.02,.04)))+
  scale_x_continuous(breaks=c(0,10,20,30,40),expand=expansion(mult=c(.01,.03)))+
  labs(x='Steps (×1,000)',y='Protein loss')+theme_fig()}
