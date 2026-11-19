# Figure 4j - cis-pQTL replication of all MR instruments: discovery vs replication effect on plasma protein (UKB-PPP
# discovery and replication sets); final-target instruments coloured. Data: exposure_replication_20261116/all_selected_instruments.csv.gz.
PANEL<-list(id='4j_pqtl_replication',blocks=c('1x1'),what='cis-pQTL replication of instruments')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 d<-read_src(file.path(MR,'exposure_replication_20261116/all_selected_instruments.csv.gz'),cols=c('gene','beta_x','beta_x_replication'))
 d<-d[is.finite(d$beta_x)&is.finite(d$beta_x_replication),];ft<-d[d$gene%in%names(FT_COL),]
 ggplot(d,aes(beta_x,beta_x_replication))+geom_abline(slope=1,intercept=0,colour=GUIDE,linewidth=gglw('thin'),linetype=REF_LTY)+
  geom_point(shape=16,colour=alpha_hex(MID,.35),size=ggsz('small'))+
  geom_point(data=ft,shape=21,fill=FT_COL[ft$gene],colour='white',stroke=.3,size=ggsz('big'))+
  scale_x_continuous(labels=num_minus,n.breaks=4)+scale_y_continuous(labels=num_minus,n.breaks=4)+
  labs(x='Discovery beta',y='Replication beta')+theme_fig()}
