# Figure 2c - condition-aware Pearson r (fill) and response-sign agreement (circle size, % of genes, scaled from 40%) for
# every method and route, rows aligned with 2b. All methods seeds 42-44 (manuscript runs). Data: human_direction_metrics.csv.
PANEL<-list(id='2c_pearson_sign',blocks=c('2x1'),what='Pearson and sign agreement by method and route')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/_bench.R'))
make_panel<-function(b){r<-bench_long('pearson');s<-bench_long('sign_agreement')
 t<-aggregate(value~method+row+col,r,mean);u<-aggregate(value~method+row+col,s,mean);names(u)[4]<-'sign'
 t<-merge(t,u);mr<-aggregate(value~method+row,t,mean);ms<-aggregate(sign~method+row,t,mean);mm<-merge(mr,ms);mm$col<-MEANX
 t<-rbind(t[,c('row','col','value','sign')],mm[,c('row','col','value','sign')]);t$sign<-t$sign*100
 lim<-c(-.3,.5);dmm<-function(v).4+(v-40)/20*3.2
 p<-ggplot(t,aes(col,row))+
  geom_tile(aes(fill=value),width=.92,height=.9,colour=NA)+
  geom_point(aes(size=sign),shape=21,fill=NA,colour=INK,stroke=.3)+
  col_header(names=FALSE)+
  scale_fill_gradientn(colours=div_cols(lim),limits=lim,oob=scales::squish)+
  scale_size_identity()+
  scale_x_continuous(expand=expansion(add=c(.5,.5)))+
  scale_y_continuous(breaks=seq_along(M_ORDER),labels=NULL,expand=expansion(add=Y_EXP))+
  bench_coord()+labs(x=NULL,y=NULL)+bench_theme()
 p$layers[[2]]$mapping<-aes(size=ggsz(dmm(sign)))
 W<-b$w;H<-b$h;y1<-H-KEY_BOT+2.2;y2<-H-1.8;x<-.6
 circ<-function(v){d<-dmm(v);list(w=d,f=function(x,yc)circleGrob(unit(x+d/2,'mm'),unit(yc,'mm'),r=unit(d/2,'mm'),gp=gpar(fill=NA,col=INK,lwd=glwd('thin'))))}
 xx<-x+sw('Sign (%)',TRUE)+1.2
 k<-list(key_title('Pearson',x,y1,H));k<-c(k,key_row(list(ki(kg_grad(div_cols(lim),num_minus_acc(.1)(lim[1]),num_minus_acc(.1)(lim[2]),w=10))),xx,y1,h=H))
 k<-c(k,list(key_title('Sign (%)',x,y2,H)),key_row(list(ki(circ(45),'45'),ki(circ(55),'55')),xx,y2,h=H))
 gTree(children=gList(fix_null_rows(gg_grob(p)),do.call(gList,k)))}
