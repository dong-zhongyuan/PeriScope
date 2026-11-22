# Figure 2b - human disease-direction leaderboard: methods x 8 blood-to-brain routes plus the equal-route mean. Tile fill
# = seed-mean Spearman; dots in the tile foot = individual seeds (position = value within the tile). All methods seeds 42-44.
PANEL<-list(id='2b_leaderboard',blocks=c('2x2'),what='Human disease-direction leaderboard')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/_bench.R'))
make_panel<-function(b){d<-bench_long('spearman');d<-d[d$seed%in%42:44,];stopifnot(length(unique(d$seed[d$method=='periscope']))==3)
 t<-aggregate(value~method+label+row+col+blood+brain,d,mean)
 sm<-aggregate(value~method+label+row+seed,d,mean);m<-aggregate(value~method+label+row,sm,mean);m$col<-MEANX
 sd_<-sm;sd_$col<-MEANX
 lim<-c(-.3,.45);pos<-function(v,x)x-.38+(pmin(pmax(v,lim[1]),lim[2])-lim[1])/diff(lim)*.76
 seeds<-rbind(data.frame(row=d$row,x=pos(d$value,d$col)),data.frame(row=sd_$row,x=pos(sd_$value,sd_$col)))
 P<-m$value[m$method=='periscope'];m$gain<-ifelse(m$method=='periscope',NA,(P/m$value-1)*100)
 tiles<-rbind(t[,c('row','col','value')],m[,c('row','col','value')])
 p<-ggplot()+
  geom_tile(data=tiles,aes(col,row,fill=value),width=.92,height=.9,colour=NA)+
  geom_point(data=seeds,aes(x,row),shape=16,size=ggsz(.8),colour=alpha_hex(INK,.6))+
  col_header()+
  scale_fill_gradientn(colours=div_cols(lim),limits=lim,oob=scales::squish)+
  scale_colour_manual(values=c(`TRUE`='white',`FALSE`=INK),guide='none')+
  scale_x_continuous(expand=expansion(add=c(.5,.5)))+
  scale_y_continuous(breaks=seq_along(M_ORDER),labels=rev(METHOD_LAB[M_ORDER]),expand=expansion(add=Y_EXP))+
  bench_coord()+labs(x=NULL,y=NULL)+bench_theme()
 W<-b$w;H<-b$h;y1<-H-KEY_BOT/2;x<-24
 k<-list(key_title('Spearman',x,y1,H));x<-x+sw('Spearman',TRUE)+1.2
 k<-c(k,key_row(list(ki(kg_grad(div_cols(lim),num_minus_acc(.1)(lim[1]),num_minus_acc(.01)(lim[2]),w=16))),x,y1,h=H))
 k<-c(k,key_row(list(ki(kg_dot(fill=alpha_hex(INK,.6),col=NA,d=.8),'Seed')),x+31,y1,h=H))
 gTree(children=gList(fix_null_rows(gg_grob(p)),do.call(gList,k)))}
