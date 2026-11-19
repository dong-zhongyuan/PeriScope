# Figure 4e - BIG40 replication of MR hit pairs. Left: discovery z vs replication z signed to the discovery direction
# (hit pairs coloured; matched non-hit pairs grey, random subset). Right: share of pairs replicating in the same
# direction with 95% Wilson intervals, by hit definition, for colocalised regions and for the two null sets.
# Data: big40_repro_20261119/repro_pairs.csv.gz, repro_summary.csv (read-only).
PANEL<-list(id='4e_big40',blocks=c('1x2'),what='BIG40 replication of MR hits')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 p<-read_src(file.path(MR,'big40_repro_20261119/repro_pairs.csv.gz'),cols=c('pair_type','analysable','z_disc_used','z_rep_signed'))
 p<-p[p$analysable%in%c(TRUE,'True')&is.finite(p$z_rep_signed),];h<-p[p$pair_type=='hit',];set.seed(9);n0<-p[p$pair_type=='null_matched',];n0<-n0[sample(nrow(n0),min(2500,nrow(n0))),]
 sc<-ggplot()+geom_hline(yintercept=0,colour=ZERO_COL,linewidth=gglw('thin'))+
  geom_point(data=n0,aes(abs(z_disc_used),z_rep_signed),shape=16,colour=alpha_hex(LIGHT,.5),size=ggsz('unit'))+
  geom_point(data=h,aes(abs(z_disc_used),z_rep_signed),shape=16,colour=alpha_hex(PERISCOPE,.45),size=ggsz('small'))+
  scale_x_continuous(limits=c(2,10.5),breaks=c(4,6,8,10))+scale_y_continuous(labels=num_minus)+
  labs(x='Discovery |z|',y='Replication z')+theme_fig()
 s<-read_src(file.path(MR,'big40_repro_20261119/repro_summary.csv'),cols=c('group','scope','n_analysable','same_sign_n'));s<-s[s$scope=='overall',]
 g<-c(`hit__drop:q_by_cis_family`='Strictest',`hit__drop:q_bh_global`='BH global',hit_ANY_definition='All hits',`coloc_H4ge0.8__instrument`='Coloc regions',null_matched='Matched null',null_random_unselected='Random pairs')
 s<-s[match(names(g),s$group),];stopifnot(!anyNA(s$group));s$lab<-g[s$group];s$y<-rev(seq_len(nrow(s)))
 w<-function(k,n,z=1.96){p<-k/n;d<-1+z^2/n;c((p+z^2/(2*n)-z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d,(p+z^2/(2*n)+z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d)}
 ci<-t(mapply(w,s$same_sign_n,s$n_analysable));s$f<-s$same_sign_n/s$n_analysable*100;s$lo<-ci[,1]*100;s$hi<-ci[,2]*100
 s$col<-ifelse(grepl('null',s$group),MID,PERISCOPE)
 la<-ggplot(s,aes(f,y))+geom_vline(xintercept=50,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_segment(aes(x=lo,xend=hi,yend=y),colour=s$col,linewidth=gglw('emph'))+
  geom_point(shape=21,fill=s$col,colour='white',stroke=.3,size=ggsz('big'))+
  scale_y_continuous(breaks=s$y,labels=s$lab,expand=expansion(add=.5))+scale_x_continuous(limits=c(45,100),breaks=c(50,75,100))+
  labs(x='Same direction (%)',y=NULL)+theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(LETTER_ZONE[['h']],2.4,.6,.6,'mm'))
 W<-b$w;H<-b$h;w1<-38
 gTree(children=gList(grobTree(gg_grob(sc),vp=viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(w1,'mm'),height=unit(H,'mm'),just=c(0,0))),
  grobTree(gg_grob(la),vp=viewport(x=unit(w1+1,'mm'),y=unit(0,'mm'),width=unit(W-w1-1,'mm'),height=unit(H,'mm'),just=c(0,0)))))}
