make_G <- function() {
  d<-read_pd('G_delta_se');stopifnot(nrow(d)==68175,all(d$se_increase_percent>=0))
  ggplot(d,aes(se_increase_percent,colour=modality))+stat_ecdf(geom='step',linewidth=.5/.pt)+
    scale_colour_manual(values=MOD,name=NULL)+
    scale_x_continuous(trans=scales::pseudo_log_trans(sigma=.1),breaks=c(0,.1,1,10),labels=num)+
    scale_y_continuous(breaks=c(0,.5,1),labels=num,limits=c(0,1),expand=c(.01,0))+
    labs(x='Delta SE increase (%)',y='Cumulative fraction')+
    theme(legend.position='top',legend.justification='left',legend.direction='horizontal')
}
