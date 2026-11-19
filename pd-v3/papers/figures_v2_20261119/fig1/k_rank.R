# Figure 1k - per-antigen agreement between predicted and measured mean ADT abundance across the 13 held-out
# donor x subtype groups: antigens ranked by Spearman, median dashed, final targets marked.
# Data: manuscript_revision_20261017/diag/protein_across_group_accuracy.csv (read-only).
PANEL<-list(id='1k_antigen_consistency',blocks=c('1x1'),what='Per-antigen across-group agreement')
make_panel<-function(b){
 d<-read_src('/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/diag/protein_across_group_accuracy.csv',cols=c('target','groups','spearman_across_groups'))
 d<-d[order(-d$spearman_across_groups),];d$rank<-seq_len(nrow(d));med<-stats::median(d$spearman_across_groups)
 t<-d[d$target%in%names(TARGET_COL),];t$col<-TARGET_COL[t$target]
 ggplot(d,aes(rank,spearman_across_groups))+
  geom_hline(yintercept=med,colour=REF_COL,linewidth=gglw('ref'),linetype=REF_LTY)+
  geom_hline(yintercept=0,colour=ZERO_COL,linewidth=gglw('thin'))+
  geom_point(shape=16,colour=LIGHT,size=ggsz('small'))+
  geom_point(data=t,shape=21,fill=t$col,colour='white',size=ggsz('big'),stroke=.3)+
  geom_text(data=t,aes(x=rank+ifelse(target=='CD71',-10,14),y=spearman_across_groups+ifelse(target=='CD71',-.13,-.07),label=target),family=FONT,size=GT,hjust=ifelse(t$target=='CD71',1,0),colour=t$col)+
  scale_y_continuous(limits=c(-.4,1),breaks=c(-.25,0,.25,.5,.75,1),labels=num_minus,expand=expansion(mult=c(.01,.02)))+
  scale_x_continuous(breaks=c(1,100,200),expand=expansion(mult=c(.02,.03)))+
  labs(x='Antigen rank',y='Spearman across groups')+theme_fig()}
