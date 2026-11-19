# Figure 4i - MR estimates for every MR + coloc pair of the final targets (CD71/TFRC, CD119/IFNGR1): beta with 95% CI per
# SD of plasma protein; filled = pair-level BH q < 0.05, open = selected by the best-trait/assay rule; glyph colour =
# target, row strip = imaging modality. Data: target_funnel_final_20261119/final_target_mr_coloc_pairs.csv.
PANEL<-list(id='4i_final_target_forest',blocks=c('1x2'),what='Final-target MR forest')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig4/_mr.R'))
make_panel<-function(b){
 f<-read_src(file.path(TF,'final_target_mr_coloc_pairs.csv'),cols=c('antigen','UKB_gene','phenotype_description','family','beta','se','q_bh_global'))
 f$lab<-gsub('_',' ',f$phenotype_description);f$lab<-sub('^QSM ','',f$lab);f$lab<-sub('rfMRI connectivity \\(ICA100 edge ([0-9]+)\\)','rfMRI edge \\1',f$lab)
 f$lab<-sub('^Median T2star in (left|right) (\\w+).*$','T2* \\1 \\2',f$lab);f$lab<-sub('^Volume of grey matter in (Left|Right) (\\w+).*$','GM volume \\L\\1 \\2',f$lab,perl=TRUE);f$lab<-vapply(strsplit(f$lab,' '),function(w)paste(head(w,4),collapse=' '),'');f$lab<-paste0(f$UKB_gene,'  ',f$lab);f<-f[order(f$UKB_gene,-abs(f$beta/f$se)),];f$y<-rev(seq_len(nrow(f)))
 f$mod<-ifelse(f$family=='QSM','QSM',ifelse(grepl('T2',f$phenotype_description),'T2*',ifelse(grepl('volume',f$family,ignore.case=TRUE),'Volume','Other')))
 mc<-c(MOD_COL,Other=LIGHT);f$sig<-is.finite(f$q_bh_global)&f$q_bh_global<.05;f$col<-FT_COL[f$UKB_gene]
 xr<-max(abs(c(f$beta-1.96*f$se,f$beta+1.96*f$se)))*1.05
 ggplot(f,aes(beta,y))+geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_segment(aes(x=-xr*1.08,xend=-xr*1.02,yend=y),colour=mc[f$mod],linewidth=gglw('emph')*2.4)+
  geom_segment(aes(x=beta-1.96*se,xend=beta+1.96*se,yend=y),colour=f$col,linewidth=gglw('emph'))+
  geom_point(shape=21,fill=ifelse(f$sig,f$col,'white'),colour=f$col,stroke=.45,size=ggsz('big'))+
  scale_y_continuous(breaks=f$y,labels=f$lab,expand=expansion(add=.6))+scale_x_continuous(labels=num_minus,n.breaks=5)+
  coord_cartesian(xlim=c(-xr*1.08,xr),clip='off')+labs(x='MR beta (per SD protein)',y=NULL)+theme_fig()+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(colour=unname(f$col)))}
