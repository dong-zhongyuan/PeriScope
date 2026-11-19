# Figure 1d - cell-type composition of the six GSE223138 blood donors as a mosaic: column width = donor cell count,
# segment height = cell-type share; diagnosis strip below with stage labels. Data: all_cohort_composition.csv (read-only).
PANEL<-list(id='1d_blood_mosaic',blocks=c('1x1'),what='Blood donor composition mosaic')
DES<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design'
CT_COL<-c(B_cell='#C9B48A',CD4_T='#E4DCCB',CD8_T='#B8AD95',NK='#8E8471',Monocyte_CD14=unname(BLOOD_COL['cMono']),
 Monocyte_CD16=unname(BLOOD_COL['nMono']),DC=unname(BLOOD_COL['cDC']),pDC=unname(BLOOD_COL['pDC']),unassigned=FAINT)
CT_LAB<-c(B_cell='B',CD4_T='CD4 T',CD8_T='CD8 T',NK='NK',Monocyte_CD14='cMono',Monocyte_CD16='nMono',DC='cDC',pDC='pDC',unassigned='Other')
make_panel<-function(b){
 d<-read_src(file.path(DES,'full_template_reuse/inputs/all_cohort_composition.csv'),cols=c('donor','condition','cell_type','n','cohort'))
 d<-d[d$cohort=='gse223138',];d$condition[d$condition=='Normal']<-'Control'
 don<-c('BL-N1','BL-N2','BL-E1','BL-E2','BL-L1','BL-L2');stopifnot(all(don%in%d$donor))
 tot<-tapply(d$n,d$donor,sum)[don];gap<-.025*sum(tot);x0<-cumsum(c(0,head(tot,-1)+gap));names(x0)<-don
 d$ct<-factor(d$cell_type,names(CT_COL));d<-d[order(d$donor,d$ct),]
 r<-do.call(rbind,lapply(don,function(k){x<-d[d$donor==k,];p<-x$n/sum(x$n);data.frame(row.names=NULL,donor=k,ct=x$ct,xmin=x0[k],xmax=x0[k]+tot[k],ymax=1-cumsum(c(0,head(p,-1))),ymin=1-cumsum(p))}))
 dx<-unique(d[,c('donor','condition')]);dx$xmin<-x0[dx$donor];dx$xmax<-x0[dx$donor]+tot[dx$donor]
 gp<-data.frame(lab=c('Ctrl','Early','Late'),lo=x0[c(1,3,5)],hi=x0[c(2,4,6)]+tot[c(2,4,6)])
 ggplot()+
  geom_rect(data=r,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=ct),colour='white',linewidth=gglw('hair'))+
  geom_rect(data=dx,aes(xmin=xmin,xmax=xmax,ymin=-.075,ymax=-.025),fill=DX_COL[dx$condition],colour=NA)+
  geom_text(data=gp,aes(x=(lo+hi)/2,y=-.1,label=lab),family=FONT,size=GT,colour=INK,vjust=1)+
  scale_fill_manual(values=CT_COL,labels=CT_LAB,breaks=names(CT_COL))+
  scale_x_continuous(expand=expansion(0))+
  scale_y_continuous(breaks=c(0,.5,1),labels=c('0','50','100'),expand=expansion(0))+
  coord_cartesian(ylim=c(-.075,1),clip='off')+
  labs(x=NULL,y='Cells (%)')+theme_fig(legend='right')+
  theme(axis.line.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),plot.margin=margin(LETTER_ZONE[['h']],.6,3.8,.6,'mm'),
   legend.key.size=unit(2.1,'mm'),legend.key.spacing.y=unit(.15,'mm'))+
  guides(fill=guide_legend(ncol=1))}
