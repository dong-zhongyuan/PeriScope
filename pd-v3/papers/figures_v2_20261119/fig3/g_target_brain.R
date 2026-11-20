# Figure 3g - matched-brain disease effect of the same ten programs (rows as 3f): PD-minus-control donor-score effect
# in the matched recipient state in Kamath (locked donors) and GSE157783; a grey bar spans the two cohorts of each program. Data: optimization_20261003/all_program_evidence.csv.
PANEL<-list(id='3g_target_brain',blocks=c('1x2'),what='Final-target programs: matched-brain effects')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_rows.R'))
make_panel<-function(b){fp<-prog_rows();n<-nrow(fp)
 e<-read_src(file.path(OPT,'all_program_evidence.csv'),cols=c('combo','kamath_donor_score_effect','independent_effect'))
 e<-e[match(fp$combo,e$combo),];stopifnot(!anyNA(e$combo))
 co<-c(Kamath='kamath',GSE157783='independent');sh<-c(Kamath=21,GSE157783=22)
 L<-do.call(rbind,lapply(names(co),function(k){ef<-if(co[[k]]=='independent')e$independent_effect else e[[paste0(co[[k]],'_donor_score_effect')]]
  data.frame(row=fp$row,cohort=k,v=ef)}))
 L$y<-L$row;L$col<-TARGET_COL[as.character(fp$antigen[match(L$row,fp$row)])]
 rg<-aggregate(v~row,L,range);rg<-data.frame(row=rg$row,lo=rg$v[,1],hi=rg$v[,2]);xr<-max(abs(L$v),na.rm=TRUE)*1.1
 brk<-min(fp$row[fp$antigen=='CD71'])-.5
 p<-ggplot()+
  geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  annotate('segment',x=-xr,xend=xr,y=brk,yend=brk,colour=FAINT,linewidth=gglw('thin'))+
  geom_segment(data=rg,aes(x=lo,xend=hi,y=row,yend=row),colour=FAINT,linewidth=gglw('emph')*1.6)+
  geom_point(data=L,aes(v,y,shape=cohort),fill=L$col,colour=L$col,stroke=.35,size=ggsz('point'))+
  scale_shape_manual(values=sh,guide='none')+
  scale_x_continuous(labels=num_minus,n.breaks=5)+scale_y_continuous(breaks=fp$row,labels=NULL)+
  row_coord(n,c(-xr,xr))+labs(x='PD − control',y=NULL)+row_theme()
 W<-b$w;H<-b$h;y1<-KEY_Y()
 kd<-function(pch,fill=MID)list(w=2,f=function(x,yc)pointsGrob(unit(x+1,'mm'),unit(yc,'mm'),pch=pch,size=unit(1.7,'mm'),gp=gpar(col=MID,fill=fill,lwd=glwd('thin'))))
 k<-key_row(lapply(names(sh),function(z)ki(kd(sh[[z]]),z)),6,y1,h=H)
 gTree(children=gList(gg_grob(p),do.call(gList,k)))}
