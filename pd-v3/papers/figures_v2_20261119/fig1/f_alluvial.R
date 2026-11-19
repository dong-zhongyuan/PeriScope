# Figure 1f - published blood-cell annotations (left) mapped to the cell states used for modelling (right): ribbon
# width = cells; counts beside each node. Annotations under 1,000 cells that are not modelled states are pooled as
# Other. Node slots are at least one text line tall so labels never collide. Data: panel_grid/objects/1E.rds.
PANEL<-list(id='1f_annotation_alluvial',blocks=c('1x2'),what='Published annotations to modelled blood states')
DES<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design'
slots<-function(ord,v,avail,gap,mh){f<-function(k)sum(pmax(v[ord]*k,mh))+gap*(length(ord)-1)-avail
 k<-stats::uniroot(f,c(1e-12,avail/sum(v)),tol=1e-14)$root;list(k=k,h=pmax(v[ord]*k,mh))}
make_panel<-function(b){
 o<-readRDS(src_log(file.path(DES,'panel_grid/objects/1E.rds')));d<-o$data[,c('source','target','count')]
 d<-aggregate(count~source+target,d,sum);d<-d[d$count>0,]
 tl<-c(B_cell='B',CD4_T='CD4 T',CD8_T='CD8 T',NK='NK',Monocyte_CD14='cMono',Monocyte_CD16='nMono',DC='cDC',pDC='pDC',unassigned='Other',Other='Other')
 d$target<-ifelse(d$target%in%names(tl),tl[d$target],d$target)
 st0<-tapply(d$count,d$source,sum);keep<-names(st0)[st0>=1000|names(st0)%in%c('cDC','pDC','CD16 mono','CD14 mono')]
 d$source[!d$source%in%keep]<-'Other';d$source<-sub(' cells$','',d$source);d<-aggregate(count~source+target,d,sum)
 tcol<-c(B='#C9B48A',`CD4 T`='#E4DCCB',`CD8 T`='#B8AD95',NK='#8E8471',cMono=unname(BLOOD_COL['cMono']),nMono=unname(BLOOD_COL['nMono']),cDC=unname(BLOOD_COL['cDC']),pDC=unname(BLOOD_COL['pDC']),Other=FAINT)
 tord<-intersect(names(tcol),unique(d$target));stopifnot(setequal(tord,unique(d$target)))
 st<-tapply(d$count,d$source,sum);tt<-tapply(d$count,d$target,sum)
 dom<-vapply(names(st),function(s){x<-d[d$source==s,];x$target[which.max(x$count)]},'')
 sord<-names(st)[order(names(st)=='Other',match(dom,tord),-st)]
 W<-b$w;H<-b$h;top<-LETTER_ZONE[['h']]+.3;bot<-.5;gap<-.5;mh<-LINE_MM*1.08;bw<-1.8
 cw<-function(v)max(sw(fint(v)))
 lw<-max(sw(sord))+.4;rw<-bw+.6+max(sw(tord))+.2
 x0<-lw+.8;x1<-W-rw
 L<-slots(sord,st,H-top-bot,gap,mh);R<-slots(tord,tt,H-top-bot,gap,mh);k<-min(L$k,R$k)
 pos<-function(ord,hs){y<-top+c(0,cumsum(hs+gap))[seq_along(ord)];setNames(y,ord)}
 yL<-pos(sord,L$h);hL<-setNames(L$h,sord);yR<-pos(tord,R$h);hR<-setNames(R$h,tord)
 # ribbons start at the slot centre region scaled by k so the ribbon width is proportional to cells everywhere
 cL<-yL+hL/2-st[sord]*k/2;cR<-yR+hR/2-tt[tord]*k/2;names(cL)<-sord;names(cR)<-tord
 offL<-cL;offR<-cR;g<-list()
 d<-d[order(match(d$source,sord),match(d$target,tord)),]
 for(i in seq_len(nrow(d))){s<-d$source[i];t<-d$target[i];h<-d$count[i]*k
  a0<-offL[s];offL[s]<-a0+h;b0<-offR[t];offR[t]<-b0+h
  g<-c(g,list(dribbon(x0+bw,a0,a0+h,x1,b0,b0+h,fill=alpha_hex(tcol[t],.6))))}
 for(s in sord){yc<-yL[s]+hL[s]/2;g<-c(g,list(drect(x0,cL[s],bw,st[s]*k,fill=MID,col=NA),
  dtxt(s,x0-.8,yc,just=c('right','centre'))))}
 for(t in tord){yc<-yR[t]+hR[t]/2;g<-c(g,list(drect(x1,cR[t],bw,tt[t]*k,fill=tcol[t],col=NA),
  dtxt(t,x1+bw+.6,yc,just=c('left','centre'))))}
 dgroup(g)}
