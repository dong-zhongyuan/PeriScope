# Figure 1c - resource passports: one row per resource (six transcriptomic resources plus the CITE-seq bridge).
# Chip = cohort colour (as in 1b), filled for resources that enter training and open for evaluation-only resources;
# unit squares = donors coloured by diagnosis (open squares = held-out CITE-seq donors P6-P8); bar = cells on a log
# scale with the count. Data: all_cohort_composition.csv, input_1A.csv, brain_diagnosis_lookup.json and the CITE-seq
# bridge donors (bridge_data.npz, counted once: 8 donors, 161,764 cells).
PANEL<-list(id='1c_resource_passports',blocks=c('1x2'),what='Resource passports: donors, diagnosis, cells, role')
DES<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design'
CITE_CELLS<-161764;CITE_DONORS<-paste0('P',1:8);CITE_HELD<-c('P6','P7','P8')
make_panel<-function(b){
 comp<-read_src(file.path(DES,'full_template_reuse/inputs/all_cohort_composition.csv'),cols=c('donor','condition','n','cohort'))
 ring<-read_src(file.path(DES,'template_semantic_review/data/input_1A.csv'),cols=c('category','count'))
 lk<-jsonlite::fromJSON(src_log(file.path(DES,'template_semantic_review/inputs/brain_diagnosis_lookup.json')))
 cn<-c(gse223138='GSE223138',moquin_beaudry2025='MB2025',gse178265='GSE178265',gse184950='GSE184950',gse157783='GSE157783',gse253975='GSE253975')
 dxm<-c(Normal='Control',Control='Control',Ctrl='Control',CTRL='Control',Early='Early PD',Late='Late PD',`Early PD`='Early PD',`Late PD`='Late PD',PD='PD',PDD='PDD',MSA='MSA',PSP='PSP',LBD='LBD')
 u<-unique(comp[,c('cohort','donor','condition')]);u$res<-cn[u$cohort]
 k<-u$res=='GSE178265';u$condition[k]<-vapply(u$donor[k],function(w)if(is.null(lk[[w]]))NA_character_ else lk[[w]][1],'');stopifnot(!anyNA(u$condition))
 u$dx<-dxm[u$condition];stopifnot(!anyNA(u$dx))
 cells<-setNames(ring$count[grepl('^cohort: ',ring$category)],sub('^gse','GSE',sub('^cohort: ','',ring$category[grepl('^cohort: ',ring$category)])))
 rows<-c('GSE223138','MB2025','CITE-seq','GSE178265','GSE184950','GSE157783','GSE253975')
 train<-c(GSE223138=TRUE,MB2025=TRUE,`CITE-seq`=TRUE,GSE178265=TRUE,GSE184950=FALSE,GSE157783=FALSE,GSE253975=FALSE)
 chip<-c(COHORT_COL,`CITE-seq`=MID);cells<-c(cells,`CITE-seq`=CITE_CELLS)
 W<-b$w;H<-b$h;top<-LETTER_ZONE[['h']]+.4;keyh<-3.6;rh<-(H-top-keyh)/length(rows)
 sq<-1.45;pg<-.3;pitch<-sq+pg;nr<-2
 xn<-3.0;xd<-xn+max(sw(clab(rows)))+1.4
 maxd<-max(table(u$res),length(CITE_DONORS));ncol<-ceiling(maxd/nr);xc<-xd+ncol*pitch+.6
 xb<-xc+.4;lw<-W-xb-4;lo<-3.7;hi<-log10(max(cells))
 g<-list()
 for(i in seq_along(rows)){r<-rows[i];yc<-top+(i-.5)*rh
  g<-c(g,list(drect(.4,yc-1.1,2.2,2.2,fill=if(train[r])chip[r]else'white',col=chip[r],lwd=glwd('emph'))))
  g<-c(g,list(dtxt(clab(r),xn,yc,just=c('left','centre'))))
  if(r=='CITE-seq'){dd<-data.frame(donor=CITE_DONORS,dx='Control')}else{dd<-u[u$res==r,c('donor','dx')]}
  dd<-dd[order(match(dd$dx,names(DX_COL)),dd$donor),];nd<-nrow(dd)
  for(j in seq_len(nd)){cc<-(j-1)%/%nr;rr<-(j-1)%%nr;x<-xd+cc*pitch;y<-yc-(nr*pitch-pg)/2+rr*pitch
   held<-r=='CITE-seq'&&dd$donor[j]%in%CITE_HELD
   g<-c(g,list(drect(x,y,sq,sq,fill=if(held)'white'else DX_COL[dd$dx[j]],col=if(held)DX_COL[dd$dx[j]]else NA,lwd=glwd('data'))))}
  L<-max(.4,(log10(cells[r])-lo)/(hi-lo)*lw)
  g<-c(g,list(drect(xb,yc-.9,L,1.8,fill=chip[r],col=NA)))}
 ya<-top+length(rows)*rh+.3;g<-c(g,list(dline(xb,ya,xb+lw,ya,col=INK,lwd=glwd('axis'))))
 for(v in c(1e4,1e5)){xx<-xb+(log10(v)-lo)/(hi-lo)*lw;g<-c(g,list(dline(xx,ya,xx,ya+.8,col=INK,lwd=glwd('tick')),dtxt(fint(v),xx,ya+1.0,just=c('centre','top'))))}
 ky<-H-keyh/2;x<-.4
 g<-c(g,list(drect(x,ky-.8,1.6,1.6,fill=MID,col=MID,lwd=glwd('emph')),dtxt('Training',x+2.2,ky,just=c('left','centre'))));x<-x+2.2+sw('Training')+2.4
 g<-c(g,list(drect(x,ky-.8,1.6,1.6,fill='white',col=MID,lwd=glwd('emph')),dtxt('Evaluation',x+2.2,ky,just=c('left','centre'))));x<-x+2.2+sw('Evaluation')+2.4
 dgroup(g)}
