# Figure 3k - model-side nomination ledger: one row per gate (all antigens, then C1-C5 on the same program), one unit
# square per antigen (213). Dark = still in, colour = removed at this gate, faint = removed earlier; the two final targets
# carry their colours in every row. Data: target_funnel_final_20261119/antigen_gate_matrix.csv.
PANEL<-list(id='3k_model_nomination',blocks=c('1x2'),what='Model-side nomination ledger C1-C5')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
make_panel<-function(b){W<-b$w;H<-b$h;assign('PANEL_H',H,envir=globalenv())
 a<-read_src(file.path(TF,'antigen_gate_matrix.csv'),cols=c('antigen','gate1','gate2','gate3','gate4','gate5'))
 G<-as.matrix(a[,paste0('gate',1:5)]);G<-G=='True'|G==TRUE;out<-apply(cbind(TRUE,G),1,function(r)if(all(r))6 else which(!r)[1]-1)
 tg<-names(TARGET_COL);a$out<-out;a$pri<-ifelse(a$antigen%in%tg,0,1);a<-a[order(-a$out,a$pri,a$antigen),];n<-nrow(a);stopifnot(n==213)
 nr<-3;nc<-ceiling(n/nr);lw<-9;pitch<-(W-lw-.6)/nc;sq<-pitch*.84;top<-LETTER_ZONE[['h']]+.8;rh<-(H-top-4.2)/6
 gl<-c('All','C1','C2','C3','C4','C5');g<-list()
 for(r in 0:5){y0<-top+r*rh;g<-c(g,list(dtxt(gl[r+1],.6,y0+nr*pitch/2,just=c('left','centre'),bold=TRUE)))
  for(i in seq_len(n)){cc<-(i-1)%/%nr;rr<-(i-1)%%nr;x<-lw+cc*pitch;y<-y0+rr*pitch
   st<-if(a$out[i]>r)'in' else if(a$out[i]==r)'out' else 'gone'
   col<-if(a$antigen[i]%in%tg&&st=='in')TARGET_COL[[a$antigen[i]]] else switch(st,`in`=MID,out=PERISCOPE_LIGHT,gone=FAINT)
   if(r==0)col<-if(a$antigen[i]%in%tg)TARGET_COL[[a$antigen[i]]] else MID
   g<-c(g,list(drect(x,y,sq,sq,fill=col,col=NA)))}}
 k<-key_row(list(ki(kg_box(MID,w=1.8,ht=1.8),'Retained'),ki(kg_box(PERISCOPE_LIGHT,w=1.8,ht=1.8),'Removed'),ki(kg_box(FAINT,w=1.8,ht=1.8),'Earlier'),
  ki(kg_box(TARGET_COL[['CD71']],w=1.8,ht=1.8),'CD71'),ki(kg_box(TARGET_COL[['CD119']],w=1.8,ht=1.8),'CD119')),lw,H-1.8,gap=2,h=H)
 dgroup(c(g,k))}
