# Figure 2f - ablation design: blood (left) and brain (right) donor groups of two conditions linked within strata.
# Full: links within the matched condition; random pairing: links shuffled across conditions; coupling off: no coupling
# term. Glyphs key the arms in 2g-2i.
PANEL<-list(id='2f_ablation_design',blocks=c('1x1'),what='Ablation arms schematic')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig2/_arms.R'))
make_panel<-function(b){W<-b$w;H<-b$h;assign('PANEL_H',H,envir=globalenv())
 C<-function(x,y,d,fill,col=NA)circleGrob(unit(x,'mm'),unit(H-y,'mm'),r=unit(d/2,'mm'),gp=gpar(fill=fill,col=col,lwd=glwd('thin')))
 cc<-c(DX_COL[['Control']],DX_COL[['Control']],DX_COL[['PD']],DX_COL[['PD']])
 top<-LETTER_ZONE[['h']]+.6;rh<-(H-top-3.4)/3;g<-list();arms<-names(ARM_SHAPE)
 xb<-26;xr<-36.6
 for(i in 1:3){y0<-top+(i-1)*rh;yy<-y0+rh*c(.2,.4,.6,.8);a<-arms[i]
  g<-c(g,dmark(1.9,y0+rh/2,switch(a,Full='periscope',`Random pairing`='cellot',`Coupling off`='ridge'),d=2.4))
  g<-c(g,list(dtxt(a,3.6,y0+rh/2,just=c('left','centre'),bold=TRUE)))
  links<-switch(a,Full=cbind(1:4,1:4),`Random pairing`=cbind(1:4,c(3,4,1,2)),`Coupling off`=NULL)
  if(!is.null(links))for(k in seq_len(nrow(links)))g<-c(g,list(dline(xb+1,yy[links[k,1]],xr-1,yy[links[k,2]],col=ARM_COL[[a]],lwd=glwd('data'))))
  for(k in 1:4)g<-c(g,list(C(xb,yy[k],1.6,cc[k]),C(xr,yy[k],1.6,cc[k])))
  g<-c(g,list(drect(xb-1.4,y0+rh*.12,2.8,rh*.8,fill=NA,col=BLOOD_COL[['cMono']],lwd=glwd('thin'),r=.6),drect(xr-1.4,y0+rh*.12,2.8,rh*.8,fill=NA,col=BRAIN_COL[['Astro']],lwd=glwd('thin'),r=.6)))
  if(i<3)g<-c(g,list(dline(1,y0+rh,W-1,y0+rh,col=FAINT,lwd=glwd('hair'))))}
 g<-c(g,list(dtxt('Blood',xb,H-.2,just=c('centre','bottom'),col=MID),dtxt('Brain',xr,H-.2,just=c('centre','bottom'),col=MID)))
 kk<-key_row(list(ki(kg_dot(DX_COL[['Control']],NA,1.6),'Control'),ki(kg_dot(DX_COL[['PD']],NA,1.6),'PD')),1,H-1.4,gap=1.4,h=H)
 dgroup(c(g,kk))}
