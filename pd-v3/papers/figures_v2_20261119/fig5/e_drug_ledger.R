# Figure 5e - every drug record of the two final targets (one square per record): rows = route (direct, ligand route,
# complex partner, gene silencing), column groups = action relative to the model restorative direction (matched,
# opposite, unknown or cell-killing); fill = modality. Data: drug_evidence_expanded_20261119/drug_records.csv.gz.
PANEL<-list(id='5e_drug_ledger',blocks=c('1x1'),what='Drug record ledger of the final targets')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig5/_fun.R'))
make_panel<-function(b){W<-b$w;H<-b$h;assign('PANEL_H',H,envir=globalenv())
 d<-read_src(file.path(DRUG,'drug_records.csv.gz'),cols=c('antigen','route_class','drug_type','direction'));d<-d[d$antigen%in%c('CD71','CD119'),]
 md<-c(CD71='down',CD119='down');d$rel<-ifelse(d$direction==md[d$antigen],'Same',ifelse(d$direction%in%c('up','down'),'Reverse','Other'))
 rt<-c(direct='Direct',ligand_route='Ligand',complex_partner='Partner',gene_silencing='Silencing');d<-d[d$route_class%in%names(rt),]
 mc<-c(Antibody=PERISCOPE,`Small molecule`=MID,Protein='#D98C2B',Oligonucleotide='#2E9C8E',Unknown=LIGHT);d$drug_type[!d$drug_type%in%names(mc)]<-'Unknown'
 rt<-rt[names(rt)%in%d$route_class];rel<-c('Same','Reverse','Other');top<-LETTER_ZONE[['h']]+3.2;lw<-10.5;cw<-(W-lw-1)/3;rh<-(H-top-7)/length(rt);sq<-1.55;pg<-.35;g<-list()
 for(j in seq_along(rel))g<-c(g,list(dtxt(rel[j],lw+(j-.5)*cw,top-1.6,col=MID)))
 for(i in seq_along(rt)){y0<-top+(i-1)*rh;g<-c(g,list(dtxt(rt[[i]],.6,y0+rh/2,just=c('left','centre'))))
  for(j in seq_along(rel)){z<-d[d$route_class==names(rt)[i]&d$rel==rel[j],];z<-z[order(match(z$drug_type,names(mc))),];if(!nrow(z))next
   nc<-floor((cw-1)/(sq+pg));for(k in seq_len(nrow(z))){cc<-(k-1)%%nc;rr<-(k-1)%/%nc
    g<-c(g,list(drect(lw+(j-1)*cw+.5+cc*(sq+pg),y0+.4+rr*(sq+pg),sq,sq,fill=mc[[z$drug_type[k]]],col=NA)))}}
  if(i<length(rt))g<-c(g,list(dline(.6,y0+rh,W-.6,y0+rh,col=FAINT,lwd=glwd('hair'))))}
 used<-names(mc)[names(mc)%in%d$drug_type];lb<-c(Antibody='Antibody',`Small molecule`='Small',Protein='Protein',Oligonucleotide='Oligo',Unknown='Other')
 h1<-head(used,2);h2<-tail(used,-2)
 k<-c(key_row(lapply(h1,function(u)ki(kg_box(mc[[u]],w=1.6,ht=1.6),lb[[u]])),.6,H-4.6,gap=1.6,h=H),if(length(h2))key_row(lapply(h2,function(u)ki(kg_box(mc[[u]],w=1.6,ht=1.6),lb[[u]])),.6,H-1.6,gap=1.6,h=H))
 dgroup(c(g,k))}
