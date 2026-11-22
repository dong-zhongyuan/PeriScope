# shared helpers for Figure 5: funnel stages per antigen and final-target program genes
TF<-'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/target_funnel_g8rel_20261122'
DRUG<-'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/target_funnel_g8rel_20261122'
OPT<-'/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003'
GATES<-c('All','C1','C2','C3','C4','C5','MR','Coloc','Direction','Drug')
GCOL<-c('gate1','gate2','gate3','gate4','gate5','gate6','gate7','gate8a_direction','gate8_expanded')
tf<-function(x)x%in%c(TRUE,'True','TRUE')
antigen_stage<-function(){a<-read_src(file.path(TF,'antigen_gate_matrix.csv'))
 G<-sapply(GCOL,function(k)tf(a[[k]]));a$out<-apply(cbind(TRUE,G),1,function(r)if(all(r))length(r) else which(!r)[1]-1);a}
DIR_COL<-c(inc='#D98C2B',dec='#3F6E9A',unr=LIGHT)
key_title<-function(t,x,y,h)textGrob(t,unit(x,'mm'),unit(h-y,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,fontface='bold',col=INK))
final_prog_genes<-function(){t<-read_src(file.path(TF,'final_targets.csv'),cols=c('antigen','programs_gate1to5'))
 fp<-do.call(rbind,lapply(seq_len(nrow(t)),function(i){p<-trimws(strsplit(t$programs_gate1to5[i],';')[[1]]);data.frame(antigen=t$antigen[i],combo=sub('\\[.*$','',sub('^raw::','',p)),cls=sub('^.*\\[(.*)\\]$','\\1',p))}))
 e<-read_src(file.path(Sys.getenv('FIGV2_ROOT'),'fig5/data/final_program_genes.csv'),cols=c('combo','n_genes','member_genes'));e<-e[match(fp$combo,e$combo),];stopifnot(!anyNA(e$combo));fp$genes<-e$member_genes;fp$n_genes<-e$n_genes
 B<-c(cDC='cDC',classical_mono='cMono',nonclassical_mono='nMono',pDC='pDC');S<-c(astro='Astro',microglia_mhc2='MG-AP')
 p<-strsplit(fp$combo,'__');fp$blood<-B[vapply(p,`[`,'',2)];fp$brain<-S[vapply(p,`[`,'',4)];fp$dir<-vapply(p,`[`,'',5)
 fp$lab<-paste0(fp$blood,' → ',fp$brain,' ',ifelse(fp$dir=='up','↑','↓'));fp$antigen<-factor(fp$antigen,c('CD71','CD119'))
 fp<-fp[order(fp$antigen,fp$brain,match(fp$blood,B)),];fp}
