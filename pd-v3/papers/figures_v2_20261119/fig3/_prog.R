# shared helpers for Figure 3: data roots and program-name parsing
FRD<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design/manuscript_review/figure_redesign_20261004/data'
OPT<-'/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003'
TF<-'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/target_funnel_final_20261119'
B_LAB<-c(cDC='cDC',classical_mono='cMono',nonclassical_mono='nMono',pDC='pDC');S_LAB<-c(astro='Astro',microglia_mhc2='MG-AP',microglia_homeostatic='MG-H')
prog_parse<-function(combo){x<-sub('^raw::','',combo);p<-strsplit(x,'__');data.frame(combo=combo,antigen=vapply(p,`[`,'',1),blood=B_LAB[vapply(p,`[`,'',2)],
 brain=S_LAB[vapply(p,`[`,'',4)],dir=vapply(p,`[`,'',5),stringsAsFactors=FALSE)}
arrow_of<-function(d)ifelse(d=='up','↑','↓')
final_programs<-function(){t<-read_src(file.path(TF,'final_targets.csv'),cols=c('antigen','programs_gate1to5'))
 do.call(rbind,lapply(seq_len(nrow(t)),function(i){p<-trimws(strsplit(t$programs_gate1to5[i],';')[[1]]);data.frame(antigen=t$antigen[i],combo=sub('\\[.*$','',sub('^raw::','',p)),
  cls=sub('^.*\\[(.*)\\]$','\\1',p),stringsAsFactors=FALSE)}))}
