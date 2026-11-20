# Figure 3j - agreement of matched-brain PD effects between the locked Kamath donors (x) and the independent GSE157783
# donors (y): every program evaluable in both cohorts (grey) and the ten final-target programs (colour).
# Data: optimization_20261003/all_program_evidence.csv; final programs as in 3f.
PANEL<-list(id='3j_donor_scores',blocks=c('1x1'),what='Cross-cohort agreement of program PD effects')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
make_panel<-function(b){
 e<-read_src(file.path(OPT,'all_program_evidence.csv'),cols=c('combo','kamath_donor_score_effect','independent_effect'))
 e<-e[is.finite(e$kamath_donor_score_effect)&is.finite(e$independent_effect),];e$x<-e$kamath_donor_score_effect;e$y<-e$independent_effect
 fp<-final_programs();f<-e[match(fp$combo,e$combo),];stopifnot(!anyNA(f$combo));f$antigen<-fp$antigen
 P<-ggplot()+
  geom_hline(yintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+geom_vline(xintercept=0,colour=ZERO_COL,linewidth=gglw('ref'),linetype=ZERO_LTY)+
  geom_point(data=e,aes(x,y),shape=16,colour=MID,alpha=.22,size=ggsz('unit'))+
  geom_point(data=f,aes(x,y),shape=21,fill=TARGET_COL[f$antigen],colour='white',stroke=.3,size=ggsz('point'))+
  scale_x_continuous(labels=num_minus,n.breaks=4)+scale_y_continuous(labels=num_minus,n.breaks=4)+
  labs(x='Kamath PD − control',y='GSE157783 PD − control')+theme_fig()+theme(plot.margin=margin(LETTER_ZONE[['h']]+.8,.6,.6,.6,'mm'))
 W<-b$w;H<-b$h
 k<-key_row(list(ki(kg_dot(TARGET_COL[['CD71']],'white',1.4),'CD71'),ki(kg_dot(TARGET_COL[['CD119']],'white',1.4),'CD119')),6,LETTER_ZONE[['h']]-1.6,h=H)
 gTree(children=gList(gg_grob(P),do.call(gList,k)))}
