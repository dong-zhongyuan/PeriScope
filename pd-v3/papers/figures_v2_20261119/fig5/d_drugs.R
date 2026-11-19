# Figure 5d - direction-matched drugs of the final targets: x = highest clinical phase, rows grouped by target; shape =
# route (circle direct on the target, diamond via the ligand), fill = antibody (colour) or other modality (grey).
# Data: target_funnel_final_20261119/final_target_drug_matches.csv.
PANEL<-list(id='5d_drug_matrix',blocks=c('1x1'),what='Direction-matched drugs of the final targets')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig5/_fun.R'))
make_panel<-function(b){
 d<-read_src(file.path(TF,'final_target_drug_matches.csv'),cols=c('antigen','drug','drug_type','route_class','max_phase'))
 d$antigen<-factor(d$antigen,c('CD71','CD119'));d<-d[order(d$antigen,-d$max_phase,d$drug),];d$y<-rev(seq_len(nrow(d)))
 d$lab<-ifelse(grepl('[0-9]',d$drug),d$drug,tools::toTitleCase(tolower(d$drug)));d$col<-TARGET_COL[as.character(d$antigen)];d$ab<-d$drug_type=='Antibody'
 ggplot(d,aes(max_phase,y))+geom_segment(aes(x=.5,xend=max_phase,yend=y),colour=FAINT,linewidth=gglw('emph'))+
  geom_point(shape=ifelse(d$route_class=='direct',21,23),fill=ifelse(d$ab,d$col,'white'),colour=d$col,stroke=.5,size=ggsz('big'))+
  scale_x_continuous(limits=c(.5,4.3),breaks=1:4,labels=c('I','II','III','IV'),expand=expansion(0))+
  scale_y_continuous(breaks=d$y,labels=d$lab,expand=expansion(add=.6))+
  labs(x='Clinical phase',y=NULL)+theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(colour=unname(d$col)))}
