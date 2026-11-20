# Figure 3f - the ten gate C1-C5 programs of the final targets (CD71 green, CD119 purple row labels): signed program
# slope in each of the ten seeds (discovery 42-46 | confirmation 47-51) and the ten-seed mean. Rows shared with 3g.
# Data: optimization_20261003/program_seed_scores.csv; target_funnel_final_20261119/final_targets.csv.
PANEL<-list(id='3f_target_programs',blocks=c('1x2'),what='Final-target programs: seed slopes')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_rows.R'))
make_panel<-function(b){fp<-prog_rows()
 s<-read_src(file.path(OPT,'program_seed_scores.csv'),cols=c('combo','seed','signed_program_slope'));s<-s[s$combo%in%fp$combo,]
 s$row<-fp$row[match(s$combo,fp$combo)];s$v<-s$signed_program_slope*1e3;s$x<-s$seed-41+ifelse(s$seed>=47,.4,0)
 m<-aggregate(v~row,s,mean);m$x<-11.9;lim<-stats::quantile(abs(s$v),.95);n<-nrow(fp)
 cols<-c('#C9C9C9','white',PERISCOPE);brk<-min(fp$row[fp$antigen=='CD71'])-.5
 p<-ggplot()+
  geom_tile(data=rbind(s[,c('x','row','v')],m[,c('x','row','v')]),aes(x,row,fill=v),width=.92,height=.86)+
  annotate('segment',x=c(.55,5.95),xend=c(5.45,10.85),y=n+.6,yend=n+.6,colour=c(MID,PERISCOPE),linewidth=gglw('emph')*1.4)+
  annotate('segment',x=-6,xend=12.45,y=brk,yend=brk,colour=FAINT,linewidth=gglw('thin'))+
  scale_fill_gradientn(colours=cols,limits=c(-lim,lim),oob=scales::squish,guide='none')+
  scale_x_continuous(breaks=c(1,5,6.4,10.4,11.9),labels=c(42,46,47,51,'Mean'))+
  scale_y_continuous(breaks=fp$row,labels=fp$lab)+
  row_coord(n,c(.5,12.45))+labs(x='Seed',y=NULL)+row_theme()+
  theme(axis.text.y=element_text(colour=unname(TARGET_COL[as.character(fp$antigen)])))
 W<-b$w;H<-b$h;y1<-KEY_Y()
 k<-key_row(list(ki(kg_line(MID,lwd=glwd('emph')*1.4,w=3.2),'Discovery'),ki(kg_line(PERISCOPE,lwd=glwd('emph')*1.4,w=3.2),'Confirmation')),5.5,y1,h=H)
 x<-attr(k,'end')+3;k<-c(k,list(key_title('Slope',x,y1,H)),key_row(list(ki(kg_grad(cols,'−','+',w=10))),x+sw('Slope',TRUE)+1,y1,h=H))
 gTree(children=gList(suppressWarnings(gg_grob(p)),do.call(gList,k)))}
