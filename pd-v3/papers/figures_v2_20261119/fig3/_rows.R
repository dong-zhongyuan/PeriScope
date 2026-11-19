# shared rows for Figure 3f and 3g (same row order, same vertical window so the two panels align); keys sit in a
# header strip above the rows, so both panels carry identical top and bottom structure
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
prog_rows<-function(){fp<-final_programs();pp<-prog_parse(fp$combo);fp<-cbind(fp,pp[,c('blood','brain','dir')])
 fp$antigen<-factor(fp$antigen,c('CD71','CD119'));fp<-fp[order(fp$antigen,match(fp$brain,c('Astro','MG-AP')),match(fp$blood,B_LAB),fp$dir),]
 fp$row<-rev(seq_len(nrow(fp)));fp$lab<-paste0(fp$blood,' → ',fp$brain,' ',arrow_of(fp$dir));fp}
ROW_TOP<-.8;HEAD<-2.8
row_coord<-function(n,xlim)coord_cartesian(xlim=xlim,ylim=c(.5,n+ROW_TOP),expand=FALSE,clip='off')
row_theme<-function()theme_fig()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(LETTER_ZONE[['h']]+HEAD,.6,.6,.6,'mm'))
key_title<-function(t,x,y,h)textGrob(t,unit(x,'mm'),unit(h-y,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,fontface='bold',col=INK))
KEY_Y<-function()LETTER_ZONE[['h']]+1.3
