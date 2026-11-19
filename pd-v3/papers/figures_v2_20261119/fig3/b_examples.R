# Figure 3b - five example programs. Top: dose-response curves at the five measured-dose quantiles (thin = seeds; thick =
# discovery mean, grey, and confirmation mean, colour). Bottom: the five strongest member genes of each program
# (discovery extremes in the program direction); colour = mean confirmation slope, area = fraction of positive seed
# slopes. Data: five_examples.csv, program_dose_curves.csv, gene_responses.csv (read-only).
PANEL<-list(id='3b_example_programs',blocks=c('1x3'),what='Five example programs: dose curves and lead genes')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig3/_prog.R'))
make_panel<-function(b){
 ex<-read_src(file.path(FRD,'five_examples.csv'),cols=c('combo'));pp<-prog_parse(ex$combo)
 pp$lab<-paste0(pp$antigen,' ',arrow_of(pp$dir),'\n',pp$blood,' → ',pp$brain);pp$f<-factor(pp$lab,pp$lab)
 cu<-read_src(file.path(FRD,'program_dose_curves.csv'),cols=c('combo','seed','quantile','response','phase'));cu$f<-pp$f[match(cu$combo,pp$combo)]
 cu$r<-cu$response*1e3;mn<-aggregate(r~f+quantile+phase,cu,mean)
 p1<-ggplot(cu,aes(quantile,r))+
  geom_hline(yintercept=0,colour=ZERO_COL,linewidth=gglw('thin'))+
  geom_line(aes(group=seed,colour=phase),linewidth=gglw('hair')*1.4,alpha=.6)+
  geom_line(data=mn,aes(colour=phase),linewidth=gglw('mean'))+
  facet_wrap(~f,nrow=1,scales='free_y')+
  scale_colour_manual(values=c(Discovery=MID,Confirmation=PERISCOPE),guide='none')+
  scale_x_continuous(breaks=c(.05,.5,.95),labels=c('5','50','95'))+scale_y_continuous(labels=num_minus,n.breaks=3)+
  labs(x='Dose quantile (%)',y='1,000 \u00d7 response')+theme_fig()+
  theme(strip.text=element_text(size=PT,face='plain',lineheight=.95,margin=margin(0,0,.6,0,'mm')),panel.spacing=unit(2.2,'mm'),axis.title.x=element_text(margin=margin(t=.4,unit='mm')),
   plot.margin=margin(LETTER_ZONE[['h']]+.2,.6,0,.6,'mm'))
 gr<-read_src(file.path(FRD,'gene_responses.csv'),cols=c('combo','gene','discovery','confirmation','positive_fraction','member'))
 top<-do.call(rbind,lapply(seq_len(nrow(pp)),function(i){z<-gr[gr$combo==pp$combo[i]&gr$member,];s<-if(pp$dir[i]=='up')-z$discovery else z$discovery
  z<-z[order(s),][1:5,];z$f<-pp$f[i];z$x<-1:5;z}))
 top$c<-top$confirmation*1e3;lim<-max(abs(top$c))
 p2<-ggplot(top,aes(x,1))+
  geom_point(aes(fill=c,size=positive_fraction),shape=21,colour=INK,stroke=.2)+
  facet_wrap(~f,nrow=1)+
  scale_fill_gradientn(colours=c('#3F6E9A','white',PERISCOPE),limits=c(-lim,lim),guide='none')+scale_size_area(max_size=ggsz(2.6),limits=c(0,1),guide='none')+
  scale_x_continuous(breaks=1:5,labels=NULL,expand=expansion(add=.6))+scale_y_continuous(NULL,breaks=NULL)+
  geom_text(aes(x,.55,label=gene),family=FONT,size=GT,angle=90,hjust=1,colour=INK)+
  coord_cartesian(ylim=c(-1.7,1.4),clip='off')+labs(x=NULL)+theme_fig()+
  theme(strip.text=element_blank(),panel.spacing=unit(2.2,'mm'),axis.line=element_blank(),axis.ticks=element_blank(),plot.margin=margin(0,.6,.3,.6,'mm'))
 g1<-gg_grob(p1);g2<-gg_grob(p2);W<-b$w;H<-b$h;h1<-26
 # p2 columns must sit under p1 panels: copy p1's widths so the facet panels align
 g2$widths<-g1$widths
 gTree(children=gList(grobTree(g1,vp=viewport(x=unit(0,'mm'),y=unit(H-h1,'mm'),width=unit(W,'mm'),height=unit(h1,'mm'),just=c(0,0))),
  grobTree(g2,vp=viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(W,'mm'),height=unit(H-h1,'mm'),just=c(0,0)))))}
