# common.R - PeriScope figures_v2 shared figure system (2026-11-19).
# Sourced by render_panel.R before every panel script; panel scripts never edit it.
# Binding user rules: one font (Arial) and one size (7 pt) for all in-panel text, axis titles bold; panel letters 9 pt bold, added only at
# page assembly; no panel titles, captions, notes or sentences; no bar chart with 1-3 bars; every panel is an integer
# block of the u x u grid (u = 41.25 mm, gutter g = 3 mm, page 4 x 5 units = 174 x 218.25 mm) and is laid out for that
# block by its own code (positions and sizes chosen per block; font sizes never change; nothing is scaled); pending
# content is an outlined blank block. Built from hh figures_r/_tools/common.R (helpers kept: strips, rails, waffles,
# rings, ribbons, hexbins, keys) and pd-v3 01_mr_imaging_scan/figures_r (Arial runtime, text audit, modality colours).
# Conventions: grid helpers take mm measured from the PANEL's top-left corner; h = panel height (default PANEL_H, set by
# the renderer); ggplot panels are built on theme_fig(); every input file is read with read_src() (logged to the manifest).
suppressPackageStartupMessages({library(grid);library(gtable);library(ggplot2)})
`%||%`<-function(a,b)if(is.null(a)||length(a)==0)b else a
if(!exists('SYS_DIR'))SYS_DIR<-Sys.getenv('FIGV2_SYS','/public/home/mengxl/dzy/pd-v3/papers/figures_v2_20261119/_system')
FIG_ROOT<-dirname(SYS_DIR)

# ================================================================== typography (one family, one size)
FONT<-'Arial'
PT<-7                     # every in-panel text: ticks, axis titles, legends, category and data labels
PT_LETTER<-9              # panel letters (bold), page assembly only
GT<-PT/.pt                # geom_text()/annotate('text') size
LINE_MM<-PT*25.4/72.27    # 2.46 mm: one 7-pt line; use >= 1.15 * LINE_MM between stacked labels
CAP_MM<-.716*PT*25.4/72.27   # Arial cap height at 7 pt (1.77 mm)

# ================================================================== unit grid
U_MM<-41.25;G_MM<-3;PAGE_COLS<-4L;PAGE_ROWS<-5L
span_mm<-function(n)n*U_MM+(n-1)*G_MM
PAGE_W<-span_mm(PAGE_COLS);PAGE_H<-span_mm(PAGE_ROWS)              # 174 x 218.25 mm
stopifnot(abs(PAGE_W-174)<1e-9)
BLOCKS<-c('1x1','1x2','1x3','1x4','2x2','2x3','2x4','2x1','3x1','4x1','3x2','4x2')   # rows x cols (height x width)
block<-function(b){if(is.list(b))return(b);b<-as.character(b);if(!b%in%BLOCKS)stop('block ',b,' is not one of ',paste(BLOCKS,collapse=' '),call.=FALSE)
 rc<-as.integer(strsplit(b,'x')[[1]]);list(name=b,rows=rc[1],cols=rc[2],w=span_mm(rc[2]),h=span_mm(rc[1]))}
block_of_size<-function(w,h,tol=.02){for(b in BLOCKS){x<-block(b);if(abs(x$w-w)<tol&&abs(x$h-h)<tol)return(b)};NA_character_}
cell_xy<-function(row,col)c(x=(col-1)*(U_MM+G_MM),y=(row-1)*(U_MM+G_MM))   # top-left of a grid cell on the page (mm)
LETTER_ZONE<-c(w=4.0,h=3.7)   # kept empty in every panel (top-left corner): the assembler prints the 9-pt letter there
LETTER_BASE<-2.8             # letter baseline, mm below the block top

# ================================================================== lines and marks (physical sizes)
# R lwd 1 = 1/96 inch = 0.75 pt on the cairo and ragg devices, so a width in points is pt / 0.75 in grid and
# pt / 0.75 / .pt in ggplot (linewidth). All panels use these widths only.
LW<-c(axis=.4,tick=.35,data=.5,thin=.25,hair=.15,mean=.9,ref=.4,frame=.3,link=.35,emph=.7)
glwd<-function(k)if(is.character(k))unname(LW[k])/.75 else k/.75          # grid gpar(lwd=)
gglw<-function(k)glwd(k)/.pt                                                # ggplot linewidth=
TICK_MM<-.8
D<-c(unit=.6,small=.9,point=1.3,mean=1.8,big=2.4)                            # point diameters (mm)
ggsz<-function(d,stroke=.25){d<-if(is.character(d))unname(D[d])else d;max((d*72.27/25.4/.75-stroke*.stroke/2)/.pt,.05)}   # ggplot size for a diameter
REF_LTY<-'22';ZERO_LTY<-'22'

# ================================================================== palette (one meaning per colour, whole paper)
INK<-'#1A1A1A';MID<-'#595959';GUIDE<-'#8C8C8C';LIGHT<-'#BDBDBD';FAINT<-'#E3E3E3';PAPER<-'#F2F2F2';MISSING<-'#E6E6E6'
FRAME_COL<-'#A6A6A6';ZERO_COL<-'#BDBDBD';REF_COL<-'#595959'
alpha_hex<-function(col,a)grDevices::adjustcolor(col,alpha.f=a)
tint<-function(col,t){m<-grDevices::col2rgb(col)/255;grDevices::rgb(1-t*(1-m[1,]),1-t*(1-m[2,]),1-t*(1-m[3,]))}   # opaque mix with white (t = share of col)
shade<-function(col,t){m<-grDevices::col2rgb(col)/255;grDevices::rgb(m[1,]*(1-t),m[2,]*(1-t),m[3,]*(1-t))}
PERISCOPE<-'#C2185B';PERISCOPE_LIGHT<-tint(PERISCOPE,.35);PERISCOPE_FAINT<-tint(PERISCOPE,.14)
# comparators: greys plus shapes (filled = learned comparator, open = summary baseline); PeriScope is the only coloured method
METHODS<-data.frame(key=c('periscope','cellot','ot_ridge','conditional_icnn','conditional_mlp','ridge',
  'state_mean','subtype_mean','shared_label_mean','shared_label_marginal','zero_response'),
 label=c('PeriScope','CellOT','OT–Ridge','Conditional ICNN','Conditional MLP','Ridge','State mean','Subtype mean','Label mean','Label marginal','Zero response'),
 fill=c(PERISCOPE,'#3C3C3C','#5C5C5C','#7A7A7A','#979797','#B0B0B0','white','white','white','white','white'),
 col=c(PERISCOPE,'#3C3C3C','#5C5C5C','#7A7A7A','#979797','#8C8C8C','#5C5C5C','#5C5C5C','#7A7A7A','#7A7A7A','#5C5C5C'),
 shape=c(21,24,22,23,25,22,22,23,24,25,4),stringsAsFactors=FALSE)
method_row<-function(k){i<-match(k,METHODS$key);if(anyNA(i))i[is.na(i)]<-match(k[is.na(i)],METHODS$label);if(anyNA(i))stop('unknown method ',paste(k[is.na(i)],collapse=', '),call.=FALSE);METHODS[i,]}
METHOD_LAB<-setNames(METHODS$label,METHODS$key);METHOD_FILL<-setNames(METHODS$fill,METHODS$label)
METHOD_COL<-setNames(METHODS$col,METHODS$label);METHOD_SHAPE<-setNames(METHODS$shape,METHODS$label)
# blood states (cool family) and brain states (warm family); routes = blood x brain
BLOOD_COL<-c(cDC='#2C5F8A',cMono='#4E9BD0',nMono='#9CC8E6',pDC='#2E9C8E')
BRAIN_COL<-c(Astro='#E09B2D',`MG-AP`='#B4532F')
BLOOD_KEY<-c(cDC='cDC',classical_mono='cMono',nonclassical_mono='nMono',pDC='pDC')          # result-table keys -> labels
BRAIN_KEY<-c(astro='Astro',microglia_mhc2='MG-AP')
ROUTES<-as.vector(outer(names(BLOOD_COL),names(BRAIN_COL),paste,sep=' > '))                 # 8 routes, blood-major
route_label<-function(key){p<-strsplit(key,'__x__',fixed=TRUE);vapply(p,function(z)paste(BLOOD_KEY[[z[1]]]%||%z[1],BRAIN_KEY[[z[2]]]%||%z[2],sep=' > '),'')}
# final targets (expanded drug definition; model direction from confirmation seeds 47-51)
TARGET_COL<-c(CD71='#2E8B57',CD119='#7B3F9E');TARGET_GENE<-c(CD71='TFRC',CD119='IFNGR1')
target_label<-function(a)paste0(a,'/',TARGET_GENE[a])
# observation type (inner ring), cohorts (shades of their tissue colour), diagnosis
TISSUE_COL<-c(`Blood cells`='#4F7CA8',`Brain nuclei`='#D9A441',`Spatial spots`='#8F9D5A')
COHORT_COL<-c(GSE223138='#3D6894',MB2025='#8DB0D2',GSE178265='#B9852B',GSE184950='#DDB462',GSE157783='#EED7A1',GSE253975='#8F9D5A')
DX_COL<-c(Control='#A6B4C2',`Early PD`='#D9B8D0',PD='#B07AA1',`Late PD`='#8A4F7D',PDD='#663A5E',LBD='#5F8A80',MSA='#8FB0A8',PSP='#C3D6D1')
MOD_COL<-c(Volume='#758AA2',QSM='#3A8B84',`T2*`='#C48A53')                                    # MR imaging modalities (pd-v3 standard)
PLATFORM_COL<-c(`UKB-PPP`='#4D4D4D',AGES='#A6A6A6')
GATE_COL<-c(pass='#3C3C3C',fail='#DADADA',pending='white')
# continuous scales
seq_col<-function(v,col,lo=0,hi=1,na=MISSING){r<-grDevices::colorRampPalette(c('#FFFFFF',col))(120)[16:120];x<-(v-lo)/(hi-lo)
 ifelse(is.na(v),na,r[1+round(pmin(pmax(x,0),1)*(length(r)-1))])}
DIV_LO<-'#2F5D8C';DIV_HI<-'#B8433A'
div_col<-function(v,lim=1){f<-grDevices::colorRamp(c(DIV_LO,'#FFFFFF',DIV_HI),space='Lab');x<-(pmin(pmax(v,-lim),lim)+lim)/(2*lim)
 ifelse(is.na(v),MISSING,grDevices::rgb(f(ifelse(is.na(x),.5,x)),maxColorValue=255))}
dens_col<-function(n,nmax,lo='#F0F0F0',hi='#3F3F3F'){f<-grDevices::colorRamp(c(lo,hi),space='Lab');x<-ifelse(n<=0,0,log1p(n)/log1p(nmax))
 grDevices::rgb(f(pmin(pmax(x,0),1)),maxColorValue=255)}
scale_fill_div<-function(lim=1,breaks=c(-lim,0,lim),guide='none',...)scale_fill_gradientn(colours=div_col(seq(-lim,lim,length.out=65),lim),
 limits=c(-lim,lim),breaks=breaks,labels=num_minus,na.value=MISSING,oob=scales::squish,guide=guide,...)
scale_fill_seq<-function(col,lim=c(0,1),breaks=lim,guide='none',...)scale_fill_gradientn(colours=seq_col(seq(0,1,length.out=65),col),
 limits=lim,breaks=breaks,labels=num_minus,na.value=MISSING,oob=scales::squish,guide=guide,...)

# ================================================================== numbers (typographic minus, separators)
fnum<-function(x,d=2)sub('^-','−',formatC(x,format='f',digits=d))
fsign<-function(x,d=2)ifelse(x>=0,paste0('+',formatC(x,format='f',digits=d)),fnum(x,d))
fint<-function(x)sub('^-','−',formatC(round(x),format='d',big.mark=','))
fpct<-function(p,d=1)paste0(formatC(100*p,format='f',digits=d),'%')
ffold<-function(x,d=2)paste0(formatC(x,format='f',digits=d),'×')
num_minus<-scales::label_number(style_negative='minus',big.mark='')
num_minus_acc<-function(acc)scales::label_number(accuracy=acc,style_negative='minus',big.mark='')
num_comma<-scales::label_number(big.mark=',',style_negative='minus')

# ================================================================== ggplot theme (7 pt everywhere, bold axis titles, no titles)
theme_fig<-function(legend='none'){
 theme_classic(base_family=FONT,base_size=PT)+theme(
  text=element_text(family=FONT,size=PT,colour=INK),line=element_line(colour=INK,linewidth=gglw('axis')),
  axis.text=element_text(size=PT,colour=INK),axis.title=element_text(size=PT,colour=INK,face='bold'),
  axis.title.x=element_text(margin=margin(t=1,unit='mm')),axis.title.y=element_text(margin=margin(r=1,unit='mm')),
  axis.text.x=element_text(margin=margin(t=.5,unit='mm')),axis.text.y=element_text(margin=margin(r=.5,unit='mm')),
  axis.line=element_line(linewidth=gglw('axis'),colour=INK,lineend='square'),axis.ticks=element_line(linewidth=gglw('tick'),colour=INK),
  axis.ticks.length=unit(TICK_MM,'mm'),panel.grid=element_blank(),panel.background=element_blank(),plot.background=element_blank(),
  panel.border=element_blank(),panel.spacing=unit(2,'mm'),
  strip.background=element_blank(),strip.text=element_text(size=PT,face='plain',colour=INK,margin=margin(0,0,.8,0,'mm')),strip.placement='outside',
  legend.position=legend,legend.title=element_blank(),legend.text=element_text(size=PT,colour=INK),legend.key=element_blank(),
  legend.key.size=unit(2.2,'mm'),legend.key.spacing.x=unit(1.5,'mm'),legend.key.spacing.y=unit(.4,'mm'),legend.background=element_blank(),
  legend.box.background=element_blank(),legend.margin=margin(0,0,0,0),legend.box.spacing=unit(1,'mm'),
  plot.title=element_blank(),plot.subtitle=element_blank(),plot.caption=element_blank(),plot.tag=element_blank(),
  plot.margin=margin(LETTER_ZONE[['h']],.6,.6,.6,'mm'))}
BASE<-theme_fig()
legend_inside<-function(x,y,just=c(0,1))theme(legend.position='inside',legend.position.inside=c(x,y),legend.justification=just,legend.direction='vertical')
guide_bar<-function(horizontal=TRUE,len=12)guide_colourbar(barwidth=if(horizontal)unit(len,'mm')else unit(1.6,'mm'),barheight=if(horizontal)unit(1.6,'mm')else unit(len,'mm'),
 ticks.colour=INK,ticks.linewidth=gglw('tick'),frame.colour=INK,frame.linewidth=gglw('tick'),title.position='top',direction=if(horizontal)'horizontal'else'vertical')

# ================================================================== data access (read-only; every file logged)
.SRC<-new.env()
src_reset<-function()rm(list=ls(.SRC),envir=.SRC)
src_log<-function(path){fp<-normalizePath(path,mustWork=TRUE);i<-file.info(fp)
 assign(fp,list(path=fp,bytes=unname(i$size),mtime=format(i$mtime,'%Y-%m-%dT%H:%M:%S'),md5=if(i$size<2e8)unname(tools::md5sum(fp))else NA_character_),envir=.SRC);fp}
src_list<-function()unname(mget(sort(ls(.SRC)),envir=.SRC))
read_src<-function(path,n=NULL,cols=NULL,...){fp<-src_log(path)
 d<-if(grepl('\\.json$',fp))jsonlite::fromJSON(fp,simplifyVector=TRUE)
  else if(grepl('\\.(tsv|txt)(\\.gz)?$',fp))utils::read.delim(fp,check.names=FALSE,stringsAsFactors=FALSE,...)
  else utils::read.csv(fp,check.names=FALSE,stringsAsFactors=FALSE,...)
 if(is.data.frame(d)){for(k in names(d))if(is.character(d[[k]])&&length(stats::na.omit(d[[k]]))&&all(stats::na.omit(d[[k]])%in%c('True','False')))d[[k]]<-d[[k]]=='True'
  if(!is.null(n)&&nrow(d)!=n)stop(sprintf('read_src %s: %d rows, expected %d',basename(fp),nrow(d),n),call.=FALSE)
  if(!is.null(cols)){m<-setdiff(cols,names(d));if(length(m))stop('read_src ',basename(fp),': missing columns ',paste(m,collapse=', '),call.=FALSE)}}
 d}

# ================================================================== bar-chart registry (rule: no bar chart with 1-3 bars)
.BARS<-new.env();bars_reset<-function()rm(list=ls(.BARS),envir=.BARS)
note_bars<-function(n,what='bars'){k<-sprintf('%03d',length(ls(.BARS))+1);assign(k,list(what=what,n_bars=as.integer(n)),envir=.BARS);invisible(n)}
bars_list<-function()unname(mget(sort(ls(.BARS)),envir=.BARS))
gg_bar_audit<-function(p){if(!inherits(p,'ggplot'))return(invisible())
 for(j in seq_along(p$layers)){ly<-p$layers[[j]];g<-class(ly$geom)[1]
  if(g%in%c('GeomBar','GeomCol','GeomRect')){ld<-ggplot2::layer_data(p,j);fl<-isTRUE(ld$flipped_aes[1])
   key<-if(g=='GeomRect')paste(ld$xmin,ld$xmax,ld$ymin,ld$ymax)else if(fl)ld$ymin else ld$xmin
   n<-tapply(key,ld$PANEL,function(z)length(unique(if(is.numeric(z))round(z,9) else z)))
   note_bars(min(n),sprintf('ggplot layer %d (%s), fewest bars in one facet',j,g))}}
 invisible()}

# ================================================================== text audit (grob level; the PDF level is in figv2_pdf.py)
text_words<-function(s){tk<-strsplit(trimws(s),'[[:space:]]+')[[1]];sum(grepl('[A-Za-z].*[A-Za-z]',tk))}
TEXT_RULES<-c(sentence='[a-z]{2,}[.!?]( |$)|;',caption=':$|^(Note|Notes|n =|Source)\\b',hyphen_minus='(^|[\\s(\\[])-[0-9.]')
text_audit<-function(g,letters_ok=FALSE){out<-list()
 walk<-function(x,fs,cex,fam,face){if(is.null(x)||inherits(x,'zeroGrob'))return()
  upd<-function(gp){if(is.null(gp))return();if(!is.null(gp$fontsize))fs<<-gp$fontsize;if(!is.null(gp$cex))cex<<-cex*gp$cex
   if(!is.null(gp$fontfamily))fam<<-gp$fontfamily;if(!is.null(gp$fontface))face<<-gp$fontface}
  if(inherits(x$vp,'viewport'))upd(x$vp$gp);upd(x$gp)
  if(inherits(x,'text')){lab<-x$label;if(is.expression(lab)||is.call(lab))lab<-vapply(as.list(lab),function(e)paste(deparse(e),collapse=''),'')
   lab<-as.character(lab);if(length(lab))out[[length(out)+1]]<<-data.frame(label=lab,pt=rep_len(fs*cex,length(lab)),family=fam%||%'',face=as.character(face%||%1)[1],stringsAsFactors=FALSE)}
  if(inherits(x,'gtable'))for(ch in x$grobs)walk(ch,fs,cex,fam,face)else if(!is.null(x$children))for(ch in x$children)walk(ch,fs,cex,fam,face)}
 walk(g,PT,1,FONT,1);d<-if(length(out))do.call(rbind,out)else data.frame(label=character(0),pt=numeric(0),family=character(0),face=character(0))
 d<-d[!is.na(d$label)&nzchar(trimws(d$label)),,drop=FALSE];d$words<-vapply(d$label,text_words,0);d$flag<-rep('',nrow(d))
 for(k in names(TEXT_RULES)){i<-grepl(TEXT_RULES[[k]],d$label,perl=TRUE);d$flag[i]<-paste0('FAIL ',k)}
 i<-d$words>5;d$flag[i]<-'FAIL more than five words'
 i<-!(d$family%in%c(FONT,'sans',''));d$flag[i]<-paste0('FAIL font ',d$family[i])
 isl<-letters_ok&abs(d$pt-PT_LETTER)<1e-6&grepl('^[a-z]$',d$label)&d$face%in%c('2','bold')
 i<-abs(d$pt-PT)>1e-6&!isl;d$flag[i]<-sprintf('FAIL %.2f pt',d$pt[i]);rownames(d)<-NULL;d}
fig_guard<-function(p,id='?'){if(!inherits(p,'ggplot'))return(invisible(TRUE))
 l<-tryCatch(p@labels,error=function(e)p$labels);bad<-c('title','subtitle','caption','tag')[vapply(c('title','subtitle','caption','tag'),function(k){v<-l[[k]];!is.null(v)&&!(is.character(v)&&all(!nzchar(trimws(v))))},TRUE)]
 if(length(bad))stop(sprintf('panel %s has a %s; panels carry no titles or captions',id,paste(bad,collapse='/')),call.=FALSE);invisible(TRUE)}

# ================================================================== measuring (text widths in mm on a scratch device)
meas_open<-function(){if(dev.cur()==1){ragg::agg_capture(width=100,height=100,units='mm',res=300);TRUE}else FALSE}
sw<-function(t,bold=FALSE){op<-meas_open();on.exit(if(op)dev.off())
 vapply(t,function(s)convertWidth(grobWidth(textGrob(s,gp=gpar(fontfamily=FONT,fontsize=PT,fontface=if(bold)'bold'else'plain'))),'mm',valueOnly=TRUE),0,USE.NAMES=FALSE)}
gt<-function(p){op<-meas_open();on.exit(if(op)dev.off());if(inherits(p,'ggplot'))ggplotGrob(p)else p}
gg_grob<-function(p,id=get0('PANEL',envir=globalenv())$id%||%'?'){fig_guard(p,id);gg_bar_audit(p);gt(p)}   # every ggplot goes through here (guard + bar audit)

# ================================================================== grid drawing helpers (mm from the panel top-left)
.H<-function()get0('PANEL_H',envir=globalenv(),ifnotfound=NA_real_)
.W<-function()get0('PANEL_W',envir=globalenv(),ifnotfound=NA_real_)
.gp<-function(fill=NA,col=INK,lwd=glwd('data'),lty=1,...)gpar(fill=fill,col=col,lwd=lwd,lty=lty,...)
dtxt<-function(t,x,y,bold=FALSE,col=INK,just='centre',h=.H(),rot=0,italic=FALSE)textGrob(t,x=unit(x,'mm'),y=unit(h-y,'mm'),just=just,rot=rot,
 gp=gpar(fontfamily=FONT,fontsize=PT,col=col,lineheight=1,fontface=if(bold&&italic)'bold.italic'else if(bold)'bold'else if(italic)'italic'else'plain'))
drect<-function(x,y,w,ht,fill='white',col=INK,lwd=glwd('frame'),h=.H(),lty=1,r=0)if(r>0)roundrectGrob(x=unit(x,'mm'),y=unit(h-y,'mm'),width=unit(w,'mm'),height=unit(ht,'mm'),just=c('left','top'),r=unit(r,'mm'),gp=.gp(fill,col,lwd,lty))else
 rectGrob(x=unit(x,'mm'),y=unit(h-y,'mm'),width=unit(w,'mm'),height=unit(ht,'mm'),just=c('left','top'),gp=.gp(fill,col,lwd,lty))
dline<-function(x,y,x2,y2,arrow=FALSE,h=.H(),lty=1,col=INK,lwd=glwd('data'),alen=1)segmentsGrob(unit(x,'mm'),unit(h-y,'mm'),unit(x2,'mm'),unit(h-y2,'mm'),
 gp=gpar(col=col,lwd=lwd,lty=lty,fill=col,lineend='butt'),arrow=if(arrow)grid::arrow(length=unit(alen,'mm'),type='closed')else NULL)
dpath<-function(xs,ys,col=INK,lwd=glwd('data'),lty=1,h=.H())linesGrob(unit(xs,'mm'),unit(h-ys,'mm'),gp=gpar(col=col,lwd=lwd,lty=lty))
dpoly<-function(xs,ys,fill,col=NA,lwd=glwd('frame'),lty=1,h=.H())polygonGrob(unit(xs,'mm'),unit(h-ys,'mm'),gp=.gp(fill,col,lwd,lty))
dcirc<-function(x,y,d,fill,col=INK,lwd=glwd('thin'),h=.H())circleGrob(unit(x,'mm'),unit(h-y,'mm'),r=unit(d/2,'mm'),gp=.gp(fill,col,lwd))
ddiam<-function(x,y,d,fill,col=INK,lwd=glwd('thin'),h=.H()){r<-d/2*1.15;dpoly(x+c(0,r,0,-r),y+c(-r,0,r,0),fill,col,lwd,h=h)}
dtri<-function(x,y,s,fill,col=NA,lwd=glwd('thin'),dir='up',h=.H()){a<-s*.5;b<-s*.45;switch(dir,up=dpoly(x+c(-b,b,0),y+c(a,a,-a),fill,col,lwd,h=h),down=dpoly(x+c(-b,b,0),y+c(-a,-a,a),fill,col,lwd,h=h),
 left=dpoly(x+c(a,a,-a),y+c(-b,b,0),fill,col,lwd,h=h),right=dpoly(x+c(-a,-a,a),y+c(-b,b,0),fill,col,lwd,h=h))}
dcurve<-function(x,y,x2,y2,curv=.3,col=INK,lwd=glwd('data'),lty=1,arrow=TRUE,h=.H())curveGrob(unit(x,'mm'),unit(h-y,'mm'),unit(x2,'mm'),unit(h-y2,'mm'),curvature=curv,square=FALSE,ncp=8,
 gp=gpar(col=col,lwd=lwd,lty=lty,fill=col),arrow=if(arrow)grid::arrow(length=unit(1,'mm'),type='closed')else NULL)
dmark<-function(x,y,method,d=D[['point']],h=.H()){m<-method_row(method);s<-m$shape   # method glyph (grid) matching the ggplot shapes
 if(s==4){r<-d*.42;return(list(dline(x-r,y-r,x+r,y+r,col=m$col,lwd=glwd('data'),h=h),dline(x-r,y+r,x+r,y-r,col=m$col,lwd=glwd('data'),h=h)))}
 switch(as.character(s),`21`=list(dcirc(x,y,d,m$fill,m$col,h=h)),`22`=list(drect(x-d*.44,y-d*.44,d*.88,d*.88,m$fill,m$col,glwd('thin'),h=h)),
  `23`=list(ddiam(x,y,d*.95,m$fill,m$col,h=h)),`24`=list(dtri(x,y,d*1.1,m$fill,m$col,dir='up',h=h)),`25`=list(dtri(x,y,d*1.1,m$fill,m$col,dir='down',h=h)))}
dgroup<-function(g,areas=NULL){gt<-gTree(children=do.call(gList,g));if(!is.null(areas))attr(gt,'areas')<-areas;gt}

# ---- mini-plot scales and axes (grid panels)
sc_lin<-function(v,d0,d1,m0,m1)m0+(v-d0)/(d1-d0)*(m1-m0)
sc_log<-function(v,d0,d1,m0,m1)m0+(log(v)-log(d0))/(log(d1)-log(d0))*(m1-m0)
panel_scales<-function(x0,x1,y0,y1,xlim,ylim)list(sx=function(v)sc_lin(v,xlim[1],xlim[2],x0,x1),sy=function(v)sc_lin(v,ylim[1],ylim[2],y1,y0))   # y0 = top, y1 = bottom
axis_b<-function(x0,x1,y,at,lab=NULL,tick=TICK_MM,line=TRUE,gap=.5,h=.H()){g<-if(line)list(dline(x0,y,x1,y,h=h,lwd=glwd('axis')))else list()
 for(i in seq_along(at)){g<-c(g,list(dline(at[i],y,at[i],y+tick,h=h,lwd=glwd('tick'))));if(!is.null(lab)&&nzchar(lab[i]))g<-c(g,list(dtxt(lab[i],at[i],y+tick+gap,just=c('centre','top'),h=h)))};g}
axis_l<-function(y0,y1,x,at,lab=NULL,tick=TICK_MM,line=TRUE,gap=.5,h=.H()){g<-if(line)list(dline(x,y0,x,y1,h=h,lwd=glwd('axis')))else list()
 for(i in seq_along(at)){g<-c(g,list(dline(x-tick,at[i],x,at[i],h=h,lwd=glwd('tick'))));if(!is.null(lab)&&nzchar(lab[i]))g<-c(g,list(dtxt(lab[i],x-tick-gap,at[i],just=c('right','centre'),h=h)))};g}
dband<-function(x,lo,hi,sx,sy,fill,h=.H()){k<-is.finite(x)&is.finite(lo)&is.finite(hi);if(sum(k)<2)return(list());x<-x[k];lo<-lo[k];hi<-hi[k]
 list(dpoly(sx(c(x,rev(x))),sy(c(hi,rev(lo))),fill,col=NA,h=h))}
dcurve_xy<-function(x,y,sx,sy,col,lwd=glwd('data'),lty=1,h=.H())list(dpath(sx(x),sy(y),col=col,lwd=lwd,lty=lty,h=h))
dhist<-function(lo,hi,count,sx,y_base,scale,dir='up',fill=MID,col=NA,outline=FALSE,lwd=glwd('thin'),h=.H()){k<-count>0   # bins from stored edges
 if(!outline){if(!any(k))return(list());return(list(rectGrob(unit(sx(lo[k]),'mm'),unit(h-y_base,'mm'),width=unit(sx(hi[k])-sx(lo[k]),'mm'),height=unit(count[k]*scale,'mm'),
  just=c('left',if(dir=='up')'bottom'else'top'),gp=gpar(fill=fill,col=col,lwd=lwd))))}
 s<-if(dir=='up')-1 else 1;o<-order(lo);lo<-lo[o];hi<-hi[o];cnt<-count[o];xs<-as.vector(rbind(sx(lo),sx(hi)));ys<-y_base+s*as.vector(rbind(cnt,cnt))*scale
 list(dpath(c(sx(lo[1]),xs,sx(hi[length(hi)])),c(y_base,ys,y_base),col=if(is.na(col))MID else col,lwd=lwd,h=h))}
# bar chart drawn in grid: registers its bar count (rule 3); x0/x1 = value axis, y positions = bar centres
dbars<-function(v,yc,sx,ht,fill,base=0,h=.H(),what='grid bars'){note_bars(length(v),what)
 lapply(seq_along(v),function(i)drect(min(sx(base),sx(v[i])),yc[i]-ht/2,abs(sx(v[i])-sx(base)),ht,fill[(i-1)%%length(fill)+1],col=NA,h=h))}

# ---- hatch (pending or unfavourable fills)
HATCH_BG<-'#F2F2F2';HATCH_LINE<-'#8C8C8C';HATCH_PITCH<-.8
.vseg<-function(x0,y0,x1,y1,col,lwd,lty=1,h)segmentsGrob(unit(x0,'mm'),unit(h-y0,'mm'),unit(x1,'mm'),unit(h-y1,'mm'),gp=gpar(col=col,lwd=lwd,lty=lty,lineend='butt'))
hatch_lines<-function(x,y,w,ht,pitch=HATCH_PITCH,col=HATCH_LINE,lwd=glwd('hair'),h=.H()){cc<-seq(-ht+pitch/2,w,by=pitch);v0<-pmax(0,-cc);v1<-pmin(ht,w-cc);k<-v1>v0
 if(!any(k))return(list());list(.vseg(x+v0[k]+cc[k],y+ht-v0[k],x+v1[k]+cc[k],y+ht-v1[k],col,lwd,h=h))}
hatch_rect<-function(x,y,w,ht,col=MID,pitch=HATCH_PITCH,h=.H())c(list(drect(x,y,w,ht,HATCH_BG,col=NA,h=h)),hatch_lines(x,y,w,ht,pitch,h=h),list(drect(x,y,w,ht,NA,col=col,h=h)))

# ---- rings (donut + outer tick ring + ring icon), ribbons (alluvial), waffles (unit ledgers)
dring<-function(n,cols,cx,cy,r1,r2,h=.H(),border='white',lwd=glwd('thin'),start=0){s<-start;tot<-sum(n);out<-list()   # 12 o'clock, clockwise
 for(i in seq_along(n)){if(n[i]<=0)next;e<-s+n[i]/tot*2*pi;t<-seq(s,e,length.out=max(6,ceiling(160*n[i]/tot)))
  out<-c(out,list(polygonGrob(unit(c(cx+r2*sin(t),cx+r1*sin(rev(t))),'mm'),unit(h-c(cy-r2*cos(t),cy-r1*cos(rev(t))),'mm'),gp=gpar(fill=cols[i],col=border,lwd=lwd))));s<-e};out}
ring_ticks<-function(cx,cy,r,len=.8,step=.1,major=.25,h=.H()){a<-seq(0,1-1e-9,by=step)*2*pi;mj<-abs((seq(0,1-1e-9,by=step)/major)-round(seq(0,1-1e-9,by=step)/major))<1e-9
 l<-ifelse(mj,len*1.6,len);list(segmentsGrob(unit(cx+r*sin(a),'mm'),unit(h-(cy-r*cos(a)),'mm'),unit(cx+(r+l)*sin(a),'mm'),unit(h-(cy-(r+l)*cos(a)),'mm'),gp=gpar(col=INK,lwd=glwd('tick'))))}
ring_icon<-function(cx,cy,which,d=2.6,col=INK,h=.H()){r<-d/2*c(.42,.71,1);g<-list()   # 3 concentric rings, ring `which` (1 inner .. 3 outer) solid
 for(k in 3:1)g<-c(g,list(dcirc(cx,cy,2*r[k],if(k==which)col else 'white',col=if(k==which)col else LIGHT,lwd=glwd('thin'),h=h)));g}
dribbon<-function(x0,y0a,y0b,x1,y1a,y1b,fill,col=NA,n=40,h=.H()){t<-seq(0,1,length.out=n);s<-3*t^2-2*t^3;xs<-x0+(x1-x0)*t
 dpoly(c(xs,rev(xs)),c(y0a+(y1a-y0a)*s,rev(y0b+(y1b-y0b)*s)),fill,col,h=h)}
waffle<-function(n_vec,cols,x,y,pitch=1.3,nrow=3,align='left',sq=pitch*.82,col=NA,h=.H()){stopifnot(length(n_vec)==length(cols),all(n_vec>=0))   # column-major unit squares
 n<-sum(n_vec);if(n==0){g<-list();attr(g,'w')<-0;attr(g,'ncol')<-0;return(g)};fill<-rep(cols,n_vec);nc<-ceiling(n/nrow);i<-seq_len(n)-1
 if(align=='left'){j<-i%/%nrow;r<-i%%nrow;x0<-x}else{n1<-n-(nc-1)*nrow;j<-ifelse(i<n1,0,1+(i-n1)%/%nrow);r<-ifelse(i<n1,i,(i-n1)%%nrow);x0<-x-nc*pitch}
 xc<-x0+j*pitch+pitch/2;yc<-y+r*pitch+pitch/2
 g<-list(rectGrob(unit(xc,'mm'),unit(h-yc,'mm'),width=unit(sq,'mm'),height=unit(sq,'mm'),gp=gpar(fill=fill,col=col,lwd=glwd('hair'))));attr(g,'w')<-nc*pitch;attr(g,'ncol')<-nc;g}
waffle_runs<-function(colseq){r<-rle(as.character(colseq));list(n=r$lengths,cols=r$values)}

# ---- bootstrap / seed draws as a density strip + CI rail (one row of a forest)
draw_strip<-function(d,lo,hi,x0,x1,yc,ht=2.2,nb=60,col=PERISCOPE_LIGHT,col_neg=LIGHT,h=.H()){d<-d[is.finite(d)];stopifnot(length(d)>0)
 sx<-function(v)sc_lin(v,lo,hi,x0,x1);br<-seq(lo,hi,length.out=nb+1);inw<-d>=lo&d<=hi
 cnt<-graphics::hist(d[inw],breaks=br,plot=FALSE,include.lowest=TRUE,right=FALSE)$counts;k<-cnt>0;g<-list()
 if(any(k)){hh<-ht*cnt[k]/max(cnt);xl<-sx(br[-length(br)][k]);xr<-sx(br[-1][k]);mid<-(br[-length(br)][k]+br[-1][k])/2
  g<-list(rectGrob(unit(xl,'mm'),unit(h-yc,'mm'),width=unit(xr-xl,'mm'),height=unit(hh,'mm'),just=c('left','centre'),gp=gpar(fill=ifelse(mid<0,col_neg,col),col=NA)))}
 nlo<-sum(d<lo);nhi<-sum(d>hi);if(nlo)g<-c(g,list(dtri(x0+.5,yc,1,MID,dir='left',h=h)));if(nhi)g<-c(g,list(dtri(x1-.5,yc,1,MID,dir='right',h=h)))
 attr(g,'n_out')<-c(lo=nlo,hi=nhi);g}
ci_rail<-function(est,lo,hi,dlo,dhi,x0,x1,yc,col=INK,fill=PERISCOPE,d=D[['mean']],h=.H()){sx<-function(v)sc_lin(pmin(pmax(v,dlo),dhi),dlo,dhi,x0,x1)
 list(.vseg(sx(lo),yc,sx(hi),yc,col,glwd('data'),h=h),dcirc(sx(est),yc,d,fill,col='white',lwd=glwd('thin'),h=h))}

# ---- hexagonal binning (display only; pointy-top cells of flat-to-flat width bw)
hex_bins<-function(x,y,bw){dy<-bw*sqrt(3)/2;i1<-round(x/bw);j1<-round(y/(2*dy));i2<-round(x/bw-.5);j2<-round(y/(2*dy)-.5)
 c1x<-i1*bw;c1y<-j1*2*dy;c2x<-(i2+.5)*bw;c2y<-(j2+.5)*2*dy;use1<-(x-c1x)^2+(y-c1y)^2<=(x-c2x)^2+(y-c2y)^2
 cx<-ifelse(use1,c1x,c2x);cy<-ifelse(use1,c1y,c2y);key<-paste(round(cx/bw*2),round(cy/dy));t<-table(key);u<-!duplicated(key)
 d<-data.frame(key=key[u],cx=cx[u],cy=cy[u]);d$n<-as.integer(t[d$key]);d$id<-seq_len(nrow(d));stopifnot(sum(d$n)==length(x));d[,c('id','cx','cy','n')]}
hex_polys<-function(bins,bw,aspect=1){R<-bw/sqrt(3);a<-pi/6+(0:5)*pi/3
 data.frame(id=rep(bins$id,each=6),n=rep(bins$n,each=6),x=rep(bins$cx,each=6)+R*cos(a),y=rep(bins$cy,each=6)+R*sin(a)*aspect)}

# ---- unit placement
swarm_offsets<-function(y,dy,max_lanes=12){o<-order(y);off<-numeric(length(y));placed<-integer(0);cand<-c(0,rbind(seq_len(max_lanes),-seq_len(max_lanes)))
 for(i in o){for(c in cand)if(!any(abs(off[placed]-c)<1&abs(y[placed]-y[i])<dy)){off[i]<-c;break};placed<-c(placed,i)};off}
spread_labels<-function(y,gap){o<-order(y);z<-y[o];for(i in seq_along(z)[-1])if(z[i]-z[i-1]<gap)z[i]<-z[i-1]+gap;z<-z-(mean(z)-mean(y[o]))
 for(i in seq_along(z)[-1])if(z[i]-z[i-1]<gap)z[i]<-z[i-1]+gap;r<-y;r[o]<-z;r}

# ---- keys: a glyph = list(w = width mm, f = function(x, yc) grob), yc from the viewport bottom
kg_dot<-function(fill=MID,col='white',d=D[['point']])list(w=d,f=function(x,yc)circleGrob(unit(x+d/2,'mm'),unit(yc,'mm'),r=unit(d/2,'mm'),gp=gpar(fill=fill,col=col,lwd=glwd('thin'))))
kg_box<-function(fill=PAPER,col=NA,w=2.2,ht=2.2,lty=1)list(w=w,f=function(x,yc)rectGrob(unit(x,'mm'),unit(yc,'mm'),width=unit(w,'mm'),height=unit(ht,'mm'),just=c('left','centre'),gp=gpar(fill=fill,col=col,lty=lty,lwd=glwd('thin'))))
kg_line<-function(col=INK,lty=1,lwd=glwd('data'),w=4)list(w=w,f=function(x,yc)segmentsGrob(unit(x,'mm'),unit(yc,'mm'),unit(x+w,'mm'),unit(yc,'mm'),gp=gpar(col=col,lty=lty,lwd=lwd)))
kg_method<-function(method,d=D[['point']])list(w=d*1.1,f=function(x,yc)gTree(children=do.call(gList,dmark(x+d*.55,0,method,d,h=yc))))
kg_ring<-function(which,d=2.6)list(w=d,f=function(x,yc)gTree(children=do.call(gList,ring_icon(x+d/2,0,which,d,h=yc))))
kg_grad<-function(cols,lo,hi,w=10,ht=1.6){a<-sw(lo);b<-sw(hi);list(w=a+w+b+1.4,f=function(x,yc)gList(textGrob(lo,unit(x,'mm'),unit(yc,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,col=INK)),
 rasterGrob(matrix(cols,nrow=1),x=unit(x+a+.7,'mm'),y=unit(yc,'mm'),width=unit(w,'mm'),height=unit(ht,'mm'),just=c('left','centre'),interpolate=TRUE),
 rectGrob(unit(x+a+.7,'mm'),unit(yc,'mm'),width=unit(w,'mm'),height=unit(ht,'mm'),just=c('left','centre'),gp=gpar(fill=NA,col=INK,lwd=glwd('hair'))),
 textGrob(hi,unit(x+a+w+1.4,'mm'),unit(yc,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,col=INK))))}
ki<-function(glyph,lab='')list(g=glyph,lab=lab)
key_row<-function(items,x0,y,gap=1.8,h=.H()){yc<-h-y;x<-x0;out<-list()   # one row of keys from x0 at y (mm from panel top); attr 'end'
 for(it in items){out<-c(out,list(it$g$f(x,yc)));x<-x+it$g$w
  if(nzchar(it$lab)){out<-c(out,list(textGrob(it$lab,unit(x+.6,'mm'),unit(yc,'mm'),just=c('left','centre'),gp=gpar(fontfamily=FONT,fontsize=PT,col=INK))));x<-x+.6+sw(it$lab)}
  x<-x+gap};attr(out,'end')<-x-gap;out}
key_width<-function(items,gap=1.8)sum(vapply(items,function(it)it$g$w+if(nzchar(it$lab)).6+sw(it$lab)else 0,0))+gap*(length(items)-1)

# ---- ggplot -> gtable helpers (sizes in mm)
.mmw<-function(u)convertWidth(u,'mm',valueOnly=TRUE);.mmh<-function(u)convertHeight(u,'mm',valueOnly=TRUE)
panel_rc<-function(g){pl<-g$layout[grepl('^panel',g$layout$name),];list(c=sort(unique(pl$l)),r=sort(unique(pl$t)))}
nonpanel_mm<-function(g){op<-meas_open();on.exit(if(op)dev.off());rc<-panel_rc(g)
 c(w=sum(.mmw(g$widths[-rc$c])),h=sum(.mmh(g$heights[-rc$r])),nc=length(rc$c),nr=length(rc$r))}
fix_panels<-function(g,w,h,hjust=.5,vjust=.5){rc<-panel_rc(g);for(c in rc$c)g$widths[c]<-unit(w,'mm');for(r in rc$r)g$heights[r]<-unit(h,'mm');g$respect<-FALSE
 g<-gtable_add_cols(g,unit(hjust,'null'),0);g<-gtable_add_cols(g,unit(1-hjust,'null'),-1);g<-gtable_add_rows(g,unit(1-vjust,'null'),0);gtable_add_rows(g,unit(vjust,'null'),-1)}
fit_panels<-function(g,W=.W(),H=.H(),square=FALSE,hjust=.5,vjust=.5,max_mm=Inf){np<-nonpanel_mm(g)   # largest facets that fit the box
 pw<-(W-np[['w']])/np[['nc']];ph<-(H-np[['h']])/np[['nr']];if(square)pw<-ph<-min(pw,ph);pw<-floor(min(pw,max_mm)*10)/10;ph<-floor(min(ph,max_mm)*10)/10
 stopifnot(pw>3,ph>3);g<-fix_panels(g,pw,ph,hjust,vjust);attr(g,'panel_mm')<-c(w=pw,h=ph);g}
place<-function(g,x,y,w,ht,h=.H())gTree(children=gList(g),vp=viewport(x=unit(x,'mm'),y=unit(h-y,'mm'),width=unit(w,'mm'),height=unit(ht,'mm'),just=c('left','top')))   # a grob in a sub-box (mm)
strip_cols<-function(g,cols,pat='^strip-t'){i<-grep(pat,g$layout$name);i<-i[order(g$layout$t[i],g$layout$l[i])];stopifnot(length(i)==length(cols))
 rec<-function(gr,col){if(inherits(gr,'text')){gr$gp<-modifyList(gr$gp%||%gpar(),list(col=col));return(gr)};if(inherits(gr,'gtable')){gr$grobs<-lapply(gr$grobs,rec,col=col);return(gr)}
  if(length(gr$children))for(n in names(gr$children))gr$children[[n]]<-rec(gr$children[[n]],col);gr}
 for(j in seq_along(i))g$grobs[[i[j]]]<-rec(g$grobs[[i[j]]],cols[j]);g}
add_to_panel<-function(g,panel,grob,name=paste0('add_',panel,'_',sample.int(1e6,1))){lay<-g$layout[g$layout$name==panel,];stopifnot(nrow(lay)==1)
 gtable_add_grob(g,grob,t=lay$t,l=lay$l,b=lay$b,r=lay$r,z=Inf,clip='off',name=name)}

# ================================================================== devices (exact block size)
# cairo_pdf truncates the page to whole points, so the device is opened at ceil(pt) + 0.25 and the panel is drawn in an
# exact-size viewport anchored top-left; figv2_pdf.py crop then sets the MediaBox to the block (no scaling anywhere).
.dev_in<-function(mm)(ceiling(mm/25.4*72)+.25)/72
.push_block<-function(w,h)pushViewport(viewport(x=0,y=unit(1,'npc'),width=unit(w,'mm'),height=unit(h,'mm'),just=c('left','top'),
 gp=gpar(fontfamily=FONT,fontsize=PT,col=INK,lineheight=1),name='block'))
open_pdf_mm<-function(file,w,h,bg='white'){grDevices::cairo_pdf(file,width=.dev_in(w),height=.dev_in(h),family=FONT,bg=bg);.push_block(w,h)}
open_png_mm<-function(file,w,h,res=600){ragg::agg_png(file,width=w,height=h,units='mm',res=res,background='white');.push_block(w,h)}
open_block_pdf<-function(file,b){b<-block(b);open_pdf_mm(file,b$w,b$h)}
open_block_png<-function(file,b,res=600){b<-block(b);open_png_mm(file,b$w,b$h,res)}
pdf_crop<-function(raw,final,w,h){st<-system2(Sys.getenv('FIGV2_PY'),c(shQuote(file.path(SYS_DIR,'figv2_pdf.py')),'crop',shQuote(raw),shQuote(final),w,h),stdout=TRUE,stderr=TRUE)
 if(!is.null(attr(st,'status'))&&attr(st,'status')!=0)stop('pdf crop failed: ',paste(st,collapse='\n'),call.=FALSE);invisible(final)}

# ================================================================== blank (pending) block: thin outline, nothing else
blank_grob<-function(w=.W(),h=.H())rectGrob(x=unit(LW[['frame']]/2*25.4/72,'mm'),y=unit(LW[['frame']]/2*25.4/72,'mm'),
 width=unit(w-LW[['frame']]*25.4/72,'mm'),height=unit(h-LW[['frame']]*25.4/72,'mm'),just=c('left','bottom'),gp=gpar(fill=NA,col=FRAME_COL,lwd=glwd('frame')))
