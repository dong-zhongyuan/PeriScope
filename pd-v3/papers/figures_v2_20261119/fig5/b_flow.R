# Figure 5b - program-level flow: every program (left, by the last model gate C1-C5 it passes on its own) flows to the
# last gate its antigen passes (right); ribbon width = programs. Data: antigen_gate_matrix.csv (n_programs,
# n_prog_gate1toK, gate columns).
PANEL<-list(id='5b_program_flow',blocks=c('1x2'),what='Program flow from model gates to antigen outcome')
source(file.path(Sys.getenv('FIGV2_ROOT'),'fig5/_fun.R'))
slots<-function(v,avail,gap,mh){f<-function(k)sum(pmax(v*k,mh))+gap*(length(v)-1)-avail;k<-stats::uniroot(f,c(1e-12,avail/sum(v)),tol=1e-14)$root;list(k=k,h=pmax(v*k,mh))}
make_panel<-function(b){W<-b$w;H<-b$h;assign('PANEL_H',H,envir=globalenv())
 a<-antigen_stage();np<-cbind(a$n_programs,a$n_prog_gate1to1,a$n_prog_gate1to2,a$n_prog_gate1to3,a$n_prog_gate1to4,a$n_prog_gate1to5);np[is.na(np)]<-0
 st<-np[,1:5]-np[,2:6];st<-cbind(st,np[,6]);colnames(st)<-c('None','C1','C2','C3','C4','C5')
 rord<-c('None','C1','C2','C3','C4','C5','MR','Coloc','Direction','Final');a$rlab<-rord[a$out];stopifnot(!anyNA(a$rlab))
F<-do.call(rbind,lapply(colnames(st),function(l)aggregate(list(n=st[,l]),list(r=a$rlab),sum)|>transform(l=l)));F<-F[F$n>0,]
 lt<-tapply(F$n,F$l,sum)[colnames(st)];lt<-lt[!is.na(lt)];rt<-tapply(F$n,F$r,sum);rt<-rt[intersect(rord,names(rt))]
 top<-LETTER_ZONE[['h']]+.3;bot<-.5;gap<-.5;mh<-LINE_MM*1.05;bw<-1.8
 x0<-max(sw(names(lt)))+1.4;x1<-W-bw-1.4-max(sw(names(rt)))
 L<-slots(lt,H-top-bot,gap,mh);Rr<-slots(rt,H-top-bot,gap,mh);k<-min(L$k,Rr$k)
 pos<-function(h)top+c(0,cumsum(h+gap))[seq_along(h)]
 yL<-setNames(pos(L$h),names(lt));hL<-setNames(L$h,names(lt));yR<-setNames(pos(Rr$h),names(rt));hR<-setNames(Rr$h,names(rt))
 cL<-yL+hL/2-lt*k/2;cR<-yR+hR/2-rt*k/2;oL<-cL;oR<-cR;g<-list()
 rc<-setNames(c(FAINT,rep(LIGHT,5),'#9A9A9A','#7A7A7A','#5C5C5C',PERISCOPE),rord)
 F<-F[order(match(F$l,names(lt)),match(F$r,names(rt))),]
 for(i in seq_len(nrow(F))){h<-F$n[i]*k;l<-F$l[i];r<-F$r[i];a0<-oL[l];oL[l]<-a0+h;b0<-oR[r];oR[r]<-b0+h
  g<-c(g,list(dribbon(x0+bw,a0,a0+h,x1,b0,b0+h,fill=alpha_hex(rc[r],.55))))}
 for(l in names(lt))g<-c(g,list(drect(x0,cL[l],bw,lt[l]*k,fill=MID,col=NA),dtxt(l,x0-.7,yL[l]+hL[l]/2,just=c('right','centre'))))
 for(r in names(rt))g<-c(g,list(drect(x1,cR[r],bw,rt[r]*k,fill=rc[r],col=NA),dtxt(r,x1+bw+.7,yR[r]+hR[r]/2,just=c('left','centre'),bold=r=='Final',col=if(r=='Final')PERISCOPE else INK)))
 dgroup(g)}
