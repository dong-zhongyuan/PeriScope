# Figure 1b - resource composition: nested rings (inner observation type, middle cohort, outer diagnosis within cohort),
# total in the hole, 10% tick ring, keyed by ring icons with each category's share. Data: original 1B inputs (read-only).
# Layouts: 2x1 (tall: rings on top, keys stacked below) and 1x2 (wide: rings left, two key stacks right). The ring radius
# and key arrangement are solved from the block, never scaled; every text is 7 pt. Key columns take their own widths.
PANEL<-list(id='1b_resources',blocks=c('2x1','1x2'),what='Resource composition rings')
DESIGN<-'/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design'

resource_data<-function(){
 comp<-read_src(file.path(DESIGN,'full_template_reuse/inputs/all_cohort_composition.csv'),cols=c('donor','condition','cell_type','n','cohort'))
 rings<-read_src(file.path(DESIGN,'template_semantic_review/data/input_1A.csv'),cols=c('category','count','ring'))
 cohort<-c(gse223138='GSE223138',moquin_beaudry2025='MB2025',gse178265='GSE178265',gse184950='GSE184950',gse157783='GSE157783',gse253975='GSE253975')
 tissue<-c(GSE223138='Blood cells',MB2025='Blood cells',GSE178265='Brain nuclei',GSE184950='Brain nuclei',GSE157783='Brain nuclei',GSE253975='Spatial spots')
 dx<-c(Control='Control',Ctrl='Control',CTRL='Control',Normal='Control',PD='PD',`Early PD`='Early PD',`Late PD`='Late PD',PDD='PDD',MSA='MSA',PSP='PSP',Disease='Disease')
 stopifnot(all(comp$cohort%in%names(cohort)),all(comp$condition%in%names(dx)))
 comp$cohort<-cohort[comp$cohort];comp$dx<-dx[comp$condition]
 a<-aggregate(n~cohort+dx,comp,sum)
 # GSE178265 'Disease' = PD + LBD; LBD occurs only in GSE178265, so its LBD count is the ring total for LBD
 lbd<-rings$count[rings$category=='condition: LBD'];k<-a$cohort=='GSE178265'&a$dx=='Disease';stopifnot(sum(k)==1,length(lbd)==1)
 a<-rbind(a[!k,],data.frame(cohort='GSE178265',dx=c('PD','LBD'),n=c(a$n[k]-lbd,lbd)))
 a$cohort<-factor(a$cohort,names(COHORT_COL));a$dx<-factor(a$dx,names(DX_COL));a<-a[order(a$cohort,a$dx),]
 a$tissue<-factor(tissue[as.character(a$cohort)],names(TISSUE_COL))
 # the rebuilt rings must reproduce the original ring totals
 chk<-function(lev,ring,pre){got<-tapply(a$n,a[[lev]],sum);want<-setNames(rings$count[rings$ring==ring],sub(pre,'',rings$category[rings$ring==ring]))
  names(want)<-sub('^gse','GSE',names(want));names(want)[names(want)=='Blood']<-'Blood cells';stopifnot(all(abs(got[names(want)]-want)<.5))}
 chk('tissue','tissue','^tissue: ');chk('dx','condition','^condition: ');chk('cohort','cohort','^cohort: ')
 a}

key_groups<-function(a){tot<-sum(a$n)
 grp<-function(v,cols,ring){s<-tapply(a$n,v,sum);s<-s[!is.na(s)&s>0];list(ring=ring,lab=names(s),pct=rep('',length(s)),col=unname(cols[names(s)]))}
 gc<-grp(a$cohort,COHORT_COL,2);gc$lab<-unname(clab(gc$lab))
 list(grp(a$tissue,TISSUE_COL,1),gc,grp(a$dx,DX_COL,3))}

# key geometry: a group = ring icon row + entries in ncol sub-columns; a stack = groups one above another
SW<-2.0;ICON<-2.6
col_ws<-function(g,ncol){nr<-ceiling(length(g$lab)/ncol);vapply(seq_len(ncol)-1,function(cc){i<-(cc*nr+1):min((cc+1)*nr,length(g$lab));SW+.6+max(sw(g$lab[i]))+.4},0)}
group_w<-function(g,ncol)sum(col_ws(g,ncol))+(ncol-1)*1.4
group_rows<-function(g,ncol)1+ceiling(length(g$lab)/ncol)
draw_group<-function(g,ncol,x,y,pitch){ws<-col_ws(g,ncol);out<-ring_icon(x+ICON/2,y+pitch/2,g$ring,ICON)
 nr<-ceiling(length(g$lab)/ncol);labw<-max(sw(g$lab))
 for(i in seq_along(g$lab)){cc<-(i-1)%/%nr;rr<-(i-1)%%nr;xx<-x+sum(head(ws,cc))+cc*1.4;ew<-ws[cc+1];yy<-y+(rr+1.5)*pitch
  out<-c(out,list(drect(xx,yy-SW/2,SW,SW,g$col[i],col=NA),dtxt(g$lab[i],xx+SW+.6,yy,just=c('left','centre'))))};out}
# stacks: list of list(g = group index, ncol); returns grobs filling box (x0, y0, w, h); pitch solved from the height
draw_keys<-function(G,stacks,x0,y0,w,h,gap_grp=1.4,gap_stack=2.4,max_pitch=3.1){
 sw_<-vapply(stacks,function(s)max(vapply(s,function(e)group_w(G[[e$g]],e$ncol),0)),0)
 if(sum(sw_)+gap_stack*(length(stacks)-1)>w+1e-6)stop(sprintf('1b keys need %.1f mm, box has %.1f mm',sum(sw_)+gap_stack*(length(stacks)-1),w),call.=FALSE)
 rows<-vapply(stacks,function(s)sum(vapply(s,function(e)group_rows(G[[e$g]],e$ncol),0)),0)
 gaps<-vapply(stacks,function(s)(length(s)-1)*gap_grp,0)
 pitch<-min(max_pitch,min((h-gaps)/rows));if(pitch<1.15*LINE_MM)stop(sprintf('1b keys: row pitch %.2f mm is below one line',pitch),call.=FALSE)
 extra<-(w-sum(sw_)-gap_stack*(length(stacks)-1))/max(length(stacks)-1,1)   # spread stacks across the box
 out<-list();x<-x0
 for(k in seq_along(stacks)){y<-y0;for(e in stacks[[k]]){out<-c(out,draw_group(G[[e$g]],e$ncol,x,y,pitch));y<-y+group_rows(G[[e$g]],e$ncol)*pitch+gap_grp}
  x<-x+sw_[k]+gap_stack+if(length(stacks)>1)extra else 0}
 attr(out,'h')<-max(rows*pitch+gaps);out}

draw_rings<-function(a,cx,cy,R){r<-R*c(.40,.60,.80,1);tot<-sum(a$n)
 tis<-tapply(a$n,a$tissue,sum);tis<-tis[!is.na(tis)];coh<-tapply(a$n,a$cohort,sum);coh<-coh[!is.na(coh)]
 c(dring(tis,TISSUE_COL[names(tis)],cx,cy,r[1],r[2]),dring(coh,COHORT_COL[names(coh)],cx,cy,r[2],r[3]),
   dring(a$n,DX_COL[as.character(a$dx)],cx,cy,r[3],r[4]),ring_ticks(cx,cy,R+.35,len=.7,step=.1,major=.25))}

make_panel<-function(b){a<-resource_data();G<-key_groups(a);W<-b$w;H<-b$h;pad<-1.0;tk<-1.5   # tk = tick ring allowance
 if(b$name=='2x1'){                              # tall: rings across the top, one key stack below (diagnosis in 2 columns)
  stacks<-list(list(list(g=1,ncol=1),list(g=2,ncol=1),list(g=3,ncol=2)))
  rows<-sum(group_rows(G[[1]],1),group_rows(G[[2]],1),group_rows(G[[3]],2));key_h<-rows*2.9+2*1.4
  R<-min((W-2*pad)/2-tk,(H-key_h-3-2*pad)/2-tk);cx<-W/2;cy<-pad+tk+R
  ky<-cy+R+tk+2.4;keys<-draw_keys(G,stacks,pad,ky,W-2*pad,H-ky-pad)
 }else if(b$name=='1x2'){                       # wide: rings at the left, two key stacks at the right
  stacks<-list(list(list(g=1,ncol=1),list(g=2,ncol=1)),list(list(g=3,ncol=1)))
  kw<-sum(vapply(stacks,function(s)max(vapply(s,function(e)group_w(G[[e$g]],e$ncol),0)),0))+2.4
  R<-min((H-2*pad)/2-tk,(W-kw-3-2*pad)/2-tk);cx<-pad+tk+R;cy<-H/2
  kx<-cx+R+tk+3;keys<-draw_keys(G,stacks,kx,pad+.3,W-kx-pad,H-2*pad-.3)
 }
 dgroup(c(draw_rings(a,cx,cy,R),keys))}
