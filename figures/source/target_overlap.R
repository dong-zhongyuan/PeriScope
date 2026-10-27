# New UpSet view: exact intersections of priority-program genes, not association edges.
target_overlap<-function(){
 d<-rd('priority_gene_membership');targets<-c('CD22','CD40','CD71')
 d$key<-apply(d[,targets],1,function(z)paste(targets[as.logical(z)],collapse=' + '))
 z<-as.data.frame(table(d$key));names(z)<-c('key','count');z<-z[order(-z$count),];z$i<-seq_len(nrow(z));z$kind<-ifelse(grepl(' [+] ',z$key),'Shared','Target-specific')
 q<-ggplot(z,aes(i,count,fill=kind))+geom_col(width=.66,color='#596C6C',linewidth=.3)+geom_text(aes(label=count),vjust=-.5,size=12/ggplot2::.pt)+scale_fill_manual(values=c(Shared='#A38AB6',`Target-specific`='#70A897'),name=NULL)+scale_x_continuous(limits=c(.4,nrow(z)+.6),expand=c(0,0))+scale_y_continuous(expand=expansion(mult=c(0,.16)))+labs(x=NULL,y='Program genes')+common+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position='top',plot.margin=margin(1,8,0,8))
 dots<-expand.grid(i=z$i,target=targets);dots$on<-mapply(function(i,t)t%in%strsplit(as.character(z$key[i]),' + ',fixed=TRUE)[[1]],dots$i,dots$target);dots$y<-match(dots$target,rev(targets))
 links<-do.call(rbind,lapply(split(dots[dots$on,],dots$i[dots$on]),function(a)data.frame(i=a$i[1],lo=min(a$y),hi=max(a$y))))
 mat<-ggplot(dots,aes(i,y))+geom_hline(yintercept=1:3,color='#F0F2F1',linewidth=7)+geom_segment(data=links,aes(x=i,xend=i,y=lo,yend=hi),inherit.aes=FALSE,linewidth=.85,color='#445D60')+geom_point(aes(fill=on),shape=21,size=3.5,color='white',stroke=.3)+scale_fill_manual(values=c(`TRUE`='#445D60',`FALSE`='#DDE4E1'),guide='none')+scale_x_continuous(limits=c(.4,nrow(z)+.6),expand=c(0,0))+scale_y_continuous(breaks=1:3,labels=rev(targets),limits=c(.5,3.5),expand=c(0,0))+theme_void(base_size=12)+theme(axis.text.y=element_text(size=12,color='#202020',face='bold'),plot.margin=margin(3,8,5,8))
 plot_grid(q,mat,ncol=1,rel_heights=c(.72,.28),align='v',axis='lr')
}
