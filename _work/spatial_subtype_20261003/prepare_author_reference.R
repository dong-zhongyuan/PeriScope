library(Matrix)
set.seed(20261003)
P <- '/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003'
W <- '/public/home/mengxl/dzy/pd_product/_work/spatial_subtype_20261003'
out <- file.path(P,'spatial_subtype_author_20261003'); inp<-file.path(out,'inputs')
dir.create(inp,recursive=TRUE,showWarnings=FALSE)
stopifnot(!dir.exists(file.path(out,'samples')))
base<-file.path(P,'spatial_subtype_20261003','inputs')
genes<-readLines(file.path(base,'genes.txt'))
obj<-readRDS('/public/home/mengxl/dzy/pd_product_assets/raw/ma_published_reference_20261003/processed_sn_obj.rds')
m<-attr(obj,'meta.data'); x<-attr(attr(obj,'assays')[['RNA']],'counts');rm(obj);gc()
stopifnot(identical(colnames(x),rownames(m)),all(x@x==round(x@x)))
cat('GENES',head(rownames(x)), '\n')
full_nUMI<-Matrix::colSums(x)
h<-read.delim('/public/home/mengxl/dzy/pd_product/_work/optimization_20261003/reference/hgnc_complete_set.txt',stringsAsFactors=FALSE,quote='',fill=TRUE)
rename<-list()
for(i in seq_len(nrow(h))){
 if(is.na(h$prev_symbol[i]))next
 for(old in strsplit(h$prev_symbol[i],'|',fixed=TRUE)[[1]])if(nzchar(old)&&!old%in%h$symbol)rename[[old]]<-unique(c(rename[[old]],h$symbol[i]))
}
rename<-rename[lengths(rename)==1];sy<-rownames(x);ii<-which(sy%in%names(rename));sy[ii]<-unlist(rename[sy[ii]])
keep<-which(sy%in%genes)
stopifnot(length(unique(sy[keep]))>20000)
collapse<-sparseMatrix(i=match(sy[keep],genes),j=keep,x=1,dims=c(length(genes),nrow(x)))
x<-collapse%*%x;rownames(x)<-genes
write.csv(data.frame(gene=genes,in_author_reference=genes%in%sy),file.path(out,'gene_alignment.csv'),row.names=FALSE)
m$state<-as.character(m$broad_lineages)
map<-c(Baseline='microglia_homeostatic',Activated='microglia_mhc2',Monocyte_Derived='BAM',Fibrous='astro_fibrous',Protoplasmic='astro_protoplasmic',Dopaminergic='neuron_dopaminergic',GABEergic='neuron_GABAergic',MEG3_Interneurons='neuron_MEG3')
i<-which(as.character(m$sublineage)%in%names(map));m$state[i]<-unname(map[as.character(m$sublineage[i])])
cross<-read.csv(file.path(W,'ma_reference_sources/GEO_author_donor_crosswalk.csv'))
cross$donor_id<-sub(', snRNA','',cross$title,fixed=TRUE)
sm<-read.csv(file.path(base,'spatial_metadata.csv'))
norm<-function(s)sub('^T-?0*','T',s)
cross$spatial_overlap<-norm(cross$donor_id)%in%norm(sm$donor)
cross$unresolved_id<-!grepl('^T-',cross$donor_id)
# One donor per condition/sex held out, fixed from inventory before outcomes.
qc_donors<-c('OJ11','OJ13','OJ9','OJ22')
allowed<-cross$author_donor[!cross$spatial_overlap & !cross$unresolved_id]
train<-setdiff(intersect(allowed,unique(m$Donor)),qc_donors)
cross$role<-ifelse(cross$spatial_overlap,'excluded_spatial_overlap',ifelse(cross$unresolved_id,'excluded_unresolved_identity',ifelse(cross$author_donor%in%qc_donors,'QC_holdout','reference')))
cross$role[!cross$author_donor%in%m$Donor]<-'absent_from_author_object'
write.csv(cross,file.path(out,'donor_roles.csv'),row.names=FALSE)
m$nUMI<-full_nUMI;m$donor<-m$Donor;m$condition<-m$Condition
support<-as.data.frame(table(m$Donor,m$state));names(support)<-c('donor','state','n');write.csv(support,file.path(out,'reference_donor_state_inventory.csv'),row.names=FALSE)
types<-sort(unique(m$state));profiles<-list();selected<-integer()
for(st in types){
 means<-list()
 for(d in train){
  ix<-which(m$Donor==d & m$state==st & m$nUMI>0)
  if(length(ix)<10)next
  if(length(ix)>250)ix<-sample(ix,250)
  selected<-c(selected,ix)
  means[[d]]<-Matrix::rowMeans(t(t(x[,ix,drop=FALSE])/m$nUMI[ix]))
 }
 if(length(means)>=2)profiles[[st]]<-rowMeans(do.call(cbind,means))
}
profiles<-do.call(cbind,profiles);rownames(profiles)<-genes;types<-colnames(profiles)
selected<-selected[m$state[selected]%in%types]
write.csv(profiles,gzfile(file.path(inp,'reference_profiles.csv.gz')))
write.csv(m[selected,],gzfile(file.path(inp,'reference_selected_cells.csv.gz')))
# Marker choice uses reference identities only, not spatial disease contrasts.
exclude<-grepl('^(MT-|RPL|RPS|TTTY|RBMY|TSPY|DAZ|BPY2|HSFY|CDY|VCY)',genes)|genes%in%c('XIST','TSIX','UTY','KDM5D','DDX3Y','EIF1AY','USP9Y','ZFY','TMSB4Y','NLGN4Y','LINC00278','TXLNGY','PRKY','AMELY','PCDH11Y','SRY')
markers<-list()
for(st in types){
 other<-apply(profiles[,setdiff(types,st),drop=FALSE],1,max)
 fc<-log2((profiles[,st]+1e-6)/(other+1e-6))
 ix<-which(!exclude & profiles[,st]>=5e-5 & fc>0)
 ix<-head(ix[order(fc[ix],decreasing=TRUE)],40)
 markers[[st]]<-data.frame(gene=genes[ix],state=st,log2FC=fc[ix])
}
markers<-do.call(rbind,markers);write.csv(markers,file.path(out,'reference_marker_selection.csv'),row.names=FALSE)
canonical<-c('GFAP','AQP4','SLC1A3','ALDH1L1','P2RY12','CX3CR1','TMEM119','SALL1','P2RY13','GPR34','HLA-DRA','HLA-DRB1','CD74','HLA-DPA1','HLA-DPB1','C3','SPP1','MRC1','CD163','CD3D','CD3E','AIF1','TYROBP','LYZ','LST1','C1QA','C1QB','C1QC')
reserved<-sort(intersect(genes,unique(c(markers$gene,canonical))))
eval<-setdiff(readLines(file.path(base,'evaluation_genes.txt')),reserved)
writeLines(reserved,file.path(inp,'deconvolution_genes.txt'));writeLines(eval,file.path(inp,'evaluation_genes.txt'))
write.csv(profiles[intersect(canonical,genes),],file.path(out,'biological_marker_profiles.csv'))
write_mtx<-function(z,name){tmp<-file.path(inp,paste0(name,'.mtx'));writeMM(z,tmp);system2('gzip',c('-f',shQuote(tmp)))}
thin<-unlist(lapply(types,function(st)head(selected[m$state[selected]==st],30)))
write.csv(m[thin,],file.path(inp,'reference_metadata.csv'));write_mtx(x[,thin],'reference_counts')
for(f in c('genes.txt','spatial_metadata.csv','spatial_counts.mtx.gz'))file.symlink(file.path(base,f),file.path(inp,f))
# Holdout nuclei generate technical mixtures only, never biological outcomes.
cal<-which(m$Donor%in%qc_donors & m$state%in%types & m$nUMI>0)
cal_x<-x[,cal,drop=FALSE]
cal_position<-setNames(seq_along(cal),cal)
present<-types[sapply(types,function(st)sum(m$state[cal]==st)>=5)]
targets<-c('astro_fibrous','astro_protoplasmic','microglia_homeostatic','microglia_mhc2')
truth<-matrix(0,400,length(types),dimnames=list(sprintf('QC_%04d',0:399),types));qc<-matrix(0,length(genes),400)
for(k in 1:400){
 if(k%%50==0)cat('QC_MIXTURE',k,'\n')
 picked<-sample(present,min(5,length(present)))
 target<-targets[(k-1)%%4+1];if(k<=300 && target%in%present && !target%in%picked)picked[1]<-target
 frac<-rgamma(length(picked),.7);frac<-frac/sum(frac);prob<-numeric(length(genes))
 for(j in seq_along(picked)){
  ix<-sample(cal[m$state[cal]==picked[j]],5)
  v<-Matrix::rowMeans(t(t(cal_x[,cal_position[as.character(ix)],drop=FALSE])/m$nUMI[ix]));v<-v/sum(v)
  prob<-prob+frac[j]*v;truth[k,picked[j]]<-frac[j]
 }
 qc[,k]<-rmultinom(1,1500,prob)
}
colnames(qc)<-rownames(truth);rownames(qc)<-genes
write_mtx(as(qc,'sparseMatrix'),'qc_counts');write.csv(truth,file.path(inp,'qc_truth.csv'))
qmd<-data.frame(x=(0:399)%%20,y=(0:399)%/%20,nUMI=1500,condition='QC_synthetic_only',donor='heldout_mixture',row.names=rownames(truth));write.csv(qmd,file.path(inp,'qc_metadata.csv'))
writeLines(c(paste('train',paste(train,collapse=',')),paste('holdout',paste(qc_donors,collapse=',')),paste('types',paste(types,collapse=',')),paste('markers',length(reserved)),paste('evaluation',length(eval))),file.path(out,'input_report.txt'))
writeLines('complete',file.path(out,'input_complete.txt'))
