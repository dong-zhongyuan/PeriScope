library(Matrix)
base <- '/public/home/mengxl/dzy/pd_product_assets/raw/ma_published_reference_20261003'
out <- '/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/spatial_reference_inventory_20261003'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
obj <- readRDS(file.path(base,'processed_sn_obj.rds'))
cat('CLASS',class(obj),'\nSLOTS',names(attributes(obj)),'\n')
md <- attr(obj,'meta.data')
stopifnot(is.data.frame(md),nrow(md)>0,!anyDuplicated(rownames(md)))
write.csv(md,gzfile(file.path(out,'author_cell_metadata.csv.gz')))
cat('METADATA',paste(colnames(md),collapse=', '),'\n')
for(c in intersect(c('Donor','donor','broad_lineages','broad_lineage','sublineage','Condition','condition','gender'),colnames(md))) {
 cat(c,'\n');print(table(md[[c]],useNA='ifany'))
}
assays <- attr(obj,'assays')
for(a in names(assays)) {
 cat('ASSAY',a,'CLASS',class(assays[[a]]),'SLOTS',names(attributes(assays[[a]])),'\n')
 for(s in intersect(c('counts','data','scale.data','layers'),names(attributes(assays[[a]])))) {
  x<-attr(assays[[a]],s);cat('SLOT',s,'CLASS',class(x),'DIM',dim(x),'\n')
  if(s=='layers')print(lapply(x,dim))
 }
}
saveRDS(md,file.path(out,'metadata.rds'))
if(!is.null(md$broad_lineages)&&!is.null(md$sublineage)) write.csv(as.data.frame(table(md$broad_lineages,md$sublineage)),file.path(out,'author_lineage_crosswalk.csv'),row.names=FALSE)
writeLines('complete',file.path(out,'inventory_complete.txt'))
