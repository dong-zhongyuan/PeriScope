work <- '/public/home/mengxl/dzy/pd_product/_work/spatial_subtype_20261003'
lib <- file.path(work,'R_library')
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
Sys.setenv(PATH=paste('/public/home/mengxl/dzy/envs/rcoloc/bin','/public/home/mengxl/dzy/envs/af3/bin',Sys.getenv('PATH'),sep=':'))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=300,Ncpus=6)
Sys.setenv(https_proxy='http://127.0.0.1:18792',http_proxy='http://127.0.0.1:18792')
desc <- read.dcf(file.path(work,'spacexr','DESCRIPTION'))
needed <- trimws(strsplit(desc[1,'Imports'],',')[[1]])
needed <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(needed)) install.packages(needed,lib=lib)
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop(paste('Missing dependencies:',paste(missing,collapse=', ')))
install.packages(file.path(work,'spacexr'),repos=NULL,type='source',lib=lib)
library(spacexr)
writeLines(capture.output(sessionInfo()),file.path(work,'R_sessionInfo.txt'))
cat('SPACEXR_READY\n')
