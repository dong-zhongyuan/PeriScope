#!/usr/bin/env Rscript
# Preserve the published pipeline's source layer: author RNA@data, not RNA@counts.
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: Rscript export_mb2025.R INPUT.rds OUTPUT_DIRECTORY')
if (dir.exists(args[2])) stop('Output directory already exists: ', args[2])
x <- readRDS(args[1])
rna <- x@assays$RNA
if (!'data' %in% slotNames(rna)) stop('Expected the frozen RNA Assay data slot')
m <- rna@data
stopifnot(identical(dim(m), c(36601L, 105081L)))
barcodes <- unname(colnames(m))
stopifnot(length(barcodes) == ncol(m), !anyDuplicated(barcodes),
          all(barcodes %in% rownames(x@meta.data)))
meta <- x@meta.data[barcodes, , drop = FALSE]
stopifnot(identical(rownames(meta), barcodes))
dir.create(args[2], recursive = TRUE)
writeMM(m, file.path(args[2], 'matrix.mtx'))
writeLines(rownames(m), file.path(args[2], 'features.tsv'))
writeLines(barcodes, file.path(args[2], 'barcodes.tsv'))
write.csv(meta, file.path(args[2], 'metadata.csv'), row.names = TRUE)
print(list(matrix = dim(m), seurat_object = as.character(packageVersion('SeuratObject')),
           matrix_package = as.character(packageVersion('Matrix'))))
