#!/usr/bin/env Rscript



### Script to extract data from Seurat object that can be used 
# to construct an equivalent Anndata object



### User defined variables

# unpack variables passed from parent shell script
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

for (e in args) {
  argname <- e[1]
  argval <- e[2]
  # regular expression to delete initial \" and trailing \"
  argval <- gsub("(^\\\"|\\\"$)", "", argval)
  assign(argname, argval)
}



### Packages

if (!require("Seurat")) install.packages("Seurat")
if (!require("patchwork")) install.packages("patchwork")
if (!require("Matrix")) install.packages("Matrix")


### Functions

## function to export data from Seurat object that can be used to generate Anndata object
velo_input <- function(seurat_obj, out_dir, project_name) {
  
  # create output directory if not existing
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  # save metadata
  seurat_obj$barcode <- colnames(seurat_obj)
  seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
  seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
  
  write.csv(seurat_obj[[]], file = paste(out_dir, "/velocity_metadata_", project_name, ".csv", sep=""), 
            quote=FALSE, row.names = FALSE)
  
  # save count matrix
  count_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
  writeMM(count_matrix, file=paste(out_dir, "/velocity_counts_", project_name, ".mtx", sep=""))
  
  # write gene names
  write.table(data.frame("gene"=rownames(count_matrix)), 
              file=paste(out_dir, "/velocity_genes_", project_name, ".csv", sep=""),
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  # write dimensionality reduction matrix
  write.csv(seurat_obj@reductions$pca@cell.embeddings, file=paste(out_dir, "/velocity_pca_", 
            project_name, ".csv", sep=""), quote=FALSE, row.names=FALSE)
}



### Analysis

# read in Seurat object
seurat_obj <- readRDS(in_file)


# export data for velocity analysis
velo_input(seurat_obj=seurat_obj, out_dir=out_dir, project_name=project_name)






