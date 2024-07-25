#!/usr/bin/env Rscript


### Run this script to generate Seurat object, remove doublets + identify cluster markers



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

#library("patchwork", lib="/ceph/home/m/mweinber/R_packages")
if (!require("patchwork")) install.packages("patchwork")
if (!require("Seurat")) install.packages("Seurat")
if (!require("pheatmap")) install.packages("pheatmap")
#if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("viridis")) install.packages("viridis")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
#if (!require("tidyr")) install.packages("tidyr")
#if (!require("DESeq2")) install.packages("DESeq2")
if (!require("writexl")) install.packages("writexl")

if (!require(DoubletFinder)) remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
if (!require(harmony)) remotes::install_github("immunogenomics/harmony")



### Functions

## function to create Seurat object and perform quality control analysis
# data_dir: directory containing cellranger output barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
# metadata_file (optional): path to file with cell metadata, cell barcodes should be in column named "barcode"
create_seurat <- function(data_dir, out_dir, project_name, metadata_file=NULL, jack_straw=FALSE,
                          min_genes=200, min_cells=3, max_genes=2500, pct_mt_cutoff=5) {
  
  print("Creating data set")
  
  # read in the 10X data
  seurat_data <- Read10X(data.dir=data_dir)
  
  if (!is.null(metadata_file)) {
    # read in metadata created using cellranger bash script
    meta_data <- read.table(metadata_file, header=TRUE)
    rownames(meta_data) <- meta_data$barcode
    meta_data$barcode <- NULL
    print(head(meta_data))
    
    # create a new Seurat object
    seurat_obj <- CreateSeuratObject(counts = seurat_data, project = project_name, 
                                     min.cells = min_cells, min.features = min_genes, meta.data=meta_data)    
  } else {
    # create a new Seurat object
    seurat_obj <- CreateSeuratObject(counts = seurat_data, project = project_name, 
                                     min.cells = min_cells, min.features = min_genes)
  }
  
  # add project name to cell barcodes
  #seurat_obj <- RenameCells(seurat_obj, new.names = paste(names(Idents(seurat_obj)), 
  #                                                        "_", project_name, sep=""))
 
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  print(seurat_obj)
  print(head(seurat_obj[[]]))
  
  # Visualize QC metrics as a violin plot
  gg <- VlnPlot(object=seurat_obj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  ggsave(gg, filename=paste(out_dir, "/", project_name, "_qc_feature_violin.pdf", sep=""), 
         device="pdf", dpi=300, useDingbats=FALSE)
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  pdf(file = paste(out_dir, "/", project_name, "_qc_feature_scatter.pdf", sep=""))
    print(plot1 + plot2)
  dev.off()
  
  print("Subsetting and normalizing data set")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < pct_mt_cutoff)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seurat_obj), 10)
  plot1 <- VariableFeaturePlot(seurat_obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  pdf(file = paste(out_dir, "/", project_name, "_qc_variable_features.pdf", sep=""))
    print(plot2)
  dev.off()
  
  # scale data and generate PCA plot
  print("Scaling data")
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  
  print("Identifying principal components")
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  # Examine and visualize PCA results a few different ways
  print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
  
  # determine number of PCs to be used in downstream analysis
  # NOTE: This process can take a long time for big datasets, comment out for expediency. More
  # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
  # computation time
  if (jack_straw == TRUE) {
    print("Computing JackStraw score")
    seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
    seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:40)
    
    pdf(file = paste(out_dir, "/", project_name, "_jack_straw.pdf", sep=""))
    print(JackStrawPlot(seurat_obj, dims = 1:40))
    dev.off()
  }
  
  pdf(file = paste(out_dir, "/", project_name, "_elbow_plot.pdf", sep=""))
    print(ElbowPlot(seurat_obj, ndims = 40, reduction = "pca"))
  dev.off()
  
  return(seurat_obj)
}




## function to cluster Seurat object and run UMAP
# n_dims: number of principal components (dimensions) to use
# resolution: resolution of the clustering (lower resolution gives less clusters)
cluster_seurat <- function(seurat_obj, out_dir, project_name, n_dims=40, harmony_var=NULL, 
                           resolution=0.4) {
  
  print("Clustering cells")
  
  if (!is.null(harmony_var)) {
    print("Peforming harmonisation")
  
    pdf(file = paste(out_dir, "/", project_name, "_harmony_run.pdf", sep=""))
      options(repr.plot.height = 2.5, repr.plot.width = 6)
      seurat_obj <- seurat_obj %>% RunHarmony(harmony_var, plot_convergence = TRUE)
    dev.off()
  
    # compare PCs before and after harmonisation
    pdf(file = paste(out_dir, "/", project_name, "_PC_without_harmony.pdf", sep=""))
      options(repr.plot.height = 5, repr.plot.width = 12)
      p1 <- DimPlot(object = seurat_obj, reduction = "pca", pt.size = .1, group.by = harmony_var)
      p2 <- VlnPlot(object = seurat_obj, features = "PC_1", group.by = harmony_var, pt.size = .1)
      print(cowplot::plot_grid(p1,p2))
    dev.off()
    
    pdf(file = paste(out_dir, "/", project_name, "_PC_with_harmony.pdf", sep=""))
      options(repr.plot.height = 5, repr.plot.width = 12)
      p1 <- DimPlot(object = seurat_obj, reduction = "harmony", pt.size = .1, group.by = harmony_var)
      p2 <- VlnPlot(object = seurat_obj, features = "harmony_1", group.by = harmony_var, pt.size = .1)
      print(cowplot::plot_grid(p1,p2))
    dev.off()
  
    # plot elbow plot again
    pdf(file = paste(out_dir, "/", project_name, "_elbow_plot_harmony.pdf", sep=""))
      print(ElbowPlot(seurat_obj, ndims = 40, reduction = "pca"))
    dev.off()
    
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:n_dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
     
    # Run non-linear dimensional reduction (UMAP/tSNE)
    # If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
    set.seed(42)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", assay = "SCT", dims = 1:n_dims)
    
  } else {
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
    
    # Run non-linear dimensional reduction (UMAP/tSNE)
    # If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
    set.seed(42)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims)
  }
  
  return(seurat_obj)
}




## function to identify doublets
# metadata_var: variable in Seurat object metadata according to which object is split before
#               identifying doublets, eg. sample
# n_dims: number of PCs to use for doublet identification
# pK: the PC neighborhood size used to compute each cell's 
#     proportion of artificial k nearest neighbors (pANN), 
#     expressed as a proportion of the merged real-artificial data
# dfr: doublet formation rate (3% for appr. 6000 cells loaded)
doublet_finder <- function(seurat_obj, out_dir, project_name, metadata_var, 
                           n_dims=20, pK=0.09, dfr=0.03) {
  
  Idents(seurat_obj) <- seurat_obj[[metadata_var]]
  
  # initiate dataframe for doublet results
  doublets <- data.frame(cell_id=str, metadata_var=str, seurat_clusters=str)
  
  for (i in unique(seurat_obj[[metadata_var]][[1]])) {
    
    print(paste("Identifying doublets in", metadata_var, i))
    
    # subset seurat object to cells annotated to metadata_var
    seurat_tmp <- subset(x = seurat_obj, idents = eval(i))
    
    # pK Identification (no ground-truth)
    #sweep.res.list <- paramSweep_v3(seurat_tmp, PCs = 1:n_dims, sct = FALSE)
    #sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    #bcmvn <- find.pK(sweep.stats)
    
    #pdf(file = paste(out_dir, "/", project_name, "_", i, "_doubletfinder_pk_bar.pdf", sep=""))
    #print(barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2))
    #dev.off()
    
    nExp_poi <- round(dfr*nrow(seurat_tmp[[]]))
    
    # Homotypic Doublet Proportion Estimate 
    #homotypic.prop <- modelHomotypic(seurat_tmp[["seurat_clusters"]])
    #nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    # Run DoubletFinder with varying classification stringencies
    seurat_tmp <- doubletFinder_v3(seurat_tmp, PCs = 1:n_dims, pN = 0.25, pK = pK, 
                                   nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    
    # name of the DF prediction can change, so extract the correct column name.
    DF.name = colnames(seurat_tmp[[]])[grepl("DF.classification", colnames(seurat_tmp[[]]))]
    
    pdf(file = paste(out_dir, "/", project_name, "_", i, "_doubletfinder_umap.pdf", sep=""))
    print(cowplot::plot_grid(ncol = 2, DimPlot(seurat_tmp, group.by = "seurat_clusters") + NoAxes(),
                             DimPlot(seurat_tmp, group.by = "tissue") + NoAxes(),
                             DimPlot(seurat_tmp, group.by = "condition") + NoAxes(), 
                             DimPlot(seurat_tmp, group.by = DF.name) + NoAxes()))
    dev.off()
    
    # extract IDs of doublets
    tmp <- as.data.frame(seurat_tmp[[]])
    tmp$cell_id <- rownames(tmp)
    tmp <- tmp[tmp[,DF.name]=="Doublet",c("cell_id", metadata_var, "seurat_clusters")]
    doublets <- rbind(doublets, tmp)
    
    print(paste("Identified", length(seurat_tmp[[DF.name]][seurat_tmp[[DF.name]]=="Doublet"]), 
                "doublets"))
  }
  
  return(doublets)
}




## function to plot UMAP
# idents_list: list of layers to overlay onto the UMAP clustering
# resolution: clustering resolution used, not the resolution of the saved images 
# if labels=TRUE, text labels will be displayed on UMAP and plot will be saved as pdf file
# if labels=FALSE, text labels will not be displayed and plot will be saved as png file
# file_type must be compatible with ggsave: "pdf", "png", etc.
umap_seurat <- function(seurat_obj, out_dir, project_name, idents_list, labels=TRUE, 
                        resolution=0.4, file_type="pdf") {
  
  for (ident in idents_list) {
    
    if (labels == TRUE) {    
      if (ident == "seurat_clusters") {
        resolution <- resolution * 100
        gg <- DimPlot(seurat_obj, reduction = "umap", group.by="seurat_clusters", 
                      label = TRUE, pt.size = 0.4) + NoLegend()
        ggsave(gg, filename=paste(out_dir, "/", project_name, "_umap_res_0_", resolution, ".", 
                                  file_type, sep=""), device=file_type, dpi=300)
      } else { 
        if (is.character(seurat_obj[[ident]][1])) {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident,
                        label = TRUE, pt.size = 0.4) + NoLegend()
          ggsave(gg, filename=paste(out_dir, "/", project_name, "_umap_", ident, ".", 
                                    file_type, sep=""), device=file_type, dpi=300)
        } else if (is.numeric(seurat_obj[[ident]][1])) {
          gg <- FeaturePlot(seurat_obj, reduction = "umap", features=ident,
                            label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=paste(out_dir, "/", project_name, "_umap_", ident, ".", 
                                    file_type, sep=""), device=file_type, dpi=300)
        }
      }
    } else if (labels == FALSE) {
      if (ident == "seurat_clusters") {
        resolution <- resolution * 100
        gg <- DimPlot(seurat_obj, reduction = "umap", group.by="seurat_clusters",
                      label = FALSE, pt.size = 0.4)
        ggsave(gg, filename=paste(out_dir, "/", project_name, "_umap_res_0_", resolution, 
                                  "_no_labels.", file_type, sep=""), device=file_type, dpi=300, width=11)
      } else { 
        if (is.character(seurat_obj[[ident]][1]))  {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident,
                        label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=paste(out_dir, "/", project_name, "_umap_", ident, "_no_labels.", 
                                    file_type, sep=""), device=file_type, dpi=300, width=11)
        } else if (is.numeric(seurat_obj[[ident]][1])) {
          gg <- FeaturePlot(seurat_obj, reduction = "umap", features=ident,
                            label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=paste(out_dir, "/", project_name, "_umap_", ident, "_no_labels.", 
                                    file_type, sep=""), device=file_type, dpi=300, width=11)
        }
      }
    }
  }

  return(seurat_obj)
}




## function to identify cluster markers
# idents indicates on which level marker genes should be identified
# n_markers indicate how many top markers per cluster to save in csv file
markers_seurat <- function(seurat_obj, out_dir, project_name, idents="seurat_clusters", 
                           n_markers=200) {
  
  Idents(seurat_obj) <- seurat_obj[[idents]]
  
  seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, 
                                       logfc.threshold = 0.25)
  
  # save the top markers of each cluster
  top_markers <- as.data.frame(seurat_obj.markers %>% group_by(cluster) %>% top_n(n = n_markers, wt = avg_log2FC))
  top_markers_df_list <- split(top_markers, f=top_markers$cluster)
  writexl::write_xlsx(top_markers_df_list, 
                      path=paste(out_dir, "/", project_name, "_markers.xlsx", sep=""),
                      col_names = TRUE,
                      format_headers = TRUE,
                      use_zip64 = FALSE)
  
  # draw heatmap of top marker expression
  top10 <- seurat_obj.markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
  
  gg <- DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
  ggsave(gg, filename=paste(out_dir, "/", project_name, "_heatmap_markers.pdf", sep=""), 
         device="pdf", dpi=300, useDingbats=FALSE)
  
  # create custom heatmap of top marker expression
  top10$cluster <- as.integer(top10$cluster)
  top10 <- top10[order(top10$cluster),]
  
  col_anno <- as.data.frame(seurat_obj[[idents]])
  colnames(col_anno) <- "cluster"
  col_anno$cluster <- as.integer(col_anno$cluster)
  col_anno$cell <- rownames(col_anno)
  col_anno <- col_anno[order(col_anno$cluster),]
  
  scaled_counts <- as.data.frame(seurat_obj@assays$RNA@scale.data)
  scaled_counts$gene <- rownames(scaled_counts)
  scaled_counts <- scaled_counts[scaled_counts$gene %in% unique(top10$gene),]
  
  scaled_counts <- merge(scaled_counts, top10[,c("gene","cluster")], by="gene")
  scaled_counts <- scaled_counts[order(scaled_counts$cluster),]
  rownames(scaled_counts) <- make.names(scaled_counts$gene, unique=TRUE)
  scaled_counts$gene <- NULL
  scaled_counts$cluster <- NULL
  
  scaled_counts <- as.data.frame(t(scaled_counts))
  scaled_counts$cell <- rownames(scaled_counts)
  scaled_counts <- merge(scaled_counts, col_anno, by="cell")
  scaled_counts <- scaled_counts[order(scaled_counts$cluster),]
  rownames(scaled_counts) <- scaled_counts$cell
  scaled_counts$cell <- NULL
  scaled_counts$cluster <- NULL
  scaled_counts <- as.data.frame(t(scaled_counts))
  
  # cap counts
  scaled_counts[scaled_counts > 2.5] <- 2.5
  scaled_counts[scaled_counts < -2.5] <- -2.5
  scaled_counts <- distinct(scaled_counts)
  col_anno$cell <- NULL
  
  # determine where to put column gaps
  cell_count <- col_anno %>% group_by(cluster) %>% summarise(count=n())
  gap_col <- cell_count$count[1]
  for (i in 2:(nrow(cell_count)-1)) {
    tmp <- cell_count$count[i] + gap_col[i-1]
    gap_col <- c(gap_col, tmp)
  }
  
  col_anno$cluster <- as.character(col_anno$cluster)
  
  pdf(file = paste(out_dir, "/", project_name, "_heatmap_markers_custom.pdf", sep=""))
  pheatmap(scaled_counts, 
    cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=FALSE,
    annotation_col = col_anno, gaps_col=gap_col,
    fontsize = 4, color = viridis(100))
  dev.off()
  
  return(seurat_obj)
}





### Analysis

min_genes <- as.numeric(min_genes)
max_genes <- as.numeric(max_genes)
max_perc_mt <- as.numeric(max_perc_mt)
min_cells <- as.numeric(min_cells)
n_pcs <- as.numeric(n_pcs)
leiden_res <- as.numeric(leiden_res)

# ensure missing harmony_var is received as NULL in function cluster_seurat
meta_data <- read.table(metadata_file, header=TRUE)

if (!harmony_var %in% colnames(meta_data)) {
  harmony_var = NULL
}
print(harmony_var)



# create output directory if it does not exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}



if (!file.exists(paste(out_dir, "/", project_name, "_scRNAseq_with_doublets.rds", sep=""))) {

  ## create Seurat object and pre-process
  seurat_obj <- create_seurat(data_dir=data_dir, out_dir=out_dir, project_name=project_name, 
                              metadata_file=metadata_file, jack_straw=FALSE,
                              min_genes=min_genes, min_cells=min_cells, max_genes=max_genes, 
                              pct_mt_cutoff=max_perc_mt)

  seurat_obj <- cluster_seurat(seurat_obj=seurat_obj, out_dir=out_dir, project_name=project_name, 
                               n_dims=n_pcs, resolution=leiden_res, harmony_var=NULL)
  print(seurat_obj)

  saveRDS(seurat_obj, paste(out_dir, "/", project_name, "_scRNAseq_with_doublets.rds", sep=""))
} else {
  seurat_obj <- readRDS(paste(out_dir, "/", project_name, "_scRNAseq_with_doublets.rds", sep=""))
}



if (!file.exists(paste(out_dir, "/", project_name, "_scRNAseq_analysed_no_doublets.rds", sep=""))) {

  ## remove doublets
  doublets <- doublet_finder(seurat_obj=seurat_obj, out_dir=out_dir, project_name=project_name, 
                             metadata_var="sample_id", n_dims=30, pK=0.09, dfr=0.03)
  write.csv(doublets, paste(out_dir, "/", project_name, "_doublets.csv", sep=""))

  seurat_obj[["cell_id"]] <- rownames(seurat_obj[[]])
  seurat_obj <- subset(x = seurat_obj, subset = cell_id %in% doublets$cell_id, invert=TRUE)
  
  
  ## cluster cells
  seurat_obj <- cluster_seurat(seurat_obj=seurat_obj, out_dir=out_dir, project_name=project_name, 
                               n_dims=n_pcs, resolution=leiden_res, harmony_var=harmony_var)

  # For some reason, seurat cluster levels are 0-based, but values are 1-based
  # -> adjust levels to be 1-based
  seurat_obj$seurat_clusters <- as.factor(as.numeric(seurat_obj$seurat_clusters)) 
  

  ## plot UMAPs with metadata overlay
  idents_list = colnames(seurat_obj[[]])
  idents_list = idents_list[idents_list != "barcode"]
  seurat_obj <- umap_seurat(seurat_obj=seurat_obj, out_dir=out_dir, project_name=project_name, 
                            resolution=leiden_res, idents_list=idents_list, labels=FALSE, 
                            file_type="png")
  
  
  ## Identify cluster markers
  seurat_obj <- markers_seurat(seurat_obj=seurat_obj, out_dir=out_dir, project_name=project_name, 
                               idents="seurat_clusters", n_markers=200)


  ## export metadata as csv
  meta_data_export <- seurat_obj[[]]
  write.csv(meta_data_export, paste(out_dir, "/", project_name, "_metadata.csv", sep=""))

  ## save Seurat object
  saveRDS(seurat_obj, paste(out_dir, "/", project_name, "_scRNAseq_analysed_no_doublets.rds", sep=""))
}




