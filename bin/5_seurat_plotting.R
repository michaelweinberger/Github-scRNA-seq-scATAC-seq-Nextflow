#!/usr/bin/env Rscript



### Run this script to add cell type annotation to Seurat object



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
if (!require("pheatmap")) install.packages("pheatmap")
#if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("viridis")) install.packages("viridis")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
#if (!require("tidyr")) install.packages("tidyr")
#if (!require("DESeq2")) install.packages("DESeq2")
#if (!require("xlsx")) install.packages("xlsx")
#if (!require("Matrix")) install.packages("Matrix")


# used for sctype
#if (!require("HGNChelper")) install.packages("HGNChelper")
#if (!require("ggraph")) install.packages("ggraph")
#if (!require("igraph")) install.packages("igraph")
#if (!require("tidyverse")) install.packages("tidyverse")
#if (!require("data.tree")) install.packages("data.tree")





### Functions

## function to set cell type names as label identities
# mode: should cell type identities be determined by using the
#  scType code (Ianevski et al. 2022, https://doi.org/10.1038/s41467-022-28803-w)
#  or a custom list of cell type identities
# annotation_df: dataframe with annotation, must contain column named "cluster" containing cell cluster IDs
#  (mode="custom")
# sctype_tissue: tissue cell types to query with scType, can be one of:
#  All, Immune system, Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart,
#  Intestine, Muscle, Placenta, Spleen, Stomach, Thymus
#  (mode="sctype")
# annotation_name: name under which the new annotation should be added (mode="sctype")
annotate_seurat <- function(seurat_obj, out_dir, project_name,  
                            mode=c("sctype","custom"), annotation_df=NULL, 
                            sctype_tissue="All", annotation_name) {
  
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
  
  # Add annotation
  if (mode == "custom") {
    seurat_meta <- as.data.frame(seurat_obj[[]])
    seurat_meta$seurat_clusters <- as.numeric(seurat_meta$seurat_clusters)
    tmp <- as.data.frame(seurat_meta %>% 
             left_join(annotation_df, by = join_by(seurat_clusters == cluster)))
    if (identical(names(seurat_obj$cell_id), tmp$cell_id)) {
      print("Annotation order matches Seurat metadata order")
    }
    rownames(tmp) <- tmp$cell_id
    seurat_obj@meta.data <- tmp
    seurat_obj$seurat_clusters <- as.factor(seurat_obj$seurat_clusters)
    Idents(seurat_obj) <- seurat_obj$cell_type
  
  } else if (mode == "sctype") {
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    
    # guess a tissue type
    # if scaled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. 
    # If you just created a Seurat object, without any scaling and normalization, 
    # set scaled = FALSE, seuratObject[[assay]]@counts will be used
    tissue_guess <- auto_detect_tissue_type(path_to_db_file = db_, seuratObject = seurat_obj, 
                                           scaled = TRUE, assay = "RNA")
    
    # plot tissue summary scores
    tissue_guess$tissue <- factor(tissue_guess$tissue, levels = rev(tissue_guess$tissue))
    y_max <- max(tissue_guess$score) + 100
    y_break_step <- round((y_max / 4) / 100) * 100
    gg <- ggplot(tissue_guess, aes(x=tissue,y=score)) +
      geom_bar(stat="identity", width=0.8, fill="#69b3a2") +
      coord_flip() +
      theme_minimal(base_size = 18) +
      theme(
        axis.text.x = element_text(size = 18, colour = "black",face = "plain", angle=30, hjust=1),
        axis.text.y = element_text(size = 18, colour = "black",face = "plain"),
        axis.title.y = element_text(size = 18, face = "plain", vjust = 0),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 1.2), 
        legend.position = "right",
        legend.text = element_text(size = 18, colour = "black"),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_text(size=16, face="plain"),
        strip.background = element_blank()) +
      scale_y_continuous(name="Summary score", limits=c(0, y_max), breaks=seq(0,y_max,by=y_break_step)) +
      labs(title="")
  
    ggsave(gg, filename=paste(out_dir, "/", project_name, "_sctype_tissue_guess.pdf", sep=""),
           device="pdf", dpi=300, useDingbats=FALSE)
    
    # if requested, make scType query all tissue/cell types by duplicating cell type database 
    # and adding cell type "All_tissues"
    if (sctype_tissue == "All") {
      db <-  openxlsx::read.xlsx(db_)
      row_count <- nrow(db)
      db_2 <- rbind(db, db)
      db_2$tissueType[(row_count+1):nrow(db_2)] <- "All"
      openxlsx::write.xlsx(db_2, "scType_db_full_modified.xlsx")
      db_ <- "scType_db_full_modified.xlsx"
    }
    
    # prepare gene sets
    gs_list <- gene_sets_prepare(db_, sctype_tissue)
    
    # get cell-type by cell matrix
    es.max <- sctype_score(scRNAseqData = seurat_obj[["RNA"]]@scale.data, scaled = TRUE, 
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
    
    # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
    # In case Seurat is used, it is either seurat_obj[["RNA"]]@scale.data (default), 
    # seurat_obj[["SCT"]]@scale.data, in case sctransform is used for normalization,
    # or seurat_obj[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
    
    # merge by cluster
    cL_results <- do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl) {
      es.max.cl <- sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$seurat_clusters==cl)), 10)
    }))
    sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
    
    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    
    # correct gamma delta T-cell name
    sctype_scores[sctype_scores$type == "γδ-T cells","type"] <- "Gamma delta T cells"
    
    # add cell type identities into metadata
    seurat_obj[[annotation_name]] <- ""
    for(j in unique(sctype_scores$cluster)){
      cl_type <- sctype_scores[sctype_scores$cluster==j,]
      seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters == j,annotation_name] = as.character(cl_type$type[1])
    }
    
    # create plot of cell type identity candidates per cluster
    # prepare edges
    cL_results <- cL_results[order(cL_results$cluster),]
    edges <- cL_results
    edges$type <- paste0(edges$type,"_",edges$cluster)
    edges$cluster <- paste0("cluster ", edges$cluster)
    edges <- edges[,c("cluster", "type")]
    colnames(edges) <- c("from", "to")
    rownames(edges) <- NULL
    
    # prepare nodes
    nodes_lvl1 <- sctype_scores[,c("cluster", "ncells")]
    nodes_lvl1$cluster <- paste0("cluster ", nodes_lvl1$cluster)
    nodes_lvl1$Colour <- "#f1f1ef"
    nodes_lvl1$ord <- 1
    nodes_lvl1$realname <- nodes_lvl1$cluster
    nodes_lvl1 <- as.data.frame(nodes_lvl1)
    nodes_lvl2 <- c()
    
    ccolss <- c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680",
                "#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c",
                "#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
    for (i in 1:length(unique(cL_results$cluster))){
      dt_tmp <- cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]
      nodes_lvl2 <- rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), 
                                                ncells = dt_tmp$scores, Colour = ccolss[i], 
                                                ord = 2, realname = dt_tmp$type))
    }
    nodes <- rbind(nodes_lvl1, nodes_lvl2)
    nodes$ncells[nodes$ncells<1] <- 1
    files_db <- openxlsx::read.xlsx(db_)[,c("cellName","shortName")]
    files_db <- unique(files_db)
    nodes <- merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", 
                   sort = F)
    nodes$shortName[is.na(nodes$shortName)] <- nodes$realname[is.na(nodes$shortName)]
    nodes <- nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
    
    mygraph <- graph_from_data_frame(edges, vertices=nodes)
    
    # Make the graph
    gg <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
      geom_node_circle(aes(filter=ord==1, fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + 
      geom_node_circle(aes(filter=ord==2, fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
      theme_void() + 
      geom_node_text(aes(filter=ord==2, label=shortName, colour="black", fill="black", 
                     repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ 
      geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(2), 
                      fill="white", parse = T), repel = !0, segment.linetype="dotted")
    
    ggsave(gg, filename=paste(out_dir, "/", project_name, "_sctype_identity_circles_", 
                              sctype_tissue, ".pdf", sep=""), 
           device="pdf", dpi=300, useDingbats=FALSE)
  }
  
  return(seurat_obj)
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
        if (is.character(seurat_obj[[ident]][,]) || is.factor(seurat_obj[[ident]][,])) {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident,
                        label = TRUE, pt.size = 0.4) + NoLegend()
          ggsave(gg, filename=paste(out_dir, "/", project_name, "_umap_", ident, ".", 
                                    file_type, sep=""), device=file_type, dpi=300)
        } else {
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
        if (is.character(seurat_obj[[ident]][,]) || is.factor(seurat_obj[[ident]][,]))  {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident,
                        label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=paste(out_dir, "/", project_name, "_umap_", ident, "_no_labels.", 
                                    file_type, sep=""), device=file_type, dpi=300, width=11)
        } else {
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





### Analysis

print(paste("Adding annotation to", in_file))

leiden_res <- as.numeric(leiden_res)


# create output directory if it does not exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}


# read in Seurat object
seurat_obj <- readRDS(in_file)


# read in cell type annotation
annotation_df <- read.csv(annotation, header=TRUE)


# add annotation to Seurat object
seurat_obj <- annotate_seurat(seurat_obj=seurat_obj, out_dir=out_dir, project_name=project_name, 
                              mode="custom", annotation_df=annotation_df)


# re-level cell type annotation
if ("order" %in% colnames(annotation_df)) {
  annotation_df <- annotation_df[order(annotation_df$order),]
  seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = unique(annotation_df$cell_type))
  Idents(seurat_obj) <- seurat_obj$cell_type
}


# plot UMAP with metadata overlays
idents_list = colnames(seurat_obj[[]])
idents_list = idents_list[!idents_list %in% c("barcode", "cell_id")]
seurat_obj <- umap_seurat(seurat_obj=seurat_obj, out_dir=out_dir, project_name=project_name, 
                          resolution=leiden_res, idents_list=idents_list, labels=FALSE, 
                          file_type="png")


# export metadata as csv
meta_data_export <- seurat_obj[[]]
write.csv(meta_data_export, paste(out_dir, "/", project_name, "_metadata.csv", sep=""))


# save final object
print(head(seurat_obj[[]]))

saveRDS(seurat_obj, paste(out_dir, "/", project_name, 
                          "_scRNAseq_no_doublets_annotated.rds", sep=""))






