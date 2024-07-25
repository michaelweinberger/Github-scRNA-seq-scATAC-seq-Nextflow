#!/usr/bin/env python



# -*- coding: utf-8 -*-
"""
Construct Anndata object from count matrix, gene and cell metadata and PCA results

"""


import scvelo as scv
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import pandas as pd
import os
import argparse



## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='Construct Anndata object from count matrix, gene and cell metadata and PCA results')
parser.add_argument('-mtx','--count_matrix', help='Path to Seurat object expression count matrix .mtx file', required=True)
parser.add_argument('-md','--meta_data', help='Path toSeurat object metadata .csv file', required=True)
parser.add_argument('-g','--genes', help='Path to Seurat object gene name .csv file', required=True)
parser.add_argument('-p','--pca', help='Path to Seurat object PCA embedding .csv file', required=True)
parser.add_argument('-o','--output_dir', help='Path to directory for saving output files', required=True)
parser.add_argument('-n','--project_name', help='Name used as prefix for plot files, etc', required=True)

args = parser.parse_args()
counts = args.count_matrix
metadata = args.meta_data
gene_names= args.genes
pca_embedding = args.pca
out_dir = args.output_dir
project_name = args.project_name


## global settings
scv.settings.verbosity = 3
scv.settings.set_figure_params('scanpy', facecolor='white', dpi=300, frameon=False)
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.autosave = True	          # do not show plots, save them in figdir
sc.settings.set_figure_params(dpi_save=300, facecolor='white')


### Functions

def anndata_from_seurat(count_matrix, metadata, genes, pca, out_dir, project_name):
    
    """
    Function to create AnnData object from Seurat-exported data
    
    Input:
    count_matrix: Path to .mtx matrix file containing scRNA-seq read counts,
                  may be normalised and log-transformed
    metadata: Path to .csv file containing cell metadata, 
              should contain a column named "barcode" with cell barcodes
    genes: Path to .csv file containing column of gene names (one name per row, no column name)
    pca: Path to .csv file containing cell x PC PCA values, PC names as column names
    out_dir: Directory to write output to
    project_name: Name of dataset, will be used in plot and output file names
    
    """
                    
    # load sparse matrix:
    counts = io.mmread(count_matrix)
    
    # create anndata object
    counts = np.transpose(counts).tocsr()
    adata = anndata.AnnData(X=counts)
    
    # load cell metadata:
    cell_meta = pd.read_csv(metadata)
    
    # load gene names:
    with open(genes, 'r') as f:
        gene_names = f.read().splitlines()
    
    # set anndata observations and index obs by barcodes, var by gene names
    adata.obs = cell_meta
    adata.obs.index = adata.obs['barcode']
    adata.var.index = gene_names
    
    # load dimensional reduction:
    pca = pd.read_csv(pca)
    pca.index = adata.obs.index
    
    # set pca and umap
    adata.obsm['X_pca'] = pca.to_numpy()
    adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
    
    # adjust plotting order
    if 'order' in adata.obs.columns:
        tmp = adata.obs.sort_values(by = ['order'])
        adata.obs['cell_type'] = pd.Categorical(
            values=adata.obs.cell_type, 
            categories=tmp['cell_type'].unique(), 
            ordered=True
            )
        
    # plot a UMAP colored by cell type to test:
    #sc.pl.umap(adata, color=['cell_type'], frameon=False, save=project_name + '_celltypes.pdf',
    #           title=project_name)
    
    # save dataset as anndata format
    adata.write(f"{out_dir}/{project_name}_from_matrix_metadata_pca.h5ad")
    
    return(adata)



### Analysis

print('Constructing Scanpy object')

adata = anndata_from_seurat(count_matrix=counts, metadata=metadata, genes=gene_names, 
                            pca=pca_embedding, out_dir=out_dir, project_name=project_name)





