#!/usr/bin/env python



# -*- coding: utf-8 -*-
"""
Preprocess and cluster scRNA-seq data with scanpy

"""


### Modules
#import matplotlib
#from matplotlib import rcParams
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import os
import random
import doubletdetection
import argparse

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.autosave = True	          # do not show plots, save them in figdir
sc.settings.set_figure_params(dpi=300, facecolor='white')

random.seed(0)



## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='Preprocess and cluster scRNA-seq data with scanpy')
parser.add_argument('-i','--input_dir', help='Directory containing cellranger output files \
                    (default: /outs/count/filtered_feature_bc_matrix)', required=True)
parser.add_argument('-m','--metadata', help='.tsv file containing "barcode" column with cell barcodes, \
                    "sample" column with sample names and additional metadata columns', required=True)
parser.add_argument('-o','--output_dir', help='Directory for output files', required=True)
parser.add_argument('-n','--output_name', help='Prefix for output file names', required=True)
parser.add_argument('-mig','--min_genes', help='Minimum number of genes detected for a cell to be kept in the dataset', required=True)
parser.add_argument('-mag','--max_genes', help='Maximum number of genes detected for a cell to be kept in the dataset', required=True)
parser.add_argument('-mam','--max_perc_mt', help='Maximum percentage of mitochondrial gene counts detected for a cell to be kept in the dataset', required=True)
parser.add_argument('-mic','--min_cells', help='Minimum number of cells in which a gene needs to be detected for it to be kept in the dataset', required=True)
parser.add_argument('-npc','--n_pcs', help='Number of principal components to use for neighbourhood graph', required=True)
parser.add_argument('-hv','--harmony_var', help='Name of metadata column to use for data integration', required=True)
parser.add_argument('-r','--leiden_res', help='Resolution of Leiden clustering', required=False)

args = parser.parse_args()
in_dir = args.input_dir
metadata = args.metadata
out_dir = args.output_dir
out_name = args.output_name
min_genes =  int(args.min_genes)
max_genes = int(args.max_genes)
max_perc_mt = int(args.max_perc_mt)
min_cells = int(args.min_cells)
n_pcs = int(args.n_pcs)
harmony_var = args.harmony_var
leiden_res = float(args.leiden_res)



### Functions

def scanpy_create_10X(cellranger_dir, out_dir, out_name, metadata=None, min_genes=200, min_cells=3):
    
    """
    Function to create an Anndata object from 10X cellranger output
    
    Input:
    cellranger_dir: Directory containing cellranger 'matrix.mtx.gz' output matrix,
                    usually .../outs/count/filtered_feature_bc_matrix/
    out_dir: Directory to write Anndata object and plots to
    out_name: Name of dataset, will be used in plot and output file names
    metadata: name of optional file with cell metadata, should be tab-separated and
              contain a column with cell barcodes named 'barcode'
    
    """
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    # set plot file suffix
    sc.settings.plot_suffix = f"_scanpy_{out_name}"
    sc.settings.figdir = out_dir
    
    # read in data
    # check how many 'matrix.mtx.gz' files there are in the directory
    matrix_file_list = []
    with os.scandir(cellranger_dir) as it:
        for entry in it:
            if entry.name.endswith('matrix.mtx.gz') and entry.is_file():
                matrix_file_list.append(entry.name)
    
    if len(matrix_file_list) == 0:
        print(f"Error: Did not find a 'matrix.mtx.gz' file in {cellranger_dir}")
    elif len(matrix_file_list) == 1:
        adata = sc.read_10x_mtx(
            cellranger_dir,
            var_names='gene_symbols',                
            cache=True)
    elif len(matrix_file_list) > 1:
        samples = []
        for entry in matrix_file_list:
            sample_name = entry.replace('_matrix.mtx.gz', '')
            s = sc.read(
                filename=f"{cellranger_dir}/{sample_name}_matrix.mtx.gz",
                cache=False
                ).T
            genes = pd.read_csv(f"{cellranger_dir}/{sample_name}_genes.tsv.gz", header=None, sep='\t')
            s.var_names = genes[0]
            s.var['gene_symbols'] = genes[1].values
            s.obs_names = pd.read_csv(f"{cellranger_dir}/{sample_name}_barcodes.tsv.gz", 
                                      header=None)[0]
            s.obs['sample_id'] = sample_name
            samples.append(s)
        
        adata = ad.concat(samples, join='outer', index_unique="_", merge='first')
        adata.var['gene_ids'] = adata.var_names
        adata.var_names = adata.var['gene_symbols']
        adata.var.index.name = None
    
    adata.obs['barcode'] = adata.obs.index
    adata.var_names_make_unique()    
    
    # merge in metadata   
    if metadata!=None:
        metadata_df = pd.read_table(metadata)
        tmp_df = adata.obs.merge(metadata_df, on='barcode', how='left', sort=False)
        tmp_df.index = tmp_df['barcode']
        #check_order = adata.obs.index.equals(tmp_df.index)
        #print(f".obs order after including metadata matches: {check_order}")
        adata.obs = tmp_df
    
    # make cell barcodes unique to dataset
    #adata.obs['barcode'] = adata.obs['barcode'] + "_" + out_name
    adata.obs.index = adata.obs['barcode']
    adata.obs.index.name = None
    
    # filter data
    sc.pl.highest_expr_genes(adata, n_top=20, save='.pdf')
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # save final scanpy output file
    adata.write(f"{out_dir}/{out_name}.h5ad", compression='gzip')
    
    return(adata)



def doublet_detection(adata, out_dir, out_name, var=None):
    
    """
    Function to analyse an Anndata object with Scanpy
    
    Input:
    adata: Anndata object
    out_dir: Directory to write Anndata object and plots to
    out_name: Name of dataset, will be used in plot and output file names
    var: Metadata variable according to which adata will be split for individual doublet detection,
        should be set if cellranger aggregation has been run on multiple samples
        
    """
    
    adata.obs['barcode'] = adata.obs.index
    adata.obs.index.name = None
    
    clf = doubletdetection.BoostClassifier(
        n_iters=10,
        clustering_algorithm="leiden",
        standard_scaling=True,
        pseudocount=0.1,
        n_jobs=-1,
        )
    
    if var==None:
        doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
        doublet_score = clf.doublet_score()
        adata.obs["doublet"] = doublets
        adata.obs["doublet_score"] = doublet_score
    else:
        barcode_list=[]
        doublets_list=[]
        doublet_score_list=[]
        
        for observ in pd.unique(adata.obs[var]):
            bdata = adata[adata.obs[var] == observ]
            doublets = clf.fit(bdata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
            doublet_score = clf.doublet_score()  
            
            barcode_list.append(bdata.obs['barcode'])
            doublets_list.append(doublets)
            doublet_score_list.append(doublet_score)
        
        barcode_list_flat = [item for sublist in barcode_list for item in sublist]
        doublets_list_flat = [item for sublist in doublets_list for item in sublist]
        doublet_score_list_flat = [item for sublist in doublet_score_list for item in sublist]

        doublet_df = pd.DataFrame(barcode_list_flat, columns=['barcode'])
        doublet_df['doublets'] = doublets_list_flat
        doublet_df['doublet_score'] = doublet_score_list_flat
        
        # exclude duplicated barcodes, keep the one with the higher 'doublets' entry
        doublet_df.sort_values(by=['doublets'], ascending=False, inplace=True)
        doublet_df.drop_duplicates(subset='barcode', keep='first', inplace=True)
        
        tmp_df = adata.obs.merge(doublet_df, on="barcode", how="left", suffixes=(None,"_1"),
                                 sort=False)
        tmp_df.index = tmp_df['barcode']
        tmp_df.to_csv(f"{out_dir}/{out_name}_doublets.csv")
        
        #check_order = adata.obs.index.equals(tmp_df.index)
        #print(f".obs order after excluding doublets matches: {check_order}")
        
        adata.obs = tmp_df
        adata.obs.index.name = None
        
    doubletdetection.plot.convergence(clf, save=f"{out_dir}/{out_name}_convergence_test.pdf", show=True, 
                                      p_thresh=1e-16, voter_thresh=0.5)
    
    # save final scanpy output file
    adata.write(f"{out_dir}/{out_name}_doublets_detected.h5ad", compression='gzip')
    
    return(adata)



def scanpy_analysis(adata, out_dir, out_name, leiden_res=0.4, min_genes=200, min_cells=3, 
                    n_genes_cutoff=2500, pct_mt_cutoff=5,
                    neighbors=10, pcs=30, 
                    batch_var=None, harmony_var=None):
    
    """
    Function to analyse an Anndata object with Scanpy
    
    Input:
    adata: Anndata object
    out_dir: Directory to write Anndata object and plots to
    out_name: Name of dataset, will be used in plot and output file names
    batch_var: Name of obs column to use as a batch key for highly variable gene detection
    leiden_res: Resolution to use for leiden clustering
    harmony_var: Name of obs column to use as a key for harmony integration
    
    """
    
    random.seed(0)
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    # set plot file suffix
    sc.settings.plot_suffix = f"_scanpy_{out_name}"
    sc.settings.figdir = out_dir
    
    # restore raw counts
    if not adata.raw is None:
        adata = adata.raw.to_adata()
    
    # filter data
    sc.pl.highest_expr_genes(adata, n_top=20, save='.png')
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # quality control
    # analyse mitochondrial gene content
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    sc.pl.violin(
        adata, 
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        jitter=0.4, 
        multi_panel=True, 
        save='_qc_violin.png'
        )
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_qc_scatter_1.png')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_qc_scatter_2.png')
    
    adata = adata[adata.obs.n_genes_by_counts < n_genes_cutoff, :]
    adata = adata[adata.obs.pct_counts_mt < pct_mt_cutoff, :]
    
    # normalise data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # identify highly variable genes
    sc.pp.highly_variable_genes(
        adata, 
        min_mean=0.0125, 
        max_mean=3, 
        min_disp=0.5, 
        batch_key=batch_var
        )
    sc.pl.highly_variable_genes(adata, save='.png')
    
    adata.raw = adata
    
    # subset dataset
    adata = adata[:, adata.var.highly_variable]
    
    # eliminate unwanted variability
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    
    # scale gene expression to unit variance. Clip values exceeding standard deviation 10.
    sc.pp.scale(adata, max_value=10)
    
    # compute PCs
    sc.tl.pca(
        adata, 
        svd_solver='arpack', 
        random_state=0
        )
    sc.pl.pca_variance_ratio(
        adata, 
        log=True, 
        save='.pdf'
        )
    
    # optionally integrate data batches with harmonypy
    if harmony_var==None:
        # compute neighbourhood graph
        sc.pp.neighbors(
            adata, 
            n_neighbors=neighbors, 
            n_pcs=pcs,
            random_state=0
            )
        
    else:
        sce.pp.harmony_integrate(adata, harmony_var)
        
        # compute neighbourhood graph
        sc.pp.neighbors(
            adata, 
            n_neighbors=neighbors, 
            n_pcs=pcs, 
            use_rep='X_pca_harmony',
            random_state=0
            )
    
    # determine clusters
    sc.tl.leiden(
        adata, 
        resolution=leiden_res,
        random_state=0
        )
    
    # compute embeddings
        
    sc.tl.paga(adata)
    sc.pl.paga(
        adata, 
        plot=True, 
        edge_width_scale=2, 
        node_size_scale=2, 
        save='.pdf'
        )
    
    sc.tl.umap(
        adata, 
        init_pos='paga', 
        random_state=0
        )
    sc.pl.umap(
        adata, 
        color=['leiden'], 
        legend_loc='right margin', 
        title='', 
        frameon=False, 
        ncols=1, 
        save=f"_leiden_{leiden_res}.png"
        )
    
    # identify marker genes
    sc.tl.rank_genes_groups(
        adata, 
        'leiden', 
        method='wilcoxon'
        )
    sc.pl.rank_genes_groups(
        adata, 
        n_genes=25, 
        sharey=False, 
        save='.pdf'
        )
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    markers = pd.DataFrame({group + '_' + key[:1]: result[key][group] 
                            for group in groups for key in ['names', 'pvals']})
    markers.to_csv(f"{out_dir}/{out_name}_markers.csv")
    
    # save final scanpy output file
    adata.write(f"{out_dir}/{out_name}_analysed.h5ad", compression='gzip')
    
    #print(adata)
    return(adata)





### Analysis
    
# create Anndata object
if not os.path.isfile(f"{out_dir}/{out_name}.h5ad"):
    adata = scanpy_create_10X(cellranger_dir=in_dir, 
                              out_dir=out_dir, 
                              out_name=out_name,
                              metadata=metadata,
                              min_genes=200, min_cells=3)
else:
    adata = sc.read_h5ad(f"{out_dir}/{out_name}.h5ad")



# detect doublets
if not os.path.isfile(f"{out_dir}/{out_name}_doublets_detected.h5ad"):
    adata = doublet_detection(adata=adata, 
                              out_dir=out_dir, 
                              out_name=out_name, 
                              var="sample_id")
else:
    adata = sc.read_h5ad(f"{out_dir}/{out_name}_doublets_detected.h5ad")



# scanpy analysis
if not os.path.isfile(f"{out_dir}/{out_name}_scRNAseq_analysed_no_doublets.h5ad"):
    adata = scanpy_analysis(adata=adata, 
                            out_dir=out_dir, 
                            out_name=out_name,
                            batch_var="sample_id",
                            harmony_var=harmony_var,
                            min_genes=min_genes, min_cells=min_cells, 
                            n_genes_cutoff=max_genes, pct_mt_cutoff=max_perc_mt,
                            neighbors=10, pcs=n_pcs, leiden_res=leiden_res)
    
    # plot doublets
    sc.pl.umap(adata, color=['doublets', 'doublet_score'], 
               legend_loc='right margin', title='Doublets', frameon=False, ncols=1, 
               save='_doublets.png')
    
    # remove doublets
    adata = adata[adata.obs['doublets'] == 0]

    # plot
    sc.pl.umap(adata, color=['sample_id'], 
               legend_loc='right margin', title='', frameon=False, ncols=1, 
               save='_sample_id.png')
    
    sc.pl.umap(adata, color=['leiden'], 
               legend_loc='right margin', title='', frameon=False, ncols=1, 
               save=f'_leiden_{leiden_res}.png')

    # save final scanpy output file
    adata.obs['dataset'] = out_name
    adata.obs['sample_id'] = adata.obs['sample_id'].cat.rename_categories(lambda x: x + "_" + out_name)

    print(adata)
    adata.obs.to_csv(f"{out_dir}/{out_name}_obs.csv")
    adata.write(f"{out_dir}/{out_name}_scRNAseq_analysed_no_doublets.h5ad", compression='gzip')

