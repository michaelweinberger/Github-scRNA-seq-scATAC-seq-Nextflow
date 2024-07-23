#!/usr/bin/env python



# -*- coding: utf-8 -*-
"""
Perform mRNA velocity analysis using scVelo

"""


import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import os
import argparse



## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='Perform mRNA velocity analysis using scVelo')
parser.add_argument('-i','--input_file', help='Path to input Seurat (suffix ".rds") or Anndata (suffix ".h5ad") object', required=True)
parser.add_argument('-l','--loom_file', help='Path to loom file containing spliced/unspliced counts', required=True)
parser.add_argument('-o','--output_dir', help='Path to directory for saving output files', required=True)
parser.add_argument('-n','--project_name', help='Name used as prefix for plot files, etc', required=True)

args = parser.parse_args()
in_file = args.input_file
loom_file = args.loom_file
out_dir = args.output_dir
project_name = args.project_name


## global settings
scv.settings.verbosity = 3
scv.settings.set_figure_params('scanpy', facecolor='white', dpi=300, frameon=False)
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.autosave = True	          # do not show plots, save them in figdir
sc.settings.set_figure_params(dpi_save=300, facecolor='white')


### Functions

def scvelo_analysis(adata, ldata, out_dir, project_name, color_palette=None):
    
    """
    Function to perform scRNA-seq mRNA velocity analysis with scVelo
    
    Input:
    adata: scRNA-seq AnnData object, count data may be (but do not need to be) normalised and log-transformed,
           adata.obsm['X_umap'] and adata.uns['umap'] fields should be present,
           adata.obs['cell_type'] with cell type labels should be present
    ldata: AnnData object containing matrices of spliced and unspliced read counts
    out_dir: Directory to write output to
    project_name: Name of dataset, will be used in plot and output file names
    color_palette: List of colors for cell types in plots
    
    """
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    # set plot file suffix
    scv.settings.figdir = out_dir
    sc.settings.plot_suffix = f"_scanpy_{project_name}"
    sc.settings.figdir = out_dir
    
    # set obs index
    #adata.obs.index = adata.obs['barcode']
          
    # merge matrices into the original adata object
    adata1 = scv.utils.merge(adata, ldata)
    
    # plot umap to check dataset
    sc.pl.umap(adata1, color=['cell_type'], frameon=False, title=project_name, 
               palette=color_palette, save='_celltypes.pdf')
    
    # show proportion of spliced versus unspliced reads
    scv.pl.proportions(adata1, groupby='cell_type', dpi=300,
                       save=project_name + '_spliced_unspliced_proportions.pdf')
    
    # pre-process
    print('Filtering and normalizing')
    """
    From the scvelo docs: 
        "Filtering and normalization is applied in the same vein to 
        spliced/unspliced counts and X. Logarithmizing is only applied to X. 
        If X is already preprocessed from former analysis, it will not be touched."
    """
    scv.pp.filter_and_normalize(adata1, n_top_genes=2000)
    print('Computing moments')
    scv.pp.moments(adata1)
    
    # compute velocity
    print('Recovering dynamics')
    scv.tl.recover_dynamics(adata1)
    print('Computing velocities')
    scv.tl.velocity(adata1, mode='dynamical')
    print('Generating velocity graph')
    scv.tl.velocity_graph(adata1)
    
    # calculate latent time (only in dynamical mode)
    print('Calculating latent time')
    scv.tl.latent_time(adata1)
    scv.pl.scatter(adata1, color='latent_time', color_map='gnuplot', size=100, dpi=300,
                   title=project_name + '_latent_time', save=project_name + '_latenttime_umap.png')
        
    top_genes = adata1.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata1, var_names=top_genes, sortby='latent_time', col_color='cell_type', 
                   n_convolve=100, save=project_name + '_latenttime_genes_heat.png')
    
    # calculate differential kinetics
    adata1.write(f"{out_dir}/{project_name}_seurat_anndata_scvelo_tmp.h5ad", compression='gzip')
    
    print('Computing differential kinetics')
    scv.tl.differential_kinetic_test(adata1, groupby='cell_type')
    print(adata1.var.head())
    scv.pl.scatter(adata1, basis=top_genes[:15], ncols=5, frameon=True, add_outline='fit_diff_kinetics',  
                   title=project_name + '_latenttime_genes_scatter',
                   save=project_name + '_latenttime_genes_scatter.png')
        
    # include differential kinetic analysis results in velocity plots
    print('Re-computing velocities')
    scv.tl.velocity(adata1, mode='dynamical', diff_kinetics=True)
    print('Re-generating velocity graph')
    scv.tl.velocity_graph(adata1)
        
    # compute PAGA
    adata1.uns['neighbors']['distances'] = adata1.obsp['distances']
    adata1.uns['neighbors']['connectivities'] = adata1.obsp['connectivities']
    print('Generating PAGA')
    scv.tl.paga(adata1, groups='cell_type')
    
    # plot graphs + embeddings
    scv.pl.velocity_graph(adata1, threshold=.2, color=['cell_type'], 
                          title=project_name + '_velocity_graph',
                          save=project_name + '_velocity_graph.png')
    
    scv.pl.paga(adata1, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.2, dpi=300, ncols=1,
                fontsize=5, arrowsize=10, title=project_name + '_paga', 
                save=project_name + '_paga.pdf')
    
    scv.pl.velocity_embedding(adata1, basis='umap', arrow_length=3, arrow_size=2, dpi=300,
                              title=project_name, palette=color_palette,
                              frameon=False, save=project_name +'_embedding.png')
    
    scv.pl.velocity_embedding_grid(adata1, basis='umap', color='cell_type', dpi=300,
                                   save=project_name + '_embedding_grid.png', 
                                   legend_loc='right margin', palette=color_palette,
                                   title=project_name + '_grid', scale=0.25, arrow_length=1)
    
    scv.pl.velocity_embedding_stream(adata1, basis='umap', color='cell_type', dpi=300,
                                     save=project_name +'_embedding_stream.png',
                                     legend_loc='right margin', size=50, alpha=0.5,
                                     title=project_name + '_stream', palette=color_palette)
    
    # identify genes driving cell state transitions
    scv.tl.rank_dynamical_genes(adata1, groupby='cell_type')
    scvelo_genes = scv.get_df(adata1, 'rank_dynamical_genes/names')
    scvelo_genes.to_csv(project_name + '_dynamical_genes.csv', sep="\t", index=False)
    
    # calculate pseudotime
    print('Calculating pseudotime')
    scv.tl.velocity_pseudotime(adata1)
    scv.pl.scatter(adata1, color='velocity_pseudotime', cmap='gnuplot', dpi=300, 
                   title=project_name + '_pseudotime',
                   save=project_name + '_pseudotime_umap.png')
    
    # plot velocity confidence
    scv.tl.velocity_confidence(adata1)
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata1, c=keys, cmap='coolwarm', perc=[5, 95], dpi=300,
                   save=project_name + '_confid_umap.png')
    
    # save final object
    print('Saving final object')
    adata1.write(f"{out_dir}/{project_name}_scvelo.h5ad", compression='gzip')
    
    return(adata1)



### Analysis

# load loom files for spliced/unspliced matrices for each sample:
ldata = scv.read(loom_file, cache=False)
ldata.obs.index = ldata.obs['obs_names']
ldata.obs.index.name = None

# run scvelo
adata = sc.read_h5ad(in_file)
adata1 = scvelo_analysis(adata=adata, ldata=ldata, out_dir=out_dir, project_name=project_name, color_palette=None)







