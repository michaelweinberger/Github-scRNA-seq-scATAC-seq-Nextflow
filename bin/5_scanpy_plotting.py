#!/usr/bin/env python



# -*- coding: utf-8 -*-
"""
Plot cell type compositions across scanpy objects

Done using Python/3.10.7

"""



### Modules
import matplotlib
#from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import random
import seaborn as sns
import argparse

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.autosave = True	          # do not show plots, save them in figdir
sc.settings.set_figure_params(dpi_save=300, facecolor='white')

np.random.seed(0)



## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='Plot cell type compositions across scanpy objects')
parser.add_argument('-i','--input_file', help='Filepath to input Anndata object', required=True)
parser.add_argument('-a','--annotation', help='Filepath to .csv file containing "leiden" column with leiden cluster numbers and \
                    "cell_type" column with cell type annotation', required=True)
parser.add_argument('-o','--output_dir', help='Directory for output files', required=True)
parser.add_argument('-n','--output_name', help='Prefix for output file names', required=True)
parser.add_argument('-r','--leiden_res', help='Resolution of Leiden clustering', required=False)

args = parser.parse_args()
in_file = args.input_file
annotation_file = args.annotation
out_dir = args.output_dir
out_name = args.output_name
leiden_res = float(args.leiden_res)



### Functions

def obs_bar(adata, column_1, column_2, out_dir, out_name,
            x_axis_levels=None, fill_levels=None):
    
    # count cell type labels by leiden cluster
    plot_df = adata.obs.groupby([column_1, column_2]).agg({column_2: 'size'})
    plot_df = plot_df.rename(columns={column_2:'count'})
    plot_df = plot_df.reset_index()
    
    # generate dataframe with percentage counts for each entry in column_2 in separate column
    plot_df_1 = plot_df.pivot(index=column_1, columns=column_2, values='count')
    plot_df_1[pd.isna(plot_df_1)] = 0
    plot_df_1.to_csv(f"{out_dir}/{out_name}_counts_{column_1}_vs_{column_2}.csv")
    plot_df_1 = plot_df_1.div(plot_df_1.sum(axis=1), axis=0)
    
    # adjust x-axis plotting order
    if x_axis_levels != None:
        plot_df_1 = plot_df_1.reindex(labels=x_axis_levels, axis=0)
        
    # adjust stacked bar plotting order
    if fill_levels != None:
        plot_df_1 = plot_df_1[reversed(fill_levels)]
    
    # generate stacked bar plot
    values = plot_df_1.to_dict(orient='list')
    labels = plot_df_1.index
    bottom = np.zeros(len(labels))
    width = 0.8
    #fig_width = len(labels) * 0.7
    
    fig, ax = plt.subplots(figsize=(14,9), frameon=False)
    
    for key, value in values.items():
        p = ax.bar(labels, value, width, label=key, bottom=bottom)
        bottom += value
        #ax.bar_label(p, label_type='center')
        
    ax.set_axisbelow(True)
    ax.set_title("")
    ax.set_xlabel(column_1, fontsize=13)
    ax.set_ylabel("Relative count", fontsize=13)
    ax.tick_params(axis='x', which='major', labelsize=13, labelrotation=30)
    ax.tick_params(axis='y', which='major', labelsize=13)
    ax.tick_params(axis='both', which='minor', labelsize=13)
    for tick in ax.xaxis.get_majorticklabels():
        tick.set_horizontalalignment("right")
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title=column_2, 
              title_fontsize=12, fontsize=13, frameon=False)
    
    ax.figure.savefig(f"{out_dir}/{out_name}_{column_1}_vs_{column_2}_stacked_bar.pdf",
                      dpi=300)
    
    plt.close()
    return(plot_df_1)





### Analysis
    
sc.settings.plot_suffix = f"_scanpy_{out_name}"
sc.settings.figdir = out_dir


if not os.path.isfile(f"{out_dir}/{out_name}_scRNAseq_no_doublets_annotated.h5ad"):

    print(f"Adding annotation to {in_file}")
    
    # read in clustered Anndata object
    adata = sc.read_h5ad(in_file)

    # add cell type annotation
    annotation_df = pd.read_csv(annotation_file)
    annotation_df['cluster'] = annotation_df['cluster'].astype(str)

    tmp_df = adata.obs.merge(annotation_df, left_on='leiden', right_on='cluster', how='left', sort=False)
    tmp_df.index = tmp_df['barcode']
    #check_order = adata.obs.index.equals(tmp_df.index)
    #print(f".obs order after including metadata matches: {check_order}")
    adata.obs = tmp_df
    adata.obs.index.name = None
    
    # adjust plotting order
    if 'order' in annotation_df.columns:
        annotation_df = annotation_df.sort_values(by = ['order'])
        adata.obs['cell_type'] = pd.Categorical(values=adata.obs.cell_type, categories=annotation_df['cell_type'].unique(), 
                                            ordered=True)
    # drop unnamed columns from metadata
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.contains('^Unnamed')]

    # plot umap with metadata overlay
    for col in adata.obs.columns:
        if col != 'barcode':
            if col != 'cluster':
                sc.pl.umap(adata, color=[col], 
                           legend_loc='right margin', title='', frameon=False, ncols=1, 
                           save=f"_{col}.png")

    # save final dataset
    adata.obs.to_csv(f"{out_dir}/{out_name}_obs.csv")
    adata.write(f"{out_dir}/{out_name}_scRNAseq_no_doublets_annotated.h5ad", compression='gzip')
else:
    adata = sc.read_h5ad(f"{out_dir}/{out_name}_scRNAseq_no_doublets_annotated.h5ad")


# save object with entire count matrix for scVI integration
if not os.path.isfile(f"{out_dir}/{out_name}_scRNAseq_no_doublets_annotated_scVI.h5ad"):

    # read in dataset with count matrix containing all genes
    adata_1 = sc.read_h5ad(f"{out_dir}/{out_name}.h5ad")
    adata_1.obs = adata_1.obs[['barcode']]
    
    # add metadata of analysed dataset
    tmp_df = adata_1.obs.merge(adata.obs, on='barcode', how='left', sort=False)
    tmp_df.index = tmp_df['barcode']
    check_order = adata_1.obs.index.equals(tmp_df.index)
    print(f"{out_name} .obs order after including metadata matches: {check_order}")
    adata_1.obs = tmp_df
    adata_1.obs.index.name = None
    adata_1 = adata_1[adata_1.obs['cell_type'].notna()]
    
    # save final dataset
    adata_1.write(f"{out_dir}/{out_name}_scRNAseq_no_doublets_annotated_scVI.h5ad", compression='gzip')
    
    
    