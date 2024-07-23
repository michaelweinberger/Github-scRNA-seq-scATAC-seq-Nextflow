#!/usr/bin/env python



# -*- coding: utf-8 -*-
"""
Merge single sample loom files containing spliced and unspliced read counts

"""


import scvelo as scv
import pandas as pd
import argparse


## unpack arguments imported from bash parent script

# Define a custom argument type for a list of strings
def list_of_strings(arg):
    return arg.split(' ')

parser = argparse.ArgumentParser(description='Merge single sample loom files containing spliced and unspliced read counts')
parser.add_argument('-l','--loom_files', nargs='+', type=list_of_strings, help='List of filepaths of loom files containing spliced and unspliced scRNA-seq read counts', required=True)
parser.add_argument('-m','--metadata', help='Path to "cellranger_aggr_cell_metadata.tsv"', required=True)
parser.add_argument('-o','--out_dir', help='Path to output directory', required=True)

args = parser.parse_args()
loom_files = [item for sublist in args.loom_files for item in sublist]
metadata = args.metadata
out_dir = args.out_dir



### Functions

def merge_loom(loom_files, barcode_info_path, out_dir):
    
    """
    Function to merge single sample loom files containing spliced and unspliced read counts
    
    Input:
    loom files: List of filepaths of loom files to merge
    meta_dir: Path to directory containing "cellranger_aggr_cell_metadata.tsv" file containing cell metadata,
              should contain columns named "sample_id" with sample identifiers and "barcode" with cell barcodes
    out_dir: Path to directory to write output to
    
    """
    
    # read in cell metadata
    barcode_info = pd.read_table(barcode_info_path, header=0)
    #print(barcode_info.head())
       
    count = 1
    for loom in loom_files:
        
        ldata = scv.read(loom, cache=True)
        
        # rename cell barcodes
        sample_name = loom.split('/')[-1].split('.')[0]
        print('Processing: ' + sample_name)
        sample_barcode = barcode_info.loc[barcode_info['sample_id']==sample_name,'barcode'].iloc[0]
        sample_number = sample_barcode.split('-')[-1]
        print('Sample number: ' + sample_number)
        
        barcodes = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
        barcodes = [bc[0:len(bc)-1] + '-' + sample_number for bc in barcodes]
        ldata.obs.index = barcodes
        
        # make variable names unique
        ldata.var_names_make_unique()
        
        if count == 1:
            ldata_merge = ldata
        else:
            ldata_merge = ldata_merge.concatenate([ldata])
        
        count += 1
    
    # reformat cell barcodes in final file
    index_part_1 = [bc.split('-')[0] for bc in ldata_merge.obs.index.tolist()]
    index_part_2 = [bc.split('-')[1] for bc in ldata_merge.obs.index.tolist()]
    index_new = [i + '-' + j for i, j in zip(index_part_1, index_part_2)]
    ldata_merge.obs.index = index_new
    print(ldata_merge.obs)
    
    ldata_merge.write_loom(out_dir + '/merged.loom')
    
    return(ldata_merge)
    



### Analysis

loom_merged = merge_loom(loom_files=loom_files, barcode_info_path=metadata, out_dir=out_dir)



