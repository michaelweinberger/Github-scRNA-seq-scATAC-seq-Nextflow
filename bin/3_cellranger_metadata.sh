#!/bin/bash



#########################################################
### This script:
# creates a gene expression count table using Cellranger
# input: 10x scRNAseq fastq files
#########################################################


#####################################################################################################################################################
################################################################ USER-DEFINED VARIABLES #############################################################


####################################################################################################################################################



set -e
set -u



### declare input arguments
while getopts ":i:" flag
do
	case $flag in
        	i) 
			cellranger_aggr_input_csv="$OPTARG"
			echo "Option -i with argument: $cellranger_aggr_input_csv"
			;;
	esac
done



# create header for cell metadata file
sed 1q $cellranger_aggr_input_csv | awk -F, '{for(i=3;i<=NF;i++) printf $i"\t"; print ""}' > cellranger_aggr_cell_metadata.tsv
sed -i "s/^/barcode\tsample_id\t/" cellranger_aggr_cell_metadata.tsv

# populate metadata file
sample_count=1

for j in $(sed '1d' $cellranger_aggr_input_csv | awk -F, '{print $1}')
do
	# extract path to cellranger count output directory
        h5_path=$(grep $j $cellranger_aggr_input_csv | awk -F, '{print $2}')
        pattern="/outs/molecule_info.h5"
        cellranger_count_outdir=${h5_path/${pattern}/}

        # get cell barcodes from cellranger count output
	gunzip ${cellranger_count_outdir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
	cat ${cellranger_count_outdir}/outs/filtered_feature_bc_matrix/barcodes.tsv | awk -F\\t '{print $(NF-0)}' > cellranger_aggr_cell_metadata_${j}.tsv

	# rename barcodes by replacing last letter "1" with the sample count (similar to what cellranger aggr does)
	sed "s/\(.\)$/${sample_count}/" cellranger_aggr_cell_metadata_${j}.tsv > cellranger_aggr_cell_metadata_${j}_2.tsv

	# add sample ID
	sample_id=$(grep $j $cellranger_aggr_input_csv | awk -F, '{print $1}')
	sed -i "s/$/\t${sample_id}/" cellranger_aggr_cell_metadata_${j}_2.tsv

	# add sample metadata
	# get metadata from aggregator input sheet, tab-delimited, trailing tab deleted
	metadata=$(grep $j $cellranger_aggr_input_csv | awk -F, '{for(i=3;i<=NF;i++) printf $i"\t"; print ""}' | sed 's/.$//')
	sed -i "s/$/\t${metadata}/" cellranger_aggr_cell_metadata_${j}_2.tsv

	# append everything to combined output file
	cat cellranger_aggr_cell_metadata_${j}_2.tsv >> cellranger_aggr_cell_metadata.tsv

	gzip ${cellranger_count_outdir}/outs/filtered_feature_bc_matrix/barcodes.tsv

	sample_count=`expr ${sample_count} + 1`
done



echo All done!