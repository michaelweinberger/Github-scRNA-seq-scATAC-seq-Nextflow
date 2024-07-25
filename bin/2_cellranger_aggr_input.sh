#!/bin/bash



#########################################################
### This script:
# generates an input .csv file for cellranger aggr
#########################################################


#####################################################################################################################################################
################################################################ USER-DEFINED VARIABLES #############################################################


####################################################################################################################################################



set -e
set -u



### declare input arguments
while getopts ":s:i:" flag
do
	case $flag in
        	s) 
			sample_sheet="$OPTARG"
			echo "Option -s with argument: $sample_sheet"
			;;
        	i) 
			cellranger_count_info="$OPTARG"
			echo "Option -i with argument: $cellranger_count_info"
			;;
	esac
done



# create header for cellranger aggregator input file
sed 1q $sample_sheet | awk -F\\t -vOFS=, '{$1=$2=""; print $0}' | sed 's/^,//' > cellranger_aggr_input.csv
sed -i "s/^/sample_id,molecule_h5/" cellranger_aggr_input.csv

# loop over sample IDs & collect molecule info h5 paths + metadata
for j in $(sed '1d' $sample_sheet | awk -F\\t '{print $1}')
do
    echo $j

    # adjust sample id
    sample=$(grep $j $sample_sheet | awk -F\\t '{print $1}' | tr , _)

    # extract additional metadata columns and paste together with comma as separator
    metadata=$(grep $j $sample_sheet | awk -F\\t -vOFS=, '{$1=$2=""; print $0}' | sed 's/^,//')

    # create molecule_info .h5 file path
    cellranger_count_outdir=$(grep $j $cellranger_count_info | awk -F\\t '{print $2}')
    h5_path="${cellranger_count_outdir}/outs/molecule_info.h5"

    # append everything to .csv file
    echo "${sample},${h5_path}${metadata}" >> cellranger_aggr_input.csv
done



echo All done!