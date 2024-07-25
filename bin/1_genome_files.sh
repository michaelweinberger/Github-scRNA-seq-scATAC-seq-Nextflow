#!/bin/bash



#############################################################
####                                                    #####
####              User defined variables                #####
####                                                    #####
#############################################################



############################################################



set -e
set -u



### declare input arguments
while getopts ":s:g:u:l:e:o:" flag
do
	case $flag in
        	s) 
			species="$OPTARG"
			echo "Option -s with argument: $species"
			;;
        	g) 
			genome="$OPTARG"
			echo "Option -g with argument: $genome"
			;;
        	u) 
			genome_ucsc="$OPTARG"
			echo "Option -u with argument: $genome_ucsc"
			;;
        	l) 
			species_latin="$OPTARG"
			echo "Option -l with argument: $species_latin"
			;;
        	e) 
			ensembl_version="$OPTARG"
			echo "Option -e with argument: $ensembl_version"
			;;
        	o) 
			outdir="$OPTARG"
			echo "Option -o with argument: $outdir"
			;;
	esac
done



# Capitalise first letter
species_latin_2=${species_latin^}



### download or generate cellranger reference data
if [ $species = "human" ] ; then
	cellranger_ref="refdata-gex-GRCh38-2020-A.tar.gz"
elif [ $species = "mouse" ] ; then
	cellranger_ref="refdata-gex-mm10-2020-A.tar.gz"
else
	cellranger_ref="refdata-cellranger-${genome}"
fi



if [ $species = "human" ] || [ $species = "mouse" ] ; then
	wget https://cf.10xgenomics.com/supp/cell-exp/${cellranger_ref}
	tar -zxvf ${cellranger_ref}
	rm ${cellranger_ref}
else
	rsync -avzP rsync://ftp.ebi.ac.uk/ensemblorg/pub/release-${ensembl_version}/fasta/${species_latin}/dna/${species_latin_2}.${genome}.dna.primary_assembly.fa.gz .
	mv ./${species_latin_2}.${genome}.dna.primary_assembly.fa.gz ./${genome}.fa.gz
	gunzip ./${genome}.fa.gz

	rsync -avzP rsync://ftp.ebi.ac.uk/ensemblorg/pub/release-${ensembl_version}/gtf/${species_latin}/${species_latin_2}.${genome}.${ensembl_version}.gtf.gz .
	mv ./${species_latin_2}.${genome}.${ensembl_version}.gtf.gz ./${genome}.${ensembl_version}.gtf.gz
	gunzip ./${genome}.${ensembl_version}.gtf.gz
fi




### download file containing positions of repetitive elements in genome
#wget -L http://hgdownload.soe.ucsc.edu/goldenPath/${genome_ucsc}/database/rmsk.txt.gz .
#gunzip ./rmsk.txt.gz
#mv ./rmsk.txt ./${genome}_rmsk.txt



echo All done!
