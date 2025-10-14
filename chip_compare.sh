#!/usr/bin/env bash

set -u
set -e

# @description Script to compare two ChIP-seq experiments based on peak counts
# @description Requires a count table (featureCounts output) and a sample table
# @description Generates an HTML report with QC and differential analysis
# @params s sampleTable = samples_experiment.txt # Tab-delimited file with columns: sample_id, sample_name, count_column, condition
# @params c conditions_colname = condition # Column name in sampleTable that defines conditions
# @params l control_condition = exp1 # Name of the control condition in conditions_colname
# @params e experiment_condition = exp2 # Name of the experimental condition in conditions_colname
# @params t countTable = featureCount_counts.txt # Path to the count table (featureCounts output)
# @params a countTableAnnot = featureCount_counts.annot.txt # Path to the annotated count table
# @params o odir = . # Output directory for the report
# @params p downsample_pairs = 1000 # Number of pairs for downsampling analysis
# @params m downsample_ma = 0 # Number of reads for MA downsampling
# @params i session = http://localhost:60151/goto? # RStudio session URL
# @params r rawcount_quantile_cutoff = 0.02 # Quantile cutoff for raw count filtering (0 to disable)
# @params h help = 0 # Show this help message

function print_usage() {
	echo ""
	echo " Usage: chip_compare.sh options"
	echo ""
	echo " Script to compare two ChIP-seq experiments based on peak counts"
	echo " Requires a count table (featureCounts output) and a sample table"
	echo " Generates an HTML report with QC and differential analysis"
	echo ""
	echo "-s (sampleTable): Tab-delimited file with columns: sample_id, sample_name, count_column, condition [samples_experiment.txt]"
	echo "-c (conditions_colname): Column name in sampleTable that defines conditions [condition]"
	echo "-l (control_condition): Name of the control condition in conditions_colname [exp1]"
	echo "-e (experiment_condition): Name of the experimental condition in conditions_colname [exp2]"
	echo "-t (countTable): Path to the count table (featureCounts output) [featureCount_counts.txt]"
	echo "-a (countTableAnnot): Path to the annotated count table [featureCount_counts.annot.txt]"
	echo "-o (odir): Output directory for the report [.]"
	echo "-p (downsample_pairs): Number of pairs for downsampling analysis [1000]"
	echo "-m (downsample_ma): Number of reads for MA downsampling [0]"
	echo "-i (session): RStudio session URL [http://localhost:60151/goto?]"
	echo "-r (rawcount_quantile_cutoff): Quantile cutoff for raw count filtering (0 to disable) [0.02]"
	echo "-h Show this help message"

	exit 1
}


######
# MAIN
######

if [[ $# == 0 ]]; then
    print_usage
    exit 1        # Exit and explain usage, if no argument(s) given.
fi

while getopts "s:c:l:e:t:a:o:p:m:i:r:h" option; do
	case $option in
		s ) sampleTable=$OPTARG;;
		c ) conditions_colname=$OPTARG;;
		l ) control_condition=$OPTARG;;
		e ) experiment_condition=$OPTARG;;
		t ) countTable=$OPTARG;;
		a ) countTableAnnot=$OPTARG;;
		o ) odir=$OPTARG;;
		p ) downsample_pairs=$OPTARG;;
		m ) downsample_ma=$OPTARG;;
		i ) session=$OPTARG;;
		r ) rawcount_quantile_cutoff=$OPTARG;;
		h ) print_usage;;
		* ) echo -e "\nERROR: Unimplemented option chosen!\n"
			print_usage;;
	esac
done

# set defaults
downsample_ma=${downsample_ma:=0}
downsample_pairs=${downsample_pairs:=1000}
rawcount_quantile_cutoff=${rawcount_quantile_cutoff:=0.02}
odir=${odir:="."}
session=${session:="http://localhost:60151/goto?"}

# print set parameters
>&2 echo "-s (sampleTable)=$sampleTable"
>&2 echo "-c (conditions_colname)=$conditions_colname"
>&2 echo "-l (control_condition)=$control_condition"
>&2 echo "-e (experiment_condition)=$experiment_condition"
>&2 echo "-t (countTable)=$countTable"
>&2 echo "-a (countTableAnnot)=$countTableAnnot"
>&2 echo "-o (odir)=$odir"
>&2 echo "-p (downsample_pairs)=$downsample_pairs"
>&2 echo "-m (downsample_ma)=$downsample_ma"
>&2 echo "-i (session)=$session"
>&2 echo "-r (rawcount_quantile_cutoff)=$rawcount_quantile_cutoff"

echo "loading modules ..."
module load build-env/f2022
module load r-bundle-bioconductor/3.18-foss-2023a-r-4.3.2
module load pandoc/2.18
echo "loading modules done."

echo "running Rmarkdown ..."
Rscript -e "rmarkdown::render(

	input='/groups/gaidt/bioinf/software/scripts/deseq_pw/chip_compare.Rmd',
	output_file='chip_compare.html',

	params=list(

		# sampleTable:
		# sample_id	sample_name	count_column condition
		# condition: labels of two conditions to be comapred
		sampleTable='$sampleTable',

		conditions_colname='$conditions_colname',
		control_condition='$control_condition',
		experiment_condition='$experiment_condition',

		countTable='$countTable',
		countTableAnnot='$countTableAnnot',

		odir='$odir',
		downsample_pairs=$downsample_pairs,
		downsample_ma=$downsample_ma,
		rawcount_quantile_cutoff=$rawcount_quantile_cutoff,
		session='$session'
	))"

