#!/bin/bash
set -o nounset  
set -o errexit  
set -o errtrace 

#
# example script on how to prepare the data
#
# please adjust: 
# tmp_dir = working directory, somewhere on scratch-cbe
# exp1_path = source path of experiment 1
# exp2_path = source path of experiment 2
# names and pathes to peak bedfiles 
# names and pathes to bam files
# names and pathes to bigwig files

#
# (1) create results direcoty on scratch
#
tmp_dir="/scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b"
mkdir -p $tmp_dir/{data/bed,data/bam,data/bigwig,results}

#
# (2) copy relevant bed files to a single, new directory
#

# first experiment
exp1_path="/groups/gaidt/user/luisa.krumwiede/01-BioInfo_Storage/ChIP/ChIP_Experiment_#3_ChIPSeq1_LK/Analysis1_241119/results_nf-core-pipeline-output/bwa/mergedLibrary"
cp -u $exp1_path/macs/broadPeak/sgAAVS1_DMSO_R1_peaks.broadPeak $tmp_dir/data/bed/
cp -u $exp1_path/macs/broadPeak/sgAAVS1_DMSO_R2_peaks.broadPeak $tmp_dir/data/bed/

# second experiment (same number of "_" in all filenames)
exp2_path="/groups/gaidt/user/luisa.krumwiede/01-BioInfo_Storage/ChIP/ChIP_R19299_ChIPSeq2_LK/Analysis_MORC3-DAXX_250820/results/bwa/merged_library"
cp -u $exp2_path/macs3/broad_peak/MORC3-V5_DMSO_M3AID_REP1_peaks.broadPeak $tmp_dir/data/bed/MORC3-V5_DMSO-M3AID_REP1_peaks.broadPeak
cp -u $exp2_path/macs3/broad_peak/MORC3-V5_DMSO_M3AID_REP2_peaks.broadPeak $tmp_dir/data/bed/MORC3-V5_DMSO-M3AID_REP2_peaks.broadPeak

# blacklist and Danbing VNTRs
cp -u /groups/gaidt/bioinf/genomes/GRCh38/hg38-blacklist.v2.bed $tmp_dir/data/
cp -u /groups/gaidt/bioinf/software/danbing-tk/2025May15/QC.20250415a.bed $tmp_dir/data/

# check
wc -l $tmp_dir/data/bed/*.broadPeak
wc -l $tmp_dir/data/*.bed

#
# (3) create bed_samples.txt file that contains the sample names in one column (without path or extension)
#

# - remove all characters after first dot in filename (removes extension)
# - do not copy everytime again, otherwise snakemake will rerun everytime
if [ ! -e $tmp_dir/data/bed_samples.txt ]
then
	ls -1 $tmp_dir/data/bed | sed 's/\..*//g' > $tmp_dir/data/bed_samples.txt
fi

# check
cat $tmp_dir/data/bed_samples.txt

#
# (4) cp relevant bam files to a single, new directory
#

# cp .bam, .bam.bai (-u : only if updated/newer)
cp -u $exp1_path/sgAAVS1_DMSO_R1.mLb.clN.sorted.bam* $tmp_dir/data/bam/
cp -u $exp1_path/sgAAVS1_DMSO_R2.mLb.clN.sorted.bam* $tmp_dir/data/bam/
cp -u $exp2_path/MORC3-V5_DMSO_M3AID_REP1.mLb.clN.sorted.bam $tmp_dir/data/bam/MORC3-V5_DMSO-M3AID_REP1.mLb.clN.sorted.bam
cp -u $exp2_path/MORC3-V5_DMSO_M3AID_REP2.mLb.clN.sorted.bam $tmp_dir/data/bam/MORC3-V5_DMSO-M3AID_REP2.mLb.clN.sorted.bam
cp -u $exp2_path/MORC3-V5_DMSO_M3AID_REP1.mLb.clN.sorted.bam.bai $tmp_dir/data/bam/MORC3-V5_DMSO-M3AID_REP1.mLb.clN.sorted.bam.bai
cp -u $exp2_path/MORC3-V5_DMSO_M3AID_REP2.mLb.clN.sorted.bam.bai $tmp_dir/data/bam/MORC3-V5_DMSO-M3AID_REP2.mLb.clN.sorted.bam.bai

# check
ls -lh $tmp_dir/data/bam/

#
# (5) create samples.txt file that contains the sample names in one column (without path or extension)
#

# - remove all characters after first dot in filename (removes extension)
if [ ! -e $tmp_dir/data/bam_samples.txt ]
then
	ls -1 $tmp_dir/data/bam | grep ".bam$" | sed 's/\..*//g' > $tmp_dir/data/bam_samples.txt
fi

# check
cat $tmp_dir/data/bam_samples.txt

#
# (6) cp associated bigwig files to results directory for visualization
#
cp -u $exp1_path/bigwig/sgAAVS1_DMSO_R1.mLb.clN.bigWig $tmp_dir/data/bigwig/
cp -u $exp1_path/bigwig/sgAAVS1_DMSO_R2.mLb.clN.bigWig $tmp_dir/data/bigwig/
cp -u $exp2_path/bigwig/MORC3-V5_DMSO_M3AID_REP1.mLb.clN.bigWig $tmp_dir/data/bigwig/MORC3-V5_DMSO-M3AID_REP1.mLb.clN.bigWig
cp -u $exp2_path/bigwig/MORC3-V5_DMSO_M3AID_REP2.mLb.clN.bigWig $tmp_dir/data/bigwig/MORC3-V5_DMSO-M3AID_REP2.mLb.clN.bigWig

# check
ls -lh $tmp_dir/data/bigwig/
