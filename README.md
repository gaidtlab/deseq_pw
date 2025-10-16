*deseq_pw* is a protocol to perform counting of reads and subsequence differential analysis.

## Overview

```
1. count_beds_setup.sh               # stage the data, cp all needed input fules to scratch-cbe
2. count_beds.yaml and count_beds.sh # count the reads in the bed files, using snakemake workflow
3. rsync                             # copy results back to group drive using rsync
4. samples_experiment.txt and chip_compare.sh # will run chip_compare.Rmd and performs DESEq2 analysis
```

Prerequisite:

```
ssh -X cbe.vbc.ac.at
srun -n 1 --mem=20g --time=8:00:00 --reservation=interactive --pty bash
export PATH=$PATH:/groups/gaidt/bioinf/software/scripts/deseq_pw

# Maybe some R libraries need to be installed.
# The easiest ist to start RStudio (version 4.3.2, https://jupyterhub.vbc.ac.at/)
# and install apeglm, EDASeq, DGEobj.utils (all BioConductor).
```

# (1) Prepare (stage) the data

We need to achieve the following directory structure, optimally somewhere on scratch-cbe:

- data/bam
- data/bed
- data/bigwig

and fill it with the associated peak bedfiles, aligned read bam files, and bigwig files.
In addition, we need two text files listing the sample names (without extensions):

- data/bam_samples.txt 
- data/bed_samples.txt

and the blacklist file and all VNTRs (for viewing in IGV lateron):

- data/hg38-blacklist.v2.bed
- data/QC.20250415a.bed

## Example (edit according to your needs:) bash script:

```
count_beds_setup.sh
```

results of this example are in /scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/ :

```
cat ./data/bam_samples.txt :

MORC3-V5_DMSO-M3AID_REP1  # changed from MORC3-V5_DMSO_M3AID_REP1 to MORC3-V5_DMSO-M3AID_REP1
MORC3-V5_DMSO-M3AID_REP2  # changed	from MORC3-V5_DMSO_M3AID_REP2 to MORC3-V5_DMSO-M3AID_REP2
sgAAVS1_DMSO_R1
sgAAVS1_DMSO_R2

cat ./data/bed_samples.txt :

MORC3-V5_DMSO-M3AID_REP1_peaks  # changed
MORC3-V5_DMSO-M3AID_REP2_peaks  # changed
sgAAVS1_DMSO_R1_peaks
sgAAVS1_DMSO_R2_peaks

./data/bed/sgAAVS1_DMSO_R1_peaks.broadPeak
./data/bed/sgAAVS1_DMSO_R2_peaks.broadPeak
./data/bed/MORC3-V5_DMSO-M3AID_REP1_peaks.broadPeak
./data/bed/MORC3-V5_DMSO-M3AID_REP2_peaks.broadPeak

./data/hg38-blacklist.v2.bed
./data/QC.20250415a.bed

./data/bam/MORC3-V5_DMSO-M3AID_REP2.mLb.clN.sorted.bam
./data/bam/MORC3-V5_DMSO-M3AID_REP2.mLb.clN.sorted.bam.bai
./data/bam/MORC3-V5_DMSO-M3AID_REP1.mLb.clN.sorted.bam
./data/bam/MORC3-V5_DMSO-M3AID_REP1.mLb.clN.sorted.bam.bai
./data/bam/sgAAVS1_DMSO_R1.mLb.clN.sorted.bam
./data/bam/sgAAVS1_DMSO_R1.mLb.clN.sorted.bam.bai
./data/bam/sgAAVS1_DMSO_R2.mLb.clN.sorted.bam
./data/bam/sgAAVS1_DMSO_R2.mLb.clN.sorted.bam.bai

./data/bigwig/sgAAVS1_DMSO_R2.mLb.clN.bigWig
./data/bigwig/sgAAVS1_DMSO_R1.mLb.clN.bigWig
./data/bigwig/MORC3-V5_DMSO-M3AID_REP1.mLb.clN.bigWig
./data/bigwig/MORC3-V5_DMSO-M3AID_REP2.mLb.clN.bigWig
```

# (2) Count the reads

Edit paths in this script according to your needs:

```
cat count_beds.yaml:

bed_indir: "/scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/data/bed"
bam_indir: "/scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/data/bam"
outdir: "/scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/results"
bed_samples: "/scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/data/bed_samples.txt"
bam_samples: "/scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/data/bam_samples.txt"
bam_ext: ".mLb.clN.sorted.bam"
blacklist: "/scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/data/hg38-blacklist.v2.bed"
threads: 4
```

go to your final results directory and launch:

```
count_beds.sh -h
Usage: ./count_beds.sh [options]
Options:
  -c, --config FILE     Specify the configuration file in YAML format (default: count_beds.yaml)
  -a, --smk_args ARGS   Additional snakemake arguments (default: ''),  --dry-run, --forceall, etc.
  -p, --processes N     Number of processes/cores to use (default: 4)
  -l, --local           Run locally without cluster support
  -h, --help            Show this help message and exit
```

```
count_beds.sh -c count_beds.yaml -a --dry-run -p 4 -l # local, but dry run (not executed, just steps are shown)
count_beds.sh -c count_beds.yaml -a --dry-run         # cluster, but dry run (not executed, just steps are shown)
count_beds.sh -c count_beds.yaml                      # cluster, executed
```

results of this example are in /scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/ :

```
./results/merged_bed/merged.bed
./results/merged_bed/merged.benchmark.txt
./results/merged_bed/merged.gtf
./results/merged_bed/merged.noblacklist.gtf
./results/merged_bed/merged.noblacklist.bed

./results/featureCount/sgAAVS1_DMSO_R1_featureCount.txt.summary
./results/featureCount/sgAAVS1_DMSO_R1_featureCount.counts.txt
./results/featureCount/sgAAVS1_DMSO_R1_featureCount.txt
./results/featureCount/sgAAVS1_DMSO_R2_featureCount.txt.summary
./results/featureCount/sgAAVS1_DMSO_R2_featureCount.counts.txt
./results/featureCount/sgAAVS1_DMSO_R2_featureCount.txt
./results/featureCount/MORC3-V5_DMSO_M3AID_REP1_featureCount.txt
./results/featureCount/MORC3-V5_DMSO_M3AID_REP1_featureCount.counts.txt
./results/featureCount/MORC3-V5_DMSO_M3AID_REP1_featureCount.txt.summary
./results/featureCount/MORC3-V5_DMSO_M3AID_REP2_featureCount.txt
./results/featureCount/MORC3-V5_DMSO_M3AID_REP2_featureCount.counts.txt
./results/featureCount/MORC3-V5_DMSO_M3AID_REP2_featureCount.txt.summary

./results/featureCount_concat/featureCount_counts.txt
./results/featureCount_concat/featureCount_counts.annot.txt

```

# (3) Copy the results to the group drive

```
# rsync -auvrn --exclude "*.bam" [results_dir]/* [target_dir_on_group_drive]/

# dry run:
#rsync -auvrn --exclude "*.bam" /scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/* /groups/gaidt/user/markus.jaritz/projects/Markus/2025Oct14/
rsync -auvrn --exclude "*.bam" /scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/* /groups/gaidt/user/markus.jaritz/projects/Markus/2025Oct14b/

# execute
rsync -auvr --exclude "*.bam" /scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/* /groups/gaidt/user/markus.jaritz/projects/Markus/2025Oct14/
rsync -auvr --exclude "*.bam" /scratch-cbe/users/markus.jaritz/Gaidt/Moritz/2025Oct14b/* /groups/gaidt/user/markus.jaritz/projects/Markus/2025Oct14b/
```

# (4) Run DESeq2 differential analysis

We need a table that holds information about the samples/conditions to be compared.

```
cat samples_experiment.txt:

sample_id	sample_name						count_column							condition
1			MORC3_V5_DMSO-M3AID_REP1_exp2	MORC3-V5_DMSO-M3AID_REP1_featureCount	exp2
2			MORC3_V5_DMSO-M3AID_REP2_exp2	MORC3-V5_DMSO-M3AID_REP2_featureCount	exp2
3			sgAAVS1_DMSO_R1_exp1			sgAAVS1_DMSO_R1_featureCount			exp1
4			sgAAVS1_DMSO_R2_exp1			sgAAVS1_DMSO_R2_featureCount			exp1
```

```
chip_compare.sh -h

Usage: chip_compare.sh options

 Script to compare two ChIP-seq experiments based on peak counts
 Requires a count table (featureCounts output) and a sample table
 Generates an HTML report with QC and differential analysis

-s (sampleTable): Tab-delimited file with columns: sample_id, sample_name, count_column, condition [samples_experiment.txt]
-c (conditions_colname): Column name in sampleTable that defines conditions [condition]
-l (control_condition): Name of the control condition in conditions_colname [exp1]
-e (experiment_condition): Name of the experimental condition in conditions_colname [exp2]
-t (countTable): Path to the count table (featureCounts output) [featureCount_counts.txt]
-a (countTableAnnot): Path to the annotated count table [featureCount_counts.annot.txt]
-o (odir): Output directory for the report [.]
-p (downsample_pairs): Number of pairs for downsampling analysis [1000]
-m (downsample_ma): Number of reads for MA downsampling [0]
-i (session): RStudio session URL [http://localhost:60151/goto?]
-h: Show this help message
```

```
chip_compare.sh -s /groups/gaidt/bioinf/software/scripts/deseq_pw/samples_experiment.txt -c condition -l exp1 -e exp2 -t /groups/gaidt/user/markus.jaritz/projects/Markus/2025Oct14b/results/featureCount_concat/featureCount_counts.txt -a /groups/gaidt/user/markus.jaritz/projects/Markus/2025Oct14b/results/featureCount_concat/featureCount_counts.annot.txt -o /groups/gaidt/user/markus.jaritz/projects/Markus/2025Oct14b/doc
```
# (5) Explore results

- see chip_compare.html
- launch IGV, load bigwigs and bedfiles (merged.noblacklist.bed, blacklist, vntrs)
- use hyperlinked differential peak lists for browsing