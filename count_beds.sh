#!/bin/bash
set -o nounset  
set -o errexit  
set -o errtrace 

module load build-env/f2022
module load snakemake/9.5.1-foss-2023b

args="--configfile count_beds.yaml --cores 4 --show-failed-logs --snakefile ./count_beds.smk --rerun-triggers mtime --rerun-incomplete"
cluster_args="--executor slurm --default-resources --cluster-generic-submit-cmd sbatch" 

# -R [rule]:  forced restart from certain rule
# -n : no run, only print what would be done
# --forceall :  force all rules to be executed, even if the output files are present and up to date
#
# no cluster
#
#cmd="snakemake -n --jobs 42 $args" # -R merge_beds" #  -R featurecount " #--allowed-rules concat_featurecount --keep-incomplete --keep-going
#cmd="snakemake --jobs 42 $args" # -R merge_beds" # featurecount " #--allowed-rules concat_featurecount all --keep-incomplete --keep-going"

#
# cluster
#
#cmd="snakemake -n --jobs 42 $args $cluster_args"
cmd="snakemake --jobs 42 $args $cluster_args --forceall"

echo $cmd; $cmd

