#!/usr/bin/env bash

set -u
set -e

# @description: Run snakemake workflow to count reads in bed files
# @ params c configFile = YAML config file
# @ params a smk_args = additional snakemake arguments, e.g. --dry-run, --forceall
# @ params p processes = number of processes/cores to use (local)
# @ params l local = run locally (no cluster)
# @ params h help = print this help message

function print_usage() {
	echo "Usage: $0 [options]"
	echo "Options:"
	echo "  -c, --config FILE     Specify the configuration file in YAML format (default: count_beds.yaml)"
	echo "  -a, --smk_args ARGS   Additional snakemake arguments (default: ''),  --dry-run, --forceall, etc."
	echo "  -p, --processes N     Number of processes/cores to use (default: 4)"
	echo "  -l, --local           Run locally without cluster support"
	echo "  -h, --help            Show this help message and exit"
	echo ""

	exit 1
}


######
# MAIN
######

if [[ $# == 0 ]]; then
    print_usage
    exit 1        # Exit and explain usage, if no argument(s) given.
fi

while getopts 'c:a:p:hl' option; do
	case $option in
		c ) config=$OPTARG ;;
		a ) smk_args=$OPTARG ;;
		p ) processes=$OPTARG ;;
		l ) local=1 ;;
		h ) print_usage ;;
		* ) echo -e "\nERROR: Unimplemented option chosen!\n"
			print_usage;;
	esac
done

# set defaults
configFile=${config:=count_beds.yaml}
smk_args=${smk_args:=""}
processes=${processes:=4}

args="--configfile $config --cores $processes --show-failed-logs --snakefile /groups/gaidt/bioinf/software/scripts/deseq_pw/count_beds.smk --rerun-triggers mtime --rerun-incomplete $smk_args"

# print set parameters
>&2 echo "-c configFile=$configFile"
>&2 echo "-a smk_args=$smk_args"
>&2 echo "-p processes=$processes"
if [ -z ${local+x} ]; then
	>&2 echo "-l local=0"
	cluster_args="--executor slurm --default-resources --cluster-generic-submit-cmd sbatch" 
else
	>&2 echo "-l local=1"
	cluster_args=""
fi


echo "loading modules ..."
module load build-env/f2022
module load snakemake/9.5.1-foss-2023b
echo "loading modules done."

echo "running snakemake ..."
cmd="snakemake --jobs 42 $args $cluster_args"
echo $cmd
$cmd

