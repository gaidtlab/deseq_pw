import os
import re

#
# Config parameters
#

# all input bedfiles must be in one input directory
BED_INDIR = config.get("bed_indir", "bed_input")
BAM_INDIR = config.get("bam_indir", "bam_input")
BAM_EXT = config.get("bam_ext", ".sorted.bam")
BLACKLIST = config.get("blacklist","blacklist_input")
OUTDIR = config.get("outdir", "output")
THREADS = config.get("threads", 4)

BED_SAMPLES = config.get("bed_samples","bed_samples.txt")
bed_sample_table = [line.strip().split() for line in open(BED_SAMPLES) if line.strip() and not line.startswith("#")]
BED_SAMPLES = [row[0] for row in bed_sample_table]

BAM_SAMPLES = config.get("bam_samples","bam_samples.txt")
bam_sample_table = [line.strip().split() for line in open(BAM_SAMPLES) if line.strip() and not line.startswith("#")]
BAM_SAMPLES = [row[0] for row in bam_sample_table]

wildcard_constraints:
	bed_indir = BED_INDIR,
	bam_indir = BAM_INDIR,
	bam_ext = re.escape(BAM_EXT),
	outdir = OUTDIR,
	bed_sample = '|'.join([re.escape(x) for x in BED_SAMPLES]), # only SAMPLE_LIST samples 
	bam_sample = '|'.join([re.escape(x) for x in BAM_SAMPLES]), # only SAMPLE_LIST samples 

rule all:
	input:
		i1 = expand("{outdir}/featureCount_concat/featureCount_counts.txt",outdir=OUTDIR),	

rule concat_featurecount:
	input:
		i1 = expand("{outdir}/featureCount/{bam_sample}_featureCount.txt",outdir=OUTDIR,bam_sample=BAM_SAMPLES),
	output:
		o1 = "{outdir}/featureCount_concat/featureCount_counts.txt",
		o2 = "{outdir}/featureCount_concat/featureCount_counts.annot.txt",
	benchmark:
			"{outdir}/featureCount_concat/featureCount_counts.benchmark.txt",
	threads: 1
	resources:
		mem_mb = 10000,
		runtime = 10,
		cpus_per_task = 1,
		tasks = 1,
	log:
		log =  "{outdir}/featureCount_concat/featureCount_counts.log",
	shell:
		"""
		# get annotation from first file
		tail -n +2 {input.i1[0]} | cut -f 1-6 > {output.o2}
  
		# get the merged peaks names
		tail -n +3 {input.i1[0]} | cut -f 1 > {output.o1}.ids

		# extract only counts
		for r in {input.i1} 
		do
			dn=$(dirname $r)
			bn=$(basename $r .txt)
			tail -n +3 $r | cut -f 7 > $dn/$bn.counts.txt
		done

		# column names (sites+samples)
		fnames_bn=$(for r in {input.i1}; do bn=$(basename $r .txt); echo $bn; done)
		echo "sites "$fnames_bn | sed 's/ /\t/g' > {output.o1}.header
		
		# paste body (sites + counts)
		fnames=$(for r in {input.i1}; do dn=$(dirname $r); bn=$(basename $r .txt); echo $dn/$bn.counts.txt; done)
		paste {output.o1}.ids $fnames > {output.o1}.body
  
		cat {output.o1}.header {output.o1}.body > {output.o1}
		rm {output.o1}.ids {output.o1}.body {output.o1}.header
		"""

rule featurecount:
	input:
		i1 = expand("{bam_indir}/{{bam_sample}}{bam_ext}",bam_indir=BAM_INDIR,bam_ext=BAM_EXT),
		i2 = expand("{outdir}/merged_bed/merged.noblacklist.gtf",outdir=OUTDIR),
	output:
		o1 = "{outdir}/featureCount/{bam_sample}_featureCount.txt",
	benchmark:
			"{outdir}/featureCount/{bam_sample}_featureCount.benchmark.txt",
	threads: THREADS
	resources:
		mem_mb = 20000,
		runtime = 30,
		cpus_per_task = THREADS,
		tasks = 1,
	log:
		log =  "{outdir}/featureCount/{bam_sample}_featureCount.log",
	shell:
		"""
		module load build-env/.f2021 build-env/f2021
		module load subread/2.0.2-gcc-10.2.0

		featureCounts -a {input.i2} -o {output.o1} -F GTF -t exon -g gene_id -M --fraction -p -B -T {threads} {input.i1} 2>> {log}
		"""

rule bed_to_gtf:
	input:
		i1 = expand("{outdir}/merged_bed/merged.noblacklist.bed",outdir=OUTDIR),
	output:
		o1 = "{outdir}/merged_bed/merged.noblacklist.gtf",
	benchmark:
			"{outdir}/merged_bed/merged_gtf.benchmark.txt",
	resources:
		mem_mb = 10000,
		runtime = 10,
		cpus_per_task = 1,
		tasks = 1,
	log:
		log =  "{outdir}/merged_bed/merged_gtf.log",
	shell:
		"""
		awk 'OFS="\t"{{ print $1,"merged_peaks","exon",$2+1,$3,"0","+",".","gene_id \\"site_"++i"\\";"}}' {input.i1} > {output.o1} 2>> {log}
		"""

rule merge_beds:
	input:
		i1 = expand("{bed_indir}/{bed_sample}.broadPeak",bed_indir=BED_INDIR,bed_sample=BED_SAMPLES),
	output:
		o1 = "{outdir}/merged_bed/merged.bed",
		o2 = "{outdir}/merged_bed/merged.noblacklist.bed",
	benchmark:
			"{outdir}/merged_bed/merged.benchmark.txt",
	params:
		blacklist = BLACKLIST,
	resources:
		mem_mb = 10000,
		runtime = 10,
		cpus_per_task = 1,
		tasks = 1,
	log:
		log =  "{outdir}/merged_bed/merged.log",
	shell:
		"""
		module load build-env/f2022
		module load bedtools/2.31.1-gcc-13.2.0
		
		cat {input.i1} | sort -k1,1 -k2,2n | bedtools merge -i - > {output.o1} 2>> {log}
		bedtools intersect -v -a {output.o1} -b {params.blacklist} > {output.o2} 2>> {log}
		echo "Merged $(cat {input.i1} | wc -l) bed peaks into $(cat {output.o2} | wc -l) merged peaks" >> {log}
		"""