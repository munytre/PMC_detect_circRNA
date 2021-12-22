#!/bin/bash

#SBATCH -t 72:00:00
#SBATCH --mem=96G
#SBATCH -c 16
#SBATCH --array=2-24%5               # eg 1-56%10 (job 1 to 56, with 10 at a time)

# Change these settings for different runs! Change SBATCH --array to the number
# of samples you have

# Error report
set -euo pipefail

# Working directory
wd="$1"

# Data directory
dd="$2"

# Set resources
ref_gtf=#Location of reference gtf
ref_STAR=#Location of STAR reference

# Load required tools (Detection)
module load python/3.6.1
module load Java/1.8.0_60
module load trimgalore/0.6.6
module load FastQC/0.11.9
module load cutadapt/3.4
module load samtools/1.12
module load STAR/2.7.8a

# Assign names for arrays
cd "${wd}/raw"
names=($(cat jobs))
selected_sample=${names[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected sample: ${selected_sample}"

# Create processed dir with underlying dirs for sample
cd "${wd}/processed/"
mkdir -p "${selected_sample}"

# Run TrimGalore on the fastq
echo -e "\n`date` Filtering and trimming ${selected_sample} ..."
trim_galore --cores 8 \
	--trim-n \
	--fastqc --fastqc_args "--outdir ${wd}/processed/${selected_sample}/trimgalore/ -t 8" \
	--gzip \
	--paired \
	"${dd}/${selected_sample}_1.fastq.gz" \
	"${dd}/${selected_sample}_2.fastq.gz" \
	-o "${wd}/processed/${selected_sample}/trimgalore"

#Run STAR
echo -e "\n`date` Mapping ${selected_sample} with STAR"
STAR --genomeDir ${ref_STAR} \
    --runThreadN 16 \
    --twopassMode Basic \
    --readFilesIn \
    "${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_1_val_1.fq.gz" \
    "${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_2_val_2.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${wd}/processed/${selected_sample}/STAR/${selected_sample}" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --outReadsUnmapped Fastx \
    --limitBAMsortRAM 150000000000

# Index by samtools
echo -e "\n`date` Indexing ${selected_sample} with samtools"
samtools index -@ 16 "${wd}/processed/${selected_sample}/STAR/${selected_sample}Aligned.sortedByCoord.out.bam"

echo -e "\n`date` Finished!"