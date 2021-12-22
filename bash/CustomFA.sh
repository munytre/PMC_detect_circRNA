#!/bin/bash

#SBATCH -t 72:00:00
#SBATCH --mem=96G
#SBATCH -c 16
#SBATCH --array=1-24%6               # eg 1-56%10 (job 1 to 56, with 10 at a time)

# Change these settings for different runs! Change SBATCH --array to the number
# of samples you have

# Error report
set -euo pipefail

# Working directory
wd="$1"

# Set resources
ref_STAR=#Location of STAR reference

# Load required tools (Detection)
module load python/3.6.1
module load Java/1.8.0_60
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

# Create analysis dir with underlying dirs for sample
cd "${wd}/analysis/"
mkdir -p "${selected_sample}"

### MATE 1 ###
#Run STAR on unmapped reads
echo -e "\n`date` Mapping ${selected_sample}mate1 with STAR"
STAR --genomeDir ${ref_STAR} \
    --runThreadN 16 \
    --twopassMode None \
    --readFilesIn \
    "${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_1_val_1.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${wd}/analysis/${selected_sample}/${selected_sample}mate1" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --alignEndsType EndToEnd \
    --outFilterMatchNmin 70 \
    --limitBAMsortRAM 150000000000

# Index by samtools and retrieval of counts
echo -e "\n`date` Indexing ${selected_sample}mate1 with samtools"
samtools index -@ 16 "${wd}/analysis/${selected_sample}/${selected_sample}mate1Aligned.sortedByCoord.out.bam"
echo -e "\n`date` Filtering unique reads from ${selected_sample}mate1 with samtools"
samtools view -q 255 -@ 16 -b "${wd}/analysis/${selected_sample}/${selected_sample}mate1Aligned.sortedByCoord.out.bam" > "${wd}/analysis/${selected_sample}/${selected_sample}mate1Unique.sortedByCoord.out.bam"
echo -e "\n`date` Indexing ${selected_sample}mate1 with samtools"
samtools index -@ 16 "${wd}/analysis/${selected_sample}/${selected_sample}mate1Unique.sortedByCoord.out.bam"

### MATE 2 ###
#Run STAR on unmapped reads
echo -e "\n`date` Mapping ${selected_sample}mate2 with STAR"
STAR --genomeDir ${ref_STAR} \
    --runThreadN 16 \
    --twopassMode None \
    --readFilesIn \
    "${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_2_val_2.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${wd}/analysis/${selected_sample}/${selected_sample}mate2" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --alignEndsType EndToEnd \
    --outFilterMatchNmin 70 \
    --limitBAMsortRAM 150000000000

# Index by samtools and retrieval of counts
echo -e "\n`date` Indexing ${selected_sample}mate2 with samtools"
samtools index -@ 16 "${wd}/analysis/${selected_sample}/${selected_sample}mate2Aligned.sortedByCoord.out.bam"
echo -e "\n`date` Filtering unique reads from ${selected_sample}mate2 with samtools"
samtools view -q 255 -@ 16 -b "${wd}/analysis/${selected_sample}/${selected_sample}mate2Aligned.sortedByCoord.out.bam" > "${wd}/analysis/${selected_sample}/${selected_sample}mate2Unique.sortedByCoord.out.bam"
echo -e "\n`date` Indexing ${selected_sample}mate2 with samtools"
samtools index -@ 16 "${wd}/analysis/${selected_sample}/${selected_sample}mate2Unique.sortedByCoord.out.bam"

echo -e "\n`date` Finished!"