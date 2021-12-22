#!/bin/bash

#SBATCH -t 72:00:00
#SBATCH --mem=96G
#SBATCH -c 16
#SBATCH --array=210-210%1

# Error report
set -euo pipefail

# Working directory
wd=$1

# Data directory
dd="$2"

# Set resources
ref_gtf=#Location of reference gtf
ref_genome=#Location of reference genome
ref_transcriptome=#Location of reference transcriptome
salmon_index=#Location of salmon index

# Load required tools (Detection)
module load salmon/1.4.0
module load trimgalore/0.6.6
module load cutadapt/3.4

# Assign names for arrays
cd "${wd}/raw"
names=($(cat jobs))
selected_sample=${names[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected sample: ${selected_sample}"

### ONLY IF GENTROME NEEDS TO BE GENERATED ###
# # Create "gentrome" for salmon indexing
# echo -e "\n`date` Creating Gentrome"
# gentrome="${wd}/processed/annotation/gentrome.fa"
# cat $ref_transcriptome ${ref_genome} > ${gentrome}

# # Grep seq names from fa and create decoys for salmon
# echo -e "\n`date` Creating decoy file"
# decoy="${wd}/processed/annotation/decoys.txt"
# grep "^>" ${ref_genome} | cut -d " " -f 1 > ${decoy}
# sed -i.bak -e 's/>//g' ${decoy}

# # Making the index for Salmon
# echo -e "\n`date` Running Salmon Index"
# salmon index \
#   --transcripts ${gentrome} \
#   --decoys ${decoy} \
#   --index "/hpc/pmc_vanheesch/shared_resources/GENOMES/Homo_sapiens.GRCh38/102/salmon/1.4.0/salmon_index" \
#   --threads 16
### ONLY IF GENTROME NEEDS TO BE GENERATED ###

# Create processed dir with underlying dirs for sample
cd "${wd}/processed/"
mkdir -p "${selected_sample}"/trimgalore
cd "${wd}/analysis/"
mkdir -p "${selected_sample}"/salmon/

# Trimming selected sample
echo -e "\n`date` Trimming ${selected_sample} ..."
trim_galore --cores 8 \
	--trim-n \
	--gzip \
	--paired \
	"${dd}/${selected_sample}_R1.fastq.gz" \
	"${dd}/${selected_sample}_R2.fastq.gz" \
	-o "${wd}/processed/${selected_sample}/trimgalore"

# Run salmon for transcript counts
echo -e "\n`date` Running Salmon quant on ${selected_sample} ..."
salmon quant \
  --index ${salmon_index} \
  --libType A \
  --mates1 "${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_R1_val_1.fq.gz" \
  --mates2 "${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_R2_val_2.fq.gz" \
  --threads 16 \
  --gcBias \
  --validateMappings \
  --output "${wd}/analysis/${selected_sample}/salmon/"

# Remove intermediate files
echo -e "\n`date` Removing intermediate files: ${selected_sample} ..."
rm -r "${wd}/processed/${selected_sample}"

echo -e "\n`date` Finished!"