#!/bin/bash

#SBATCH -t 72:00:00
#SBATCH --mem=96G
#SBATCH -c 16
#SBATCH --array=1-1%1               # eg 1-56%10 (job 1 to 56, with 10 at a time)

# Change these settings for different runs! Change SBATCH --array to the number
# of samples you have

# Error report
set -euo pipefail

# Working directory
wd="$1"

# Set resources
custom_bed=#Location of custom bed
ref_gtf=#Location of reference gtf
ref_genome=#Location of reference genome

# Assign names for arrays
cd "${wd}/raw"
names=($(cat jobs))
selected_sample=${names[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected sample: ${selected_sample}"

# Create analysis dir and underlying dirs for sample
cd "${wd}/analysis/"
mkdir -p "${selected_sample}"/CIRIquant

# Set resources (Quantification)
CIRIquant_config=/hpc/local/CentOS7/pmc_vanheesch/pipelines/detect_circrna/CIRIquant_config.yaml

# Load required tools (Quantification)
module load CIRIquant/1.1.2

# Run CIRIquant
echo -e  "\n`date` Running CIRIquant for ${selected_sample}"
CIRIquant -t 16 \
	-1 "${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_1_val_1.fq.gz" \
    -2 "${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_2_val_2.fq.gz" \
    --config ${CIRIquant_config} \
    -o "${wd}/analysis/${selected_sample}/CIRIquant" \
    --bed ${custom_bed}

echo -e "\n`date` Finished!"