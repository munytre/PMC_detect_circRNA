#!/bin/bash

#SBATCH -t 72:00:00
#SBATCH --mem=16G
#SBATCH -c 16

# Error report
set -euo pipefail

# Working directory
wd=$1

# Set resources
ref_gtf=#Location of reference gtf

# Load required tools (Detection)
module load subread/2.0.2

# Counting
featureCounts -v
cd "${wd}/processed"
echo -e "\n`date` Counting all samples (STAR) with featureCounts"
featureCounts -a ${ref_gtf} \
    -o STAR_counts_38 \
    ${wd}/processed/*/STAR/*Aligned.sortedByCoord.out.bam \
    -p \
    -T 16

echo -e "\n`date` Finished!"