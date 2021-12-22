#!/bin/bash

#SBATCH -t 6:00:00

# Error report
set -euo pipefail

# Starting time
echo -e "\n`date` Start"

# Working directory
wd="$1"
cd "${wd}"

# Unzipping
echo -e "\n`date` Unzipping"
gzip -dc neuroblastoma/SAMPLE_1.fastq.gz > neuroblastoma/SAMPLE_1.fastq
gzip -dc neuroblastoma/SAMPLE_2.fastq.gz > neuroblastoma/SAMPLE_2.fastq

# Subsampling with seqtk
echo -e "\n`date` Subsampling"
/hpc/local/CentOS7/pmc_vanheesch/bin/seqtk/seqtk sample -s2021 neuroblastoma/SAMPLE_1.fastq 1000000 > Subsample_1.fastq
/hpc/local/CentOS7/pmc_vanheesch/bin/seqtk/seqtk sample -s2021 neuroblastoma/SAMPLE_2.fastq 1000000 > Subsample_2.fastq

# Deleting raw.fastq
echo -e "\n`date` Deleting unzipped files"
rm neuroblastoma/SAMPLE_1.fastq
rm neuroblastoma/SAMPLE_2.fastq

# Gzipping subsamples
echo -e "\n`date` Gzipping subsamples"
gzip Subsample_1.fastq 
gzip Subsample_2.fastq

echo -e "\n`date` Done with subsampling"