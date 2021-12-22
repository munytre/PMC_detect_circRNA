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

# Data directory
dd="$2"

# Set resources
ref_gtf=#Location of reference gtf
ref_genome=#Location of reference genome

# Load required tools (Detection)
module load python/3.6.1
module load Java/1.8.0_60
module load trimgalore/0.6.6
module load bwa/0.7.17
module load FastQC/0.11.9
module load cutadapt/3.4
CIRI2=/hpc/local/CentOS7/pmc_vanheesch/software/CIRI2-2.0.6/CIRI2.pl

# Assign names for arrays
cd "${wd}/raw"
names=($(cat jobs))
selected_sample=${names[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected sample: ${selected_sample}"

# Create processed dir with underlying dirs for sample
cd "${wd}/processed/"
mkdir -p "${selected_sample}"/{trimgalore,bwa}

# Create analysis dir and underlying dirs for sample
cd "${wd}/analysis/"
mkdir -p "${selected_sample}"/{CIRI,CIRIquant}

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

# Run bwa index to prepare for bwa mem if not done yet
cd "${wd}"
if [[ -f ${ref_genome}.amb ]] && [[ -f ${ref_genome}.ann ]] && [[ -f ${ref_genome}.bwt ]] && [[ -f ${ref_genome}.pac ]] && [[ -f ${ref_genome}.sa ]]
then
	echo -e "\n`date` Indexing with bwa already done, skipping to Mapping!"
else
	echo -e "\n`date` Indexing reference genome ${ref_genome}"
	bwa index -a bwtsw ${ref_genome}
fi

# Run bwa mem to prepare files for CIRI2
echo -e "\n`date` Mapping ${selected_sample} to reference genome ${ref_genome}"
bwa mem \
	-T 19 \
	-t 16 \
	${ref_genome} \
	"${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_1_val_1.fq.gz" \
	"${wd}/processed/${selected_sample}/trimgalore/${selected_sample}_2_val_2.fq.gz" \
	> "${wd}/processed/${selected_sample}/bwa/${selected_sample}_PE.sam"

# Run CIRI2
echo -e  "\n`date` Running CIRI2 for ${selected_sample}"
perl ${CIRI2} -T 12\
	-I "${wd}/processed/${selected_sample}/bwa/${selected_sample}_PE.sam" \
	-A ${ref_gtf} \
	-F ${ref_genome} \
	-O "${wd}/analysis/${selected_sample}/CIRI/${selected_sample}"

# New settings!
echo -e  "\n`date` Settings are changed for CIRIquant!"

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
    --circ "${wd}/analysis/${selected_sample}/CIRI/${selected_sample}" \
    --tool CIRI2

# Load samtools
module load samtools/1.12

# Converting sam to bam (Save Storage)
echo -e "\n`date` Converting sam to bam to save space!"
cd ${wd}/processed/${selected_sample}/bwa/
samtools view -S -b ${selected_sample}_PE.sam -o ${selected_sample}_PE.bam
rm ${selected_sample}_PE.sam

echo -e "\n`date` Finished!"