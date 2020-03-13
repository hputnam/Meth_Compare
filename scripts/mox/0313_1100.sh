#!/bin/bash
## Job Name
#SBATCH --job-name=tg-bismark
##  This script is meant to align trim galore data to both genomes
##  And generate files for downstream analyses
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes (We only get 1, so this is fixed)
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=5-00:00:00
## Memory per node
#SBATCH --mem=500G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr320@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/sr320/$$$$$$


# Directories and programs
bismark_dir="/gscratch/srlab/programs/Bismark-0.21.0"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
samtools="/gscratch/srlab/programs/samtools-1.9/samtools"

source /gscratch/srlab/programs/scripts/paths.sh


# ${bismark_dir}/bismark_genome_preparation \
# --verbose \
# --parallel 28 \
# --path_to_aligner ${bowtie2_dir} \
# ${genome_folder}



# Aligning to genomes. Will need trim reads in distinct directories
# Need to confirm filenames jive with code.
# AND set output directory

genome_folder="/gscratch/srlab/sr320/data/froger/Mcap_Genome/"
reads_dir="/gscratch/srlab/sr320/data/froger/trim/Mc/"


find ${reads_dir}*2020*_R1_001.fastq.gz \
| xargs basename -s _R1_001.fastq.gz | xargs -I{} ${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
-genome ${genome_folder} \
-p 4 \
-score_min L,0,-0.6 \
--non_directional \
-1 ${reads_dir}{}_R1_001.fastq.gz \
-2 ${reads_dir}{}_R2_001.fastq.gz \
-o Mcap_tg


genome_folder="/gscratch/srlab/sr320/data/froger/Pact_Genome/"
reads_dir="/gscratch/srlab/sr320/data/froger/trim/Pa/"


find ${reads_dir}*2020*_R1_001.fastq.gz \
| xargs basename -s _R1_001.fastq.gz | xargs -I{} ${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
-genome ${genome_folder} \
-p 4 \
-score_min L,0,-0.6 \
--non_directional \
-1 ${reads_dir}{}_R1_001.fastq.gz \
-2 ${reads_dir}{}_R2_001.fastq.gz \
-o Pact_tg



reads_dir="/gscratch/srlab/sr320/data/froger/trim/Pa/"
genome_folder="/gscratch/srlab/sr320/data/froger/Pdam_Genome/"



find ${reads_dir}*2020*_R1_001.fastq.gz \
| xargs basename -s _R1_001.fastq.gz | xargs -I{} ${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
-genome ${genome_folder} \
-p 4 \
-score_min L,0,-0.6 \
-u 10000000 \
--non_directional \
-1 ${reads_dir}{}_R1_001.fastq.gz \
-2 ${reads_dir}{}_R2_001.fastq.gz \
-o Pdam_full-u1M



reads_dir="/gscratch/scrubbed/samwhite/outputs/20200305_methcompare_fastp_trimming/"
genome_folder="/gscratch/srlab/sr320/data/lambda/"


find ${reads_dir}*2020*_R1_001.fastq.gz \
| xargs basename -s _R1_001.fastq.gz | xargs -I{} ${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
-genome ${genome_folder} \
-p 4 \
-score_min L,0,-0.6 \
-u 10000000 \
--non_directional \
-1 ${reads_dir}{}_R1_001.fastq.gz \
-2 ${reads_dir}{}_R2_001.fastq.gz \
-o lambda_full-u1M
