#!/bin/bash
## Job Name
#SBATCH --job-name=SubsamplePicardPact
##  This script is meant to align trim galore data to both genomes
##  And generate files for downstream analyses
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes (We only get 1, so this is fixed)
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-00:00:00
## Memory per node
#SBATCH --mem=500G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/strigg/analyses/20200513


#Pact dedup alignments

find /gscratch/scrubbed/strigg/analyses/20200507/Pact_tg/dedup/*.sorted.bam|\
xargs basename -s .sorted.bam \
| xargs -I{} java -jar /gscratch/srlab/programs/picard-2.9.1/picard.jar CollectWgsMetrics \
I=/gscratch/scrubbed/strigg/analyses/20200507/Pact_tg/dedup/{}.sorted.bam \
O={}.collect_wgs_metrics.txt \
R=/gscratch/srlab/sr320/data/froger/Pact_Genome/Pocillopora_acuta_genome_v1.fasta

#Pact nodedup alignments

find /gscratch/scrubbed/strigg/analyses/20200507/Pact_tg/nodedup/*.sorted.bam|\
xargs basename -s .sorted.bam \
| xargs -I{} java -jar /gscratch/srlab/programs/picard-2.9.1/picard.jar CollectWgsMetrics \
I=/gscratch/scrubbed/strigg/analyses/20200507/Pact_tg/nodedup/{}.sorted.bam \
O={}.collect_wgs_metrics.txt \
R=/gscratch/srlab/sr320/data/froger/Pact_Genome/Pocillopora_acuta_genome_v1.fasta