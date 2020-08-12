#!/bin/bash
## Job Name
#SBATCH --job-name=mQC
##  Running multiqc on https://github.com/hputnam/Meth_Compare/issues/17
## Allocation Definition
#SBATCH --account=coenv
#SBATCH --partition=coenv
## Resources
## Nodes (We only get 1, so this is fixed)
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-00:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr320@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/sr320/032320-mqc



source /gscratch/srlab/programs/scripts/paths.sh



/gscratch/srlab/programs/anaconda3/bin/multiqc \
-o Mcap-multiqc \
/gscratch/scrubbed/sr320/031520-TG-100g/Mcap_tg/



/gscratch/srlab/programs/anaconda3/bin/multiqc \
-o Pact-multiqc \
/gscratch/scrubbed/sr320/031520-TG-100g/Pact_tg/
