# Methylation Quantification pipeline


[NF Core Methylseq](https://github.com/nf-core/methylseq/)  
NEXTFLOW  ~  version 20.04.1  
nf-core/methylseq v1.5  

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=400GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/hputnam/Meth_Compare
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

module load Nextflow

# run nextflow methylseq

nextflow run nf-core/methylseq -profile singularity \
--aligner bismark \
--fasta /data/putnamlab/REFS/Mcap/Mcap.genome_assembly.fa \
--save_reference \
--reads '/data/putnamlab/hputnam/Meth_Compare/raw/*_R{1,2}_001.fastq.gz' \
--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--bismark_align_cpu_per_multicore \
--outdir Mcap_Meth_Results \
-name Mcap_Meth_Results
```



# P. acuta

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=400GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/hputnam/Meth_Compare
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

module load Nextflow

# run nextflow methylseq

nextflow run nf-core/methylseq -profile singularity \
--aligner bismark \
--fasta /data/putnamlab/REFS/Pact/Pocillopora_acuta_genome_v1.fasta \
--save_reference \
--reads '/data/putnamlab/hputnam/Meth_Compare/Pact_raw/*_R{1,2}_001.fastq.gz' \
--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--bismark_align_cpu_per_multicore \
--outdir Pact_Meth_Results \
-name Pact_Meth_Results
```
