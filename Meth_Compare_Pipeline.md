# Methods
## RRBS Library Prep
https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/RRBS-Meth-Comp/

## MBDBS Enrichment
https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/MBD-Meth-Comp/

## WGBS and MBDBS Library Prep
https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/WGBS-MC-1/
https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/meth-comp-PMS/
https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Redos-PMS/

# Coral WGBS MBD-BS and RRBS data

Bismark Bisulfite Mapper VX “ map bisulfite treated sequencing reads to a genome of interest and perform methylation calls in a single step”

Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications

requires bowtie2, samtools, perl, trimmomatic

* Bowtie 2 version X by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
* Program: samtools (Tools for alignments in the SAM format) Version: X
* perl X 
* Bismark Bisulfite Mapper VX 
* Trimmomatic
* multiqc

## Obtain Genome files and Run Bismark Genome preparation
``mkdir GENOME``  
``cd GENOME``
``mkdir Mcap_Genome`` 
``cd Mcap_Genome`` 
``wget http://cyanophora.rutgers.edu/montipora/Mcap.genome_assembly.fa.gz`` 
``mkdir Pact_Genome`` 
``cd Pact_Genome`` 
``wget http://ihpe.univ-perp.fr/telechargement/Data_to_downoload.rar ``

https://www.ncbi.nlm.nih.gov/assembly/GCF_003627195.1
``mkdir Lambda_Genome`` 
``cd Lambda_Genome`` 

``scp -P 2292  /Users/hputnam/Desktop/20190622/20190503/Pacuta_genome/Pocillopora_acuta_genome_v1.fasta hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/GENOME/Lambda_Genome
``

# Pdamicornis genome

http://pdam.reefgenomics.org/download/

``mkdir Pdam_Genome`` 
``cd Pdam_Genome`` 

```wget http://pdam.reefgenomics.org/download/pdam_scaffolds.fasta.gz```

#### Mcapitata 
``bismark_genome_preparation Mcap_Genome`` 

#### Pacuta 
``bismark_genome_preparation Pact_Genome`` 

#### Lambda 
``bismark_genome_preparation Lambda_Genome`` 

## Checking correct transfer of files from sequencer

from Genewiz
63a1cc9e23dcedd71a66dac00f36467d  Meth10_R1_001.fastq.gz
ac966619cb5a2483870f973b696a38a1  Meth10_R2_001.fastq.gz
06278c40d97f7a3efdae47218921ac6a  Meth11_R1_001.fastq.gz
07e71cb470c1622ced336478df70204a  Meth11_R2_001.fastq.gz
692992366c6709c1c50419bde55397cb  Meth12_R1_001.fastq.gz
27803c9f39f0456ebc56aaa07620b7c5  Meth12_R2_001.fastq.gz
9988a9eefd526231f6608a4549d2cfa7  Meth13_R1_001.fastq.gz
04db1696a84174d2a146fdd2e2b8a1ca  Meth13_R2_001.fastq.gz
131cb0c436de334bf17c375d49e02977  Meth14_R1_001.fastq.gz
fc7f543911fffce636899e05ec29dee3  Meth14_R2_001.fastq.gz
448f06b9bf1e19ee535ab74ce0f744fa  Meth15_R1_001.fastq.gz
6bd07498df6b72f060b8e2c7ce624e01  Meth15_R2_001.fastq.gz
46a49f6745c7478346f25ed7f509a704  Meth16_R1_001.fastq.gz
a6a1dd2355ce0e6e625c24f23373dc86  Meth16_R2_001.fastq.gz
f0a894a39a8a46d8729cd953c2531035  Meth17_R1_001.fastq.gz
bd0c3257b16d91e181b688d64d7d0b54  Meth17_R2_001.fastq.gz
d1b6b17e76c2634af2d9918f9ef7d9e9  Meth18_R1_001.fastq.gz
e9ad11546163f7e22ed7aa427fa7d0c1  Meth18_R2_001.fastq.gz
2560348651ae8a94fb351d0e54eff54e  Meth1_R1_001.fastq.gz
45faf1e9cbf5a664e5c6818817e85ef1  Meth1_R2_001.fastq.gz
a4963c92a3b6f62f2480f4f4bd72ec8f  Meth2_R1_001.fastq.gz
beac9fbd4b2e460bc7f8375e3e6caaa8  Meth2_R2_001.fastq.gz
1ae956f0fb4023165aaa4c58df9d437c  Meth3_R1_001.fastq.gz
6c6a24173eb24f2db81290b21edf17c0  Meth3_R2_001.fastq.gz
2ad73235560c58d8e68c3f75ada90840  Meth4_R1_001.fastq.gz
ad106065fbdbcf9030cd804296290e18  Meth4_R2_001.fastq.gz
8b966a3413ad40a087dcbf68204e84d6  Meth5_R1_001.fastq.gz
779e2fee610f06fe006c2ef8a4535114  Meth5_R2_001.fastq.gz
e481ea3fd24c09dafdfad82dec5d2d0a  Meth6_R1_001.fastq.gz
3d876fd2cd176450c8a63bf739ac5857  Meth6_R2_001.fastq.gz
d5199e48ed87773f67c4e3d8d9db3fce  Meth7_R1_001.fastq.gz
5d8467eb30d47bb00d8133fa4a743acf  Meth7_R2_001.fastq.gz
9b59ce33ab17158f2fea64f2d0d0ade7  Meth8_R1_001.fastq.gz
37c495ced527f3e30ae2f3388a9146bf  Meth8_R2_001.fastq.gz
c8ec0b73721a067729124c043904dc0e  Meth9_R1_001.fastq.gz
b896f80209fbaaee4ca5316f435839ff  Meth9_R2_001.fastq.gz


URI download
63a1cc9e23dcedd71a66dac00f36467d  Meth10_R1_001.fastq.gz
ac966619cb5a2483870f973b696a38a1  Meth10_R2_001.fastq.gz
06278c40d97f7a3efdae47218921ac6a  Meth11_R1_001.fastq.gz
07e71cb470c1622ced336478df70204a  Meth11_R2_001.fastq.gz
692992366c6709c1c50419bde55397cb  Meth12_R1_001.fastq.gz
27803c9f39f0456ebc56aaa07620b7c5  Meth12_R2_001.fastq.gz
9988a9eefd526231f6608a4549d2cfa7  Meth13_R1_001.fastq.gz
04db1696a84174d2a146fdd2e2b8a1ca  Meth13_R2_001.fastq.gz
131cb0c436de334bf17c375d49e02977  Meth14_R1_001.fastq.gz
fc7f543911fffce636899e05ec29dee3  Meth14_R2_001.fastq.gz
448f06b9bf1e19ee535ab74ce0f744fa  Meth15_R1_001.fastq.gz
6bd07498df6b72f060b8e2c7ce624e01  Meth15_R2_001.fastq.gz
46a49f6745c7478346f25ed7f509a704  Meth16_R1_001.fastq.gz
a6a1dd2355ce0e6e625c24f23373dc86  Meth16_R2_001.fastq.gz
f0a894a39a8a46d8729cd953c2531035  Meth17_R1_001.fastq.gz
bd0c3257b16d91e181b688d64d7d0b54  Meth17_R2_001.fastq.gz
d1b6b17e76c2634af2d9918f9ef7d9e9  Meth18_R1_001.fastq.gz
e9ad11546163f7e22ed7aa427fa7d0c1  Meth18_R2_001.fastq.gz
2560348651ae8a94fb351d0e54eff54e  Meth1_R1_001.fastq.gz
45faf1e9cbf5a664e5c6818817e85ef1  Meth1_R2_001.fastq.gz
a4963c92a3b6f62f2480f4f4bd72ec8f  Meth2_R1_001.fastq.gz
beac9fbd4b2e460bc7f8375e3e6caaa8  Meth2_R2_001.fastq.gz
1ae956f0fb4023165aaa4c58df9d437c  Meth3_R1_001.fastq.gz
6c6a24173eb24f2db81290b21edf17c0  Meth3_R2_001.fastq.gz
2ad73235560c58d8e68c3f75ada90840  Meth4_R1_001.fastq.gz
ad106065fbdbcf9030cd804296290e18  Meth4_R2_001.fastq.gz
8b966a3413ad40a087dcbf68204e84d6  Meth5_R1_001.fastq.gz
779e2fee610f06fe006c2ef8a4535114  Meth5_R2_001.fastq.gz
e481ea3fd24c09dafdfad82dec5d2d0a  Meth6_R1_001.fastq.gz
3d876fd2cd176450c8a63bf739ac5857  Meth6_R2_001.fastq.gz
d5199e48ed87773f67c4e3d8d9db3fce  Meth7_R1_001.fastq.gz
5d8467eb30d47bb00d8133fa4a743acf  Meth7_R2_001.fastq.gz
9b59ce33ab17158f2fea64f2d0d0ade7  Meth8_R1_001.fastq.gz
37c495ced527f3e30ae2f3388a9146bf  Meth8_R2_001.fastq.gz
c8ec0b73721a067729124c043904dc0e  Meth9_R1_001.fastq.gz
b896f80209fbaaee4ca5316f435839ff  Meth9_R2_001.fastq.gz


# Checking Sequence Quality
nohup fastqc /home/hputnam/Meth_Compare/RAW/*fastq.gz

scp -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/QC/raw_qc/multiqc_report.html /Users/hputnam/MyProjects/Meth_Compare

## Trimming and cleaning reads with TrimGalore
https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

mkdir RRBS_cleaned

### Trim RRBS
```
nohup sh -c 'for file in "Meth4" "Meth5" "Meth6"  "Meth13" "Meth14" "Meth15"
do
/home/hputnam/prog/TrimGalore-0.6.4/trim_galore \
--rrbs \
--non_directional \
--length 50 \
--quality 20 \
--phred33 \
--output_dir /home/hputnam/Meth_Compare/RRBS_cleaned \
--paired \
/home/hputnam/Meth_Compare/RAW/${file}_R1_001.fastq.gz \
/home/hputnam/Meth_Compare/RAW/${file}_R2_001.fastq.gz
done'
```

### Trim WGBS and MBDBS
```
nohup sh -c 'for file in "Meth1" "Meth2" "Meth3"  "Meth7" "Meth8" "Meth9" "Meth10" "Meth11" "Meth12" "Meth16" "Meth17" "Meth18"
do
/home/hputnam/prog/TrimGalore-0.6.4/trim_galore \
--length 20 \
--quality 20 \
--phred33 \
--output_dir /home/hputnam/Meth_Compare/cleaned_reads \
--paired \
/home/hputnam/Meth_Compare/RAW/${file}_R1_001.fastq.gz \
/home/hputnam/Meth_Compare/RAW/${file}_R2_001.fastq.gz
done'
```

### Count Clean reads

```zgrep -c "@A00387" *.fastq.gz```

### Quality control of cleaned sequencing reads (FASTQC)
https://github.com/s-andrews/FastQC

```
nohup fastqc Meth*fq.gz
```

### Using Multiqc to summarize cleaned FastQC results
https://github.com/ewels/MultiQC

```
multiqc .

scp -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/cleaned_reads/multiqc_report.html /Users/hputnam/MyProjects/Meth_Compare
```

Rationale for internal conversion efficiency
"Because only 5mC in a CpG site should still remain as cytosine following bisulfite treatment, the percentage of remaining non-CpG cytosines typically is inversely proportional to the conversion efficiency. Therefore, examining any remaining C at these non-CpG sites is known as a bisulfite treatment control."

# Mapping to determine conversion efficiency

#### MBD_BS LAMBDA
* Set for non-directional
* minimum alignment score function L,0,-0.6
``nohup bismark --genome /home/hputnam/Meth_Compare/GENOME/Lambda_Genome  --bowtie2 / --non_directional --score_min L,0,-0.6 -1 /home/hputnam/Meth_Compare/cleaned_reads/X.fq.gz -2 /home/hputnam/Meth_Compare/cleaned_reads/X.fq.gz`` 


# Deduplication

### MBD_BS Deduplicating
 ``deduplicate_bismark /home/hputnam/Meth_Compare/Mapped/X.bam
 ``
 
 ### WGBS Deduplication

# Extract Methylation
``bismark_methylation_extractor --gzip -p --ignore_r2 2 --bedGraph --zero_based --no_overlap --multicore 20 --buffer_size 20G --cytosine_report --report --genome_folder /home/hputnam/Meth_Compare/GENOME/Lambda_Genome  /home/hputnam/Meth_Compare/DeDup/``



