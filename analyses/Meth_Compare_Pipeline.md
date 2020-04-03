# Library Prep Methods
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
`mkdir GENOME`

`cd GENOME`

`mkdir Mcap_Genome`

`cd Mcap_Genome`

`wget http://cyanophora.rutgers.edu/montipora/Mcap.genome_assembly.fa.gz`

`mkdir Pact_Genome`

`cd Pact_Genome`

`wget http://ihpe.univ-perp.fr/telechargement/Data_to_downoload.rar `

### Lambda genome
https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161521

`mkdir Lambda_Genome`

`cd Lambda_Genome`

`scp -P 2292  /Users/hputnam/Desktop/20190622/20190503/Pacuta_genome/Pocillopora_acuta_genome_v1.fasta hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/GENOME/Lambda_Genome
`

# Pdamicornis genome

http://pdam.reefgenomics.org/download/

`mkdir Pdam_Genome`

`cd Pdam_Genome`

`wget http://pdam.reefgenomics.org/download/pdam_scaffolds.fasta.gz`

#### Mcapitata
``bismark_genome_preparation Mcap_Genome``

#### Pacuta
``bismark_genome_preparation Pact_Genome``

#### Lambda
``bismark_genome_preparation Lambda_Genome``

## Checking correct transfer of files from sequencer

from Genewiz

```
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
```

URI download

```
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
```

Copy to Gannet server
```
63a1cc9e23dcedd71a66dac00f36467d  ./Meth10_R1_001.fastq.gz
ac966619cb5a2483870f973b696a38a1  ./Meth10_R2_001.fastq.gz
06278c40d97f7a3efdae47218921ac6a  ./Meth11_R1_001.fastq.gz
07e71cb470c1622ced336478df70204a  ./Meth11_R2_001.fastq.gz
692992366c6709c1c50419bde55397cb  ./Meth12_R1_001.fastq.gz
27803c9f39f0456ebc56aaa07620b7c5  ./Meth12_R2_001.fastq.gz
9988a9eefd526231f6608a4549d2cfa7  ./Meth13_R1_001.fastq.gz
04db1696a84174d2a146fdd2e2b8a1ca  ./Meth13_R2_001.fastq.gz
131cb0c436de334bf17c375d49e02977  ./Meth14_R1_001.fastq.gz
fc7f543911fffce636899e05ec29dee3  ./Meth14_R2_001.fastq.gz
448f06b9bf1e19ee535ab74ce0f744fa  ./Meth15_R1_001.fastq.gz
6bd07498df6b72f060b8e2c7ce624e01  ./Meth15_R2_001.fastq.gz
46a49f6745c7478346f25ed7f509a704  ./Meth16_R1_001.fastq.gz
a6a1dd2355ce0e6e625c24f23373dc86  ./Meth16_R2_001.fastq.gz
f0a894a39a8a46d8729cd953c2531035  ./Meth17_R1_001.fastq.gz
bd0c3257b16d91e181b688d64d7d0b54  ./Meth17_R2_001.fastq.gz
d1b6b17e76c2634af2d9918f9ef7d9e9  ./Meth18_R1_001.fastq.gz
e9ad11546163f7e22ed7aa427fa7d0c1  ./Meth18_R2_001.fastq.gz
2560348651ae8a94fb351d0e54eff54e  ./Meth1_R1_001.fastq.gz
45faf1e9cbf5a664e5c6818817e85ef1  ./Meth1_R2_001.fastq.gz
a4963c92a3b6f62f2480f4f4bd72ec8f  ./Meth2_R1_001.fastq.gz
beac9fbd4b2e460bc7f8375e3e6caaa8  ./Meth2_R2_001.fastq.gz
1ae956f0fb4023165aaa4c58df9d437c  ./Meth3_R1_001.fastq.gz
6c6a24173eb24f2db81290b21edf17c0  ./Meth3_R2_001.fastq.gz
2ad73235560c58d8e68c3f75ada90840  ./Meth4_R1_001.fastq.gz
ad106065fbdbcf9030cd804296290e18  ./Meth4_R2_001.fastq.gz
8b966a3413ad40a087dcbf68204e84d6  ./Meth5_R1_001.fastq.gz
779e2fee610f06fe006c2ef8a4535114  ./Meth5_R2_001.fastq.gz
e481ea3fd24c09dafdfad82dec5d2d0a  ./Meth6_R1_001.fastq.gz
3d876fd2cd176450c8a63bf739ac5857  ./Meth6_R2_001.fastq.gz
d5199e48ed87773f67c4e3d8d9db3fce  ./Meth7_R1_001.fastq.gz
5d8467eb30d47bb00d8133fa4a743acf  ./Meth7_R2_001.fastq.gz
9b59ce33ab17158f2fea64f2d0d0ade7  ./Meth8_R1_001.fastq.gz
37c495ced527f3e30ae2f3388a9146bf  ./Meth8_R2_001.fastq.gz
c8ec0b73721a067729124c043904dc0e  ./Meth9_R1_001.fastq.gz
b896f80209fbaaee4ca5316f435839ff  ./Meth9_R2_001.fastq.gz
```

Copy to Mox
```
63a1cc9e23dcedd71a66dac00f36467d  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth10_R1_001.fastq.gz
ac966619cb5a2483870f973b696a38a1  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth10_R2_001.fastq.gz
06278c40d97f7a3efdae47218921ac6a  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth11_R1_001.fastq.gz
07e71cb470c1622ced336478df70204a  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth11_R2_001.fastq.gz
692992366c6709c1c50419bde55397cb  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth12_R1_001.fastq.gz
27803c9f39f0456ebc56aaa07620b7c5  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth12_R2_001.fastq.gz
9988a9eefd526231f6608a4549d2cfa7  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth13_R1_001.fastq.gz
04db1696a84174d2a146fdd2e2b8a1ca  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth13_R2_001.fastq.gz
131cb0c436de334bf17c375d49e02977  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth14_R1_001.fastq.gz
fc7f543911fffce636899e05ec29dee3  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth14_R2_001.fastq.gz
448f06b9bf1e19ee535ab74ce0f744fa  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth15_R1_001.fastq.gz
6bd07498df6b72f060b8e2c7ce624e01  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth15_R2_001.fastq.gz
46a49f6745c7478346f25ed7f509a704  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth16_R1_001.fastq.gz
a6a1dd2355ce0e6e625c24f23373dc86  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth16_R2_001.fastq.gz
f0a894a39a8a46d8729cd953c2531035  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth17_R1_001.fastq.gz
bd0c3257b16d91e181b688d64d7d0b54  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth17_R2_001.fastq.gz
d1b6b17e76c2634af2d9918f9ef7d9e9  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth18_R1_001.fastq.gz
e9ad11546163f7e22ed7aa427fa7d0c1  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth18_R2_001.fastq.gz
2560348651ae8a94fb351d0e54eff54e  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth1_R1_001.fastq.gz
45faf1e9cbf5a664e5c6818817e85ef1  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth1_R2_001.fastq.gz
a4963c92a3b6f62f2480f4f4bd72ec8f  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth2_R1_001.fastq.gz
beac9fbd4b2e460bc7f8375e3e6caaa8  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth2_R2_001.fastq.gz
1ae956f0fb4023165aaa4c58df9d437c  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth3_R1_001.fastq.gz
6c6a24173eb24f2db81290b21edf17c0  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth3_R2_001.fastq.gz
2ad73235560c58d8e68c3f75ada90840  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth4_R1_001.fastq.gz
ad106065fbdbcf9030cd804296290e18  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth4_R2_001.fastq.gz
8b966a3413ad40a087dcbf68204e84d6  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth5_R1_001.fastq.gz
779e2fee610f06fe006c2ef8a4535114  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth5_R2_001.fastq.gz
e481ea3fd24c09dafdfad82dec5d2d0a  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth6_R1_001.fastq.gz
3d876fd2cd176450c8a63bf739ac5857  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth6_R2_001.fastq.gz
d5199e48ed87773f67c4e3d8d9db3fce  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth7_R1_001.fastq.gz
5d8467eb30d47bb00d8133fa4a743acf  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth7_R2_001.fastq.gz
9b59ce33ab17158f2fea64f2d0d0ade7  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth8_R1_001.fastq.gz
37c495ced527f3e30ae2f3388a9146bf  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth8_R2_001.fastq.gz
c8ec0b73721a067729124c043904dc0e  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth9_R1_001.fastq.gz
b896f80209fbaaee4ca5316f435839ff  /gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth9_R2_001.fastq.gz

#ran the code below and it generated no output indicating files on Gannet and Mox are identical
diff <(cut -d' ' -f1 /gscratch/scrubbed/sr320/froger-raw/00_fastq/md5sum_list.txt) <(cut -d' ' -f1 /gscratch/scrubbed/strigg/analyses/20200311/md5sum.txt)
```


# Checking Sequence Quality and Trimming
fastp

```
[sr320@mox1 20200305_methcompare_fastp_trimming]$ cat *.sh
#!/bin/bash
## Job Name
#SBATCH --job-name=pgen_fastp_trimming_EPI
## Allocation Definition
#SBATCH --account=coenv
#SBATCH --partition=coenv
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=1-00:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samwhite@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/samwhite/outputs/20200305_methcompare_fastp_trimming


### WGBS and RRBS trimming using fastp.
### FastQ files were provide by Hollie Putnam.
### See this GitHub repo for more info:
### https://github.com/hputnam/Meth_Compare

# Exit script if any command fails
# set -e

# Load Python Mox module for Python module availability

module load intel-python3_2017

# Document programs in PATH (primarily for program version ID)

{
date
echo ""
echo "System PATH for $SLURM_JOB_ID"
echo ""
printf "%0.s-" {1..10}
echo "${PATH}" | tr : \\n
} >> system_path.log

# Set number of CPUs to use
threads=27

# Paths to programs
fastp=/gscratch/srlab/programs/fastp-0.20.0/fastp
multiqc=/gscratch/srlab/programs/anaconda3/bin/multiqc

# Programs array
programs_array=("${fastp}" "${multiqc}")


# Capture program options
for program in "${!programs_array[@]}"
do
	echo "Program options for ${programs_array[program]}: "
	echo ""
	${programs_array[program]} -h &>> program_options.log
	echo ""
	echo ""
	echo "----------------------------------------------"
	echo ""
	echo ""
done


# Input/output files
trimmed_checksums=trimmed_fastq_checksums.md5

# Inititalize arrays
# These were provided by Hollie Putnam
# See https://github.com/hputnam/Meth_Compare/blob/master/Meth_Compare_Pipeline.md
rrbs_array=(Meth4 Meth5 Meth6 Meth13 Meth14 Meth15)
wgbs_array=(Meth1 Meth2 Meth3  Meth7 Meth8 Meth9 Meth10 Meth11 Meth12 Meth16 Meth17 Meth18)

# Assign file suffixes to variables
read1="_R1_001.fastq.gz"
read2="_R2_001.fastq.gz"

# Create list of fastq files used in analysis
for fastq in *.gz
do
  echo "${fastq}" >> fastq.list.txt
done

# Run fastp on RRBS files
# Specifies removal of first 2bp from 3' end of read1 and
# removes 2bp from 5' end of read2, per Bismark instructions for RRBS
# https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
for index in "${!rrbs_array[@]}"
do
	timestamp=$(date +%Y%m%d%M%S)
	${fastp} \
	--in1 "${rrbs_array[index]}${read1}" \
	--in2 "${rrbs_array[index]}${read2}" \
	--detect_adapter_for_pe \
	--trim_tail1 2 \
	--trim_front2 2 \
	--thread ${threads} \
	--html "${rrbs_array[index]}.fastp-trim.${timestamp}.report.html" \
	--json "${rrbs_array[index]}.fastp-trim.${timestamp}.report.json" \
	--out1 "${rrbs_array[index]}.fastp-trim.${timestamp}${read1}" \
	--out2 "${rrbs_array[index]}.fastp-trim.${timestamp}${read2}"

	# Generate md5 checksums for newly trimmed files
	{
		md5sum "${rrbs_array[index]}.fastp-trim.${timestamp}${read1}"
		md5sum "${rrbs_array[index]}.fastp-trim.${timestamp}${read2}"
	} >> "${trimmed_checksums}"

done

# Run fastp on WGBS files
# Specifies removal of first 10bp from 5' and 3' end of all reads
# per Bismark instructions for WGBS Zymo/Swift library kits
# https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
for index in "${!wgbs_array[@]}"
do
	timestamp=$(date +%Y%m%d%M%S)
	${fastp} \
	--in1 "${wgbs_array[index]}${read1}" \
	--in2 "${wgbs_array[index]}${read2}" \
	--detect_adapter_for_pe \
	--trim_front1 10 \
	--trim_tail1 10 \
	--trim_front2 10 \
	--trim_tail2 10 \
	--thread ${threads} \
	--html "${wgbs_array[index]}.fastp-trim.${timestamp}.report.html" \
	--json "${wgbs_array[index]}.fastp-trim.${timestamp}.report.json" \
	--out1 "${wgbs_array[index]}.fastp-trim.${timestamp}${read1}" \
	--out2 "${wgbs_array[index]}.fastp-trim.${timestamp}${read2}"

	# Generate md5 checksums for newly trimmed files
	{
		md5sum "${wgbs_array[index]}.fastp-trim.${timestamp}${read1}"
		md5sum "${wgbs_array[index]}.fastp-trim.${timestamp}${read2}"
	} >> "${trimmed_checksums}"

done

# Run multiqc
${multiqc} .
```


Here is the contents of the output directory

```
[sr320@mox1 20200305_methcompare_fastp_trimming]$ ls
20200305_methcompare_fastp_trimming.sh	Meth3_R1_001.fastq.gz
fastq.list.txt				Meth3_R2_001.fastq.gz
Meth10_R1_001.fastq.gz			Meth4.fastp-trim.202003055221_R1_001.fastq.gz
Meth10_R2_001.fastq.gz			Meth4.fastp-trim.202003055221_R2_001.fastq.gz
Meth11_R1_001.fastq.gz			Meth4.fastp-trim.202003055221.report.html
Meth11_R2_001.fastq.gz			Meth4.fastp-trim.202003055221.report.json
Meth12_R1_001.fastq.gz			Meth4_R1_001.fastq.gz
Meth12_R2_001.fastq.gz			Meth4_R2_001.fastq.gz
Meth13_R1_001.fastq.gz			Meth5.fastp-trim.202003051832_R1_001.fastq.gz
Meth13_R2_001.fastq.gz			Meth5.fastp-trim.202003051832_R2_001.fastq.gz
Meth14_R1_001.fastq.gz			Meth5_R1_001.fastq.gz
Meth14_R2_001.fastq.gz			Meth5_R2_001.fastq.gz
Meth15_R1_001.fastq.gz			Meth6_R1_001.fastq.gz
Meth15_R2_001.fastq.gz			Meth6_R2_001.fastq.gz
Meth16_R1_001.fastq.gz			Meth7_R1_001.fastq.gz
Meth16_R2_001.fastq.gz			Meth7_R2_001.fastq.gz
Meth17_R1_001.fastq.gz			Meth8_R1_001.fastq.gz
Meth17_R2_001.fastq.gz			Meth8_R2_001.fastq.gz
Meth18_R1_001.fastq.gz			Meth9_R1_001.fastq.gz
Meth18_R2_001.fastq.gz			Meth9_R2_001.fastq.gz
Meth1_R1_001.fastq.gz			program_options.log
Meth1_R2_001.fastq.gz			slurm-2015259.out
Meth2_R1_001.fastq.gz			system_path.log
Meth2_R2_001.fastq.gz			trimmed_fastq_checksums.md5
```

### TrimGalore
```
#!/bin/bash
## Job Name
#SBATCH --job-name=TrimGfroger
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=4-23:30:00
## Memory per node
#SBATCH --mem=100G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/strigg/analyses/20200311

%%bash

#run md5 on fastqs to be sure the files in the directory (/gscratch/scrubbed/sr320/froger-raw/00_fastq/*.gz) match what was on Gannet (/gscratch/scrubbed/sr320/froger-raw/00_fastq/md5sum_list.txt).

md5sum /gscratch/scrubbed/sr320/froger-raw/00_fastq/*.gz > md5sum.txt

#TrimGalore on WGBS and MBD-BS
/gscratch/srlab/programs/TrimGalore-0.4.5/trim_galore \
--output_dir /gscratch/scrubbed/strigg/analyses/20200311/WGBS_MBD \
--paired \
--fastqc_args \
"--outdir /gscratch/scrubbed/strigg/analyses/20200311/WGBS_MBD/FASTQC \
--threads 28" \
--illumina \
--clip_R1 10 \
--clip_R2 10 \
--three_prime_clip_R1 10 \
--three_prime_clip_R2 10 \
--path_to_cutadapt /gscratch/srlab/programs/miniconda3/bin/cutadapt \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth10_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth10_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth11_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth11_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth12_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth12_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth1_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth1_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth2_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth2_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth3_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth3_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth7_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth7_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth8_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth8_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth9_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth9_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth16_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth16_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth17_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth17_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth18_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth18_R2_001.fastq.gz

#TrimGalore on RRBS
/gscratch/srlab/programs/TrimGalore-0.4.5/trim_galore \
--output_dir /gscratch/scrubbed/strigg/analyses/20200311/RRBS \
--fastqc_args \
"--outdir /gscratch/scrubbed/strigg/analyses/20200311/RRBS/FASTQC \
--threads 28" \
--non_directional \
--rrbs \
--paired \
--path_to_cutadapt /gscratch/srlab/programs/miniconda3/bin/cutadapt \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth13_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth13_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth14_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth14_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth15_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth15_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth4_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth4_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth5_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth5_R2_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth6_R1_001.fastq.gz \
/gscratch/scrubbed/sr320/froger-raw/00_fastq/Meth6_R2_001.fastq.gz

#run multiqc for WGBS and MBD samples

/gscratch/srlab/strigg/bin/anaconda3/bin/multiqc \
/gscratch/scrubbed/strigg/analyses/20200311/WGBS_MBD/FASTQC/.

#run multiqc for RRBS samples

/gscratch/srlab/strigg/bin/anaconda3/bin/multiqc \
/gscratch/scrubbed/strigg/analyses/20200311/RRBS/FASTQC/.
```
---

# Genome preparation

```
#!/bin/bash
## Job Name
#SBATCH --job-name=fr-01
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes (We only get 1, so this is fixed)
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=6-00:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr320@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/sr320/030520-fr01/

# Prepping 3 coral genomes

# Directories and programs
bismark_dir="/gscratch/srlab/programs/Bismark-0.21.0"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
samtools="/gscratch/srlab/programs/samtools-1.9/samtools"
reads_dir="/gscratch/srlab/sr320/data/olurida-bs/decomp/"
#genome_folder="/gscratch/srlab/sr320/data/olurida-genomes/v081/"

source /gscratch/srlab/programs/scripts/paths.sh


${bismark_dir}/bismark_genome_preparation \
--verbose \
--parallel 28 \
--path_to_aligner ${bowtie2_dir} \
/gscratch/srlab/sr320/data/froger/Mcap_Genome/


${bismark_dir}/bismark_genome_preparation \
--verbose \
--parallel 28 \
--path_to_aligner ${bowtie2_dir} \
/gscratch/srlab/sr320/data/froger/Pact_Genome/
```


# Mapping





# Methylation Extract
