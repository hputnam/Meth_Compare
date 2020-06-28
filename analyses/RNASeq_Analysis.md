# Required Programs

fastqc 
multiqc  
fastp  
samtools  
hisat2 
stringtie 

# RNASeq Files

#### File location RiboDeplete ```/RAID_STORAGE2/hputnam/20191010_HoloInt_Compare```
Sample1_R1.fastq.gz:16487768
Sample1_R2.fastq.gz:16487768
Sample2_R1.fastq.gz:10755630
Sample2_R2.fastq.gz:10755630
Sample3_R1.fastq.gz:13534820
Sample3_R2.fastq.gz:13534820
Sample4_R1.fastq.gz:10275396
Sample4_R2.fastq.gz:10275396
Sample5_R1.fastq.gz:16601774
Sample5_R2.fastq.gz:16601774
Sample6_R1.fastq.gz:13440121
Sample6_R2.fastq.gz:13440121

#### File location PolyA 
```/RAID_STORAGE2/hputnam/20191104_HoloInt/30-274849628```
1041_R1_001.fastq.gz
1041_R2_001.fastq.gz
1471_R1_001.fastq.gz
1471_R2_001.fastq.gz
1637_R1_001.fastq.gz
1637_R2_001.fastq.gz

```/RAID_STORAGE2/hputnam/20200404_HoloInt_Batch3```
1101_R1_001.fastq.gz
1101_R2_001.fastq.gz
1548_R1_001.fastq.gz
1548_R2_001.fastq.gz
1628_R1_001.fastq.gz
1628_R2_001.fastq.gz


#### Sample Plug IDs
PA_1041
PA_1471
PA_1637
MC_1101
MC_1548
MC_1628

# Main Folder

`mkdir RiboDep_RNASeq`
`mkdir PolyA_RNASeq`

`mkdir ref`
`cd ref`
`wget http://cyanophora.rutgers.edu/montipora/Mcap.genome_assembly.fa.gz`
`wget http://cyanophora.rutgers.edu/montipora/Mcap.GFFannotation.gff`

`scp -P 2292  /Users/hputnam/Desktop/20190622/20190503/Pacuta_genome/Pocillopora_acuta_genome_v1.fasta hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/ref`

`scp -P 2292  /Users/hputnam/Desktop/20190622/20190503/Pacuta_genome/Structural_annotation_abintio.gff hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/ref`

`mv Structural_annotation_abintio.gff Pact_Structural_annotation_abintio.gff`

`sed 's/cov/_cov/' Pact_Structural_annotation_abintio.gff > Pact.GFFannotation.gff`

# Link Raw from RAID2
#### File location RiboDeplete ```/RAID_STORAGE2/hputnam/20191010_HoloInt_Compare```

#### File location PolyA 
```/RAID_STORAGE2/hputnam/20191104_HoloInt/30-274849628```

# QC raw files

`mkdir fastqc_raw`

`cd fastqc_raw`

`fastqc /RAID_STORAGE2/hputnam/20191010_HoloInt_Compare/*.fastq.gz  -o /home/hputnam/Meth_Compare/RiboDep_RNASeq/fastqc_raw`

`fastqc /home/hputnam/Meth_Compare/PolyA_RNASeq/RAW/1101_R*_001.fastq.gz  -o /home/hputnam/Meth_Compare/PolyA_RNASeq/fastqc_raw`
`fastqc /home/hputnam/Meth_Compare/PolyA_RNASeq/RAW/1548_R*_001.fastq.gz  -o /home/hputnam/Meth_Compare/PolyA_RNASeq/fastqc_raw`
`fastqc /home/hputnam/Meth_Compare/PolyA_RNASeq/RAW/1628_R*_001.fastq.gz  -o /home/hputnam/Meth_Compare/PolyA_RNASeq/fastqc_raw`


`multiqc .`

# Trimming

`mkdir cleaned_reads`


```
sh -c 'for file in "Sample1" "Sample2" "Sample3" "Sample4" "Sample5" "Sample6" 
do
fastp \
--in1 /RAID_STORAGE2/hputnam/20191010_HoloInt_Compare/${file}_R1.fastq.gz \
--in2 /RAID_STORAGE2/hputnam/20191010_HoloInt_Compare/${file}_R2.fastq.gz \
--out1 ${file}_R1_clean.fastq.gz \
--out2 ${file}_R2_clean.fastq.gz \
--failed_out cleaned_reads/${file}_failed.txt \
--qualified_quality_phred 20 \
--unqualified_percent_limit 10 \
--length_required 50 detect_adapter_for_pe \
--cut_right cut_right_window_size 5 cut_right_mean_quality 20
done'
```
"1041" "1101" "1471" "1548" "1628" "1637" 

```
sh -c 'for file in "1637" 
do
fastp \
--in1 /home/hputnam/Meth_Compare/PolyA_RNASeq/RAW/${file}_R1_001.fastq.gz \
--in2 /home/hputnam/Meth_Compare/PolyA_RNASeq/RAW/${file}_R2_001.fastq.gz \
--out1 ${file}_R1_clean.fastq.gz \
--out2 ${file}_R2_clean.fastq.gz \
--failed_out cleaned_reads/${file}_failed.txt \
--qualified_quality_phred 20 \
--unqualified_percent_limit 10 \
--length_required 50 detect_adapter_for_pe \
--cut_right cut_right_window_size 5 cut_right_mean_quality 20
done'
```


# QC trimmed files
mkdir fastqc_cleaned

`fastqc /home/hputnam/Meth_Compare/RiboDep_RNASeq/fastqc_cleaned/*.fastq.gz  -o /home/hputnam/Meth_Compare/RiboDep_RNASeq/fastqc_cleaned`

`fastqc /home/hputnam/Meth_Compare/PolyA_RNASeq/cleaned_reads/*.fastq.gz  -o /home/hputnam/Meth_Compare/PolyA_RNASeq/fastqc_cleaned`

`multiqc fastqc_cleaned`


# Download and view files

scp -r -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/fastqc_raw/multiqc_report.html /Users/hputnam/MyProjects/Meth_Compare/RNASeq

scp -r -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/cleaned_reads/multiqc_report.html /Users/hputnam/MyProjects/Meth_Compare/RNASeq

scp -r -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/PolyA_RNASeq/fastqc_raw/multiqc_report.html /Users/hputnam/MyProjects/Meth_Compare/RNASeq/PolyA

scp -r -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/PolyA_RNASeq/fastqc_cleaned/multiqc_report.html /Users/hputnam/MyProjects/Meth_Compare/RNASeq/PolyA


# Align Reads to Reference Genome
*HISAT2 is a fast and sensitive alignment program for mapping next-generation DNA and RNA sequencing reads to a reference genome.*

## Index the reference genome

Index the reference genome in the reference directory.

++HISAT2-build Alignment Arguments Used++:  
- <reference_in> - name of reference files  
- <gt2_base> -  basename of index files to write  
- -f -  reference file is a FASTA file

```
cd /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref
hisat2-build Mcap.genome_assembly.fa Mcap_ref
hisat2-build Pocillopora_acuta_genome_v1.fasta Pact_ref
```

# Alignment of clean reads to the reference genome

### Aligning paired end reads

rna strandedness - XX for Zymo RiboDep
Libraries are minus-stranded: The Read 1 sequence will be antisense to the RNA transcript from which it originates.

mkdir mapped
cd mapped

## Mcap Mapping for RiboDep
`
sh -c 'for i in "Sample4" "Sample5" "Sample6"
do
hisat2 -p 8 --dta -q -x /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap_ref \
-1 ${i}_R2_clean.fastq.gz \
-2 $(echo ${i}_R2_clean.fastq.gz|sed s/_R1/_R2/) -S ${i}.sam 
done'
`

## Pact Mapping for RiboDep
`
sh -c 'for i in "Sample1" "Sample2" "Sample3"
do
hisat2 -p 8 --dta -q -x /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact_ref \
-1 ${i}_R2_clean.fastq.gz \
-2 $(echo ${i}_R2_clean.fastq.gz|sed s/_R1/_R2/) -S ${i}.sam 
done'
`

rna strandedness - RF for Illumina
https://github.com/griffithlab/rnaseq_tutorial/blob/master/manuscript/supplementary_tables/supplementary_table_5.md


## Mcap Mapping for PolyA
`
sh -c 'for i in "1101" "1548" "1628"
do
hisat2 -p 8 --dta -q -x /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap_ref \
--rna-strandness RF \
-1 /home/hputnam/Meth_Compare/PolyA_RNASeq/cleaned_reads/${i}_R1_clean.fastq.gz \
-2 /home/hputnam/Meth_Compare/PolyA_RNASeq/cleaned_reads/${i}_R2_clean.fastq.gz -S ${i}.sam 
done'
`

## Pact Mapping for PolyA
`
sh -c 'for i in "1041" "1471" "1637"
do
hisat2 -p 8 --dta -q -x /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact_ref \
--rna-strandness RF \
-1 /home/hputnam/Meth_Compare/PolyA_RNASeq/cleaned_reads/${i}_R1_clean.fastq.gz \
-2 /home/hputnam/Meth_Compare/PolyA_RNASeq/cleaned_reads/${i}_R2_clean.fastq.gz -S ${i}.sam 
done'
`

# Convert and sort Sam to BAM
## RiboDep
`
sh -c 'for i in "Sample1" "Sample2" "Sample3" "Sample4" "Sample5" "Sample6"
do
samtools sort -@ 8 -o ${i}.bam ${i}.sam
done'
`
## PolyA
`
sh -c 'for i in "1041" "1471" "1637" "1101" "1548" "1628"
do
samtools sort -@ 8 -o ${i}.bam ${i}.sam
done'
`


# Remove Sam
`
rm *.sam

`

# Assemble aligned reads and quantify transcripts 

---

*StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts.*

- Reference-guided assembly with novel transcript discovery
- Merge output GTF files and assess the assembly performance
- Merged-GTF guided assembly without novel transcript discovery
- Compilation of GTF-files into gene and transcript count matrices

For our assembly, we will have to run StringTie. This run is necessary in order to compile our GTF files into gene count matrices that we will need for our differential expression analysis, because the StringTie script that compiles the GTF files ```prepDE.py``` only runs if the ```-e``` option is "on" during the previous assembly.

# Reference-guided assembly to genes only

++StringTie Arguments Used++:  
- -p - Specify number of processers
- -G - Specify annotation file
- -o - Name of output file
- -e - only estimate the abundance of given reference transcripts (requires -G)
- -A gene abundance estimation output file

mkdir stie

# Ignore Novel

#Mcap RiboDep Ignore Novel
`
stringtie Sample4.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o Sample4.gtf 
stringtie Sample5.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o Sample5.gtf 
stringtie Sample6.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o Sample6.gtf 
`

#Pact RiboDep Ignore Novel
`
stringtie Sample1.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -o Sample1.gtf 
stringtie Sample2.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -o Sample2.gtf 
stringtie Sample3.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -o Sample3.gtf 
`

#Mcap PolyA Ignore Novel
`
stringtie 1101.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.fixed.gff -o 1101.gtf 
stringtie 1548.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o 1548.gtf 
stringtie 1628.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o 1628.gtf 
`
#Pact PolyA Ignore Novel
`
stringtie 1041.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -o 1041.gtf 
stringtie 1471.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -o 1471.gtf 
stringtie 1637.bam -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -o 1637.gtf 
`


#### Assess the performance of the assembly

*Gffcompare is a tool that can compare, merge, annotate and estimate accuracy of GFF/GTF files when compared with a reference annotation*

Using the StringTie merge mode, merge the assembly-generated GTF files to assess how well the predicted transcripts track to the reference annotation file. This step requires the TXT file,  (mergelist.txt). This file lists all of the file names to be merged. *Make sure ```mergelist.txt``` is in the StringTie program directory*.

++StringTie Arguments Used++:  
- --merge - Distinct from the assembly usage mode used above, in the merge mode, StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts.
- -p - Specify number of processers
- -G - Specify reference annotation file. With this option, StringTie assembles the transfrags from the input GTF files with the reference sequences
- -o - Name of output file
- <mergelist.txt> - File listing all filenames to be merged. Include full path.

#Mcap RiboDep

`
nano Mcap_mergelist.txt

Sample4.gtf
Sample5.gtf
Sample6.gtf
`

`
stringtie --merge -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o Mcap_stringtie_merged.gtf Mcap_mergelist.txt
`

#Pact RiboDep
`
nano Pact_mergelist.txt

Sample1.gtf
Sample2.gtf
Sample3.gtf
`

`
stringtie --merge -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -o Pact_stringtie_merged.gtf Pact_mergelist.txt
`

#Mcap PolyA

`
nano Mcap_mergelist.txt

1101.gtf  
1548.gtf
1628.gtf 
`

`
stringtie --merge -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -o Mcap_stringtie_merged.gtf Mcap_mergelist.txt
`

#Pact PolyA
`
nano Pact_mergelist.txt

1041.gtf
1471.gtf
1637.gtf
`

`
stringtie --merge -p 8 -e -G /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -o Pact_stringtie_merged.gtf Pact_mergelist.txt
`


Now we can use the program gffcompare to compare the merged GTF to our reference genome.

++Gffcompare Arguments Used++:  
- -r - Specify reference annotation file
- -G - Compare all the transcripts in our input file ```stringtie_merged.gtf```
- -o - Prefix of all output files

#Mcap RiboDep
`
gffcompare -r /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -G -o Mcap_merged Mcap_stringtie_merged.gtf
`

#Pact RiboDep
`
gffcompare -r /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -G -o Pact_merged Pact_stringtie_merged.gtf
`

#Mcap PolyA
`
gffcompare -r /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Mcap.GFFannotation.gff -G -o Mcap_merged Mcap_stringtie_merged.gtf
`

#Pact PolyA
`
gffcompare -r /home/hputnam/Meth_Compare/RiboDep_RNASeq/ref/Pact.GFFannotation.gff -G -o Pact_merged Pact_stringtie_merged.gtf
`

Some of the output files you will see are... 
- merged.stats
- merged.tracking
- merged.annotated.gtf
- merged.stringtie_merged.gtf.refmap
- merged.loci
- merged.stringtie_merged.gtf.tmap

We are most interested in the files ```merged.annotated.gtf``` and ```merged.stats```. The file ```merged.annotated.gtf``` tells you how well the predicted transcripts track to the reference annotation file and the file ```merged.stats``` file shows the sensitivity and precision statistics and total number for different features (genes, exons, transcripts).  



#### Compilation of GTF-files into gene and transcript count matrices

The StringTie program includes a script, ```prepDE.py``` that compiles your assembly files into gene and transcript count matrices. This script requires as input the list of sample names and their full file paths, (sample_list.txt).

. Run ```prepDE.py``` to merge assembled files together into a DESeq2-friendly version.

++StringTie prepDE.py Arguments Used++:  
- -i - Specify that input is a TXT file
- -g - Require output gene count file, default name is ```gene_count_matrix.csv```
- -t - Require output transcript count gene count file, default name is ```transcript_count_matrix.csv```

# Mcap RiboDep
nano Mcap_sample_list.txt

1101	Sample4.gtf
1548	Sample5.gtf
1628	Sample6.gtf

`
prepDE.py -i Mcap_sample_list.txt -g Mcap_gene_count_matrix.csv
`

# Pact RiboDep

nano Pact_sample_list.txt

1041	Sample1.gtf
1471	Sample2.gtf
1637	Sample3.gtf

`
prepDE.py -i Pact_sample_list.txt -g Pact_gene_count_matrix.csv
`

# Mcap PolyA
nano Mcap_sample_list.txt

1101	1101.gtf
1548	1548.gtf
1628	1628.gtf

`
prepDE.py -i Mcap_sample_list.txt -g Mcap_gene_count_matrix.csv
`

# Pact PolyA

nano Pact_sample_list.txt

1041	1041.gtf
1471	1471.gtf
1637	1637.gtf

`
prepDE.py -i Pact_sample_list.txt -g Pact_gene_count_matrix.csv
`

## remove bam file before downloading
rm *.bam

#download mapping and counts matrices

scp -r -P 2292  hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/PolyA_RNASeq/mapped/ /Users/hputnam/MyProjects/Meth_Compare/RNASeq/PolyA

scp -r -P 2292  hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/mapped/ /Users/hputnam/MyProjects/Meth_Compare/RNASeq/RiboDep


