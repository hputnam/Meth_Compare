# Required Programs

fastqc 
multiqc  
fastp  
samtools  
hisat2 
stringtie 



# Main Folder

`mkdir RiboDep_RNASeq`

`mkdir ref`
`cd ref`
`wget http://cyanophora.rutgers.edu/montipora/Mcap.genome_assembly.fa.gz`
`wget http://cyanophora.rutgers.edu/montipora/Mcap.GFFannotation.gff`

scp -P 2292  /Users/hputnam/Desktop/20190622/20190503/Pacuta_genome/Pocillopora_acuta_genome_v1.fasta hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/ref

scp -P 2292  /Users/hputnam/Desktop/20190622/20190503/Pacuta_genome/Structural_annotation_abintio.gff hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/ref

mv Structural_annotation_abintio.gff Pact_Structural_annotation_abintio.gff

sed 's/cov/_cov/' Pact_Structural_annotation_abintio.gff > Pact.GFFannotation.gff

# Link Raw from RAID2
#### File location ```/RAID_STORAGE2/hputnam/20191010_HoloInt_Compare```


# QC raw files

`mkdir fastqc_raw`

`cd fastqc_raw`

`fastqc /RAID_STORAGE2/hputnam/20191010_HoloInt_Compare/*.fastq.gz  -o /home/hputnam/Meth_Compare/RiboDep_RNASeq/fastqc_raw`

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


# QC trimmed files
mkdir fastqc_cleaned

`fastqc /home/hputnam/Meth_Compare/RiboDep_RNASeq/fastqc_cleaned/*.fastq.gz  -o /home/hputnam/Meth_Compare/RiboDep_RNASeq/fastqc_cleaned`

`multiqc fastqc_cleaned`


# Download and view files

scp -r -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/fastqc_raw/multiqc_report.html /Users/hputnam/MyProjects/Meth_Compare/RNASeq

scp -r -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/cleaned_reads/multiqc_report.html /Users/hputnam/MyProjects/Meth_Compare/RNASeq

*HISAT2 is a fast and sensitive alignment program for mapping next-generation DNA and RNA sequencing reads to a reference genome.*

The version of StringTie available on Bioconda is not the most recent version (v2.1.0). The version installed in conda (v2.0) has errors when running with the '-e' option that we need for this next step in StringTie. We will have to install StringTie outside of the conda environment. The following commands will install the latest version and test the binary. This only took about 3 min to run.


- Index the reference genome
- Alignment of clean reads to the reference genome

Create a subdirectory within data for HISAT2


`mkdir hisat2`  



#### Index the reference genome

Index the reference genome in the reference directory.

++HISAT2-build Alignment Arguments Used++:  
- <reference_in> - name of reference files  
- <gt2_base> -  basename of index files to write  
- -f -  reference file is a FASTA file

```
cd ref
hisat2-build Mcap.genome_assembly.fa Mcap_ref
hisat2-build Pocillopora_acuta_genome_v1.fasta Pact_ref
```

# Alignment of clean reads to the reference genome

### Aligning paired end reads
#### Has the R1 in array1 because the sed in the for loop changes it to an R2. SAM files are of both forward and reverse reads


## Mcap Mapping
`
sh -c 'for i in "Sample4" "Sample5" "Sample6"
do
hisat2 -p 8 --dta -q -x ../ref/Mcap_ref \
-1 ${i}_R2_clean.fastq.gz \
-2 $(echo ${i}_R2_clean.fastq.gz|sed s/_R1/_R2/) -S ${i}.sam 
done'
`

## Pact Mapping
`
sh -c 'for i in "Sample1" "Sample2" "Sample3"
do
hisat2 -p 8 --dta -q -x ../ref/Pact_ref \
-1 ${i}_R2_clean.fastq.gz \
-2 $(echo ${i}_R2_clean.fastq.gz|sed s/_R1/_R2/) -S ${i}.sam 
done'
`

# Convert and sort Sam to BAM

`
sh -c 'for i in "Sample1" "Sample2" "Sample3" "Sample4" "Sample5" "Sample6"
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

For our assembly, we will have to run StringTie twice. The first run will be a reference guided assembly that will allow for discovery of novel transcripts (by leaving out the -e option). Then we will merge the output GTF files and examine the sensitivity of our assembly. We will use the merged-GTF from our first assembly to guide our second StringTie run (including the ```-e``` option). This second run is necessary in order to compile our GTF files into gene and transcript count matrices that we will need for our differential expression analysis, because the StringTie script that compiles the GTF files ```prepDE.py``` only runs if the ```-e``` option is "on" during the previous assembly.

#### Reference-guided assembly with novel transcript discovery


Create the StringTie reference-guided assembly script, ```McapStringTie-assembly-to-ref.sh``` *inside of the StringTie program directory.*  

++StringTie Arguments Used++:  
- -p - Specify number of processers
- -G - Specify annotation file
- -o - Name of output file

#Mcap 
`
stringtie Sample4.bam -p 8 -G ../ref/Mcap.GFFannotation.gff -o Sample4.gtf 
stringtie Sample5.bam -p 8 -G ../ref/Mcap.GFFannotation.gff -o Sample5.gtf 
stringtie Sample6.bam -p 8 -G ../ref/Mcap.GFFannotation.gff -o Sample6.gtf 
`

#Pact 
`
stringtie Sample1.bam -p 8 -G ../ref/Pact.GFFannotation.gff -o Sample1.gtf 
stringtie Sample2.bam -p 8 -G ../ref/Pact.GFFannotation.gff -o Sample2.gtf 
stringtie Sample3.bam -p 8 -G ../ref/Pact.GFFannotation.gff -o Sample3.gtf 
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

#Mcap

`nano Mcap_mergelist.txt

Sample4.gtf
Sample5.gtf
Sample6.gtf
`

`
stringtie --merge -p 8 -G ../ref/Mcap.GFFannotation.gff -o Mcap_stringtie_merged.gtf Mcap_mergelist.txt
`

#Pact
`
nano Pact_mergelist.txt

Sample1.gtf
Sample2.gtf
Sample3.gtf
`

`
stringtie --merge -p 8 -G ../ref/Pact.GFFannotation.gff -o Pact_stringtie_merged.gtf Pact_mergelist.txt
`


Now we can use the program gffcompare to compare the merged GTF to our reference genome.

++Gffcompare Arguments Used++:  
- -r - Specify reference annotation file
- -G - Compare all the transcripts in our input file ```stringtie_merged.gtf```
- -o - Prefix of all output files

#Mcap
`
gffcompare -r ../ref/Mcap.GFFannotation.gff -G -o Mcap_merged Mcap_stringtie_merged.gtf
`

#Pact
`
gffcompare -r ../ref/Pact.GFFannotation.gff -G -o Pact_merged Pact_stringtie_merged.gtf
`

Some of the output files you will see are... 
- merged.stats
- merged.tracking
- merged.annotated.gtf
- merged.stringtie_merged.gtf.refmap
- merged.loci
- merged.stringtie_merged.gtf.tmap

Move all of the gffcompare output files to the output directory. We are most interested in the files ```merged.annotation.gtf``` and ```merged.stats```. The file ```merged.annotated.gtf``` tells you how well the predicted transcripts track to the reference annotation file and the file ```merged.stats``` file shows the sensitivity and precision statistics and total number for different features (genes, exons, transcripts).  


#### Merged-GTF guided assembly without novel transcript discovery

Now run StringTie again, but this time the assembly guided by the merged GTF file ```stringtie_merged.gtf```. Here, StringTie skips over novel sequences because we include the ```-e``` option. This is okay now because we identified novel transcripts in the previous StringTie run. The ```-e``` option is the important part of this second run because the ```prepDE.py``` script used in the next step to compile the output GTFs only works if the ```-e``` option is included here.

++StringTie Arguments Used++:  
- -p - Specify number of processers
- -e - Limit the estimation and output of transcripts to only those that match the reference (in this case, our merged GTF)
- -G - Specify annotation file
- -o - Name of output file


#Mcap 
`
stringtie Sample4.bam -p 8 -e -G Mcap_stringtie_merged.gtf -o Sample4_merged2.gtf 
stringtie Sample5.bam -p 8 -e -G Mcap_stringtie_merged.gtf -o Sample5_merged2.gtf 
stringtie Sample6.bam -p 8 -e -G Mcap_stringtie_merged.gtf -o Sample6_merged2.gtf  
`

#Pact 
`
stringtie Sample1.bam -p 8 -e -G Pact_stringtie_merged.gtf -o Sample1_merged2.gtf 
stringtie Sample2.bam -p 8 -e -G Pact_stringtie_merged.gtf -o Sample2_merged2.gtf 
stringtie Sample3.bam -p 8 -e -G Pact_stringtie_merged.gtf -o Sample3_merged2.gtf 
`


#### Compilation of GTF-files into gene and transcript count matrices

The StringTie program includes a script, ```prepDE.py``` that compiles your assembly files into gene and transcript count matrices. This script requires as input the list of sample names and their full file paths, (sample_list.txt).

. Run ```prepDE.py``` to merge assembled files together into a DESeq2-friendly version.

++StringTie prepDE.py Arguments Used++:  
- -i - Specify that input is a TXT file
- -g - Require output gene count file, default name is ```gene_count_matrix.csv```
- -t - Require output transcript count gene count file, default name is ```transcript_count_matrix.csv```

# Mcap
nano Mcap_sample_list.txt

1101	Sample4_merged2.gtf
1548	Sample5_merged2.gtf
1628	Sample6_merged2.gtf

`
prepDE.py -i Mcap_sample_list.txt -g Mcap_gene_count_matrix.csv -t Mcap_transcript_count_matrix.csv
`

# Pact

nano Pact_sample_list.txt

1041	Sample1_merged2.gtf
1471	Sample2_merged2.gtf
1637	Sample3_merged2.gtf

`
prepDE.py -i Pact_sample_list.txt -g Pact_gene_count_matrix.csv -t Pact_transcript_count_matrix.csv
`

scp -r -P 2292  hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare/RiboDep_RNASeq/cleaned_reads/ /Users/hputnam/MyProjects/Meth_Compare/RNASeq/



