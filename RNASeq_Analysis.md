#Required Programs
fastqc 
multiqc  
fastp  
samtools  
hisat2 
stringtie 


`git clone https://github.com/gpertea/stringtie`  
`cd stringtie`  
`make release`  
`make test`  

#Load raw sequence files onto server
`mkdir raw`  


#QC raw files
`mkdir fastqc_raw`  
`cd fastqc_raw`  
`fastqc raw/*.fastq.gz`  


#Trimming
`mkdir cleaned_reads'  


sh -c 'for file in "Sample1" "Sample2" "Sample3" "Sample4" "Sample5" "Sample6" 
do
fastp --in1 ${file}_R1.fastq.gz --in2 ${file}_R2.fastq.gz --out1 ../cleaned_reads/${file}_R1_clean.fastq.gz --out2 ../cleaned_reads/${file}_R2_clean.fastq.gz --failed_out ../cleaned_reads/${file}_failed.txt --qualified_quality_phred 20 --unqualified_percent_limit 10 --length_required 100 detect_adapter_for_pe --cut_right cut_right_window_size 5 cut_right_mean_quality 20
done'


#QC trimmed files

`fastqc cleaned_reads/*.fastq.gz`  



*HISAT2 is a fast and sensitive alignment program for mapping next-generation DNA and RNA sequencing reads to a reference genome.*

- Index the reference genome
- Alignment of clean reads to the reference genome

Create a subdirectory within data for HISAT2
```
mkdir hisat2
cd hisat2
```

#### Index the reference genome

Index the reference genome in the reference directory.

++HISAT2-build Alignment Arguments Used++:  
- <reference_in> - name of reference files  
- <gt2_base> -  basename of index files to write  
- -f -  reference file is a FASTA file

```
cd ref
hisat2-build -f ../ref/Mcap.genome_assembly.fa.gz ./Mcap_ref
```

#### Alignment of clean reads to the reference genome

Align your reads to the index files. We will do this by writing a script we will call ```McapHISAT2.sh```. This script will also take the output SAM files from our HISAT2 alignment and covert them into the sorted BAM files that are the necessary input for our assembly tool, StringTie. We do this by calling SAMtools in our script.

++HISAT2 Alignment Arguments Used++:   
- -x <hisat2-idx> - Basename of index files to read  
- -1 <m1> - List of forward sequence files  
- -2 <m1> - List of reverse sequence files  
- -S - Name of output files
- -q - Input files are in FASTQ format  
- -p - Number processors
- --dta - Adds the XS tag to indicate the genomic strand that produced the RNA from which the read was sequenced. As noted by StringTie... "be sure to run HISAT2 with the --dta option for alignment, or your results will suffer."

++SAMtools Options Arguments Used++:  
- -@ - Number threads  
- -o - Output file  

```
nano McapHISAT2.sh
```
```
##!/bin/bash

#Specify working directory
F=/home/echille/mcap2019/data/hisat2

#Aligning paired end reads
#Has the R1 in array1 because the sed in the for loop changes it to an R2. SAM files are of both forward and reverse reads
array1=($(ls $F/*_R1_001.clean.fastq.gz))

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

for i in ${array1[@]}; do
        hisat2 -p 8 --dta -q -x Mcap_ref -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -S ${i}.sam
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 PE ${i}" $(date)
done
```

Now, make the file executable by the user (you) and run the script.
```
chmod u+x McapHISAT2.sh
./McapHISAT2.sh
```
Now we've got some sorted BAM files that can be used in our assembly!!

### Assemble aligned reads and quantify transcripts 

---

*StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts.*

- Reference-guided assembly with novel transcript discovery
- Merge output GTF files and assess the assembly performance
- Merged-GTF guided assembly without novel transcript discovery
- Compilation of GTF-files into gene and transcript count matrices

For our assembly, we will have to run StringTie twice. The first run will be a reference guided assembly that will allow for discovery of novel transcripts (by leaving out the -e option). Then we will merge the output GTF files and examine the sensitivity of our assembly. We will use the merged-GTF from our first assembly to guide our second StringTie run (including the ```-e``` option). This second run is necessary in order to compile our GTF files into gene and transcript count matrices that we will need for our differential expression analysis, because the StringTie script that compiles the GTF files ```prepDE.py``` only runs if the ```-e``` option is "on" during the previous assembly.

#### Reference-guided assembly with novel transcript discovery

First, create and enter into StringTie directory. Then create a symbolic link to our reference genome and copy our BAM files to a special directory inside our stringtie directory. This is where our output GTF files will live too.
```
mkdir ../stringtie
cd stringtie
ln -s ../ref/Mcap.GFFannotation.gff ./
mkdir BAM
cd BAM
ln -s ../../hisat2/*.bam ./
cd ../
```

Create the StringTie reference-guided assembly script, ```McapStringTie-assembly-to-ref.sh``` *inside of the StringTie program directory.*  

++StringTie Arguments Used++:  
- -p - Specify number of processers
- -G - Specify annotation file
- -o - Name of output file

```
cd stringtie
nano ./McapStringTie-assembly-to-ref.sh
```
```
##!/bin/bash

#Specify working directory
F=/home/echille/mcap2019/data/stringtie

#StringTie reference-guided assembly
#Has the R1 in array1 because of the naming convention in the former script. However, these BAM files contain both forward and reverse reads.
array1=($(ls $F/BAM/*_R1_001.clean.fastq.gz.bam))

#Running without the -e option the first time to increase novel transcript discovery

for i in ${array1[@]}; do
        /stringtie/stringtie -p 8 -G Mcap.GFFannotation.gff -o ${i}.gtf ${i}
        mv /ref-guided-gtfs/${i}.gtf
        echo "StringTie-assembly-to-ref ${i}" $(date)
done
```

Now, make the file executable by the user and run the script.
```
chmod u+x McapStringTie-assembly-to-ref.sh
./McapStringTie-assembly-to-ref.sh
```

#### Assess the performance of the assembly

*Gffcompare is a tool that can compare, merge, annotate and estimate accuracy of GFF/GTF files when compared with a reference annotation*

Using the StringTie merge mode, merge the assembly-generated GTF files to assess how well the predicted transcripts track to the reference annotation file. This step requires the TXT file,  [```mergelist.txt```](https://github.com/echille/Montipora_OA_Development_Timeseries/blob/master/RNAseq_Analyses/mergelist.txt). This file lists all of the file names to be merged. *Make sure ```mergelist.txt``` is in the StringTie program directory*.

++StringTie Arguments Used++:  
- --merge - Distinct from the assembly usage mode used above, in the merge mode, StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts.
- -p - Specify number of processers
- -G - Specify reference annotation file. With this option, StringTie assembles the transfrags from the input GTF files with the reference sequences
- -o - Name of output file
- <mergelist.txt> - File listing all filenames to be merged. Include full path.

```
stringtie --merge -p 8 -G ../Mcap.GFFannotation.gff -o ../stringtie_merged.gtf mergelist.txt
```

Now we can use the program gffcompare to compare the merged GTF to our reference genome.

++Gffcompare Arguments Used++:  
- -r - Specify reference annotation file
- -G - Compare all the transcripts in our input file ```stringtie_merged.gtf```
- -o - Prefix of all output files

```
gffcompare -r ../ref/Mcap.GFFannotation.gff -G -o ../merged ../stringtie_merged.gtf
```

Some of the output files you will see are... 
- merged.stats
- merged.tracking
- merged.annotated.gtf
- merged.stringtie_merged.gtf.refmap
- merged.loci
- merged.stringtie_merged.gtf.tmap

Move all of the gffcompare output files to the output directory. We are most interested in the files ```merged.annotation.gtf``` and ```merged.stats```. The file ```merged.annotation.gtf``` tells you how well the predicted transcripts track to the reference annotation file and the file ```merged.stats``` file shows the sensitivity and precision statistics and total number for different features (genes, exons, transcripts).  Then, from the local host securely copy ```merged.stats``` to a local directory. Unfortunately, ```merged.annotation.gtf``` is too big to store locally, but we can view it remotely.

```
mv ./merged.* ../../ouput
```
```
scp -P xxxx echille@kitt.uri.edu:<path_to_output>/merged.stats /Users/user/<path_to_local_directory>
```

#### Merged-GTF guided assembly without novel transcript discovery

Now run StringTie again, but this time the assembly guided by the merged GTF file ```stringtie_merged.gtf```. Here, StringTie skips over novel sequences because we include the ```-e``` option. This is okay now because we identified novel transcripts in the previous StringTie run. The ```-e``` option is the important part of this second run because the ```prepDE.py``` script used in the next step to compile the output GTFs only works if the ```-e``` option is included here.

Create the StringTie merged GTF-guided assembly script, ```McapStringTie-assembly-to-merged.sh``` *inside of the StringTie program directory,*. We should already be there, but check first.  

++StringTie Arguments Used++:  
- -p - Specify number of processers
- -e - Limit the estimation and output of transcripts to only those that match the reference (in this case, our merged GTF)
- -G - Specify annotation file
- -o - Name of output file

```
pwd
nano ./McapStringTie-assembly-to-merged.sh
```
```
##!/bin/bash

#StringTie reference-guided assembly
#Has the R1 in array1 because of the naming convention in the former script. However, these BAM files contain both forward and reverse reads.
array1=($(ls /home/echille/mcap2019/data/stringtie/BAM/*_R1_001.clean.fastq.gz.bam))

#Running without the -e option so that the output files will be compatible with the prepDE.py script that will be used later

for i in ${array1[@]}; do
        ./stringtie -p 8 -e -G ../stringtie_merged.gtf -o ${i}.merged2.gtf ${i}
        echo "StringTie-assembly-to-merged ${i}" $(date)
done
```

Now, make the file executable by the user and run the script.
```
chmod u+x McapStringTie-assembly-to-merged.sh
./McapStringTie-assembly-to-merged.sh
```

#### Compilation of GTF-files into gene and transcript count matrices

The StringTie program includes a script, ```prepDE.py``` that compiles your assembly files into gene and transcript count matrices. This script requires as input the list of sample names and their full file paths, [```sample_list.txt```](https://github.com/echille/Montipora_OA_Development_Timeseries/blob/master/RNAseq_Analyses/sample_list.txt). This file will live in StringTie program directory.

Go back into your stringtie directory (the one I should have named assembly). Run ```prepDE.py``` to merge assembled files together into a DESeq2-friendly version.

++StringTie prepDE.py Arguments Used++:  
- -i - Specify that input is a TXT file
- -g - Require output gene count file, default name is ```gene_count_matrix.csv```
- -t - Require output transcript count gene count file, default name is ```transcript_count_matrix.csv```

```
cd ../
./stringtie/prepDE.py -g -t -i ./stringtie/sample_list.txt
```

Finally, move your count matrices into the output directory and securely copy them to your local directory, from your local host.
```
mv ./*.csv ../../output



