# Required Programs

fastqc 
multiqc  
fastp  
samtools  
hisat2 
stringtie 

# RNASeq Files


#### File location PolyA 
##### P. acuta
```/RAID_STORAGE2/hputnam/20191104_HoloInt/30-274849628```
1041_R1_001.fastq.gz
1041_R2_001.fastq.gz
1471_R1_001.fastq.gz
1471_R2_001.fastq.gz
1637_R1_001.fastq.gz
1637_R2_001.fastq.gz

##### M. capitata
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

### Obtain Mcapitata Genome and GFF
`wget http://cyanophora.rutgers.edu/montipora/Mcap.genome_assembly.fa.gz`
`wget http://cyanophora.rutgers.edu/montipora/Mcap.GFFannotation.gff`
### Use script to fix the gff format issues introduced
### from the combination of AUGUSTUS and GeMoMa genes Fix_GeMoMa.R

	#Load Mcap gene gff
	Mcap.gff <- read.csv(file="genome-feature-files/Mcap.GFFannotation.gff.1", header=FALSE, sep="\t") 

	#rename columns
	colnames(Mcap.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")


	Mcap.gff$gene <- gsub("_cds[0123456789]", ".t1", Mcap.gff$gene) #change the second part of the GeMoMa genes from cds to t.1 to match Augustus
	Mcap.gff$gene <-sub("(;[^;]+);.*", "\\1", Mcap.gff$gene) #remove everything after the second ; in the gene column
	Mcap.gff$gene <- gsub("Parent=", "  gene_id ", Mcap.gff$gene) #remove ID= from GeMoMa genes

	#If id ==CDS replace ID= with transcript_id, else replace with nothing
	Mcap.gff <- Mcap.gff %>% 
  	mutate(gene = ifelse(Gene.Predict %in% c("GeMoMa") & 
                         id == "CDS" ,  
                       gsub("ID=", "transcript_id ", gene), gsub("ID=", "", gene)))

	#If id ==gene remove everything after the ; else replace with nothing
	Mcap.gff <- Mcap.gff %>% 
  	mutate(gene = ifelse(Gene.Predict %in% c("GeMoMa") & 
                         id == "gene" ,  
                       gsub(";.*", "", gene), gsub("ID=", "", gene)))

	# sub to add quotes around the transcript name
	Mcap.gff$gene <- gsub("transcript_id ", "transcript_id \"", Mcap.gff$gene) 
	Mcap.gff$gene <- gsub(";", "\";", Mcap.gff$gene) 

	#add quotes before the gene_id
	Mcap.gff$gene <- gsub("gene_id ", "gene_id \"", Mcap.gff$gene) 

	#If id ==CDS add "; at the end else replace with nothing
	Mcap.gff <- Mcap.gff %>% 
  	mutate(gene = ifelse(id == "CDS" ,  
                         paste0(gene, "\";"), paste0(gene, "")))

	#save file
	write.table(Mcap.gff, file="/Users/hputnam/MyProjects/Meth_Compare/genome-feature-files/Mcap.GFFannotation.fixed.gff", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)

### Obtain Pacuta Genome and GFF
`scp -P 2292  /Users/hputnam/Downloads/Pacuta_genome/Pocillopora_acuta_genome_v1.fasta hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare_RNASeq`

`scp -P 2292  /Users/hputnam/Downloads/Pacuta_genome/Structural_annotation_abintio.gff hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare_RNASeq`

`mv Structural_annotation_abintio.gff Pact_Structural_annotation_abintio.gff`

`sed 's/cov/_cov/' Pact_Structural_annotation_abintio.gff > Pact.GFFannotation.gff`

# Link Raw from RAID2
#### File location PolyA 
```/RAID_STORAGE2/hputnam/20191104_HoloInt/30-274849628```

# QC raw files

`fastqc *.fastq.gz`

`multiqc .`

# Trimming and QC

"1041" "1101" "1471" "1548" "1628" "1637" 

```
sh -c 'for file in "1041" "1101" "1471" "1548" "1628" "1637" 
do
fastp \
--in1 ${file}_R1_001.fastq.gz \
--in2 ${file}_R2_001.fastq.gz \
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
	fastqc 1637*clean.fastq.gz 

	multiqc .
# Download and view files

scp -r -P 2292 hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare_RNASeq/multiqc_report.html /Users/hputnam/MyProjects/Meth_Compare/RNASeq/PolyA



# Align Reads to Reference Genome
*HISAT2 is a fast and sensitive alignment program for mapping next-generation DNA and RNA sequencing reads to a reference genome.*

## Index the reference genome

Index the reference genome in the reference directory.

++HISAT2-build Alignment Arguments Used++:  
- <reference_in> - name of reference files  
- <gt2_base> -  basename of index files to write  
- -f -  reference file is a FASTA file

```
hisat2-build Mcap.genome_assembly.fa Mcap_ref
hisat2-build Pocillopora_acuta_genome_v1.fasta Pact_ref
```

# Alignment of clean reads to the reference genome

### Aligning paired end reads

rna strandedness - RF for Illumina
https://github.com/griffithlab/rnaseq_tutorial/blob/master/manuscript/supplementary_tables/supplementary_table_5.md


## Mcap Mapping for PolyA
"1101" "1548" "1628"
`
sh -c 'for i in "1101"
do
hisat2 -p 8 --dta -q -x Mcap_ref \
--rna-strandness RF \
-1 ${i}_R1_clean.fastq.gz \
-2 ${i}_R2_clean.fastq.gz -S ${i}.sam 
done'
`

## Pact Mapping for PolyA
"1041" "1471" "1637"
`
sh -c 'for i in "1637" 
do
hisat2 -p 8 --dta -q -x Pact_ref \
--rna-strandness RF \
-1 ${i}_R1_clean.fastq.gz \
-2 ${i}_R2_clean.fastq.gz -S ${i}.sam 
done'
`

# Convert and sort Sam to BAM
## RiboDep

## PolyA
"1041" "1471" "1637" "1101" "1548" "1628"
`
sh -c 'for i in "1637"
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

- Merged-GTF guided assembly without novel transcript discovery
- Compilation of GTF-files into gene and transcript count matrices

For our assembly, we will have to run StringTie. This run is necessary in order to compile our GTF files into gene count matrices that we will need for our differential expression analysis, because the StringTie script that compiles the GTF files ```prepDE.py``` only runs if the ```-e``` option is "on" during the previous assembly.

# Reference-guided assembly to genes only

++StringTie Arguments Used++:  
- -p - Specify number of processers
- -G - Specify annotation file
- -o - Name of output file
- -e - only estimate the abundance of given reference transcripts (requires -G)


#Mcap PolyA Ignore Novel
`
stringtie 1101.bam -p 8 -A 1101_gene_abund.tab -e -G Mcap.GFFannotation.fixed.gff -o 1101.gtf 
stringtie 1548.bam -p 8 -A 1548_gene_abund.tab -e -G Mcap.GFFannotation.fixed.gff -o 1548.gtf 
stringtie 1628.bam -p 8 -A 1628_gene_abund.tab -e -G Mcap.GFFannotation.fixed.gff -o 1628.gtf 
`
#Pact PolyA Ignore Novel
`
stringtie 1041.bam -p 8 -A 1041_gene_abund.tab -e -G Pact.GFFannotation.gff -o 1041.gtf 
stringtie 1471.bam -p 8 -A 1471_gene_abund.tab -e -G Pact.GFFannotation.gff -o 1471.gtf 
stringtie 1637.bam -p 8 -A 1637_gene_abund.tab -e -G Pact.GFFannotation.gff -o 1637.gtf 
`


#### Compilation of GTF-files into gene and transcript count matrices

The StringTie program includes a script, ```prepDE.py``` that compiles your assembly files into gene and transcript count matrices. This script requires as input the list of sample names and their full file paths, (sample_list.txt).

. Run ```prepDE.py``` to merge assembled files together into a DESeq2-friendly version.

++StringTie prepDE.py Arguments Used++:  
- -i - Specify that input is a TXT file
- -g - Require output gene count file, default name is ```gene_count_matrix.csv```
```


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


#download mapping and counts matrices

scp -r -P 2292  hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare_RNASeq/*.gtf /Users/hputnam/MyProjects/Meth_Compare/RNASeq/PolyA

scp -r -P 2292  hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare_RNASeq/*.csv /Users/hputnam/MyProjects/Meth_Compare/RNASeq/PolyA

scp -r -P 2292  hputnam@kitt.uri.edu:/home/hputnam/Meth_Compare_RNASeq/*gene_abund.tab /Users/hputnam/MyProjects/Meth_Compare/RNASeq/PolyA



