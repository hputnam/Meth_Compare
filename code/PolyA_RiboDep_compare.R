#last modified 20200530
#Copyright 2020 HM Putnam
#Use of this code must be accompanied by a citation to XXXX
#Data should not be used without permission from HM Putnam
#See Readme

rm(list=ls()) # removes all prior objects

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("ggplot2")
library("gplots")
library("tidyverse")
library("spdep") 
library("adegenet")
library("UpSetR")
library("goseq")
library("GO.db")
#library("GSEABase")

###
sample.info <- read.csv(file="metadata/RNASeq_metadata.csv", header=T, sep=",", row.names=1) #load sample info
Pact.sample.info <- subset(sample.info, Species == "Pocillopora acuta")
Mcap.sample.info <- subset(sample.info, Species == "Montipora capitata")
Mcap.annot <-  read.csv(file="genome-feature-files/Mcap-GO-KO-Kegg.tab", header=FALSE, sep="\t") #load sample info
colnames(Mcap.annot) <- c("Uniprot", "X.gene_id", "eval", "Prot.ID", "Rev", "Prot.Name.Long", "Prot.Name.Short", "Taxa", "Num", "GO1", "GO2", "GO3", "GO4", "GO.IDs","KEGG", "KEGG.Path")  
Mcap.annot <- Mcap.annot %>% 
  distinct(X.gene_id, .keep_all = TRUE)
Mcap.annot$X.gene_id <- gsub("augustus.", "", Mcap.annot$X.gene_id)
Mcap.annot$X.gene_id <- gsub(".t1", "", Mcap.annot$X.gene_id)


### Pacuta
PolyA.Pact.counts <- read.csv(file="RNASeq/PolyA/mapped/Pact_gene_count_matrix.csv", header=T, sep=",") #Load expression matrix from trinity
RiboDep.Pact.counts <- read.csv(file="RNASeq/RiboDep/mapped/Pact_gene_count_matrix.csv", header=T, sep=",") #Load expression matrix from trinity

Pact.counts <- left_join(PolyA.Pact.counts,RiboDep.Pact.counts, by="gene_id")
colnames(Pact.counts) <- c("gene_id", "PolyA_1041", "PolyA_1471", "PolyA_1637", "RiboDep_1041", "RiboDep_1471", "RiboDep_1637")
rownames(Pact.counts) <- Pact.counts$gene_id
Pact.counts <- Pact.counts[,-1]
Pact.counts <- as.matrix(Pact.counts)
str(Pact.counts)

filt <- filterfun(pOverA(0.5,0)) #set filter values for PoverA, P percent of the samples have counts over A
tfil <- genefilter(Pact.counts, filt) #create filter for the counts data
keep <- Pact.counts[tfil,] #identify transcripts to keep by count filter
Pact.keep <- rownames(keep) #identify transcript list
Pact.counts <- as.matrix(Pact.counts[which(rownames(Pact.counts) %in% Pact.keep),]) #data filtered in PoverA, P percent of the samples have counts over A

Pact.data <- DESeqDataSetFromMatrix(countData = Pact.counts, colData = Pact.sample.info, design = ~ Method) #create a DESeqDataSet object
Pact.rld <- rlog(Pact.data, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size

# Expression Visualization
sampleDists <- dist(t(assay(Pact.rld))) #calculate distance matix
sampleDistMatrix <- as.matrix(sampleDists) #distance matrix
rownames(sampleDistMatrix) <- colnames(Pact.rld) #assign row names
colnames(sampleDistMatrix) <- NULL #assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
pheatmap(sampleDistMatrix, #plot matrix
         clustering_distance_rows=sampleDists, #cluster rows
         clustering_distance_cols=sampleDists, #cluster columns
         col=colors) #set colors

plotPCA(Pact.rld, intgroup = c("Method")) #plot PCA of samples with all data

# Differential Gene Expression Analysis
#Interaction Test: test of the factor of "group" with all combinations of the original factors as groups
DEG.int <- DESeq(Pact.data) #run differential expression test by group using the wald test
DEG.int.res <- results(DEG.int) #save DE results
resultsNames(DEG.int) #view DE results
DEG.int.res <- results(DEG.int, contrast=c("Method", "PolyA", "RiboDep"))
DEG.int.res
DEG.sig.num <- sum(DEG.int.res$padj <0.05 & DEG.int.res$log2FoldChange<abs(1), na.rm=T) #identify the number of significant p values with 5%FDR 
DEG.sig <- subset(DEG.int.res, padj<0.05 & DEG.int.res$log2FoldChange<abs(1)) #
DEG.sig.list <- Pact.data[which(rownames(Pact.data) %in% rownames(DEG.sig)),] #subset list of sig transcripts from original count data
write.csv(counts(DEG.sig.list), file="Output/Pact_DEG_PolyA_vs_RiboDep.csv")

DEG.rsig <- rlog(DEG.sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
PCA.plot <- plotPCA(DEG.rsig, intgroup = c("Method")) #Plot PCA of all samples for DEG only
PCA.plot #view plot

df <- as.data.frame(colData(DEG.rsig)[,c("Method")]) #make dataframe
rownames(df) <- colnames(mat)
rownames(df)

mat <- as.matrix(counts(DEG.sig.list)) #make a matrix
colnames(mat)

pdf(file="Output/Unique_Heatmap.DEG.Annotated.pdf") #save file
pheatmap(mat, annotation_col=df, scale="row",
         show_rownames =T, fontsize_row = 4, cluster_cols = FALSE,
         show_colnames =T) #plot heatmap of all DEG by group
dev.off()



### Mcapitata
PolyA.Mcap.counts <- read.csv(file="RNASeq/PolyA/mapped/Mcap_gene_count_matrix.csv", header=T, sep=",") #Load expression matrix from trinity
RiboDep.Mcap.counts <- read.csv(file="RNASeq/RiboDep/mapped/Mcap_gene_count_matrix.csv", header=T, sep=",") #Load expression matrix from trinity

Mcap.counts <- left_join(PolyA.Mcap.counts,RiboDep.Mcap.counts, by="gene_id")
colnames(Mcap.counts) <- c("gene_id", "PolyA_1101", "PolyA_1548", "PolyA_1628", "RiboDep_1101", "RiboDep_1548", "RiboDep_1628")
rownames(Mcap.counts) <- Mcap.counts$gene_id
Mcap.counts <- Mcap.counts[,-1]
Mcap.counts <- as.matrix(Mcap.counts)
str(Mcap.counts)

filt <- filterfun(pOverA(0.5,0)) #set filter values for PoverA, P percent of the samples have counts over A
tfil <- genefilter(Mcap.counts, filt) #create filter for the counts data
keep <- Mcap.counts[tfil,] #identify transcripts to keep by count filter
Mcap.keep <- rownames(keep) #identify transcript list
Mcap.counts <- as.matrix(Mcap.counts[which(rownames(Mcap.counts) %in% Mcap.keep),]) #data filtered in PoverA, P percent of the samples have counts over A

Mcap.data <- DESeqDataSetFromMatrix(countData = Mcap.counts, colData = Mcap.sample.info, design = ~ Method) #create a DESeqDataSet object
Mcap.rld <- rlog(Mcap.data, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size

# Expression Visualization
sampleDists <- dist(t(assay(Mcap.rld))) #calculate distance matix
sampleDistMatrix <- as.matrix(sampleDists) #distance matrix
rownames(sampleDistMatrix) <- colnames(Mcap.rld) #assign row names
colnames(sampleDistMatrix) <- NULL #assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
pheatmap(sampleDistMatrix, #plot matrix
         clustering_distance_rows=sampleDists, #cluster rows
         clustering_distance_cols=sampleDists, #cluster columns
         col=colors) #set colors

plotPCA(Mcap.rld, intgroup = c("Method")) #plot PCA of samples with all data

# Differential Gene Expression Analysis
#Interaction Test: test of the factor of "group" with all combinations of the original factors as groups
DEG.int <- DESeq(Mcap.data) #run differential expression test by group using the wald test
DEG.int.res <- results(DEG.int) #save DE results
resultsNames(DEG.int) #view DE results
DEG.int.res <- results(DEG.int, contrast=c("Method", "PolyA", "RiboDep"))
DEG.int.res
DEG.sig.num <- sum(DEG.int.res$padj <0.05 & DEG.int.res$log2FoldChange<abs(1), na.rm=T) #identify the number of significant p values with 5%FDR 
DEG.sig <- subset(DEG.int.res, padj<0.05 & DEG.int.res$log2FoldChange<abs(1)) #
DEG.sig.list <- Mcap.data[which(rownames(Mcap.data) %in% rownames(DEG.sig)),] #subset list of sig transcripts from original count data
write.csv(counts(DEG.sig.list), file="Output/Mcap_DEG_PolyA_vs_RiboDep.csv")

DEG.rsig <- rlog(DEG.sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
PCA.plot <- plotPCA(DEG.rsig, intgroup = c("Method")) #Plot PCA of all samples for DEG only
PCA.plot #view plot

df <- as.data.frame(colData(DEG.rsig)[,c("Method")]) #make dataframe
rownames(df) <- colnames(mat)
rownames(df)

mat <- as.matrix(counts(DEG.sig.list)) #make a matrix
colnames(mat)

pdf(file="Output/Unique_Heatmap.DEG.Annotated.pdf") #save file
pheatmap(mat, annotation_col=df, scale="row",
         show_rownames =T, fontsize_row = 4, cluster_cols = FALSE,
         show_colnames =T) #plot heatmap of all DEG by group
dev.off()

```{bash}
#Intersect CpG with genes
for f in RAnalysis/Data/OSF/Time/*5x.tab
do
intersectBed \
-wb \
-a ${f} \
-b RAnalysis/Data/Genome/Panopea-generosa-v1.0.a4.gene.gff3 \
> ${f}_gene
done
```
#merge with annotation

#identify GO terms


