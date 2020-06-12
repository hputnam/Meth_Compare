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

sample.info <- read.csv(file="metadata/RNASeq_metadata.csv", header=T, sep=",", row.names=1) #load sample info
Pact.sample.info <- subset(sample.info, Species == "Pocillopora acuta")
Mcap.sample.info <- subset(sample.info, Species == "Montipora capitata")

PolyA.Pact.counts <- read.csv(file="RNASeq/PolyA/mapped/Pact_gene_count_matrix.csv", header=T, sep=",") #Load expressin matrix from trinity
RiboDep.Pact.counts <- read.csv(file="RNASeq/RiboDep/mapped/Pact_gene_count_matrix.csv", header=T, sep=",") #Load expressin matrix from trinity
PolyA.Mcap.counts <- read.csv(file="RNASeq/PolyA/mapped/Mcap_gene_count_matrix.csv", header=T, sep=",") #Load expressin matrix from trinity
RiboDep.Mcap.counts <- read.csv(file="RNASeq/RiboDep/mapped/Mcap_gene_count_matrix.csv", header=T, sep=",") #Load expressin matrix from trinity

Pact.counts <- left_join(PolyA.Pact.counts,RiboDep.Pact.counts, by="gene_id")
colnames(Pact.counts) <- c("gene_id", "PolyA_1041", "PolyA_1471", "PolyA_1637", "RiboDep_1041", "RiboDep_1471", "RiboDep_1637")
rownames(Pact.counts) <- Pact.counts$gene_id
Pact.counts <- Pact.counts[,-1]
Pact.counts <- as.matrix(Pact.counts)
str(Pact.counts)

# filt <- filterfun(pOverA(0.5,0)) #set filter values for PoverA, P percent of the samples have counts over A
# tfil <- genefilter(Pact.counts, filt) #create filter for the counts data
# keep <- Pact.counts[tfil,] #identify transcripts to keep by count filter
# Pact.keep <- rownames(keep) #identify transcript list
# Pact.counts <- as.matrix(Pact.counts[which(rownames(Pact.counts) %in% Pact.keep),]) #data filtered in PoverA, P percent of the samples have counts over A

Pact.data <- DESeqDataSetFromMatrix(countData = Pact.counts, colData = Pact.sample.info, design = ~ Method) #create a DESeqDataSet object
Pact.trans.data <- rlog(Pact.counts, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
row.names(Pact.trans.data) <- row.names(Pact.counts)
max(Pact.trans.data)

# pdf(file="Output/Pact_counts.pdf") #save file
# out <- pheatmap(Pact.trans.data, show_rownames =F, cluster_cols=FALSE, cluster_rows=TRUE,
#          show_colnames =T) #plot heatmap of all DEG by group
# dev.off()

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
DEG.sig.num <- sum(DEG.int.res$padj <0.05, na.rm=T) #identify the number of significant p values with 10%FDR (padj<0.1)
DEG.sig <- subset(DEG.int.res, padj<0.05,) #identify signficant pvalues with 10%FDR
DEG.sig.list <- Pact.data[which(rownames(Pact.data) %in% rownames(DEG.sig)),] #subset list of sig transcripts from original count data
write.csv(counts(DEG.sig.list), file="Output/Pact_DEG_PolyA_vs_RiboDep.csv")

DEG.rsig <- rlog(DEG.sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
PCA.plot <- plotPCA(DEG.rsig, intgroup = c("Method")) #Plot PCA of all samples for DEG only
PCA.plot #view plot

df <- as.data.frame(colData(DEG.rsig)[,c("Method")]) #make dataframe
mat <- as.matrix(counts(DEG.sig.list)) #make a matrix

#pdf(file="Output/Unique_Heatmap.DEG.Annotated.pdf") #save file
pheatmap(mat, annotation_col=df, scale="row",
         show_rownames =T, fontsize_row = 4, cluster_cols = FALSE,
         show_colnames =F) #plot heatmap of all DEG by group
#dev.off()








Mcap.counts <- left_join(PolyA.Mcap.counts,RiboDep.Mcap.counts, by="gene_id")
colnames(Mcap.counts) <- c("gene_id", "PolyA_1101", "PolyA_1548", "PolyA_1628", "RiboDep_1101", "RiboDep_1548", "RiboDep_1628")
rownames(Mcap.counts) <- Mcap.counts$gene_id
Mcap.counts <- Mcap.counts[,-1]
Mcap.counts <- as.matrix(Mcap.counts)
str(Mcap.counts)






#K means clustering Elbow Method
set.seed(124) #Specify which set of random numbers to use. This will make sure that the same set of random numbers is generated each timea random set of numbers is called. 
wss <- function(k) {kmeans(trans.data, k, iter.max = 100, nstart = 25 )$tot.withins} # function to compute total within-cluster sum of square 
k.values <- 1:15# Compute and plot wss for k = 1 to k = 15
wss_values <- map_dbl(k.values, wss) # extract wss for 2-15 clusters
pdf(file="Output/Mcap_elbow.pdf")
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters (K)",
     ylab="Total within-clusters sum of squares")
dev.off()

##### Extract the clusters
Mcap.clusts <- as.data.frame(sort(cutree(out$tree_row, k=6)))
Mcap.clusts$X.gene_id <- rownames(Mcap.clusts)
Mcap.annot <- read.csv(file="Data/Mcap-GO-KO-Kegg.tab", header=FALSE, sep="\t") #Load expression matrix f
colnames(Mcap.annot) <- c("Uniprot", "X.gene_id", "eval", "Prot.ID", "Rev", "Prot.Name.Long", "Prot.Name.Short", "Taxa", "Num", "GO1", "GO2", "GO3", "GO4", "GO.IDs","KEGG", "KEGG.Path")  
Mcap.annot <- Mcap.annot %>% 
  distinct(X.gene_id, .keep_all = TRUE)
Mcap.annot$X.gene_id <- gsub("augustus.", "", Mcap.annot$X.gene_id)
Mcap.annot$X.gene_id <- gsub(".t1", "", Mcap.annot$X.gene_id)
trans.df <- as.data.frame(trans.data)
trans.df$X.gene_id <- rownames(trans.data)
Mcap.gene.annot <- left_join(trans.df, Mcap.annot)
Mcap.gene.clust.annot <- left_join(Mcap.gene.annot, Mcap.clusts)
colnames(Mcap.gene.clust.annot)[25] <- "Cluster"

Mcap.mean.clust <- Mcap.gene.clust.annot %>%
  group_by(Cluster) %>%
  summarise(rlog.exp= mean(X153:X160))

#subset highest clusters 2, 5 and 6
Mcap.clusters <- Mcap.gene.clust.annot %>%
  filter(Cluster %in% 2 | Cluster %in% 5 | Cluster %in% 6)

#Amillepora
sample.info <- read.csv(file="Data/Amil_embryo_sample.info.csv", header=T, sep=",", row.names=3) #load sample info
counts <- read.csv(file="Data/Amil_gene_count_matrix.csv", header=T, sep=",") #Load expressin matrix from trinity
counts <- counts[,1:5]
rownames(counts) <- counts$gene_id
counts <- counts[,-1]
str(counts)
#Filter reads by proportion of samples containing the cutoff value
filt <- filterfun(pOverA(1.0,10)) #set filter values for PoverA, P percent of the samples have counts over A
tfil <- genefilter(counts, filt) #create filter for the counts data
keep <- counts[tfil,] #identify transcripts to keep by count filter
gn.keep <- rownames(keep) #identify transcript list
counts.5x <- as.matrix(counts[which(rownames(counts) %in% gn.keep),]) #data filtered in PoverA, P percent of the samples have counts over A
storage.mode(counts.5x) = "integer" #store counts data as integer
data <- DESeqDataSetFromMatrix(countData = counts.5x, colData = sample.info, design = ~ 1) #create a DESeqDataSet object
trans.data <- rlog(counts.5x, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
row.names(trans.data) <- row.names(counts.5x)
max(trans.data)

pdf(file="Output/Amillepora_rlog_Embryo_counts.pdf") #save file
out <- pheatmap(trans.data, show_rownames =T, fontsize_row = 2, cluster_cols=FALSE, cluster_rows=TRUE, cutree_rows=6,
         show_colnames =T) #plot heatmap of all DEG by group
dev.off()

#K means clustering Elbow Method
set.seed(124) #Specify which set of random numbers to use. This will make sure that the same set of random numbers is generated each timea random set of numbers is called. 
wss <- function(k) {kmeans(trans.data, k, iter.max = 100, nstart = 25 )$tot.withins} # function to compute total within-cluster sum of square 
k.values <- 1:15# Compute and plot wss for k = 1 to k = 15
wss_values <- map_dbl(k.values, wss) # extract wss for 2-15 clusters
pdf(file="Ouput/Amill_elbow.pdf")
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters (K)",
     ylab="Total within-clusters sum of squares")
dev.off()

##### Extract the clusters
Amil.clusts <- as.data.frame(sort(cutree(out$tree_row, k=6)))
Amil.clusts$X.gene_id <- rownames(Amil.clusts)
Amil.annot <- read.csv(file="Data/Amillepora_trinotate_annotation_report.csv", header=T, sep=",") #Load expressin matrix from trinity
rownames(Amil.annot) <- Amil.annot$X.gene_id
trans.df <- as.data.frame(trans.data)
trans.df$X.gene_id <- rownames(trans.data)
Amil.gene.annot <- left_join(trans.df, Amil.annot)
Amil.gene.clust.annot <- left_join(Amil.gene.annot, Amil.clusts)
colnames(Amil.gene.clust.annot)[21] <- "Cluster"

Amil.mean.clust <- Amil.gene.clust.annot %>%
  group_by(Cluster) %>%
  summarise(rlog.exp= mean(T0_A:T0_D))

#subset highest clusters 4, 5 and 6
Amil.clusters <- Amil.gene.clust.annot %>%
  filter(Cluster %in% 3 | Cluster %in% 4 | Cluster %in% 5 | Cluster %in% 6)


#clean up annotation columns
Amil.clusters$gene_ontology_blast2 <- Amil.clusters$gene_ontology_blast #add another copy of the GO information
Amil.clusters.x <- separate(data = Amil.clusters, col = sprot_Top_BLASTX_hit, into = c("Protein.Name", "Full.Name"), sep = "\\^", extra="merge") #split columns to get protein ID
Amil.clusters.x <- separate(data = Amil.clusters.x, col = Full.Name, into = c( "Full.Name", "Name"), sep = "\\=", extra="merge") #split columns to get protein name
Amil.clusters.x <- separate(data = Amil.clusters.x, col = Name, into = c( "Full.Name"), sep = "\\;", extra="drop") #split columns to get protein name
Amil.clusters.x <- separate(data = Amil.clusters.x, col = gene_ontology_blast, into = c("CC", "MF"), sep = "GO\\:[0-9]*\\^molecular_function\\^", extra="merge") #split columns to get CC
Amil.clusters.x <- separate(data = Amil.clusters.x, col = gene_ontology_blast2, into = c("Funct", "BP"), sep = "GO\\:[0-9]*\\^biological_process\\^", extra="merge") #split columns to get BP
Amil.clusters.x <- separate(data = Amil.clusters.x, col = CC, into = c("CC"), sep = "GO\\:[0-9]*\\^biological_process\\^", extra="drop") #split columns to get CC
Amil.clusters.x <- separate(data = Amil.clusters.x, col = MF, into = c("MF"), sep = "GO\\:[0-9]*\\^biological_process\\^", extra="drop") #split columns to get MF

