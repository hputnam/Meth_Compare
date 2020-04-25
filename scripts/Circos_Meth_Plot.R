library(circlize) 
library(tidyverse)
library(ComplexHeatmap)

#https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html#customize-chromosome-track
#The input data for circos.genomicInitialize() is also a data frame with at least three columns. The first column is genomic category (for cytoband data, it is chromosome name), and the next two columns are positions in each genomic category.

##### POCILLOPORA #####
#load in all genes gff from a genome
Pact.cytoband <- read.table("data/Pact.GFFannotation.Genes.gff", sep = "\t",  
                       header = FALSE,colClasses = c("character","character","character", "numeric",
                       "numeric",  "numeric", "character", "character", "character") )

#remove extra coloumns from gene gff
Pact.cytoband <- Pact.cytoband[,c(1,4,5,6,9)]
str(Pact.cytoband)
head(Pact.cytoband)

n_distinct(Pact.cytoband$V1)

#examine scaffolds for the number of genes in a scaffold 
sizes <- Pact.cytoband %>% 
  group_by(V1) %>% 
  summarise(count=n()) %>%
  arrange(-count)

Pact.cytoband.df <- Pact.cytoband %>% filter(V1 == "scaffold502cov107"| V1 == "scaffold3517cov108" | 
                                                  V1 == "scaffold179cov105"| V1 == "scaffold382cov105" | 
                                                  V1 == "scaffold2581cov101" | V1 == "scaffold1396cov105" |
                                                  V1 == "scaffold150315cov104" | V1 == "scaffold2081cov104"| 
                                                  V1 == "scaffold3008cov108" | V1 == "scaffold222cov106" )


#read in the 5x coverage of CpG 
Pact.1 <- read.table("data/Pact_tab/Meth1_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Pact.1) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Pact.1)
str(Pact.1)

Pact.2 <- read.table("data/Pact_tab/Meth2_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Pact.2) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Pact.2)
str(Pact.2)

Pact.3 <- read.table("data/Pact_tab/Meth3_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Pact.3) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Pact.3)
str(Pact.3)

Pact.4 <- read.table("data/Pact_tab/Meth4_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Pact.4) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Pact.4)
str(Pact.4)

Pact.5 <- read.table("data/Pact_tab/Meth5_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Pact.5) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Pact.5)
str(Pact.5)

Pact.6 <- read.table("data/Pact_tab/Meth6_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Pact.6) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Pact.6)
str(Pact.6)

Pact.7 <- read.table("data/Pact_tab/Meth7_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Pact.7) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Pact.7)
str(Pact.7)

Pact.8 <- read.table("data/Pact_tab/Meth8_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Pact.8) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Pact.8)
str(Pact.8)

Pact.9 <- read.table("data/Pact_tab/Meth9_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Pact.9) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Pact.9)
str(Pact.9)

Pact.1.1 <- Pact.1 %>% filter(chrom == "scaffold502_cov107"| chrom == "scaffold3517_cov108" | 
                                chrom == "scaffold179_cov105"| chrom == "scaffold382_cov105" | 
                                chrom == "scaffold2581_cov101" | chrom == "scaffold1396_cov105" |
                                chrom == "scaffold150315_cov104" | chrom == "scaffold2081_cov104"| 
                                chrom == "scaffold3008_cov108" | chrom == "scaffold222_cov106")
Pact.1.1$chrom <- gsub("_", "", Pact.1.1$chrom)
Pact.1.1$avg <- Pact.1.1$per.meth/100
Pact.1.1 <- Pact.1.1[,c(1,2,3,7)]

Pact.2.1 <- Pact.2 %>% filter(chrom == "scaffold502_cov107"| chrom == "scaffold3517_cov108" | 
                                chrom == "scaffold179_cov105"| chrom == "scaffold382_cov105" | 
                                chrom == "scaffold2581_cov101" | chrom == "scaffold1396_cov105" |
                                chrom == "scaffold150315_cov104" | chrom == "scaffold2081_cov104"| 
                                chrom == "scaffold3008_cov108" | chrom == "scaffold222_cov106")
Pact.2.1$chrom <- gsub("_", "", Pact.2.1$chrom)
Pact.2.1$avg <- Pact.2.1$per.meth/100
Pact.2.1 <- Pact.2.1[,c(1,2,3,7)]

Pact.3.1 <- Pact.3 %>% filter(chrom == "scaffold502_cov107"| chrom == "scaffold3517_cov108" | 
                                chrom == "scaffold179_cov105"| chrom == "scaffold382_cov105" | 
                                chrom == "scaffold2581_cov101" | chrom == "scaffold1396_cov105" |
                                chrom == "scaffold150315_cov104" | chrom == "scaffold2081_cov104"| 
                                chrom == "scaffold3008_cov108" | chrom == "scaffold222_cov106")
Pact.3.1$chrom <- gsub("_", "", Pact.3.1$chrom)
Pact.3.1$avg <- Pact.3.1$per.meth/100
Pact.3.1 <- Pact.3.1[,c(1,2,3,7)]

Pact.4.1 <- Pact.4 %>% filter(chrom == "scaffold502_cov107"| chrom == "scaffold3517_cov108" | 
                                chrom == "scaffold179_cov105"| chrom == "scaffold382_cov105" | 
                                chrom == "scaffold2581_cov101" | chrom == "scaffold1396_cov105" |
                                chrom == "scaffold150315_cov104" | chrom == "scaffold2081_cov104"| 
                                chrom == "scaffold3008_cov108" | chrom == "scaffold222_cov106")
Pact.4.1$chrom <- gsub("_", "", Pact.4.1$chrom)
Pact.4.1$avg <- Pact.4.1$per.meth/100
Pact.4.1 <- Pact.4.1[,c(1,2,3,7)]

Pact.5.1 <- Pact.5 %>% filter(chrom == "scaffold502_cov107"| chrom == "scaffold3517_cov108" | 
                                chrom == "scaffold179_cov105"| chrom == "scaffold382_cov105" | 
                                chrom == "scaffold2581_cov101" | chrom == "scaffold1396_cov105" |
                                chrom == "scaffold150315_cov104" | chrom == "scaffold2081_cov104"| 
                                chrom == "scaffold3008_cov108" | chrom == "scaffold222_cov106")
Pact.5.1$chrom <- gsub("_", "", Pact.5.1$chrom)
Pact.5.1$avg <- Pact.5.1$per.meth/100
Pact.5.1 <- Pact.5.1[,c(1,2,3,7)]

Pact.6.1 <- Pact.6 %>% filter(chrom == "scaffold502_cov107"| chrom == "scaffold3517_cov108" | 
                                chrom == "scaffold179_cov105"| chrom == "scaffold382_cov105" | 
                                chrom == "scaffold2581_cov101" | chrom == "scaffold1396_cov105" |
                                chrom == "scaffold150315_cov104" | chrom == "scaffold2081_cov104"| 
                                chrom == "scaffold3008_cov108" | chrom == "scaffold222_cov106")
Pact.6.1$chrom <- gsub("_", "", Pact.6.1$chrom)
Pact.6.1$avg <- Pact.6.1$per.meth/100
Pact.6.1 <- Pact.6.1[,c(1,2,3,7)]

Pact.7.1 <- Pact.7 %>% filter(chrom == "scaffold502_cov107"| chrom == "scaffold3517_cov108" | 
                                chrom == "scaffold179_cov105"| chrom == "scaffold382_cov105" | 
                                chrom == "scaffold2581_cov101" | chrom == "scaffold1396_cov105" |
                                chrom == "scaffold150315_cov104" | chrom == "scaffold2081_cov104"| 
                                chrom == "scaffold3008_cov108" | chrom == "scaffold222_cov106")
Pact.7.1$chrom <- gsub("_", "", Pact.7.1$chrom)
Pact.7.1$avg <- Pact.7.1$per.meth/100
Pact.7.1 <- Pact.7.1[,c(1,2,3,7)]

Pact.8.1 <- Pact.8 %>% filter(chrom == "scaffold502_cov107"| chrom == "scaffold3517_cov108" | 
                                chrom == "scaffold179_cov105"| chrom == "scaffold382_cov105" | 
                                chrom == "scaffold2581_cov101" | chrom == "scaffold1396_cov105" |
                                chrom == "scaffold150315_cov104" | chrom == "scaffold2081_cov104"| 
                                chrom == "scaffold3008_cov108" | chrom == "scaffold222_cov106")
Pact.8.1$chrom <- gsub("_", "", Pact.8.1$chrom)
Pact.8.1$avg <- Pact.8.1$per.meth/100
Pact.8.1 <- Pact.8.1[,c(1,2,3,7)]

Pact.9.1 <- Pact.9 %>% filter(chrom == "scaffold502_cov107"| chrom == "scaffold3517_cov108" | 
                                chrom == "scaffold179_cov105"| chrom == "scaffold382_cov105" | 
                                chrom == "scaffold2581_cov101" | chrom == "scaffold1396_cov105" |
                                chrom == "scaffold150315_cov104" | chrom == "scaffold2081_cov104"| 
                                chrom == "scaffold3008_cov108" | chrom == "scaffold222_cov106")
Pact.9.1$chrom <- gsub("_", "", Pact.9.1$chrom)
Pact.9.1$avg <- Pact.9.1$per.meth/100
Pact.9.1 <- Pact.9.1[,c(1,2,3,7)]



str(Pact.1.1)
str(Pact.2.1)
str(Pact.cytoband.df)

# #scale the data from 0-1
# Pact_union.WGBS <- Pact_union %>% mutate(avg = rowMeans(.[4:6], na.rm=TRUE))
# Pact_union.WGBS <- Pact_union.WGBS[,c(1,2,3,13)]
# Pact_union.WGBS$avg <- Pact_union.WGBS$avg/100
# Pact_union.RRBS <- Pact_union %>% mutate(avg = rowMeans(.[7:9], na.rm=TRUE))
# Pact_union.RRBS <- Pact_union.RRBS[,c(1,2,3,13)]
# Pact_union.RRBS$avg <- Pact_union.RRBS$avg/100
# Pact_union.MBDBS <- Pact_union %>% mutate(avg = rowMeans(.[10:12], na.rm=TRUE))
# Pact_union.MBDBS <- Pact_union.MBDBS[,c(1,2,3,13)]
# Pact_union.MBDBS$avg <- Pact_union.MBDBS$avg/100


#plot top 10 scaffold with the highest gene numbers
pdf('Output/Pact_genes_CpG_test.pdf')
circos.clear()

circos.initializeWithIdeogram(Pact.cytoband.df, species = NULL, sort.chr = TRUE)
circos.genomicTrack(Pact.1.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="blue")
                    })
circos.genomicTrack(Pact.2.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="blue")
                    })
circos.genomicTrack(Pact.3.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="blue")
                    })
circos.genomicTrack(Pact.4.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="green")
                    })
circos.genomicTrack(Pact.5.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="green")
                    })
circos.genomicTrack(Pact.6.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="green")
                    })
circos.genomicTrack(Pact.7.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="cyan")
                    })
circos.genomicTrack(Pact.8.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="cyan")
                    })
circos.genomicTrack(Pact.9.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="cyan")
                    })
dev.off()











##### MONTIPORA #####
#load in all genes gff from a genome
Mcap.cytoband <- read.table("data/Mcap.GFFannotation.gene.gff", sep = "\t",  
                            header = FALSE)

#remove extra coloumns from gene gff
Mcap.cytoband <- Mcap.cytoband[,c(1,4,5,6,9)]
str(Mcap.cytoband)
head(Mcap.cytoband)

#examine scaffolds for the number of genes in a scaffold 
sizes <- Mcap.cytoband.df %>% 
  group_by(V1) %>% 
  summarise(count=n()) %>%
  arrange(-count)

Mcap.cytoband.df <- Mcap.cytoband %>% filter(V1 == 1 | V1 == 2 | V1 == 3 | V1 == 4| V1 == 5 |
                                               V1 == 6 | V1 == 7 | V1 == 8 | V1 == 9| V1 == 10) 

#read in the 5x coverage of CpG 
Mcap.10 <- read.table("data/Mcap_tab/Meth10_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mcap.10) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Mcap.10)
str(Mcap.10)

Mcap.11 <- read.table("data/Mcap_tab/Meth11_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mcap.11) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Mcap.11)
str(Mcap.11)

Mcap.12 <- read.table("data/Mcap_tab/Meth12_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mcap.12) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Mcap.12)
str(Mcap.12)

Mcap.13 <- read.table("data/Mcap_tab/Meth13_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mcap.13) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Mcap.13)
str(Mcap.13)

Mcap.14 <- read.table("data/Mcap_tab/Meth14_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mcap.14) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Mcap.14)
str(Mcap.14)

Mcap.15 <- read.table("data/Mcap_tab/Meth15_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mcap.15) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Mcap.15)
str(Mcap.15)

Mcap.16 <- read.table("data/Mcap_tab/Meth16_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mcap.16) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Mcap.16)
str(Mcap.16)

Mcap.17 <- read.table("data/Mcap_tab/Meth17_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mcap.17) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Mcap.17)
str(Mcap.17)

Mcap.18 <- read.table("data/Mcap_tab/Meth18_R1_001_val_1_bismark_bt2_pe._5x.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mcap.18) <- c("chrom", "start", "stop", "per.meth", "meth", "unmeth")
head(Mcap.18)
str(Mcap.18)

Mcap.10.1 <- Mcap.10 %>% filter(chrom == 1 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10)
Mcap.10.1$chrom <- gsub("_", "", Mcap.10.1$chrom)
Mcap.10.1$avg <- Mcap.10.1$per.meth/100
Mcap.10.1 <- Mcap.10.1[,c(1,2,3,7)]

Mcap.11.1 <- Mcap.11 %>% filter(chrom == 1 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10)
Mcap.11.1$chrom <- gsub("_", "", Mcap.11.1$chrom)
Mcap.11.1$avg <- Mcap.11.1$per.meth/100
Mcap.11.1 <- Mcap.11.1[,c(1,2,3,7)]

Mcap.12.1 <- Mcap.12 %>% filter(chrom == 1 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10)
Mcap.12.1$chrom <- gsub("_", "", Mcap.12.1$chrom)
Mcap.12.1$avg <- Mcap.12.1$per.meth/100
Mcap.12.1 <- Mcap.12.1[,c(1,2,3,7)]

Mcap.13.1 <- Mcap.13 %>% filter(chrom == 1 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10)
Mcap.13.1$chrom <- gsub("_", "", Mcap.13.1$chrom)
Mcap.13.1$avg <- Mcap.13.1$per.meth/100
Mcap.13.1 <- Mcap.13.1[,c(1,2,3,7)]

Mcap.14.1 <- Mcap.14 %>% filter(chrom == 1 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10)
Mcap.14.1$chrom <- gsub("_", "", Mcap.14.1$chrom)
Mcap.14.1$avg <- Mcap.14.1$per.meth/100
Mcap.14.1 <- Mcap.14.1[,c(1,2,3,7)]

Mcap.15.1 <- Mcap.15 %>% filter(chrom == 1 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10)
Mcap.15.1$chrom <- gsub("_", "", Mcap.15.1$chrom)
Mcap.15.1$avg <- Mcap.15.1$per.meth/100
Mcap.15.1 <- Mcap.15.1[,c(1,2,3,7)]

Mcap.16.1 <- Mcap.16 %>% filter(chrom == 1 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10)
Mcap.16.1$chrom <- gsub("_", "", Mcap.16.1$chrom)
Mcap.16.1$avg <- Mcap.16.1$per.meth/100
Mcap.16.1 <- Mcap.16.1[,c(1,2,3,7)]

Mcap.17.1 <- Mcap.17 %>% filter(chrom == 1 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10)
Mcap.17.1$chrom <- gsub("_", "", Mcap.17.1$chrom)
Mcap.17.1$avg <- Mcap.17.1$per.meth/100
Mcap.17.1 <- Mcap.17.1[,c(1,2,3,7)]

Mcap.18.1 <- Mcap.18 %>% filter(chrom == 1 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10)
Mcap.18.1$chrom <- gsub("_", "", Mcap.18.1$chrom)
Mcap.18.1$avg <- Mcap.18.1$per.meth/100
Mcap.18.1 <- Mcap.18.1[,c(1,2,3,7)]


#plot first 10 scaffolds
pdf('Output/Mcap_genes_CpG_test.pdf')
circos.clear()

circos.initializeWithIdeogram(Mcap.cytoband.df, species = NULL, sort.chr = TRUE)
circos.genomicTrack(Mcap.10.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="blue")
                    })
circos.genomicTrack(Mcap.11.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="blue")
                    })
circos.genomicTrack(Mcap.12.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="blue")
                    })
circos.genomicTrack(Mcap.13.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="green")
                    })
circos.genomicTrack(Mcap.14.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="green")
                    })
circos.genomicTrack(Mcap.15.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="green")
                    })
circos.genomicTrack(Mcap.16.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="cyan")
                    })
circos.genomicTrack(Mcap.17.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="cyan")
                    })
circos.genomicTrack(Mcap.18.1, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="cyan")
                    })
dev.off()



##### UNIONS #####
#Pocillopora#

#Montipora#





