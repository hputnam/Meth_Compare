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
Pact.cytoband.df <- Pact.cytoband[,c(1,4,5,6,9)]
str(Pact.cytoband.df)
head(Pact.cytoband.df)

#read in the 5x coverage bedgraph of CpG for all 18 samples, 3 reps x 3 methods x 2 species
Pact <- read.table("data/Pact_union_5x.bedgraph", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(Pact)
str(Pact)

#remove extra coloumns from gene gff
Pact.cytoband.df <- Pact.cytoband[,c(1,4,5,6,9)]
str(Pact.cytoband.df)
head(Pact.cytoband.df)

#examine scaffolds for the number of genes in a scaffold 
sizes <- Pact.cytoband.df %>% 
  group_by(V1) %>% 
  summarise(count=n()) %>%
  arrange(-count)

which(Pact.cytoband.df$V1=='scaffold502cov107')

Pact.cytoband.df <- Pact.cytoband.df %>%
  filter(V1 == "scaffold502cov107"| V1 == "scaffold3517cov108" | V1 == "scaffold179cov105"| V1 == "scaffold382cov105" | V1 == "scaffold2581cov101"
         | V1 == "scaffold1396cov105" | V1 == "scaffold150315cov104" | V1 == "scaffold2081cov104"| V1 == "scaffold3008cov108" | V1 == "scaffold222cov106" )

Pact_union <- Pact %>%
  filter(chrom== "scaffold502_cov107"| chrom== "scaffold3517_cov108" | chrom== "scaffold179_cov105"| chrom== "scaffold382_cov105" | chrom== "scaffold2581_cov101"
         | chrom== "scaffold1396_cov105" | chrom== "scaffold150315_cov104" | chrom== "scaffold2081_cov104"| chrom== "scaffold3008_cov108" | chrom== "scaffold222_cov106")

Pact_union$chrom <- gsub("_", "", Pact_union$chrom)

which(Pact_union$chrom=='scaffold502cov107')

#set all CpG data as numeric
Pact_union$X1 <-as.numeric(Pact_union$X1)
Pact_union$X2 <-as.numeric(Pact_union$X2)
Pact_union$X3 <-as.numeric(Pact_union$X3)
Pact_union$X4 <-as.numeric(Pact_union$X4)
Pact_union$X5 <-as.numeric(Pact_union$X5)
Pact_union$X6 <-as.numeric(Pact_union$X6)
Pact_union$X7 <-as.numeric(Pact_union$X7)
Pact_union$X8 <-as.numeric(Pact_union$X8)
Pact_union$X9 <-as.numeric(Pact_union$X9)

#scale the data from 0-1
Pact_union.WGBS <- Pact_union %>% mutate(avg = rowMeans(.[4:6], na.rm=TRUE))
Pact_union.WGBS <- Pact_union.WGBS[,c(1,2,3,13)]
Pact_union.WGBS$avg <- Pact_union.WGBS$avg/100
Pact_union.RRBS <- Pact_union %>% mutate(avg = rowMeans(.[7:9], na.rm=TRUE))
Pact_union.RRBS <- Pact_union.RRBS[,c(1,2,3,13)]
Pact_union.RRBS$avg <- Pact_union.RRBS$avg/100
Pact_union.MBDBS <- Pact_union %>% mutate(avg = rowMeans(.[10:12], na.rm=TRUE))
Pact_union.MBDBS <- Pact_union.MBDBS[,c(1,2,3,13)]
Pact_union.MBDBS$avg <- Pact_union.MBDBS$avg/100

#plot top 10 scaffold with the highest gene numbers
pdf('Output/Pact_genes_CpG_3methods_top10.pdf')
circos.clear()

circos.initializeWithIdeogram(Pact.cytoband.df, species = NULL, sort.chr = TRUE)
circos.genomicTrack(Pact_union.WGBS, stack=FALSE, ylim=c(0,1), track.height = 0.1,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="blue")
                    })
circos.genomicTrack(Pact_union.MBDBS, stack=FALSE, ylim=c(0,1), track.height = 0.1,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="cyan")
                    })
circos.genomicTrack(Pact_union.RRBS, stack=FALSE, ylim=c(0,1),track.height = 0.1,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="green")
                    })
dev.off()

##### MONTIPORA #####
#load in all genes gff from a genome
Mcap.cytoband <- read.table("data/Mcap.GFFannotation.gene.gff", sep = "\t",  header = FALSE)

Mcap.cytoband <- Mcap.cytoband %>%
  filter(V3 == 'gene')

#remove extra coloumns from gene gff
Mcap.cytoband.df <- Mcap.cytoband[,c(1,4,5,6,9)]
str(Mcap.cytoband.df)
head(Mcap.cytoband.df)

#read in the 5x coverage bedgraph of CpG for all 18 samples, 3 reps x 3 methods x 2 species
Mcap <- read.table("data/Mcap_union_5x.bedgraph", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(Mcap)
str(Mcap)

#examine scaffolds for the number of genes in a scaffold 
sizes <- Mcap.cytoband.df %>% 
  group_by(V1) %>% 
  summarise(count=n()) %>%
  arrange(-count)

which(Mcap.cytoband.df$V1=='2')

Mcap.cytoband.df <- Mcap.cytoband.df %>%
  filter(V1 == "2"| V1 == "19" | V1 == "28"| V1 == "10" | V1 == "1"
         | V1 == "12" | V1 == "105" | V1 == "68"| V1 == "74" | V1 == "26" )

Mcap_union <- Mcap %>%
  filter(chrom== "2"| chrom== "19" | chrom== "28"| chrom== "10" | chrom== "1"
         | chrom== "12" | chrom== "105" | chrom== "68"| chrom== "74" | chrom== "26" )

which(Mcap_union$chrom=='2')

#set all CpG data as numeric
Mcap_union$X10 <-as.numeric(Mcap_union$X10)
Mcap_union$X11 <-as.numeric(Mcap_union$X11)
Mcap_union$X12 <-as.numeric(Mcap_union$X12)
Mcap_union$X13 <-as.numeric(Mcap_union$X13)
Mcap_union$X14 <-as.numeric(Mcap_union$X14)
Mcap_union$X15 <-as.numeric(Mcap_union$X15)
Mcap_union$X16 <-as.numeric(Mcap_union$X16)
Mcap_union$X17 <-as.numeric(Mcap_union$X17)
Mcap_union$X18 <-as.numeric(Mcap_union$X18)


#scale the data from 0-1
Mcap_union.WGBS <- Mcap_union %>% mutate(avg = rowMeans(.[4:6], na.rm=TRUE))
Mcap_union.WGBS <- Mcap_union.WGBS[,c(1,2,3,13)]
Mcap_union.WGBS$avg <- Mcap_union.WGBS$avg/100
Mcap_union.RRBS <- Mcap_union %>% mutate(avg = rowMeans(.[7:9], na.rm=TRUE))
Mcap_union.RRBS <- Mcap_union.RRBS[,c(1,2,3,13)]
Mcap_union.RRBS$avg <- Mcap_union.RRBS$avg/100
Mcap_union.MBDBS <- Mcap_union %>% mutate(avg = rowMeans(.[10:12], na.rm=TRUE))
Mcap_union.MBDBS <- Mcap_union.MBDBS[,c(1,2,3,13)]
Mcap_union.MBDBS$avg <- Mcap_union.MBDBS$avg/100

#plot CpG of top 10 scaffolds with the highest numbers of genes
#outer track = WGBS in blue
#middle track = RRBS in cyan
#inner track = RRBS in green
pdf('Output/Mcap_genes_CpG_3methods_top10.pdf')
circos.clear()

circos.initializeWithIdeogram(Mcap.cytoband.df, species = NULL, sort.chr = TRUE)
circos.genomicTrack(Mcap_union.WGBS, stack=FALSE, ylim=c(0,1), track.height = 0.1,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="blue")
                    })
circos.genomicTrack(Mcap_union.MBDBS, stack=FALSE, ylim=c(0,1), track.height = 0.1,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="cyan")
                    })
circos.genomicTrack(Mcap_union.RRBS, stack=FALSE, ylim=c(0,1),track.height = 0.1,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="green")
                    })
dev.off()