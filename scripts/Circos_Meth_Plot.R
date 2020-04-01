library(circlize) 
library(tidyverse)
library(ComplexHeatmap)

#https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html#customize-chromosome-track
#The input data for circos.genomicInitialize() is also a data frame with at least three columns. The first column is genomic category (for cytoband data, it is chromosome name), and the next two columns are positions in each genomic category.

cytoband <- read.table("data/Pact_CpG.gff", sep = "\t", skip=3, header = FALSE,colClasses = c("character","character","character", "numeric",
                                                                                                 "numeric",  "numeric", "character", "character", "character"), )

Pact <- read.table("data/Pact_union_5x.bedgraph", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#, colClasses = c("character","numeric","numeric",  "numeric", "numeric",  "numeric", "numeric","numeric",  "numeric", "numeric","numeric",  "numeric"),
head(Pact)
str(Pact)

cytoband.df <- cytoband[,c(1,4,5,6,9)]
str(cytoband.df)
head(unique(cytoband.df$V1), 20)

cytoband.df <- cytoband.df %>% 
  filter(V1 == "scaffold1_cov55"| V1 == "scaffold2_cov51" |V1 == "scaffold3_cov83" | V1 == "scaffold4_cov57" | V1 == "scaffold6_cov64" | V1 == "scaffold7_cov100" | V1 == "scaffold8_cov45"
         | V1 == "scaffold9_cov118" |V1 == "scaffold10_cov103"| V1 == "scaffold11_cov60" | V1 == "scaffold12_cov67"| V1 == "scaffold13_cov99"
         |V1 == "scaffold14_cov75"| V1 == "scaffold15_cov110"| V1 == "scaffold16_cov58"| V1 == "scaffold17_cov151"  | V1 == "scaffold18_cov131"| 
           V1 == "scaffold19_cov103"|V1 == "scaffold20_cov103"| V1 == "scaffold21_cov102")
 
Pact_union <- Pact %>% 
  filter(chrom == "scaffold1_cov55"| chrom == "scaffold2_cov51" |chrom == "scaffold3_cov83" | chrom == "scaffold4_cov57" | chrom == "scaffold6_cov64" | chrom == "scaffold7_cov100" | chrom == "scaffold8_cov45"
         | chrom == "scaffold9_cov118" |chrom == "scaffold10_cov103"| chrom == "scaffold11_cov60" | chrom == "scaffold12_cov67"| chrom == "scaffold13_cov99"
         |chrom == "scaffold14_cov75"| chrom == "scaffold15_cov110"| chrom == "scaffold16_cov58"| chrom == "scaffold17_cov151"  | chrom == "scaffold18_cov131"| 
           chrom == "scaffold19_cov103"|chrom == "scaffold20_cov103"| chrom == "scaffold21_cov102")

Pact_union$X1 <-as.numeric(Pact_union$X1)
Pact_union$X2 <-as.numeric(Pact_union$X2)
Pact_union$X3 <-as.numeric(Pact_union$X3)
Pact_union$X4 <-as.numeric(Pact_union$X4)
Pact_union$X5 <-as.numeric(Pact_union$X5)
Pact_union$X6 <-as.numeric(Pact_union$X6)
Pact_union$X7 <-as.numeric(Pact_union$X7)
Pact_union$X8 <-as.numeric(Pact_union$X8)
Pact_union$X9 <-as.numeric(Pact_union$X9)

# Pact.all <- Pact_union
# 
# Pact.all$X1[Pact.all$X1 >= 0] <- 1
# Pact.all$X2[Pact.all$X2 >= 0] <- 1
# Pact.all$X3[Pact.all$X3 >= 0] <- 1
# Pact.all$X4[Pact.all$X4 >= 0] <- 1
# Pact.all$X5[Pact.all$X5 >= 0] <- 1
# Pact.all$X6[Pact.all$X6 >= 0] <- 1
# Pact.all$X7[Pact.all$X7 >= 0] <- 1
# Pact.all$X8[Pact.all$X8 >= 0] <- 1
# Pact.all$X9[Pact.all$X9 >= 0] <- 1
# Pact.all <- Pact.all %>% mutate_at(c(4:12), ~replace(., is.na(.), 0))

#Pact_union$meth <- Pact_union %>% mutate(Total = select(.,X1:X9) %>% rowSums())

Pact_union.WGBS <- Pact_union %>% mutate(avg = rowMeans(.[4:6], na.rm=TRUE))
Pact_union.WGBS <- Pact_union.WGBS[,c(1,2,3,13)]
Pact_union.WGBS$avg <- Pact_union.WGBS$avg/100
Pact_union.RRBS <- Pact_union %>% mutate(avg = rowMeans(.[7:9], na.rm=TRUE))
Pact_union.RRBS <- Pact_union.RRBS[,c(1,2,3,13)]
Pact_union.RRBS$avg <- Pact_union.RRBS$avg/100
Pact_union.MBDBS <- Pact_union %>% mutate(avg = rowMeans(.[10:12], na.rm=TRUE))
Pact_union.MBDBS <- Pact_union.MBDBS[,c(1,2,3,13)]
Pact_union.MBDBS$avg <- Pact_union.MBDBS$avg/100

#Pact_union <- Pact_union[,c(1,2,3,13)]

# circos.clear()
# circos.initializeWithIdeogram(cytoband.df, species = NULL, sort.chr = TRUE)
# 
# #bed_list <- list(Pact_union.WGBS, Pact_union.RRBS,Pact_union.MBDBS )
# circos.genomicTrack(Pact_union.WGBS, stack = TRUE, 
#                     panel.fun = function(region, value, ...) {
#                       i = getI(...)
#                       circos.genomicLines(region, value, type = "h")
#                     })
# circos.genomicTrack(Pact_union.RRBS, stack = TRUE, 
#                     panel.fun = function(region, value, ...) {
#                       i = getI(...)
#                       circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = "red", ...)
#                     })
# circos.genomicTrack(Pact_union.MBDBS, stack = TRUE, 
#                     panel.fun = function(region, value, ...) {
#                       i = getI(...)
#                       circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = "blue", ...)
#                     })
# circos.clear()
# 
# 
# circos.initializeWithIdeogram(cytoband.df, species = NULL, sort.chr = TRUE)
# bed <- Pact.all
# circos.genomicTrack(bed, stack=FALSE,
#                     panel.fun = function(region, value, ...) {
#                       circos.genomicLines(region, value, col = 4:12, ...)
#                     })
circos.clear()

circos.initializeWithIdeogram(cytoband.df, species = NULL, sort.chr = TRUE)
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

# lgd_points <-  Legend(at = c("WGBS", "MBDBS", "RRBS"), type = "points", legend_gp = gpar(col = 1:2), 
#                     title_position = "topleft", title = "Method")
# 
# grid.draw(lgd_points)
# upViewport()


