library(circlize) 
library(tidyverse)
library(ComplexHeatmap)

#https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html#customize-chromosome-track
#The input data for circos.genomicInitialize() is also a data frame with at least three columns. The first column is genomic category (for cytoband data, it is chromosome name), and the next two columns are positions in each genomic category.

cytoband <- read.table("data/Mcap_CpG.gff", sep = "\t", skip=3, header = FALSE,colClasses = c("character","character","character", "numeric",
                                                                                                 "numeric",  "numeric", "character", "character", "character"), )

Pact <- read.table("data/Mcap_union_5x.bedgraph", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#, colClasses = c("character","numeric","numeric",  "numeric", "numeric",  "numeric", "numeric","numeric",  "numeric", "numeric","numeric",  "numeric"),
head(Pact)
str(Pact)

cytoband.df <- cytoband[,c(1,4,5,6,9)]
str(cytoband.df)

cytoband.df <- cytoband.df %>% 
  filter(V1 == 1| V1 == 2|V1 == 3| V1 == 4| V1 == 5| V1 == 6 | V1 == 7| V1 == 8|V1 == 9| V1 == 10| V1 == 11| V1 == 12
         |V1 == 13| V1 == 14| V1 == 15| V1 == 16 | V1 == 17| V1 == 18|V1 == 19| V1 == 20)
 
Pact_union <- Pact %>% 
  filter(chrom == 1| chrom == 2|chrom == 3| chrom == 4| chrom == 5| chrom == 6 | chrom == 7| chrom == 8|chrom == 9| chrom == 10| chrom == 11| chrom == 12
         |chrom == 13| chrom == 14| chrom == 15| chrom == 16 | chrom == 17| chrom == 18|chrom == 19| chrom == 20)

Pact_union$X10 <-as.numeric(Pact_union$X10)
Pact_union$X11 <-as.numeric(Pact_union$X11)
Pact_union$X12 <-as.numeric(Pact_union$X12)
Pact_union$X13 <-as.numeric(Pact_union$X13)
Pact_union$X14 <-as.numeric(Pact_union$X14)
Pact_union$X15 <-as.numeric(Pact_union$X15)
Pact_union$X16 <-as.numeric(Pact_union$X16)
Pact_union$X17 <-as.numeric(Pact_union$X17)
Pact_union$X18 <-as.numeric(Pact_union$X18)

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
circos.genomicTrack(Pact_union.WGBS, stack=FALSE, ylim=c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.1, col="blue")
                    })
circos.genomicTrack(Pact_union.MBDBS, stack=FALSE, ylim=c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.1, col="cyan")
                    })
circos.genomicTrack(Pact_union.RRBS, stack=FALSE, ylim=c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.1, col="green")
                    })

# lgd_points <-  Legend(at = c("WGBS", "MBDBS", "RRBS"), type = "points", legend_gp = gpar(col = 1:2), 
#                     title_position = "topleft", title = "Method")
# 
# grid.draw(lgd_points)
# upViewport()


