
library(methylKit)
library(dplyr)
library(ggplot2)
setwd("~/Documents/Pact_Mcap_mergedcov/")

###Use methylKit package to analyze methylation across WGBS, RRBS and MBD-BS libraries###

#P.acuta
####P.acuta all methods#####
#make a list of cov files
Pact_all_TG.list=list('Meth1_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth2_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth3_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth7_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth8_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth9_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth4_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth5_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth6_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov')

#The following command reads coverage files - minimum 5x coverage
myobj_Pact_all_TG=methRead(Pact_all_TG.list, sample.id = list("WGBS 1","WGBS 2","WGBS 3","MBD 1","MBD 2","MBD 3","RRBS 1","RRBS 2","RRBS 3"), assembly = "Pact_genome", treatment = c(0,0,0,0,0,0,1,1,1), context = "CpG", pipeline = "bismarkCoverage", mincov=5)

#Identify CG covered across all samples and generate a correlation plot
#unite#
meth_Pact_all_TG<-unite(myobj_Pact_all_TG)
nrow(meth_Pact_all_TG)
#correlation plot#
jpeg("MethCompare_correlation_Pact_all.jpg", width = 1000, height = 600)
getCorrelation(meth_Pact_all_TG,plot=TRUE)
dev.off()

#P. acuta - pairwise method comparisons to assess quantiative methylation calls between methods
#WGBS v MBD####################################################################################
Pact_WvM_TG.list=list('Meth1_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth2_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth3_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth7_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth8_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth9_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov')

#The following command reads coverage files - minimum 5x coverage
myobj_Pact_WvM_TG=methRead(Pact_WvM_TG.list, sample.id = list("Pact_W_1","Pact_W_2","Pact_W_3","Pact_M_7","Pact_M_8","Pact_M_9"), 
                           assembly = "Pact_genome", treatment = c(0,0,0,1,1,1), context = "CpG", pipeline = "bismarkCoverage", mincov = 5)


###unite - perform differential methylation analysis with various cutoffs
#unite
meth_Pact_WvM_TG<-unite(myobj_Pact_WvM_TG)
nrow(meth_Pact_WvM_TG)
#calculate differential methylation
myDiff_Pact_WvM_TG<-calculateDiffMeth(meth_Pact_WvM_TG,num.cores = 8,weighted.mean=FALSE,slim=TRUE)

#get significant difference at 50% and qvalue <0.01
myDiff_Pact_WvM_TG_50p <- getMethylDiff(myDiff_Pact_WvM_TG, difference = 50, qvalue = 0.01)
nrow(myDiff_Pact_WvM_TG_50p)
length(which(myDiff_Pact_WvM_TG_50p$meth.diff > 0))

#RRBS v MBD####################################################################################

Pact_RvM_TG.list=list( 'Meth4_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth5_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth6_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth7_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth8_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth9_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov')

#The following command reads coverage files
myobj_Pact_RvM_TG=methRead(Pact_RvM_TG.list, sample.id = list("Pact_R_4","Pact_R_5","Pact_R_6","Pact_M_7","Pact_M_8","Pact_M_9"),
                           assembly = "Pact_genome", treatment = c(0,0,0,1,1,1), context = "CpG", pipeline = "bismarkCoverage", mincov = 5)

##unite - perform differential methylation analysis with various cutoffs
#unite
meth_Pact_RvM_TG<-unite(myobj_Pact_RvM_TG)
nrow(meth_Pact_RvM_TG)
#calculate differential methylation
myDiff_Pact_RvM_TG<-calculateDiffMeth(meth_Pact_RvM_TG,num.cores = 8,weighted.mean=FALSE,slim=TRUE)

#get significant difference at 50% and qvalue <0.01
myDiff_Pact_RvM_TG_50p <- getMethylDiff(myDiff_Pact_RvM_TG, difference = 50, qvalue = 0.01)
nrow(myDiff_Pact_RvM_TG_50p)
#how many are hypermethylated in MBDBS?
length(which(myDiff_Pact_RvM_TG_50p$meth.diff > 0))
#WGBS v RRBS####################################################################################

Pact_WvR_TG.list=list( 'Meth1_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth2_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth3_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth4_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth5_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth6_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov')

#The following command reads coverage files
myobj_Pact_WvR_TG=methRead(Pact_WvR_TG.list, sample.id = list("Pact_W_1","Pact_W_2","Pact_W_3","Pact_R_4","Pact_R_5","Pact_R_6"),
                           assembly = "Pact_genome", treatment = c(0,0,0,1,1,1), context = "CpG", pipeline = "bismarkCoverage", mincov = 5)

##unite - perform differential methylation analysis with various cutoffs
#unite
meth_Pact_WvR_TG<-unite(myobj_Pact_WvR_TG)
nrow(meth_Pact_WvR_TG)
#calculate differential methylation
myDiff_Pact_WvR_TG<-calculateDiffMeth(meth_Pact_WvR_TG,num.cores = 8,weighted.mean=FALSE,slim=TRUE)

#get significant difference at 50% and qvalue <0.01
myDiff_Pact_WvR_TG_50p <- getMethylDiff(myDiff_Pact_WvR_TG, difference = 50, qvalue = 0.01)
nrow(myDiff_Pact_WvR_TG_50p)
length(which(myDiff_Pact_WvR_TG_50p$meth.diff > 0))




#M. capitata ####
####M. capitata all methods#####
#make a list of cov files
Mcap_all_TG.list=list('Meth10_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth11_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth12_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth16_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth17_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth18_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth13_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth14_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth15_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov')

#The following command reads coverage files - minimum 5x coverage
myobj_Mcap_all_TG=methRead(Mcap_all_TG.list, sample.id = list("WGBS 1","WGBS 2","WGBS 3","MBD 1","MBD 2","MBD 3","RRBS 1","RRBS 2","RRBS 3"), assembly = "Mcap_genome", treatment = c(0,0,0,0,0,0,1,1,1), context = "CpG", pipeline = "bismarkCoverage", mincov = 5)

#Identify CG covered across all samples and generate a correlation plot
#unite#
meth_Mcap_all_TG<-unite(myobj_Mcap_all_TG)
nrow(meth_Mcap_all_TG)
#correlation plot#
jpeg("MethCompare_correlation_Mcap_all.jpg", width = 1000, height = 600)
getCorrelation(meth_Mcap_all_TG,plot=TRUE)
dev.off()

#M. capitata - pairwise method comparisons to assess quantiative methylation calls between methods
#WGBS v MBD####################################################################################
Mcap_WvM_TG.list=list('Meth10_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth11_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth12_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth16_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth17_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth18_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov')

#The following command reads coverage files
myobj_Mcap_WvM_TG=methRead(Mcap_WvM_TG.list, sample.id = list("Mcap_W_10","Mcap_W_11","Mcap_W_12","Mcap_M_16","Mcap_M_17","Mcap_M_18"), 
                           assembly = "Mcap_genome", treatment = c(0,0,0,1,1,1), context = "CpG", pipeline = "bismarkCoverage", mincov = 5)


###unite - perform differential methylation analysis with various cutoffs

#unite
meth_Mcap_WvM_TG<-unite(myobj_Mcap_WvM_TG)
nrow(meth_Mcap_WvM_TG)
#calculate differential methylation
myDiff_Mcap_WvM_TG<-calculateDiffMeth(meth_Mcap_WvM_TG,num.cores = 8,weighted.mean=FALSE,slim=TRUE)
#get significant difference at 50% and qvalue <0.01
myDiff_Mcap_WvM_TG_50p <- getMethylDiff(myDiff_Mcap_WvM_TG, difference = 50, qvalue = 0.01)
nrow(myDiff_Mcap_WvM_TG_50p)
length(which(myDiff_Mcap_WvM_TG_50p$meth.diff > 0))

#RRBS v MBD####################################################################################

Mcap_RvM_TG.list=list('Meth13_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth14_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth15_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth16_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth17_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                      'Meth18_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov')

#The following command reads coverage files
myobj_Mcap_RvM_TG=methRead(Mcap_RvM_TG.list, sample.id = list("Mcap_R_13","Mcap_R_14","Mcap_R_15","Mcap_M_16","Mcap_M_17","Mcap_M_18"),
                           assembly = "Mcap_genome", treatment = c(0,0,0,1,1,1), context = "CpG", pipeline = "bismarkCoverage", mincov = 5)

##unite - perform differential methylation analysis with various cutoffs
#unite
meth_Mcap_RvM_TG<-unite(myobj_Mcap_RvM_TG)
nrow(meth_Mcap_RvM_TG)
#calculate differential methylation
myDiff_Mcap_RvM_TG<-calculateDiffMeth(meth_Mcap_RvM_TG,num.cores = 8,weighted.mean=FALSE,slim=TRUE)

#get significant difference at 50% and qvalue <0.01
myDiff_Mcap_RvM_TG_50p <- getMethylDiff(myDiff_Mcap_RvM_TG, difference = 50, qvalue = 0.01)
nrow(myDiff_Mcap_RvM_TG_50p)
length(which(myDiff_Mcap_RvM_TG_50p$meth.diff > 0))

#WGBS v RRBS####################################################################################

Mcap_WvR_TG.list=list( 'Meth10_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth11_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth12_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth13_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth14_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov',
                       'Meth15_R1_001_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov')

#The following command reads coverage files
myobj_Mcap_WvR_TG=methRead(Mcap_WvR_TG.list, sample.id = list("Mcap_W_10","Mcap_W_11","Mcap_W_12","Mcap_R_13","Mcap_R_14","Mcap_R_15"),
                           assembly = "Mcap_genome", treatment = c(0,0,0,1,1,1), context = "CpG", pipeline = "bismarkCoverage", mincov = 5)

##unite - performing differential methylation analysis with various cutoffs
#unite
meth_Mcap_WvR_TG<-unite(myobj_Mcap_WvR_TG)
nrow(meth_Mcap_WvR_TG)
#calculate differential methylation
myDiff_Mcap_WvR_TG<-calculateDiffMeth(meth_Mcap_WvR_TG,num.cores = 8,weighted.mean=FALSE,slim=TRUE)

#get significant difference at 50% and qvalue <0.01
myDiff_Mcap_WvR_TG_50p <- getMethylDiff(myDiff_Mcap_WvR_TG, difference = 50, qvalue = 0.01)
nrow(myDiff_Mcap_WvR_TG_50p)
length(which(myDiff_Mcap_WvR_TG_50p$meth.diff > 0))
