sessionInfo()

#  load libraries
#```{r}
library(data.table)
library(dplyr)
library(ggplot2)
#install.packages("ggpubr")
#library(ggpubr)
library(tidyr)
#install.packages("psych")
#library(psych)
#```

#read in union bedgraph for Pact
#```{r}

#Pact_meth <- read.table("Pact_union_Meth.bedgraph", header = T, sep = "\t", stringsAsFactors = F)

#Pact_cov <- read.table("Pact_union_CpG_cov.bedgraph", header = T, sep = "\t", stringsAsFactors = F)

#Mcap_meth <- read.table("Mcap_union_Meth.bedgraph", header = T, sep = "\t", stringsAsFactors = F)

#Mcap_cov <- read.table("Mcap_union_Cov.bedgraph", header = T, sep = "\t", stringsAsFactors = F)

#```

#convert unionBed file to long format
#```{r}
#Pact
#Pact_meth_l <- tidyr::gather(data.frame(Pact_meth),"Sample","perc_meth",4:12)
#Pact_cov_l <- tidyr::gather(data.frame(Pact_cov),"Sample","reads",4:12)

#Mcap
#Mcap_meth_l <- tidyr::gather(data.frame(Mcap_meth),"Sample","perc_meth",4:12)
#Mcap_cov_l <- tidyr::gather(data.frame(Mcap_cov),"Sample","reads",4:12)

#```


#```{r}
#remove lines with NA (no data)
#Pact_cov_l[Pact_cov_l == "N/A"] <- NA
#Pact_meth_l[Pact_meth_l == "N/A"] <- NA

#Pact_cov_l <- data.frame(Pact_cov_l[complete.cases(Pact_cov_l),])
#Pact_meth_l <- data.frame(Pact_meth_l[complete.cases(Pact_meth_l),])

#edit create loci column with scaffold, start, and end as one string
#Pact_cov_l$loci <- paste0(Pact_cov_l$chrom,":", Pact_cov_l$start,"-",Pact_cov_l$end)

#Pact_meth_l$loci <- paste0(Pact_meth_l$chrom,":", Pact_meth_l$start,"-",Pact_meth_l$end)

#remove the scaffold, start, and end columns 
#Pact_cov_l <- Pact_cov_l[,-c(1:3)]
#Pact_meth_l <- Pact_meth_l[,-c(1:3)]

#remove the X in sample names
#Pact_cov_l$Sample <- gsub("X","", Pact_cov_l$Sample)
#Pact_meth_l$Sample <- gsub("X","", Pact_meth_l$Sample)

#convery sample column to numeric
#Pact_cov_l$Sample <- as.numeric(as.character(Pact_cov_l$Sample))
#Pact_meth_l$Sample <- as.numeric(as.character(Pact_meth_l$Sample))

#convert reads and methylation to numeric
#Pact_cov_l$reads <- as.numeric(as.character(Pact_cov_l$reads))
#Pact_meth_l$perc_meth <- as.numeric(as.character(Pact_meth_l$perc_meth))

#add method column
#Pact_cov_l$method <- Pact_cov_l$Sample
#Pact_cov_l$method <- gsub("1|2|3","WGBS",Pact_cov_l$method)
#Pact_cov_l$method <- gsub("4|5|6","RRBS",Pact_cov_l$method)
#Pact_cov_l$method <- gsub("7|8|9","MBD",Pact_cov_l$method)

#Pact_meth_l$method <- Pact_meth_l$Sample
#Pact_meth_l$method <- gsub("1|2|3","WGBS",Pact_meth_l$method)
#Pact_meth_l$method <- gsub("4|5|6","RRBS",Pact_meth_l$method)
#Pact_meth_l$method <- gsub("7|8|9","MBD",Pact_meth_l$method)

#preview objects
#head(Pact_cov_l)
#head(Pact_meth_l)


## Do the same for Mcap ##


#remove lines with NA (no data)
#Mcap_cov_l[Mcap_cov_l == "N/A"] <- NA
#Mcap_meth_l[Mcap_meth_l == "N/A"] <- NA

#Mcap_cov_l <- data.frame(Mcap_cov_l[complete.cases(Mcap_cov_l),])
#Mcap_meth_l <- data.frame(Mcap_meth_l[complete.cases(Mcap_meth_l),])

#edit create loci column with scaffold, start, and end as one string
#Mcap_cov_l$loci <- paste0(Mcap_cov_l$chrom,":", Mcap_cov_l$start,"-",Mcap_cov_l$end)

#Mcap_meth_l$loci <- paste0(Mcap_meth_l$chrom,":", Mcap_meth_l$start,"-",Mcap_meth_l$end)

#remove the scaffold, start, and end columns 
#Mcap_cov_l <- Mcap_cov_l[,-c(1:3)]
#Mcap_meth_l <- Mcap_meth_l[,-c(1:3)]

#remove the X in sample names
#Mcap_cov_l$Sample <- gsub("X","", Mcap_cov_l$Sample)
#Mcap_meth_l$Sample <- gsub("X","", Mcap_meth_l$Sample)

#convery sample column to numeric
#Mcap_cov_l$Sample <- as.numeric(as.character(Mcap_cov_l$Sample))
#Mcap_meth_l$Sample <- as.numeric(as.character(Mcap_meth_l$Sample))

#convert reads and methylation to numeric
#Mcap_cov_l$reads <- as.numeric(as.character(Mcap_cov_l$reads))
#Mcap_meth_l$perc_meth <- as.numeric(as.character(Mcap_meth_l$perc_meth))

#add method column
#Mcap_cov_l$method <- Mcap_cov_l$Sample
#Mcap_cov_l$method <- gsub("1|2|3","WGBS",Mcap_cov_l$method)
#Mcap_cov_l$method <- gsub("4|5|6","RRBS",Mcap_cov_l$method)
#Mcap_cov_l$method <- gsub("7|8|9","MBD",Mcap_cov_l$method)

#Mcap_meth_l$method <- Mcap_meth_l$Sample
#Mcap_meth_l$method <- gsub("1|2|3","WGBS",Mcap_meth_l$method)
#Mcap_meth_l$method <- gsub("4|5|6","RRBS",Mcap_meth_l$method)
#Mcap_meth_l$method <- gsub("7|8|9","MBD",Mcap_meth_l$method)

#preview objects
#head(Mcap_cov_l)
#head(Mcap_meth_l)

#```



#Make column to specify individual
#```{r}

#Pact_cov_perc$individual <- Pact_cov_perc$Sample
#Pact_cov_perc$individual <- gsub("1|4|7","A",Pact_cov_perc$individual)
#Pact_cov_perc$individual <- gsub("2|5|8","B",Pact_cov_perc$individual)
#Pact_cov_perc$individual <- gsub("3|6|9","C",Pact_cov_perc$individual)
#```

# bin read coverage
#```{r}
# # specify interval/bin labels
# meth_tags <- c("[0-10%]","[10-20%]", "[20-30%]", "[30-40%]", "[40-50%]", "[50-60%]", "[60-70%]", "[70-80%]","[80-90%]", "[90-100%]")
# 
# # specify interval/bin labels
# cov_tags <- c("[1-5]","[5-10]","[10-20]", "[20-30]", "[30-40]", "[40-50]", "[50-100]", "[100-500]", "[500-1000]", "[>1000]")
# 
# # bucketing values into bins
# 
# Pact_meth_tag <- as_tibble(Pact_meth_l) %>% mutate(meth_tag = case_when(perc_meth < 10 ~ meth_tags[1],perc_meth >= 10 & perc_meth < 20 ~ meth_tags[2],perc_meth >= 20 & perc_meth < 30 ~ meth_tags[3],perc_meth >= 30 & perc_meth < 40 ~ meth_tags[4],perc_meth >= 40 & perc_meth < 50 ~ meth_tags[5], perc_meth >= 50 & perc_meth < 60 ~ meth_tags[6], perc_meth >= 60 & perc_meth < 70 ~ meth_tags[7], perc_meth >= 70 & perc_meth < 80 ~ meth_tags[8], perc_meth >= 80 & perc_meth < 90 ~ meth_tags[9], perc_meth >= 90 & perc_meth <= 100 ~ meth_tags[10]))
# 
# Pact_cov_tag <- as_tibble(Pact_cov_l) %>% mutate(cov_tag = case_when(reads < 5 ~ cov_tags[1],reads >= 5 & reads < 10 ~ cov_tags[2],reads >= 10 & reads < 20 ~ cov_tags[3],reads >= 20 & reads < 30 ~ cov_tags[4],reads >= 30 & reads < 40 ~ cov_tags[5], reads >= 40 & reads < 50 ~ cov_tags[6], reads >= 50 & reads < 100 ~ cov_tags[7], reads >=100 & reads < 500 ~ cov_tags[8], reads >=500 & reads < 1000 ~ cov_tags[9], reads >=1000 ~ cov_tags[10]))
#```


#combine coverage and methylation data
#```{r}
#Pact_cov_perc <- merge(Pact_cov_l, Pact_meth_l, by = c("Sample", "loci", "method"), all = T)

#head(Pact_cov_perc)

#clean up environment
#rm(Pact_cov)
#rm(Pact_meth)
#rm(Pact_cov_l)
#rm(Pact_meth_l)
#rm(Mcap_cov)
#rm(Mcap_meth)

#Mcap_cov_perc <- merge(Mcap_cov_l, Mcap_meth_l, by = c("Sample", "loci", "method"), all = T)

#head(Mcap_cov_perc)

#clean up environment
#rm(Mcap_cov_l)
#rm(Mcap_meth_l)

#write out files
#write.csv(Mcap_cov_perc, "~/strigg/analyses/20200910/Mcap_cov_perc.txt", sep = "\t", row.names = F, quote = F)

#write.csv(Pact_cov_perc, "~/strigg/analyses/20200910/Pact_cov_perc.txt", sep = "\t", row.names = F, quote = F)

Mcap_cov_perc <- read.csv("~/strigg/analyses/20200910/Mcap_cov_perc.txt", stringsAsFactors = F, header = T)

Pact_cov_perc <- read.csv("~/strigg/analyses/20200910/Pact_cov_perc.txt", stringsAsFactors = F, header = T)

#```

#calculate median methylation and coverage
#```{r}
Pact_meth_stats <- Pact_cov_perc  %>% group_by(loci,method) %>% summarise(medn_meth = median(perc_meth),medn_cov = median(reads),indivd = n())

Mcap_meth_stats <- Mcap_cov_perc  %>% group_by(loci,method) %>% summarise(medn_meth = median(perc_meth),medn_cov = median(reads),indivd = n())

nrow(Pact_meth_stats)
str(Pact_meth_stats)
head(Pact_meth_stats)

nrow(Mcap_meth_stats)
str(Mcap_meth_stats)
head(Mcap_meth_stats)


#```

#add labels for % meth and coverage bins
#```{r}
# specify % meth interval bin labels
meth_tags <- c("[0-10%]","[10-20%]", "[20-30%]", "[30-40%]", "[40-50%]", "[50-60%]", "[60-70%]", "[70-80%]","[80-90%]", "[90-100%]")

# specify coverage bin labels
cov_tags <- c("[1-5]","[5-10]","[10-20]", "[20-30]", "[30-40]", "[40-50]", "[50-100]", "[100-500]", "[500-1000]", "[>1000]")

# bucketing values into bins

Pact_meth_stats_tag <- as_tibble(Pact_meth_stats) %>% mutate(meth_tag = case_when(medn_meth < 10 ~ meth_tags[1],medn_meth >= 10 & medn_meth < 20 ~ meth_tags[2],medn_meth >= 20 & medn_meth < 30 ~ meth_tags[3],medn_meth >= 30 & medn_meth < 40 ~ meth_tags[4],medn_meth >= 40 & medn_meth < 50 ~ meth_tags[5], medn_meth >= 50 & medn_meth < 60 ~ meth_tags[6], medn_meth >= 60 & medn_meth < 70 ~ meth_tags[7], medn_meth >= 70 & medn_meth < 80 ~ meth_tags[8], medn_meth >= 80 & medn_meth < 90 ~ meth_tags[9], medn_meth >= 90 & medn_meth <= 100 ~ meth_tags[10]), cov_tag = case_when(medn_cov < 5 ~ cov_tags[1],medn_cov >= 5 & medn_cov < 10 ~ cov_tags[2],medn_cov >= 10 & medn_cov < 20 ~ cov_tags[3],medn_cov >= 20 & medn_cov < 30 ~ cov_tags[4],medn_cov >= 30 & medn_cov < 40 ~ cov_tags[5], medn_cov >= 40 & medn_cov < 50 ~ cov_tags[6], medn_cov >= 50 & medn_cov < 100 ~ cov_tags[7], medn_cov >=100 & medn_cov < 500 ~ cov_tags[8], medn_cov >=500 & medn_cov < 1000 ~ cov_tags[9], medn_cov >=1000 ~ cov_tags[10]))


head(Pact_meth_stats_tag)
nrow(Pact_meth_stats_tag)

#Mcap
Mcap_meth_stats_tag <- as_tibble(Mcap_meth_stats) %>% mutate(meth_tag = case_when(medn_meth < 10 ~ meth_tags[1],medn_meth >= 10 & medn_meth < 20 ~ meth_tags[2],medn_meth >= 20 & medn_meth < 30 ~ meth_tags[3],medn_meth >= 30 & medn_meth < 40 ~ meth_tags[4],medn_meth >= 40 & medn_meth < 50 ~ meth_tags[5], medn_meth >= 50 & medn_meth < 60 ~ meth_tags[6], medn_meth >= 60 & medn_meth < 70 ~ meth_tags[7], medn_meth >= 70 & medn_meth < 80 ~ meth_tags[8], medn_meth >= 80 & medn_meth < 90 ~ meth_tags[9], medn_meth >= 90 & medn_meth <= 100 ~ meth_tags[10]), cov_tag = case_when(medn_cov < 5 ~ cov_tags[1],medn_cov >= 5 & medn_cov < 10 ~ cov_tags[2],medn_cov >= 10 & medn_cov < 20 ~ cov_tags[3],medn_cov >= 20 & medn_cov < 30 ~ cov_tags[4],medn_cov >= 30 & medn_cov < 40 ~ cov_tags[5], medn_cov >= 40 & medn_cov < 50 ~ cov_tags[6], medn_cov >= 50 & medn_cov < 100 ~ cov_tags[7], medn_cov >=100 & medn_cov < 500 ~ cov_tags[8], medn_cov >=500 & medn_cov < 1000 ~ cov_tags[9], medn_cov >=1000 ~ cov_tags[10]))


head(Mcap_meth_stats_tag)
nrow(Mcap_meth_stats_tag)
#```


#create df for loci where all three individuals per method show coverage
#```{r}
Pact_meth_stats_tag_3x <- Pact_meth_stats_tag[which(Pact_meth_stats_tag$indivd == 3),]

Mcap_meth_stats_tag_3x <- Mcap_meth_stats_tag[which(Mcap_meth_stats_tag$indivd == 3),]
#```


#add group column
#```{r}
#create function for converting list to string
uniqMethod.FUN <- function(x){paste(sort(unique(unlist(strsplit(x, split =", "), use.names = F))), collapse = ", ")}

#store methods that loci were detected in as a comma separated list
#Pact_cov_perc_group <- Pact_meth_stats_tag %>% group_by(loci) %>% mutate(group = toString(method))

##converting group class (list to string)
#Pact_cov_perc_group$group_simp <- sapply(Pact_cov_perc_group$group, uniqMethod.FUN, USE.NAMES = F)

#do the same for 3x
#store methods that loci were detected in as a comma separated list
Pact_cov_perc_group_3x <- Pact_meth_stats_tag_3x %>% group_by(loci) %>% mutate(group = toString(method))

##converting group class (list to string)
Pact_cov_perc_group_3x$group_simp <- sapply(Pact_cov_perc_group_3x$group, uniqMethod.FUN, USE.NAMES = F)

#check df integrity
#str(Pact_cov_perc_group)
#table(Pact_cov_perc_group$group_simp)
#nrow(Pact_cov_perc_group)

str(Pact_cov_perc_group_3x)
table(Pact_cov_perc_group_3x$group_simp)
nrow(Pact_cov_perc_group_3x)

## Mcap ##

#store methods that loci were detected in as a comma separated list

#do the same for 3x
#store methods that loci were detected in as a comma separated list
Mcap_cov_perc_group_3x <- Mcap_meth_stats_tag_3x %>% group_by(loci) %>% mutate(group = toString(method))

##converting group class (list to string)
Mcap_cov_perc_group_3x$group_simp <- sapply(Mcap_cov_perc_group_3x$group, uniqMethod.FUN, USE.NAMES = F)


str(Mcap_cov_perc_group_3x)
table(Mcap_cov_perc_group_3x$group_simp)
nrow(Mcap_cov_perc_group_3x)
#```

#clean up environment
rm(Pact_cov_perc)
rm(Mcap_cov_perc)
rm(Pact_meth_stats)
rm(Mcap_meth_stats)
rm(Pact_meth_stats_tag)
rm(Mcap_meth_stats_tag)
rm(Pact_meth_stats_tag_3x)
rm(Mcap_meth_stats_tag_3x)

#write out *cov_perc_group_3x data
write.csv(Mcap_cov_perc_group_3x, "~/strigg/analyses/20200910/Mcap_cov_perc_group_3x.csv", row.names = F, quote = F)
write.csv(Pact_cov_perc_group_3x, "~/strigg/analyses/20200910/Pact_cov_perc_group_3x.csv", row.names = F, quote = F)



#summarize data
#```{r}
#collapse loci and just count the number of loci for each methylation and coverage bin
#Pact_cov_perc_STACKED <- data.frame(Pact_cov_perc_group) %>% group_by(group_simp, method,meth_tag, cov_tag) %>% summarise(count=n())


# do the same for 3x loci
#Pact_cov_perc_3x_STACKED <- data.frame(Pact_cov_perc_group_3x) %>% group_by(group_simp, method, meth_tag, cov_tag) %>% summarise(count=n())


#Pact_cov_perc_STACKED %>% group_by(group_simp) %>% summarise(method_total_count = sum(count))

# calculate the percent CpGs for each 
#Pact_cov_perc_STACKED <- Pact_cov_perc_STACKED %>% group_by(group_simp,method) %>% mutate(perc_CpG=count/sum(count))

#Pact_cov_perc_3x_STACKED <- Pact_cov_perc_3x_STACKED %>% group_by(group_simp,method) %>% mutate(perc_CpG=count/sum(count))

#```

#create a df with MBD vs. WGBS data
#```{r}
#remove RRBS data
Pact_cov_perc_group_3x_mbd_rel <- Pact_cov_perc_group_3x[which(Pact_cov_perc_group_3x$method != "RRBS"),]
#subset WGBS groups
Pact_cov_perc_group_3x_mbd_rel <- Pact_cov_perc_group_3x_mbd_rel[grep("WGBS", Pact_cov_perc_group_3x_mbd_rel$group_simp),]


#set cov_tag to zero for CpGs not covered by MBD but covered by WGBS
Pact_cov_perc_group_3x_mbd_rel$mbd_cov_tag <- ifelse(Pact_cov_perc_group_3x_mbd_rel$method == "WGBS" & !grepl("MBD",Pact_cov_perc_group_3x_mbd_rel$group_simp),"[0]", Pact_cov_perc_group_3x_mbd_rel$cov_tag)

#set WGBS cov_tag to NA for CpGs covered by WGBS and by MBD since we're gonna use the MBD cov_tags
Pact_cov_perc_group_3x_mbd_rel$mbd_cov_tag <- ifelse(Pact_cov_perc_group_3x_mbd_rel$method == "WGBS" & grepl("MBD",Pact_cov_perc_group_3x_mbd_rel$group_simp),NA, Pact_cov_perc_group_3x_mbd_rel$mbd_cov_tag)

#Set methylation to NA for CpGs covered by MBD
Pact_cov_perc_group_3x_mbd_rel$wgbs_meth_tag <- ifelse(Pact_cov_perc_group_3x_mbd_rel$method == "WGBS",Pact_cov_perc_group_3x_mbd_rel$meth_tag,NA)


#remove all columns except loci, group, mbd_cov_tag, and wgbs_meth_tag
#convert to long format and filter out rows with NAs 
#convert to wide format so that only mbd_cov_tag and wgbs_meth_tag values remain
Pact_cov_perc_group_3x_mbd_rel_reshape <- Pact_cov_perc_group_3x_mbd_rel[,c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")


table(Pact_cov_perc_group_3x_mbd_rel_reshape$group_simp)

#summarize the data for plotting
Pact_cov_perc_group_3x_mbd_rel_sum <- Pact_cov_perc_group_3x_mbd_rel_reshape %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

jpeg("~/strigg/analyses/20200910/MBDcov_x_WGBSmeth_heatmap_Pact.jpg", width = 6, height = 5, units = "in", res = 300 )
ggplot(data = Pact_cov_perc_group_3x_mbd_rel_sum, aes(x=wgbs_meth_tag, y=factor(mbd_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from MBD data)", fill = "number of\nCpGs (log10)") + scale_fill_continuous(high = "#132B43", low = "#b5dfff")
dev.off()


## Mcap ##


#remove RRBS data
Mcap_cov_perc_group_3x_mbd_rel <- Mcap_cov_perc_group_3x[which(Mcap_cov_perc_group_3x$method != "RRBS"),]
#subset WGBS groups
Mcap_cov_perc_group_3x_mbd_rel <- Mcap_cov_perc_group_3x_mbd_rel[grep("WGBS", Mcap_cov_perc_group_3x_mbd_rel$group_simp),]


#set cov_tag to zero for CpGs not covered by MBD but covered by WGBS
Mcap_cov_perc_group_3x_mbd_rel$mbd_cov_tag <- ifelse(Mcap_cov_perc_group_3x_mbd_rel$method == "WGBS" & !grepl("MBD",Mcap_cov_perc_group_3x_mbd_rel$group_simp),"[0]", Mcap_cov_perc_group_3x_mbd_rel$cov_tag)

#set WGBS cov_tag to NA for CpGs covered by WGBS and by MBD since we're gonna use the MBD cov_tags
Mcap_cov_perc_group_3x_mbd_rel$mbd_cov_tag <- ifelse(Mcap_cov_perc_group_3x_mbd_rel$method == "WGBS" & grepl("MBD",Mcap_cov_perc_group_3x_mbd_rel$group_simp),NA, Mcap_cov_perc_group_3x_mbd_rel$mbd_cov_tag)

#Set methylation to NA for CpGs covered by MBD
Mcap_cov_perc_group_3x_mbd_rel$wgbs_meth_tag <- ifelse(Mcap_cov_perc_group_3x_mbd_rel$method == "WGBS",Mcap_cov_perc_group_3x_mbd_rel$meth_tag,NA)


#remove all columns except loci, group, mbd_cov_tag, and wgbs_meth_tag
#convert to long format and filter out rows with NAs 
#convert to wide format so that only mbd_cov_tag and wgbs_meth_tag values remain
Mcap_cov_perc_group_3x_mbd_rel_reshape <- Mcap_cov_perc_group_3x_mbd_rel[,c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")


table(Mcap_cov_perc_group_3x_mbd_rel_reshape$group_simp)

#summarize the data for plotting
Mcap_cov_perc_group_3x_mbd_rel_sum <- Mcap_cov_perc_group_3x_mbd_rel_reshape %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

jpeg("~/strigg/analyses/20200910/MBDcov_x_WGBSmeth_heatmap_Mcap.jpg", width = 6, height = 5, units = "in", res = 300 )
ggplot(data = Mcap_cov_perc_group_3x_mbd_rel_sum, aes(x=wgbs_meth_tag, y=factor(mbd_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from MBD data)", fill = "number of\nCpGs (log10)") + scale_fill_continuous(high = "#132B43", low = "#b5dfff")
dev.off()
#```

#filter for loci covered by WGBS at 5x, 10x, 50x, 100x
#```{r}
loci_lessthan5x <- Pact_cov_perc_group_3x_mbd_rel[which(Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]"),"loci"]

loci_lessthan10x <- Pact_cov_perc_group_3x_mbd_rel[which(Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[5-10]" ),"loci"]

loci_lessthan20x <- Pact_cov_perc_group_3x_mbd_rel[which(Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[5-10]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[10-20]"),"loci"]

loci_lessthan30x <- Pact_cov_perc_group_3x_mbd_rel[which(Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[5-10]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[10-20]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[20-30]"),"loci"]

loci_lessthan40x <- Pact_cov_perc_group_3x_mbd_rel[which(Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[5-10]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[10-20]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[20-30]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[30-40]"),"loci"]

loci_greaterthan50x <- Pact_cov_perc_group_3x_mbd_rel[which(Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[50-100]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[100-500]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[500-1000]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[>1000]"),"loci"]

loci_greaterthan100x <- Pact_cov_perc_group_3x_mbd_rel[which(Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[100-500]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[500-1000]" | Pact_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Pact_cov_perc_group_3x_mbd_rel$cov_tag == "[>1000]"),"loci"]


Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs5x <- Pact_cov_perc_group_3x_mbd_rel[which(!(Pact_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan5x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs10x <- Pact_cov_perc_group_3x_mbd_rel[which(!(Pact_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan10x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs20x <- Pact_cov_perc_group_3x_mbd_rel[which(!(Pact_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan20x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs30x <- Pact_cov_perc_group_3x_mbd_rel[which(!(Pact_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan30x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs40x <- Pact_cov_perc_group_3x_mbd_rel[which(!(Pact_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan40x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs50x <- Pact_cov_perc_group_3x_mbd_rel[which(Pact_cov_perc_group_3x_mbd_rel$loci %in% loci_greaterthan50x$loci),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs100x <- Pact_cov_perc_group_3x_mbd_rel[which(Pact_cov_perc_group_3x_mbd_rel$loci %in% loci_greaterthan100x$loci),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")
#summarize the data for plotting
Pact_cov_perc_group_3x_mbd_rel_sum_wgbs5x <- Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs5x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_mbd_rel_sum_wgbs10x <- Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs10x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_mbd_rel_sum_wgbs20x <- Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs20x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_mbd_rel_sum_wgbs30x <- Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs30x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_mbd_rel_sum_wgbs40x <- Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs40x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_mbd_rel_sum_wgbs50x <- Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs50x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_mbd_rel_sum_wgbs100x <- Pact_cov_perc_group_3x_mbd_rel_reshape_wgbs100x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())


#create a df with all subsets
Pact_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs <- dplyr::bind_rows(list(wgbs_1x=Pact_cov_perc_group_3x_mbd_rel_sum, wgbs_5x=Pact_cov_perc_group_3x_mbd_rel_sum_wgbs5x, wgbs_10x= Pact_cov_perc_group_3x_mbd_rel_sum_wgbs10x, wgbs_20x= Pact_cov_perc_group_3x_mbd_rel_sum_wgbs20x,wgbs_30x= Pact_cov_perc_group_3x_mbd_rel_sum_wgbs30x,wgbs_40x= Pact_cov_perc_group_3x_mbd_rel_sum_wgbs40x,wgbs_50x= Pact_cov_perc_group_3x_mbd_rel_sum_wgbs50x, wgbs_100x= Pact_cov_perc_group_3x_mbd_rel_sum_wgbs100x), .id = "id")

Pact_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs <- Pact_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs %>% group_by(id, wgbs_meth_tag) %>% mutate(perc_CpG = CpGs/sum(CpGs)*100)

#make wgbs capitalized
Pact_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs$id <- gsub("wgbs", "WGBS", Pact_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs$id)







## Mcap ##


loci_lessthan5x <- Mcap_cov_perc_group_3x_mbd_rel[which(Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]"),"loci"]

loci_lessthan10x <- Mcap_cov_perc_group_3x_mbd_rel[which(Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[5-10]" ),"loci"]

loci_lessthan20x <- Mcap_cov_perc_group_3x_mbd_rel[which(Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[5-10]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[10-20]"),"loci"]

loci_lessthan30x <- Mcap_cov_perc_group_3x_mbd_rel[which(Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[5-10]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[10-20]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[20-30]"),"loci"]

loci_lessthan40x <- Mcap_cov_perc_group_3x_mbd_rel[which(Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[1-5]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[5-10]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[10-20]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[20-30]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[30-40]"),"loci"]

loci_greaterthan50x <- Mcap_cov_perc_group_3x_mbd_rel[which(Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[50-100]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[100-500]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[500-1000]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[>1000]"),"loci"]

loci_greaterthan100x <- Mcap_cov_perc_group_3x_mbd_rel[which(Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[100-500]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[500-1000]" | Mcap_cov_perc_group_3x_mbd_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_mbd_rel$cov_tag == "[>1000]"),"loci"]


Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs5x <- Mcap_cov_perc_group_3x_mbd_rel[which(!(Mcap_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan5x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs10x <- Mcap_cov_perc_group_3x_mbd_rel[which(!(Mcap_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan10x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs20x <- Mcap_cov_perc_group_3x_mbd_rel[which(!(Mcap_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan20x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs30x <- Mcap_cov_perc_group_3x_mbd_rel[which(!(Mcap_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan30x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs40x <- Mcap_cov_perc_group_3x_mbd_rel[which(!(Mcap_cov_perc_group_3x_mbd_rel$loci %in% loci_lessthan40x$loci)),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs50x <- Mcap_cov_perc_group_3x_mbd_rel[which(Mcap_cov_perc_group_3x_mbd_rel$loci %in% loci_greaterthan50x$loci),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs100x <- Mcap_cov_perc_group_3x_mbd_rel[which(Mcap_cov_perc_group_3x_mbd_rel$loci %in% loci_greaterthan100x$loci),c("loci","group_simp", "mbd_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", mbd_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")
#summarize the data for plotting
Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs5x <- Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs5x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs10x <- Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs10x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs20x <- Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs20x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs30x <- Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs30x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs40x <- Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs40x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs50x <- Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs50x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs100x <- Mcap_cov_perc_group_3x_mbd_rel_reshape_wgbs100x %>% group_by(mbd_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())


#create a df with all subsets
Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs <- dplyr::bind_rows(list(wgbs_1x=Mcap_cov_perc_group_3x_mbd_rel_sum, wgbs_5x=Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs5x, wgbs_10x= Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs10x, wgbs_20x= Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs20x,wgbs_30x= Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs30x,wgbs_40x= Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs40x,wgbs_50x= Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs50x, wgbs_100x= Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs100x), .id = "id")

Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs <- Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs %>% group_by(id, wgbs_meth_tag) %>% mutate(perc_CpG = CpGs/sum(CpGs)*100)

#make wgbs capitalized
Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs$id <- gsub("wgbs", "WGBS", Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs$id)

#```


#plot data
#```{r}
jpeg("~/strigg/analyses/20200910/MBDcov_x_WGBSmeth_heatmap_Pact_thresholds.jpg", width = 10, height = 5, units = "in", res = 300 )
ggplot(data = Pact_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs, aes(x=wgbs_meth_tag, y=factor(mbd_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from MBD data)", fill = "number of\nCpGs (log10)") + facet_wrap(~factor(id, levels = c("WGBS_1x","WGBS_5x","WGBS_10x","WGBS_20x","WGBS_30x","WGBS_40x","WGBS_50x","WGBS_100x")), nrow = 2) + scale_fill_continuous(high = "#132B43", low = "#b5dfff")
dev.off()

#Mcap
jpeg("~/strigg/analyses/20200910/MBDcov_x_WGBSmeth_heatmap_Mcap_thresholds.jpg", width = 10, height = 5, units = "in", res = 300 )
ggplot(data = Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs, aes(x=wgbs_meth_tag, y=factor(mbd_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from MBD data)", fill = "number of\nCpGs (log10)") + facet_wrap(~factor(id, levels = c("WGBS_1x","WGBS_5x","WGBS_10x","WGBS_20x","WGBS_30x","WGBS_40x","WGBS_50x","WGBS_100x")), nrow = 2) + scale_fill_continuous(high = "#132B43", low = "#b5dfff")
dev.off()

#d <- ggplot(data = Pact_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs, aes(x=wgbs_meth_tag, y=factor(mbd_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from MBD data)", fill = "number of\nCpGs (log10)", title = "P. acuta") + facet_wrap(~factor(id, levels = c("WGBS_1x","WGBS_5x","WGBS_10x","WGBS_20x","WGBS_30x","WGBS_40x","WGBS_50x","WGBS_100x")), nrow = 2) + scale_fill_continuous(high = "#132B43", low = "#b5dfff") + theme(plot.title = element_text(face = "italic"))


#c <- ggplot(data = Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs, aes(x=wgbs_meth_tag, y=factor(mbd_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from MBD data)", fill = "number of\nCpGs (log10)", title = "M. capitata") + facet_wrap(~factor(id, levels = c("WGBS_1x","WGBS_5x","WGBS_10x","WGBS_20x","WGBS_30x","WGBS_40x","WGBS_50x","WGBS_100x")), nrow = 2) + scale_fill_continuous(high = "#132B43", low = "#b5dfff") + theme(plot.title = element_text(face = "italic"))


#jpeg("~/strigg/analyses/20200910/MBDcov_x_WGBSmeth_heatmap_thresholds.jpg", width = 10, height = 5, units = "in", res = 300 )
#ggarrange(c,d, nrow = 2,labels = "AUTO",common.legend = T)
#dev.off()
#```



#This code generates 3 heatmaps showing %methylation x # of reads for each method to see if any method is biased against detecting highly methylated reads (e.g. WGBS)
#```{r}
#subset loci for those only covered by WGBS data and only keep 

#subset loci for those only covered by WGBS data,
#only keep the loci, cov_tag, and meth_tag columns,
# and total the CpGs for each bin
Pact_cov_perc_group_3x_sum <- Pact_cov_perc_group_3x %>% group_by(method,cov_tag, meth_tag) %>% summarise(CpGs =n())

#convert column to ordered factor
Pact_cov_perc_group_3x_sum$cov_tag <- factor(Pact_cov_perc_group_3x_sum$cov_tag, levels = cov_tags)

#plot the WGBS % meth x coverage
jpeg("~/strigg/analyses/20200910/method_x_cov_x_meth_heatmap_Pact.jpg", width = 10, height = 4, units = "in", res = 300 )
ggplot(data = Pact_cov_perc_group_3x_sum, aes(x = meth_tag, y=cov_tag, fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + facet_wrap(~method)  + labs(x = "% methylation", y = "number of reads", fill = "number of\nCpGs (log10)", title = "P. acuta") + scale_fill_continuous(high = "#132B43", low = "#b5dfff") 
dev.off()



#subset loci for those only covered by WGBS data,
#only keep the loci, cov_tag, and meth_tag columns,
# and total the CpGs for each bin
Mcap_cov_perc_group_3x_sum <- Mcap_cov_perc_group_3x %>% group_by(method,cov_tag, meth_tag) %>% summarise(CpGs =n())

#convert column to ordered factor
Mcap_cov_perc_group_3x_sum$cov_tag <- factor(Mcap_cov_perc_group_3x_sum$cov_tag, levels = cov_tags)

#plot the WGBS % meth x coverage
jpeg("~/strigg/analyses/20200910/method_x_cov_x_meth_heatmap_Mcap.jpg", width = 10, height = 4, units = "in", res = 300 )
ggplot(data = Mcap_cov_perc_group_3x_sum, aes(x = meth_tag, y=cov_tag, fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + facet_wrap(~method)  + labs(x = "% methylation", y = "number of reads", fill = "number of\nCpGs (log10)", title = "M. capitata") + scale_fill_continuous(high = "#132B43", low = "#b5dfff")
dev.off()


#e <- ggplot(data = Mcap_cov_perc_group_3x_sum, aes(x = meth_tag, y=cov_tag, fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + facet_wrap(~method)  + labs(x = "% methylation", y = "number of reads", fill = "number of\nCpGs (log10)", title = "M. capitata") + scale_fill_continuous(high = "#132B43", low = "#b5dfff") + theme(plot.title = element_text(face = "italic"))


#f <- ggplot(data = Pact_cov_perc_group_3x_sum, aes(x = meth_tag, y=cov_tag, fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + facet_wrap(~method)  + labs(x = "% methylation", y = "number of reads", fill = "number of\nCpGs (log10)", title = "P. acuta") + scale_fill_continuous(high = "#132B43", low = "#b5dfff") + theme(plot.title = element_text(face = "italic"))


#jpeg("~/strigg/analyses/20200910/MBDcov_x_WGBSmeth_heatmap_thresholds.jpg", width = 10, height = 5, units = "in", res = 300 )
#ggarrange(e,f, nrow = 2,labels = "AUTO",common.legend = T)
#dev.off()
#```

#RRBS vs WBGS
#create a df with RRBS vs. WGBS data
#```{r}
#remove RRBS data
Pact_cov_perc_group_3x_RRBS_rel <- Pact_cov_perc_group_3x[which(Pact_cov_perc_group_3x$method != "MBD"),]
#subset WGBS groups
Pact_cov_perc_group_3x_RRBS_rel <- Pact_cov_perc_group_3x_RRBS_rel[grep("WGBS", Pact_cov_perc_group_3x_RRBS_rel$group_simp),]


#set cov_tag to zero for CpGs not covered by RRBS but covered by WGBS
Pact_cov_perc_group_3x_RRBS_rel$RRBS_cov_tag <- ifelse(Pact_cov_perc_group_3x_RRBS_rel$method == "WGBS" & !grepl("RRBS",Pact_cov_perc_group_3x_RRBS_rel$group_simp),"[0]", Pact_cov_perc_group_3x_RRBS_rel$cov_tag)

#set WGBS cov_tag to NA for CpGs covered by WGBS and by RRBS since we're gonna use the RRBS cov_tags
Pact_cov_perc_group_3x_RRBS_rel$RRBS_cov_tag <- ifelse(Pact_cov_perc_group_3x_RRBS_rel$method == "WGBS" & grepl("RRBS",Pact_cov_perc_group_3x_RRBS_rel$group_simp),NA, Pact_cov_perc_group_3x_RRBS_rel$RRBS_cov_tag)

#Set methylation to NA for CpGs covered by RRBS
Pact_cov_perc_group_3x_RRBS_rel$wgbs_meth_tag <- ifelse(Pact_cov_perc_group_3x_RRBS_rel$method == "WGBS",Pact_cov_perc_group_3x_RRBS_rel$meth_tag,NA)


#remove all columns except loci, group, RRBS_cov_tag, and wgbs_meth_tag
#convert to long format and filter out rows with NAs 
#convert to wide format so that only RRBS_cov_tag and wgbs_meth_tag values remain
Pact_cov_perc_group_3x_RRBS_rel_reshape <- Pact_cov_perc_group_3x_RRBS_rel[,c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")


table(Pact_cov_perc_group_3x_RRBS_rel_reshape$group_simp)

#summarize the data for plotting
Pact_cov_perc_group_3x_RRBS_rel_sum <- Pact_cov_perc_group_3x_RRBS_rel_reshape %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

jpeg("~/strigg/analyses/20200910/RRBScov_x_WGBSmeth_heatmap_Pact.jpg", width = 6, height = 5, units = "in", res = 300 )
ggplot(data = Pact_cov_perc_group_3x_RRBS_rel_sum, aes(x=wgbs_meth_tag, y=factor(RRBS_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from RRBS data)", fill = "number of\nCpGs (log10)") + scale_fill_continuous(high = "#132B43", low = "#b5dfff")
dev.off()


## Mcap ##


#remove RRBS data
Mcap_cov_perc_group_3x_RRBS_rel <- Mcap_cov_perc_group_3x[which(Mcap_cov_perc_group_3x$method != "MBD"),]
#subset WGBS groups
Mcap_cov_perc_group_3x_RRBS_rel <- Mcap_cov_perc_group_3x_RRBS_rel[grep("WGBS", Mcap_cov_perc_group_3x_RRBS_rel$group_simp),]


#set cov_tag to zero for CpGs not covered by RRBS but covered by WGBS
Mcap_cov_perc_group_3x_RRBS_rel$RRBS_cov_tag <- ifelse(Mcap_cov_perc_group_3x_RRBS_rel$method == "WGBS" & !grepl("RRBS",Mcap_cov_perc_group_3x_RRBS_rel$group_simp),"[0]", Mcap_cov_perc_group_3x_RRBS_rel$cov_tag)

#set WGBS cov_tag to NA for CpGs covered by WGBS and by RRBS since we're gonna use the RRBS cov_tags
Mcap_cov_perc_group_3x_RRBS_rel$RRBS_cov_tag <- ifelse(Mcap_cov_perc_group_3x_RRBS_rel$method == "WGBS" & grepl("RRBS",Mcap_cov_perc_group_3x_RRBS_rel$group_simp),NA, Mcap_cov_perc_group_3x_RRBS_rel$RRBS_cov_tag)

#Set methylation to NA for CpGs covered by RRBS
Mcap_cov_perc_group_3x_RRBS_rel$wgbs_meth_tag <- ifelse(Mcap_cov_perc_group_3x_RRBS_rel$method == "WGBS",Mcap_cov_perc_group_3x_RRBS_rel$meth_tag,NA)


#remove all columns except loci, group, RRBS_cov_tag, and wgbs_meth_tag
#convert to long format and filter out rows with NAs 
#convert to wide format so that only RRBS_cov_tag and wgbs_meth_tag values remain
Mcap_cov_perc_group_3x_RRBS_rel_reshape <- Mcap_cov_perc_group_3x_RRBS_rel[,c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")


table(Mcap_cov_perc_group_3x_RRBS_rel_reshape$group_simp)

#summarize the data for plotting
Mcap_cov_perc_group_3x_RRBS_rel_sum <- Mcap_cov_perc_group_3x_RRBS_rel_reshape %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

jpeg("~/strigg/analyses/20200910/RRBScov_x_WGBSmeth_heatmap_Mcap.jpg", width = 6, height = 5, units = "in", res = 300 )
ggplot(data = Mcap_cov_perc_group_3x_RRBS_rel_sum, aes(x=wgbs_meth_tag, y=factor(RRBS_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from RRBS data)", fill = "number of\nCpGs (log10)") + scale_fill_continuous(high = "#132B43", low = "#b5dfff")
dev.off()
#```

#filter for loci covered by WGBS at 5x, 10x, 50x, 100x
#```{r}
loci_lessthan5x <- Pact_cov_perc_group_3x_RRBS_rel[which(Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]"),"loci"]

loci_lessthan10x <- Pact_cov_perc_group_3x_RRBS_rel[which(Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[5-10]" ),"loci"]

loci_lessthan20x <- Pact_cov_perc_group_3x_RRBS_rel[which(Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[5-10]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[10-20]"),"loci"]

loci_lessthan30x <- Pact_cov_perc_group_3x_RRBS_rel[which(Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[5-10]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[10-20]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[20-30]"),"loci"]

loci_lessthan40x <- Pact_cov_perc_group_3x_RRBS_rel[which(Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[5-10]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[10-20]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[20-30]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[30-40]"),"loci"]

loci_greaterthan50x <- Pact_cov_perc_group_3x_RRBS_rel[which(Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[50-100]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[100-500]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[500-1000]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[>1000]"),"loci"]

loci_greaterthan100x <- Pact_cov_perc_group_3x_RRBS_rel[which(Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[100-500]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[500-1000]" | Pact_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Pact_cov_perc_group_3x_RRBS_rel$cov_tag == "[>1000]"),"loci"]


Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs5x <- Pact_cov_perc_group_3x_RRBS_rel[which(!(Pact_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan5x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs10x <- Pact_cov_perc_group_3x_RRBS_rel[which(!(Pact_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan10x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs20x <- Pact_cov_perc_group_3x_RRBS_rel[which(!(Pact_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan20x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs30x <- Pact_cov_perc_group_3x_RRBS_rel[which(!(Pact_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan30x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs40x <- Pact_cov_perc_group_3x_RRBS_rel[which(!(Pact_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan40x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs50x <- Pact_cov_perc_group_3x_RRBS_rel[which(Pact_cov_perc_group_3x_RRBS_rel$loci %in% loci_greaterthan50x$loci),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs100x <- Pact_cov_perc_group_3x_RRBS_rel[which(Pact_cov_perc_group_3x_RRBS_rel$loci %in% loci_greaterthan100x$loci),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")
#summarize the data for plotting
Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs5x <- Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs5x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs10x <- Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs10x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs20x <- Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs20x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs30x <- Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs30x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs40x <- Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs40x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs50x <- Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs50x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs100x <- Pact_cov_perc_group_3x_RRBS_rel_reshape_wgbs100x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())


#create a df with all subsets
Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs <- dplyr::bind_rows(list(wgbs_1x=Pact_cov_perc_group_3x_RRBS_rel_sum, wgbs_5x=Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs5x, wgbs_10x= Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs10x, wgbs_20x= Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs20x,wgbs_30x= Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs30x,wgbs_40x= Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs40x,wgbs_50x= Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs50x, wgbs_100x= Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs100x), .id = "id")

Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs <- Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs %>% group_by(id, wgbs_meth_tag) %>% mutate(perc_CpG = CpGs/sum(CpGs)*100)

#make wgbs capitalized
Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs$id <- gsub("wgbs", "WGBS", Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs$id)







## Mcap ##


loci_lessthan5x <- Mcap_cov_perc_group_3x_RRBS_rel[which(Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]"),"loci"]

loci_lessthan10x <- Mcap_cov_perc_group_3x_RRBS_rel[which(Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[5-10]" ),"loci"]

loci_lessthan20x <- Mcap_cov_perc_group_3x_RRBS_rel[which(Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[5-10]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[10-20]"),"loci"]

loci_lessthan30x <- Mcap_cov_perc_group_3x_RRBS_rel[which(Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[5-10]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[10-20]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[20-30]"),"loci"]

loci_lessthan40x <- Mcap_cov_perc_group_3x_RRBS_rel[which(Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[1-5]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[5-10]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[10-20]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[20-30]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[30-40]"),"loci"]

loci_greaterthan50x <- Mcap_cov_perc_group_3x_RRBS_rel[which(Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[50-100]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[100-500]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[500-1000]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[>1000]"),"loci"]

loci_greaterthan100x <- Mcap_cov_perc_group_3x_RRBS_rel[which(Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[100-500]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[500-1000]" | Mcap_cov_perc_group_3x_RRBS_rel$method=="WGBS"  & Mcap_cov_perc_group_3x_RRBS_rel$cov_tag == "[>1000]"),"loci"]


Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs5x <- Mcap_cov_perc_group_3x_RRBS_rel[which(!(Mcap_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan5x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs10x <- Mcap_cov_perc_group_3x_RRBS_rel[which(!(Mcap_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan10x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs20x <- Mcap_cov_perc_group_3x_RRBS_rel[which(!(Mcap_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan20x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs30x <- Mcap_cov_perc_group_3x_RRBS_rel[which(!(Mcap_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan30x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs40x <- Mcap_cov_perc_group_3x_RRBS_rel[which(!(Mcap_cov_perc_group_3x_RRBS_rel$loci %in% loci_lessthan40x$loci)),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs50x <- Mcap_cov_perc_group_3x_RRBS_rel[which(Mcap_cov_perc_group_3x_RRBS_rel$loci %in% loci_greaterthan50x$loci),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")

Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs100x <- Mcap_cov_perc_group_3x_RRBS_rel[which(Mcap_cov_perc_group_3x_RRBS_rel$loci %in% loci_greaterthan100x$loci),c("loci","group_simp", "RRBS_cov_tag", "wgbs_meth_tag")] %>% gather("tag", "value", RRBS_cov_tag:wgbs_meth_tag) %>% na.omit() %>% spread("tag", "value")
#summarize the data for plotting
Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs5x <- Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs5x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs10x <- Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs10x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs20x <- Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs20x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs30x <- Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs30x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs40x <- Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs40x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs50x <- Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs50x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())

Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs100x <- Mcap_cov_perc_group_3x_RRBS_rel_reshape_wgbs100x %>% group_by(RRBS_cov_tag, wgbs_meth_tag) %>% summarise(CpGs =n())


#create a df with all subsets
Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs <- dplyr::bind_rows(list(wgbs_1x=Mcap_cov_perc_group_3x_RRBS_rel_sum, wgbs_5x=Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs5x, wgbs_10x= Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs10x, wgbs_20x= Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs20x,wgbs_30x= Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs30x,wgbs_40x= Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs40x,wgbs_50x= Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs50x, wgbs_100x= Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs100x), .id = "id")

Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs <- Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs %>% group_by(id, wgbs_meth_tag) %>% mutate(perc_CpG = CpGs/sum(CpGs)*100)

#make wgbs capitalized
Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs$id <- gsub("wgbs", "WGBS", Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs$id)

#```


#plot data
#```{r}
jpeg("~/strigg/analyses/20200910/RRBScov_x_WGBSmeth_heatmap_Pact_thresholds.jpg", width = 10, height = 5, units = "in", res = 300 )
ggplot(data = Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs, aes(x=wgbs_meth_tag, y=factor(RRBS_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from RRBS data)", fill = "number of\nCpGs (log10)") + facet_wrap(~factor(id, levels = c("WGBS_1x","WGBS_5x","WGBS_10x","WGBS_20x","WGBS_30x","WGBS_40x","WGBS_50x","WGBS_100x")), nrow = 2) + scale_fill_continuous(high = "#132B43", low = "#b5dfff")
dev.off()

#Mcap
jpeg("~/strigg/analyses/20200910/RRBScov_x_WGBSmeth_heatmap_Mcap_thresholds.jpg", width = 10, height = 5, units = "in", res = 300 )
ggplot(data = Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs, aes(x=wgbs_meth_tag, y=factor(RRBS_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from RRBS data)", fill = "number of\nCpGs (log10)") + facet_wrap(~factor(id, levels = c("WGBS_1x","WGBS_5x","WGBS_10x","WGBS_20x","WGBS_30x","WGBS_40x","WGBS_50x","WGBS_100x")), nrow = 2) + scale_fill_continuous(high = "#132B43", low = "#b5dfff")
dev.off()

#d <- ggplot(data = Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs, aes(x=wgbs_meth_tag, y=factor(RRBS_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from RRBS data)", fill = "number of\nCpGs (log10)", title = "P. acuta") + facet_wrap(~factor(id, levels = c("WGBS_1x","WGBS_5x","WGBS_10x","WGBS_20x","WGBS_30x","WGBS_40x","WGBS_50x","WGBS_100x")), nrow = 2) + scale_fill_continuous(high = "#132B43", low = "#b5dfff") + theme(plot.title = element_text(face = "italic"))


#c <- ggplot(data = Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs, aes(x=wgbs_meth_tag, y=factor(RRBS_cov_tag, levels = c("[0]",cov_tags)), fill=log(CpGs,10))) + geom_tile() +  theme(axis.text.x = element_text(size = 7, angle = 60, hjust = 1)) + labs(x = "% methylation\n(from WGBS data)", y = "number of reads\n(from RRBS data)", fill = "number of\nCpGs (log10)", title = "M. capitata") + facet_wrap(~factor(id, levels = c("WGBS_1x","WGBS_5x","WGBS_10x","WGBS_20x","WGBS_30x","WGBS_40x","WGBS_50x","WGBS_100x")), nrow = 2) + scale_fill_continuous(high = "#132B43", low = "#b5dfff") + theme(plot.title = element_text(face = "italic"))


#jpeg("~/strigg/analyses/20200910/RRBScov_x_WGBSmeth_heatmap_thresholds.jpg", width = 10, height = 5, units = "in", res = 300 )
#ggarrange(c,d, nrow = 2,labels = "AUTO",common.legend = T)
#dev.off()
#```

#write out data underneath each plot

#MBD vs WGBS data with coverage cutoffs
write.csv(Pact_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs, "~/strigg/analyses/Pact_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs.csv", row.names = F, quote = F)
write.csv(Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs, "~/strigg/analyses/20200910/Mcap_cov_perc_group_3x_mbd_rel_sum_wgbs_cutoffs.csv", row.names = F, quote = F)

#RRBS vs WGBS data with coverage cutoffs
write.csv(Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs, "~/strigg/analyses/Pact_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs.csv", row.names = F, quote = F)
write.csv(Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs, "~/strigg/analyses/20200910/Mcap_cov_perc_group_3x_RRBS_rel_sum_wgbs_cutoffs.csv", row.names = F, quote = F)

#all loci from different methods binned by cov and meth
write.csv(Pact_cov_perc_group_3x_sum, "~/strigg/analyses/Pact_cov_perc_group_3x_sum.csv", row.names = F, quote = F)
write.csv(Mcap_cov_perc_group_3x_sum, "~/strigg/analyses/20200910/Mcap_cov_perc_group_3x_sum", row.names = F, quote = F)
