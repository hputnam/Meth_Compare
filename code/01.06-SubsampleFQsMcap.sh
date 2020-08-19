#!/bin/bash
## Job Name
#SBATCH --job-name=SubsampleFQs
## Allocation Definition 
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes 
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=5-23:30:00
## Memory per node
#SBATCH --mem=500G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/strigg/analyses/20200515


### This script pools sample fastqs and randomly subsamples the files 
### The point it to get an estimate of % genome covered by X reads

%%bash

#record each command that is called in slurm file

set -ex


################
###   Mcap   ###
################


#copy data over to mox

#rsync --archive --progress --verbose --files-from=:/Volumes/web/metacarcinus/FROGER_meth_compare/20200515/Mcap_trim_reads_file_paths.txt strigg@ostrich.fish.washington.edu:/ .


#define RRBS and WGBS dir

RRBS_dir="Volumes/web/metacarcinus/FROGER_meth_compare/20200311/RRBS"
WGBS_MBD_dir="Volumes/web/metacarcinus/FROGER_meth_compare/20200311/WGBS_MBD"


### 200M ###

## Mcap WGBS

###create randomly subsampled read ID file
zcat \
${WGBS_MBD_dir}/Meth10_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth11_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth12_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 200000000 \
> Mcap_WGBS_200M_ID.txt

###Pull out R1 from read IDs
zcat \
${WGBS_MBD_dir}/Meth10_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth11_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth12_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_WGBS_200M_ID.txt \
- \
> Mcap_WGBS_200M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${WGBS_MBD_dir}/Meth10_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth11_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth12_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_WGBS_200M_ID.txt \
- \
> Mcap_WGBS_200M_R2_001_val_2.fq

## Mcap RRBS

###create randomly subsampled read ID file
zcat \
${RRBS_dir}/Meth13_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth14_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth15_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 200000000 \
> Mcap_RRBS_200M_ID.txt

###Pull out R1 from read IDs
zcat \
${RRBS_dir}/Meth13_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth14_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth15_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_RRBS_200M_ID.txt \
- \
> Mcap_RRBS_200M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${RRBS_dir}/Meth13_R2_001_val_2.fq.gz \
${RRBS_dir}/Meth14_R2_001_val_2.fq.gz \
${RRBS_dir}/Meth15_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_RRBS_200M_ID.txt \
- \
> Mcap_RRBS_200M_R2_001_val_2.fq


## Mcap MBD

###create randomly subsampled read ID file
zcat \
${WGBS_MBD_dir}/Meth16_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth17_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth18_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 200000000 \
> Mcap_MBD_200M_ID.txt

###Pull out R1 from read IDs
zcat \
${WGBS_MBD_dir}/Meth16_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth17_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth18_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_MBD_200M_ID.txt \
- \
> Mcap_MBD_200M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${WGBS_MBD_dir}/Meth16_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth17_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth18_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_MBD_200M_ID.txt \
- \
> Mcap_MBD_200M_R2_001_val_2.fq

gzip *.fq



### 150M ###

## Mcap WGBS

###create randomly subsampled read ID file
zcat \
${WGBS_MBD_dir}/Meth10_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth11_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth12_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 150000000 \
> Mcap_WGBS_150M_ID.txt

###Pull out R1 from read IDs
zcat \
${WGBS_MBD_dir}/Meth10_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth11_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth12_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_WGBS_150M_ID.txt \
- \
> Mcap_WGBS_150M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${WGBS_MBD_dir}/Meth10_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth11_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth12_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_WGBS_150M_ID.txt \
- \
> Mcap_WGBS_150M_R2_001_val_2.fq

## Mcap RRBS

###create randomly subsampled read ID file
zcat \
${RRBS_dir}/Meth13_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth14_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth15_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 150000000 \
> Mcap_RRBS_150M_ID.txt

###Pull out R1 from read IDs
zcat \
${RRBS_dir}/Meth13_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth14_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth15_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_RRBS_150M_ID.txt \
- \
> Mcap_RRBS_150M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${RRBS_dir}/Meth13_R2_001_val_2.fq.gz \
${RRBS_dir}/Meth14_R2_001_val_2.fq.gz \
${RRBS_dir}/Meth15_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_RRBS_150M_ID.txt \
- \
> Mcap_RRBS_150M_R2_001_val_2.fq


## Mcap MBD

###create randomly subsampled read ID file
zcat \
${WGBS_MBD_dir}/Meth16_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth17_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth18_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 150000000 \
> Mcap_MBD_150M_ID.txt

###Pull out R1 from read IDs
zcat \
${WGBS_MBD_dir}/Meth16_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth17_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth18_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_MBD_150M_ID.txt \
- \
> Mcap_MBD_150M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${WGBS_MBD_dir}/Meth16_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth17_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth18_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_MBD_150M_ID.txt \
- \
> Mcap_MBD_150M_R2_001_val_2.fq

gzip *.fq


### 100M ###

## Mcap WGBS

###create randomly subsampled read ID file
zcat \
${WGBS_MBD_dir}/Meth10_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth11_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth12_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 100000000 \
> Mcap_WGBS_100M_ID.txt

###Pull out R1 from read IDs
zcat \
${WGBS_MBD_dir}/Meth10_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth11_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth12_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_WGBS_100M_ID.txt \
- \
> Mcap_WGBS_100M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${WGBS_MBD_dir}/Meth10_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth11_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth12_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_WGBS_100M_ID.txt \
- \
> Mcap_WGBS_100M_R2_001_val_2.fq


## Mcap RRBS

###create randomly subsampled read ID file
zcat \
${RRBS_dir}/Meth13_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth14_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth15_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 100000000 \
> Mcap_RRBS_100M_ID.txt

###Pull out R1 from read IDs
zcat \
${RRBS_dir}/Meth13_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth14_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth15_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_RRBS_100M_ID.txt \
- \
> Mcap_RRBS_100M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${RRBS_dir}/Meth13_R2_001_val_2.fq.gz \
${RRBS_dir}/Meth14_R2_001_val_2.fq.gz \
${RRBS_dir}/Meth15_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_RRBS_100M_ID.txt \
- \
> Mcap_RRBS_100M_R2_001_val_2.fq


## Mcap MBD

###create randomly subsampled read ID file
zcat \
${WGBS_MBD_dir}/Meth16_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth17_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth18_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 100000000 \
> Mcap_MBD_100M_ID.txt

###Pull out R1 from read IDs
zcat \
${WGBS_MBD_dir}/Meth16_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth17_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth18_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_MBD_100M_ID.txt \
- \
> Mcap_MBD_100M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${WGBS_MBD_dir}/Meth16_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth17_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth18_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_MBD_100M_ID.txt \
- \
> Mcap_MBD_100M_R2_001_val_2.fq

gzip *.fq

### 50M ###

## Mcap WGBS

###create randomly subsampled read ID file
zcat \
${WGBS_MBD_dir}/Meth10_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth11_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth12_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 50000000 \
> Mcap_WGBS_50M_ID.txt

###Pull out R1 from read IDs
zcat \
${WGBS_MBD_dir}/Meth10_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth11_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth12_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_WGBS_50M_ID.txt \
- \
> Mcap_WGBS_50M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${WGBS_MBD_dir}/Meth10_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth11_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth12_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_WGBS_50M_ID.txt \
- \
> Mcap_WGBS_50M_R2_001_val_2.fq


## Mcap RRBS

###create randomly subsampled read ID file
zcat \
${RRBS_dir}/Meth13_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth14_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth15_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 50000000 \
> Mcap_RRBS_50M_ID.txt

###Pull out R1 from read IDs
zcat \
${RRBS_dir}/Meth13_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth14_R1_001_val_1.fq.gz \
${RRBS_dir}/Meth15_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_RRBS_50M_ID.txt \
- \
> Mcap_RRBS_50M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${RRBS_dir}/Meth13_R2_001_val_2.fq.gz \
${RRBS_dir}/Meth14_R2_001_val_2.fq.gz \
${RRBS_dir}/Meth15_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_RRBS_50M_ID.txt \
- \
> Mcap_RRBS_50M_R2_001_val_2.fq


## Mcap MBD

###create randomly subsampled read ID file
zcat \
${WGBS_MBD_dir}/Meth16_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth17_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth18_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 50000000 \
> Mcap_MBD_50M_ID.txt

###Pull out R1 from read IDs
zcat \
${WGBS_MBD_dir}/Meth16_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth17_R1_001_val_1.fq.gz \
${WGBS_MBD_dir}/Meth18_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_MBD_50M_ID.txt \
- \
> Mcap_MBD_50M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
${WGBS_MBD_dir}/Meth16_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth17_R2_001_val_2.fq.gz \
${WGBS_MBD_dir}/Meth18_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Mcap_MBD_50M_ID.txt \
- \
> Mcap_MBD_50M_R2_001_val_2.fq

gzip *.fq