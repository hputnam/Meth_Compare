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
#SBATCH --mem=100G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/strigg/analyses/20200505


### This script pools sample fastqs and randomly subsamples the files 
### The point it to get an estimate of % genome covered by X reads

%%bash

#record each command that is called in slurm file

set -ex

################
###   Pact   ###
################

### 200M ###

## Pact WGBS

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth1_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth2_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth3_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 200000000 \
> Pact_WGBS_200M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth1_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth2_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth3_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_WGBS_200M_ID.txt \
- \
> Pact_WGBS_200M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth1_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth2_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth3_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_WGBS_200M_ID.txt \
- \
> Pact_WGBS_200M_R2_001_val_2.fq


## Pact RRBS

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 200000000 \
> Pact_RRBS_200M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_RRBS_200M_ID.txt \
- \
> Pact_RRBS_200M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_RRBS_200M_ID.txt \
- \
> Pact_RRBS_200M_R2_001_val_2.fq


## Pact MBD

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 200000000 \
> Pact_MBD_200M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_MBD_200M_ID.txt \
- \
> Pact_MBD_200M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_MBD_200M_ID.txt \
- \
> Pact_MBD_200M_R2_001_val_2.fq

gzip *.fq


### 150M ###

## Pact RRBS

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 150000000 \
> Pact_RRBS_150M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_RRBS_150M_ID.txt \
- \
> Pact_RRBS_150M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_RRBS_150M_ID.txt \
- \
> Pact_RRBS_150M_R2_001_val_2.fq


## Pact MBD

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 150000000 \
> Pact_MBD_150M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_MBD_150M_ID.txt \
- \
> Pact_MBD_150M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_MBD_150M_ID.txt \
- \
> Pact_MBD_150M_R2_001_val_2.fq

gzip *.fq


### 100M ###

## Pact WGBS

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth1_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth2_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth3_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 100000000 \
> Pact_WGBS_100M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth1_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth2_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth3_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_WGBS_100M_ID.txt \
- \
> Pact_WGBS_100M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth1_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth2_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth3_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_WGBS_100M_ID.txt \
- \
> Pact_WGBS_100M_R2_001_val_2.fq


## Pact RRBS

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 100000000 \
> Pact_RRBS_100M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_RRBS_100M_ID.txt \
- \
> Pact_RRBS_100M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_RRBS_100M_ID.txt \
- \
> Pact_RRBS_100M_R2_001_val_2.fq


## Pact MBD

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 100000000 \
> Pact_MBD_100M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_MBD_100M_ID.txt \
- \
> Pact_MBD_100M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_MBD_100M_ID.txt \
- \
> Pact_MBD_100M_R2_001_val_2.fq

gzip *.fq

### 50M ###

## Pact WGBS

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth1_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth2_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth3_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 50000000 \
> Pact_WGBS_50M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth1_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth2_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth3_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_WGBS_50M_ID.txt \
- \
> Pact_WGBS_50M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth1_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth2_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth3_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_WGBS_50M_ID.txt \
- \
> Pact_WGBS_50M_R2_001_val_2.fq


## Pact RRBS

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 50000000 \
> Pact_RRBS_50M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_RRBS_50M_ID.txt \
- \
> Pact_RRBS_50M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth4_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth5_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth6_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_RRBS_50M_ID.txt \
- \
> Pact_RRBS_50M_R2_001_val_2.fq


## Pact MBD

###create randomly subsampled read ID file
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R1_001_val_1.fq.gz |\
awk '{if($1 ~/@A00387:/)print $1}' |\
shuf -n 50000000 \
> Pact_MBD_50M_ID.txt

###Pull out R1 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R1_001_val_1.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R1_001_val_1.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_MBD_50M_ID.txt \
- \
> Pact_MBD_50M_R1_001_val_1.fq

###Pull out R2 from read IDs
zcat \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth7_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth8_R2_001_val_2.fq.gz \
/gscratch/scrubbed/sr320/041720-pan-u/Pact_trim/Meth9_R2_001_val_2.fq.gz |\
awk 'NR==FNR{a[$1]=$1;next}$1 in a{x=NR+3}(NR<=x){print}' \
Pact_MBD_50M_ID.txt \
- \
> Pact_MBD_50M_R2_001_val_2.fq

gzip *.fq