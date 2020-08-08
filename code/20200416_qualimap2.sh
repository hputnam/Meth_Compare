%%bash

set -ex


## Run Qualimap to generate individual bamqc files for each sample

#Mcap
FILES=/Users/strigg/Desktop/20200416/Mcap/*.bam

for f in $FILES
do
/Users/Shared/bioinformatics/qualimap_v2.2.1/qualimap bamqc \
--java-mem-size=44G \
-nt 28 \
-bam ${f}
done

#Pact
FILES=/Users/strigg/Desktop/20200416/Pact/*.bam

for f in $FILES
do
/Users/Shared/bioinformatics/qualimap_v2.2.1/qualimap bamqc \
--java-mem-size=44G \
-nt 28 \
-bam ${f} 
done

#C1
FILES=/Users/strigg/Desktop/20200416/Pact_C1/*.bam

for f in $FILES
do
/Users/Shared/bioinformatics/qualimap_v2.2.1/qualimap bamqc \
--java-mem-size=44G \
-nt 28 \
-bam ${f} 
done

## Run Qualimap to generate multibamqc report

Mcap
/Users/Shared/bioinformatics/qualimap_v2.2.1/qualimap multi-bamqc \
--java-mem-size=44G \
-d /Users/strigg/Desktop/20200416/Mcap/multibamQC/multiBamqc_data_paths.txt \
-outdir /Users/strigg/Desktop/20200416/Mcap/multibamQC 


#Pact
/Users/Shared/bioinformatics/qualimap_v2.2.1/qualimap multi-bamqc \
--java-mem-size=44G \
-d /Users/strigg/Desktop/20200416/Pact/multibamQC/multiBamqc_data_paths.txt \
-outdir /Users/strigg/Desktop/20200416/Pact/multibamQC


#C1
/Users/Shared/bioinformatics/qualimap_v2.2.1/qualimap multi-bamqc \
--java-mem-size=44G \
-d /Users/strigg/Desktop/20200416/Pact_C1/multibamQC/multiBamqc_data_paths.txt \
-outdir /Users/strigg/Desktop/20200416/Pact_C1/multibamQC
