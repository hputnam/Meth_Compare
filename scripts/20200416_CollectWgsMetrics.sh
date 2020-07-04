## Run picard tools CollectWgsMetrics

#Mcap alignments
find /Users/strigg/Desktop/20200416/Mcap/*.sorted.bam|\
xargs basename -s .sorted.bam \
| xargs -I{} java -jar ~/anaconda3/bin/picard.jar CollectWgsMetrics \
I=/Users/strigg/Desktop/20200416/Mcap/{}.sorted.bam \
O=/Users/strigg/Desktop/20200416/Mcap/{}.collect_wgs_metrics.txt \
R=/Volumes/web/seashell/bu-mox/data/froger/Mcap_Genome/Mcap.genome_assembly.fa

#Pact alignments
find /Users/strigg/Desktop/20200416/Pact/*.sorted.bam|\
xargs basename -s .sorted.bam \
| xargs -I{} java -jar ~/anaconda3/bin/picard.jar CollectWgsMetrics \
I=/Users/strigg/Desktop/20200416/Pact/{}.sorted.bam \
O=/Users/strigg/Desktop/20200416/Pact/{}.collect_wgs_metrics.txt \
R=/Volumes/web/seashell/bu-mox/data/froger/Pact_Genome/Pocillopora_acuta_genome_v1.fasta

#Pact_C1 alignments
find /Users/strigg/Desktop/20200416/Pact_C1/*.sorted.bam|\
xargs basename -s .sorted.bam \
| xargs -I{} java -jar ~/anaconda3/bin/picard.jar CollectWgsMetrics \
I=/Users/strigg/Desktop/20200416/Pact_C1/{}.sorted.bam \
O=/Users/strigg/Desktop/20200416/Pact_C1/{}.collect_wgs_metrics.txt \
R=/Volumes/web/seashell/bu-mox/data/froger/C1_Genome/SymbC1.Genome.Scaffolds.fasta
