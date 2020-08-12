# code

## DNA Sequencing Processing
- [00-Meth_Compare_Pipeline.md](00-Meth_Compare_Pipeline.md): Markdown file with code used corresponding toe "DNA Sequence Processing" section in Methods. From raw data to calling of methlyation state at all loci.


## Quality Control
- [MethCompare_MultiQC.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/MethCompare_MultiQC.ipynb):  Script used to run MultiQC on FastQC output for trimmed BS reads that generates the input files for [FormatMultiQC.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/analyses/FormatMultiQC/FormatMultiQC.Rmd)
- [FormatMultiQC.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/scripts/FormatMultiQC.Rmd):  This script produces Table 1 and Supplementary Tables X-X from MultiQC output files, including lambda alignments. It also calculates lambda conversion efficiency based on the ratio of the sum of all unmethylated cytosines in CHG and CHH context to the sum of methylated and unmethylated cytosines in CHG and CHH.
- [CompareConversionEfficiency.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/scripts/CompareConversionEfficiency.Rmd):  This script generates a table showing conversion efficiency based on lambda alignments and estimated conversion efficiency based on coral alignments for each sample.
- [20200416_qualimap2.sh](https://github.com/hputnam/Meth_Compare/blob/master/scripts/20200416_qualimap2.sh):  This script runs Qualimap bamqc and multi-bamQC on deduplicated WBGS and MDBBS and non-deduplicated RRBS sorted bam files for each species, and produces both individual sample Qualimap reports and multi-sample BAM QC reports (which include PCAs for each species).
- [Qualimap\_MultiBamQC\_PCA.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/analyses/Coverage_analysis/Qualimap_MultiBamQC_PCA.Rmd):  This script takes in the sample summary tables reported in the multi-sample BAM QC reports, runs PCA, and generates score plots for each species.

## Coverage analysis

#### Downsampling analysis
- [20200505_SubsampleFQs2.sh](https://github.com/hputnam/Meth_Compare/blob/master/scripts/20200505_SubsampleFQs2.sh):  This script pools trimmed reads by method for P. acuta samples and randomly downsamples read pairs at 50M, 100M, 150M, and 200M.
- [20200525_SubsampleFQsMcap.sh](https://github.com/hputnam/Meth_Compare/blob/master/scripts/20200525_SubsampleFQsMcap.sh):  This script pools trimmed reads by method for M. capitata samples and randomly downsamples read pairs at 50M, 100M, 150M, and 200M.
- [20200512_SubsampleBmrkPact.sh](https://github.com/hputnam/Meth_Compare/blob/master/scripts/20200512_SubsampleBmrkPact.sh):  This script performs Bismark alignment of P. acuta downsampled reads, methylation extraction, deduplicates WGBS and MBDBS downsampled reads, generates cytosine reports, generates sorted bam files, and generates 5x coverage files.
- [20200527_SubsampleBmrkMcap.sh](https://github.com/hputnam/Meth_Compare/blob/master/scripts/20200527_SubsampleBmrkMcap.sh):  This script performs Bismark alignment of M. capitata downsampled reads, methylation extraction, deduplicates WGBS and MBDBS downsampled reads, generates cytosine reports, generates sorted bam files, and generates 5x coverage files.
- [20200513_SubsamplePicardPact.sh](https://github.com/hputnam/Meth_Compare/blob/master/scripts/20200513_SubsamplePicardPact.sh):  This script runs Picard CollectWgsMetrics on on deduplicated WBGS and MDBBS and non-deduplicated RRBS sorted bam files from downsampled P.acuta reads.
- [20200527_SubsamplePicardMcap.sh](https://github.com/hputnam/Meth_Compare/blob/master/scripts/20200527_SubsamplePicardMcap.sh):  This script runs Picard CollectWgsMetrics on on deduplicated WBGS and MDBBS and non-deduplicated RRBS sorted bam files from downsampled P.acuta reads.

#### CpG Coverage
- [Mcap\_CpG\_coverageXdepth.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Mcap_CpG_coverageXdepth.ipynb):  This script totals CpGs at different levels of coverage for individual _M. capitata_ samples and for pooled _M. capitata_ samples based on method and downsampled to different sequencing depths.
- [Pact\_CpG\_coverageXdepth.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Pact_CpG_coverageXdepth.ipynb):  This script totals CpGs at different levels of coverage for individual _P. acuta_ samples and for pooled _P. acuta_ samples based on method and downsampled to different sequencing depths.

#### Genome Coverage
- [20200416_CollectWgsMetrics.sh](https://github.com/hputnam/Meth_Compare/blob/master/scripts/20200416_CollectWgsMetrics.sh):  This script runs Picard CollectWgsMetrics on on deduplicated WBGS and MDBBS and non-deduplicated RRBS sorted bam files for each sample.

#### Coverage plots
- [Genome\_CpG\_coverage_analysis.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/analyses/Coverage_analysis/Genome_CpG_coverage_analysis.Rmd):  This script generates genome and CpG coverage plots from CpG_coverageXdepth script outputs and CollectWgsMetrics output for individual samples and downsampling analysis output.


## Upset plots
- [Generate\_UpsetPlot\_input.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Generate_UpsetPlot_input.ipynb):  This script sorts and merges by method sample bedgraphs of CpG loci with 5x coverage for each species. It then generates a union bed file for each species from merged begraphs, and genome CpG bedgraph (contains all CpG loci in the genome).
- [GenerateUpsetPlot.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/scripts/GenerateUpsetPlot.Rmd):  This script generates an upset plot from the union bed graph for each species.

## Methylation status and genomic locations

- [Generating-Genome-Feature-Tracks.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Generating-Genome-Feature-Tracks.ipynb): Create gene, CDS, intron, flanking region, and intergenic feature tracks for *M. capitata* and *P. acuta* genomes
- [Characterizing-CpG-Methylation-5x.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Characterizing-CpG-Methylation-5x.ipynb): Characterize methylation status and genomic locations of CpGs in individual sample data for both species using `bedtools`
- [Characterizing-CpG-Methylation-5x-Union.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Characterizing-CpG-Methylation-5x-Union.ipynb): Characterize methylation status and genomic locations of CpGs from 5x union bedgraphs for both species using `bedtools`
- [Characterizing-CpG-Methylation-5x-Union-Summary-Plots.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Characterizing-CpG-Methylation-5x-Union-Summary-Plots.Rmd): Create summary tables and stacked barplots to understand methylation status and genomic locations of 5x union CpGs in both species
- [Identifying-Genomic-Locations.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Identifying-Genomic-Locations.ipynb): Characterize genomic locations of CpGs in upset plots and method-associated DMC using `bedtools`
- [Identifying-Genome-Features-Summary.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Identifying-Genome-Features-Summary.Rmd): Create summary tables for upset plot data and various method-associated DMC categories
