SupplementaryMaterial
================
Shelly Trigg
11/6/2020

## Supplemental Methods

**Bisulfite Conversion Efficiency Assessment** Trimmed sequence reads
were aligned to the genome of E. coli strain K-12 MG1655 (Riley et al.,
2006) using Bismark v0.21.0 with the –non\_directional option and
alignment stringency set by -score\_min L,0,-0.6. Conversion efficiency
was calculated as the ratio of the sum of all unmethylated cytosines in
CHG and CHH context to the sum of methylated and unmethylated cytosines
in CHG and CHH. Conversion efficiency was also estimated from coral
alignments as the ratio of the sum of unmethylated cytosines in CHG and
CHH context to the sum of methylated and unmethylated cytosines in CHG
and CHH. ANOVA was used to test for an effect of library preparation
method on conversion efficiency within each species (conversion
efficiency \~ library preparation method) for both estimated and lambda
alignment calculated conversion efficiencies. A two-sample t-test was
used to test if conversion efficiency calculated from lambda alignments
was the same as estimated conversion efficiency for each different
library preparation methods within each species (lambda conversion
efficiency \~ estimated conversion efficiency).

**Methylation Characterization by PCoA** Sequencing method differences
in the proportion of high, moderate, and low methylated CpGs were
investigated using Principal Coordinates Analysis (PCoA) and
permutational multivariate analysis of variance (PERMANOVA). The
hypothesis was that each method detected the same proportion of highly,
moderately, and lowly methylated CpGs within a species. A euclidean
dissimilarity matrix calculated using centered log-ratio transformed
proportion data was used for both the PCoA and PERMANOVA. Samples were
compared visually with a PCoA using cmdscale in the R package vegan
(Oksanen et al., 2016). A global PERMANOVA test was conducted to discern
significant differences, with pairwise PERMANOVA tests conducted for
significant global PERMANOVA results. A centroid beta dispersion model
(betadisper) and Analysis of Variance (ANOVA) determined if significant
PERMANOVA results were due to centroid or variance differences between
sequencing methods.

Significant global PERMANOVA tests were investigated further using
generalized linear mixed effect models with a beta error distribution
and logit link on untransformed proportion data as the response
variable, and the explanatory variables of sequencing method (fixed
effect) and biological replicate (random effect). Since a full model
violates independence assumptions, a separate model for each of the
three methylation statuses was constructed using glmmTMB (Brooks et al.,
2017). Post-hoc estimated marginal means were obtained using emmeans,
and a log odds ratio test was conducted to determine significant
pairwise differences (Lenth, Singmann, Love, Buerkner, & Herve, 2018). A
Benjamini-Hochberg FDR threshold was used to correct P-values for
multiple comparisons, with adjusted P-values less than 0.05 considered
significant.

## CpG Coverage

Qualimap v2.2.1 MultiBamQC (Okonechnikov, Conesa, & García-Alcalde,
2016) was run on deduplicated WBGS and MBDBS and non-deduplicated RRBS
bam files that were sorted using SAMtools v.1.9 (H. Li et al., 2009).
PCA plots produced by MultiBamQC were further formatted in R by
exporting the sample summary tables from the MultiBamQC reports,
re-running the PCA on the summary tables using the R function prcomp (R
Core Team, 2020), and plotting the output with custom colors using the R
package ggplot2 (Gómez-Rubio, 2017). To assess average genome-wide CpG
coverage, the number of Cytosines passing different read depth
thresholds were totaled from sample CpG coverage reports produced by . A
downsampling analysis was performed to estimate overall genome-wide CpG
coverage by pooling all sample reads within a method and species.
Briefly, fastq files of trimmed reads were concatenated for each method
and species then randomly subsampled to 50, 100, 150, and 200 million
reads. Next, alignment and methylation calling were carried out as
described above on each subset, and the number of cytosines passing
different read depth thresholds were totaled from CpG coverage reports
from each subset. The totaled CpGs were then imported into R and
coverage was plotted across different read depth thresholds using
ggplot2 (Gómez-Rubio, 2017). Sequencing saturation was estimated from a
Michaelis-Menten model with the ‘drm’ function from the R package drc
(Ritz, Baty, Streibig, & Gerhard, 2015) using CpG coverage reports from
subsampled data as input, and both observed CpG coverage from subsampled
data and estimated CpG coverage were plotted using the R package ggplot2
(Gómez-Rubio, 2017).

To assess overlap among CpGs in each sample, CpG coverage reports were
filtered for 5x coverage, saved in bedgraph format, and for each method
and species were sorted and merged using bedtools v2.26.0 (Quinlan &
Hall, 2010). To include overlap with CpGs in the genome, a CpG gff was
generated using fuzznuc on the Galaxy web interface (version 5.0.2)
(Jalili et al., 2020) and converted into bedgraph format. For each
species, a union bed file was generated with the function unionBedGraphs
in bedtools v2.26.0 (Quinlan & Hall, 2010) using the merged bedgraph
files for each method and the genome CpG bedgraph as input. The union
bed files were imported into R and upset plots were generated for each
species using the R package ComplexHeatmap (Gu, Eils, & Schlesner,
2016), with legend added manually.
