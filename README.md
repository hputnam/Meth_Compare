# Analysis of DNA methylation in invertebrates requires consideration of genome characteristics and methylation landscape

Shelly A Trigg\*1, Yaamini R Venkataraman\*1, Mackenzie R Gavery2, Steven B Roberts1, Debashish Bhattacharya3, Alan Downey-Wall4, Jose Eirin-Lopez5, Kevin M Johnson6, Katie E Lotterhos4, Jonathan B. Puritz7 and Hollie M Putnam7+ 


1 University of Washington, School of Aquatic and Fishery Sciences 1122 NE Boat St. Seattle, WA, 98195, USA  
2 Environmental and Fisheries Sciences Division, Northwest Fisheries Science Center, National Marine Fisheries Service, National Oceanic and Atmospheric Administration, 2725 Montlake Blvd E, Seattle, WA, 98112, USA  
3 Rutgers  
4 Department of Marine and Environmental Sciences, Northeastern University, 430 Nahant Road, Nahant, MA 01908  
5 Florida International University, Environmental Epigenetics Laboratory, Institute of Environment 3000 NE 151 St. North Miami, FL, 33181, USA  
6 Department of Biological Sciences, University of Rhode Island, Kingston, RI 02881, USA  

+ corresponding author: hputnam@uri.edu  
\* equal contribution

Keywords:  bisulfite sequencing, coral, epigenetics, marine invertebrate

[![issues](https://img.shields.io/github/issues/hputnam/Meth_Compare.svg)](https://img.shields.io/github/issues/hputnam/Meth_Compare)
[![NSF-1921149](https://img.shields.io/badge/NSF-1921149-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1921149)
[![NSF-1921149](https://img.shields.io/badge/NSF-1921465-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1921465)
[![OSF](https://img.shields.io/badge/OSF-x5waz-blueviolet.svg)](https://osf.io/x5waz/)


## Repository Structure

- [metadata](https://github.com/hputnam/Meth_Compare/tree/master/metadata): Sample metadata
- [lab-methods](https://github.com/hputnam/Meth_Compare/tree/master/lab-methods): Commercial protocols, kit information, and lab notebook entries for DNA extraction, enrichment, and library preparation
- [genome-feature-tracks](https://github.com/hputnam/Meth_Compare/tree/master/genome-feature-files): Genome feature tracks generated for *M. capitata* and *P. acuta*
- [code](https://github.com/hputnam/Meth_Compare/tree/master/code): Bash scripts, R Markdown files, and Jupyter notebooks used to analyze data. Files are organized by corresponding methods sections.
- [output](https://github.com/hputnam/Meth_Compare/tree/master/output): Individual subdirectories for each analysis along with intermediate output
    - [intermediate-files]()
         - _subdirectories_
    - [figures](https://github.com/hputnam/Meth_Compare/tree/master/output/figures): Figures in the manuscript main text
    - [supplemental-material](https://github.com/hputnam/Meth_Compare/tree/master/output/supplemental-material): supplemental material referenced in manuscript

## Abstract

Elucidation of genome to phenome pathways is essential for our understanding of biology, ecology, and evolution. The study of mechanisms underlying the emergence of phenotypic plasticity from genome-environment interactions has recently broadened to include epigenetic mechanisms, notably DNA methylation. Given the suite of anthropogenic impacts on the world’s oceans, marine organisms are under particular threat. In marine invertebrates, a broad group with diverse life-histories, including sessile organisms that need to rapidly respond to changing environmental factors, there has been a growing focus on epigenetic regulation of gene expression through DNA methylation. Unfortunately, genome-wide DNA methylation studies using bisulfite sequencing have been hampered by high costs and potential biases due to methods in non-model organisms. To improve understanding of differences in the bisulfite sequencing methods and ultimately how genome-wide DNA methylation studies are carried out, we compared three methods for quantifying DNA methylation at single base pair resolution — Whole Genome Bisulfite Sequencing (WGBS), Reduced Representation Bisulfite Sequencing (RRBS), and Methyl-CpG Binding Domain Bisulfite Sequencing (MBDBS) — using replicate individuals from two reef-building coral species with contrasting environmental sensitivity, the more resistant Montipora capitata and the more sensitive Pocillopora acuta. All bisulfite sequencing methods reveal the majority of CpG methylation in both species occurs in gene bodies and flanking regions. There is substantially greater methylation in M. capitata (11.4%) than P. acuta (2.9%). In both species, at our sequencing coverage and depth, using MBDBS generates the greatest capacity for detecting CpGs in coding regions, which is beneficial for analysis of species with predominantly gene body methylation. However MBDBS does have the potential limitations associated with intra-sample methylation heterogeneity. RRBS offers robust information for specific loci, though the proportion of the genome coverage is significantly reduced. As expected, the relative genome size strongly influences the number and location of CpGs covered by each method, providing additional nuances in cross-species comparisons. Together these findings reinforce the role and importance of DNA methylation to phenotypic plasticity in critical marine invertebrate taxa and provide a resource for informed decision making for examining the functional role of DNA methylation.


## Results

- Figures
	- [Figure 1]: Comparison of CpG site coverage at different sequencing depths across library preparation methods. (A) Fraction of CpG sites in the genome covered at sequencing coverage (X) by M. capitata samples from each library preparation (see table 1 for number of reads). (B) Fraction of CpG sites in the genome covered at each sequencing depth (X) across different numbers of read pairs (M reads) for M. capitata samples from each library preparation. (C) Fraction of CpG sites in the genome covered at sequencing coverage (X) by P. acuta samples from each library preparation. (D) Fraction of CpG sites in the genome covered at each sequencing depth (X) across different numbers of read pairs (M reads) for P. acuta samples from each library preparation. For (A,B), Loess smoothing was applied and 95% confidence intervals are shown by shaded areas. For (F,G), results obtained from downsampling analyses are shown for sequencing depths ranging from 1–50X for library sizes of 50–200 M read pairs. All samples were pooled for the downsampling analyses (*code [here]*).
	- [Figure 2]: The number of unique CpGs (vertical bars) for each intersection set (connected dots) for A) M. capitata and B) P. acuta. Color key shows intersection sets in venn diagram format with CpGs covered by all three methods (black), by any two methods (blue), by only one method (red), or CpGs not covered by any method (yellow). The CpG set size for the genome and for each method is shown in the horizontal bars (gray) (*code [here]*).
	- [Figure 3]: Percent of highly methylated (≥ 50%; darkest shade), moderately methylated (10-50%; medium shade), and lowly methylated CpGs (< 10%; lightest shade) detected by each method A) for M. capitata and B) P. acuta (*code [here]*).
	- [Figure 4]: Percent of CpGs detected by sequencing methods in coding sequences (CDS), introns, 1 kb flanking regions upstream (Upstream Flank) or downstream of genes (Downstream Flank), and intergenic regions A) for M. capitata and B) P. acuta (*code [here]*).
	- [Figure 5]: Matrix of pairwise scatter plots for CpG loci covered at > 5x across all samples for A) M. capitata (n=4,666 common loci) and B) P. acuta (n=93,714 common loci). The red lines represent linear regression fits and the green lines are polynomial regression fits. Pearson correlation coefficients for each pairwise comparison are presented in the upper right boxes. Methods are color coded on the X and Y axes (WGBS = green, MBDBS = purple, and RRBS = orange) and replicate samples are indicated on the diagonal along with histograms of % CpG methylation (*code [here]*).
	- [Figure 6]: Mean proportion (n=3 samples per method) of CpGs per gene that have at least 5x coverage in all of the one-to-one orthologous genes, as identified by OrthoFinder (Supplemental File X) for A) Montipora capitata and B) Pocillopora acuta. The heatmaps are ordered by orthogroup for the groups with one to one gene matches (*code [here]*).
