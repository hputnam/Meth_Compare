## Data and code
associated with the manuscript:


**"Invertebrate methylomes provide insight into mechanisms of environmental tolerance and reveal methodological biases"**

Shelly A Trigg\*<sup>1</sup>, Yaamini R Venkataraman\*<sup>1</sup>, Mackenzie R Gavery<sup>2</sup>, Steven B Roberts<sup>1</sup>, Debashish Bhattacharya<sup>3</sup>, Alan Downey-Wall<sup>4</sup>, Jose Eirin-Lopez<sup>5</sup>, Kevin M Johnson<sup>6,7</sup>, Katie E Lotterhos<sup>4</sup>, Jonathan B. Puritz<sup>8</sup> and Hollie M Putnam<sup>8</sup>+ 


<sup>1</sup> University of Washington, School of Aquatic and Fishery Sciences 1122 NE Boat St. Seattle, WA, 98195, USA  
<sup>2</sup> Environmental and Fisheries Sciences Division, Northwest Fisheries Science Center, National Marine Fisheries Service, National Oceanic and Atmospheric Administration, 2725 Montlake Blvd E, Seattle, WA, 98112, USA  
<sup>3</sup> Department of Biochemistry and Microbiology, Rutgers University, New Brunswick, NJ 08901 USA  
<sup>4</sup> Department of Marine and Environmental Sciences, Northeastern University, 430 Nahant Road, Nahant, MA 01908  
<sup>5</sup> Florida International University, Environmental Epigenetics Laboratory, Institute of Environment 3000 NE 151 St. North Miami, FL, 33181, USA  
<sup>6</sup> Center for Coastal Marine Sciences, California Polytechnic State University, San Luis Obispo, CA, 93407, USA  
<sup>7</sup> California Sea Grant, University of California San Diego, La Jolla, CA, 92093  
<sup>8</sup> Department of Biological Sciences, University of Rhode Island, Kingston, RI 02881, USA  

\+ corresponding author: hputnam@uri.edu  
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
    - [intermediate-files](https://github.com/hputnam/Meth_Compare/tree/master/output/intermediate-files)
    - [figures](https://github.com/hputnam/Meth_Compare/tree/master/output/figures): Figures in the manuscript main text
    - [supplemental-material](https://github.com/hputnam/Meth_Compare/tree/master/output/supplemental-material): supplemental material referenced in manuscript

[Raw Sequence Data NCBI PRJNA691891](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA691891)

## Abstract

There is a growing focus on the role of DNA methylation in the ability of marine invertebrates to rapidly respond to changing environmental factors and anthropogenic impacts. However, genome-wide DNA methylation studies in non-model organisms are currently hampered by limited understanding of methodological biases. Here we compare three methods for quantifying DNA methylation at single base pair resolution — Whole Genome Bisulfite Sequencing (WGBS), Reduced Representation Bisulfite Sequencing (RRBS), and Methyl-CpG Binding Domain Bisulfite Sequencing (MBDBS) — using multiple individuals from two reef-building coral species with contrasting environmental sensitivity. All methods reveal substantially greater methylation in Montipora capitata (11.4%) than the more sensitive Pocillopora acuta (2.9%). The majority of CpG methylation in both species occurs in gene bodies and flanking regions. In both species, MBDBS has the greatest capacity for detecting CpGs in coding regions at our sequencing depth, however MBDBS may be limited by intra-sample methylation heterogeneity. RRBS yields robust information for specific loci albeit without enrichment of any particular genome feature and with significantly reduced genome coverage. Relative genome size strongly influences the number and location of CpGs detected by each method when sequencing depth is limited, illuminating nuances in cross-species comparisons. These findings reinforce the role and importance of DNA methylation underlying environmental sensitivity in critical marine invertebrate taxa, and provide a genomic resource for investigating the functional role of DNA methylation in environmental tolerance.



## Figures

- [**Figure 1**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_1.jpg): Experimental design. Three biological replicate coral samples were obtained from both *M. capitata* and *P. acuta*. DNA was extracted from each coral sample and split for use in Whole Genome Bisulfite Sequencing (WGBS), Reduced Representation Bisulfite Sequencing (RRBS), and Methyl-CpG Binding Domain Bisulfite Sequencing (MBDBS) library preparation methods. Three libraries were generated for each of the three methods, yielding nine libraries for each species and 18 libraries total. ([code here](https://github.com/hputnam/Meth_Compare/tree/master/code))
- [**Figure 2**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_2.jpg): Circos plots showing the mean percent methylation of CpGs with 5x cover for each of the three methods on the largest scaffolds of each genome. The outer track shows the scaffold locations and dots indicate the percent methylation as indicated by the y-axes from 0-100% for each of the inner tracks. ([code here](https://github.com/hputnam/Meth_Compare/tree/master/code))
- [**Figure 3**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_3.jpg): Matrix of pairwise scatter plots for  shared CpG loci (i.e., CpG covered at > 5x across all samples) for (**A**) *M. capitata* (n = 4,666 common loci) and (**B**) *P. acuta* (n = 93,714 common loci). The red lines represent linear regression fits and the green lines are polynomial regression fits. Pearson correlation coefficients for each pairwise comparison are presented in the upper right boxes. Methods are color coded on the X and Y axes (WGBS = green, MBDBS = purple, and RRBS = orange) and replicate samples are indicated on the diagonal along with histograms of % CpG methylation. ([code here](https://github.com/hputnam/Meth_Compare/tree/master/code))
- [**Figure 4**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_4.jpg): CpG site coverage across library preparation methods. Mean fraction of CpG sites in the genome covered at different sequencing depths (read depths) by (**A**) MBDBS libraries, (**B**) RRBS libraries, and (**C**) WGBS libraries with standard deviations shown by shaded areas (see [Supplementary Table 2](https://github.com/hputnam/Meth_Compare/tree/master/output/supplemental-material) for number of reads in each sample). ([code here](https://github.com/hputnam/Meth_Compare/tree/master/code))
- [**Figure 5**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_5.jpg): Proportion of CpGs detected by sequencing methods in coding sequences (CDS), introns, 1 Kb flanking regions upstream (Upstream Flank) or downstream of genes (Downstream Flank), and intergenic regions for (**A**) *M. capitata* and (**B**) *P. acuta*. Principal Coordinate Analyses associated with perMANOVA and beta-dispersion tests related to [Supplementary Table 6](https://github.com/hputnam/Meth_Compare/tree/master/output/supplemental-material) that show differences in proportion of CpGs in various genomic locations (CDS, introns, upstream flanks, downstream flanks, and intergenic regions) for (**C**) *M. capitata* and **D**) *P. acuta*. WGBS is represented by green circles, RRBS by purple triangles, and MBDBS by orange diamonds. Percent variation explained by each PCoA axis is included in the axis label. Ellipses depict 95% confidence intervals for each sequencing method. All eigenvectors are significant at the α = 0.05 level. ([code here](https://github.com/hputnam/Meth_Compare/tree/master/code))
