## Data and code
associated with the manuscript:


**"Analysis of DNA methylation in invertebrates requires consideration of genome characteristics and methylation landscape"**

Shelly A Trigg\*1, Yaamini R Venkataraman\*1, Mackenzie R Gavery2, Steven B Roberts1, Debashish Bhattacharya3, Alan Downey-Wall4, Jose Eirin-Lopez5, Kevin M Johnson6, Katie E Lotterhos4, Jonathan B. Puritz7 and Hollie M Putnam7+ 


1 University of Washington, School of Aquatic and Fishery Sciences 1122 NE Boat St. Seattle, WA, 98195, USA
2 Environmental and Fisheries Sciences Division, Northwest Fisheries Science Center, National Marine Fisheries Service, National Oceanic and Atmospheric Administration, 2725 Montlake Blvd E, Seattle, WA, 98112, USA
3 Department of Biochemistry and Microbiology, Rutgers University, New Brunswick, NJ 08901 USA
4 Department of Marine and Environmental Sciences, Northeastern University, 430 Nahant Road, Nahant, MA 01908
5 Florida International University, Environmental Epigenetics Laboratory, Institute of Environment 3000 NE 151 St. North Miami, FL, 33181, USA
6 Center for Coastal Marine Sciences, California Polytechnic State University, San Luis Obispo, CA, 93407, USA
7 California Sea Grant, University of California San Diego, La Jolla, CA, 92093
8 Department of Biological Sciences, University of Rhode Island, Kingston, RI 02881, USA
+corresponding author: hputnam@uri.edu
* equal contribution

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
    - [intermediate-files](https://github.com/hputnam/Meth_Compare/tree/master/output/intermediate-files)
         - _subdirectories_
    - [figures](https://github.com/hputnam/Meth_Compare/tree/master/output/figures): Figures in the manuscript main text
    - [supplemental-material](https://github.com/hputnam/Meth_Compare/tree/master/output/supplemental-material): supplemental material referenced in manuscript

## Abstract

*Background*
Elucidating genome to phenome pathways is essential for our understanding of biology, ecology, and evolution. There has been a growing focus on gene expression regulation through DNA methylation in marine invertebrates that need to rapidly respond to changing environmental factors and anthropogenic impacts. However, limited understanding of potential methodological biases hampers genome-wide DNA methylation studies in non-model organisms.

*Results*
We compare three methods for quantifying DNA methylation at single base pair resolution — Whole Genome Bisulfite Sequencing (WGBS), Reduced Representation Bisulfite Sequencing (RRBS), and Methyl-CpG Binding Domain Bisulfite Sequencing (MBDBS) — using multiple individuals from two reef-building coral species with contrasting environmental sensitivity. All methods reveal substantially greater methylation in *M. capitata* (11.4%) than the more sensitive *P. acuta* (2.9%), and the majority of CpG methylation in both species occurs in gene bodies and flanking regions. In both species, MBDBS has the greatest capacity for detecting CpGs in coding regions at our sequencing depth, however MBDBS may be limited by intra-sample methylation heterogeneity. RRBS yields robust information for specific loci albeit without enrichment of any particular genome feature and with significantly reduced genome coverage.

*Conclusions*
Relative genome size strongly influences the number and location of CpGs detected by each method when sequencing depth is limited, providing nuances in cross-species comparisons. These findings reinforce the role and importance of DNA methylation underlying environmental sensitivity in critical marine invertebrate taxa, and provide a resource for important considerations when examining the functional role of DNA methylation.



## Results

- [**Figure 1**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_1.jpg): Experimental design. Three biological replicate coral samples were obtained from both *M. capitata* and *P. acuta*. DNA was extracted from each coral sample and split for use in Whole Genome Bisulfite Sequencing (WGBS), Reduced Representation Bisulfite Sequencing (RRBS), and Methyl-CpG Binding Domain Bisulfite Sequencing (MBDBS) library preparation methods. Three libraries were generated for each of the three methods, yielding nine libraries for each species and 18 libraries total. ([code here]())
- [**Figure 2**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_2.jpg): Circos plots showing the mean percent methylation of CpGs with 5x cover for each of the three methods on the largest scaffolds of each genome. The outer track shows the scaffold locations and dots indicate the percent methylation as indicated by the y-axes from 0-100% for each of the inner tracks. ([code here]())
- [**Figure 3**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_3.jpg): Matrix of pairwise scatter plots for  shared CpG loci (i.e., CpG covered at > 5x across all samples) for (**A**) *M. capitata* (n = 4,666 common loci) and (**B**) *P. acuta* (n = 93,714 common loci). The red lines represent linear regression fits and the green lines are polynomial regression fits. Pearson correlation coefficients for each pairwise comparison are presented in the upper right boxes. Methods are color coded on the X and Y axes (WGBS = green, MBDBS = purple, and RRBS = orange) and replicate samples are indicated on the diagonal along with histograms of % CpG methylation. ([code here]())
- [**Figure 4**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_4.jpg): CpG site coverage across library preparation methods. Mean fraction of CpG sites in the genome covered at different sequencing depths (read depths) by (**A**) MBDBS libraries, (**B**) RRBS libraries, and (**C**) WGBS libraries with standard deviations shown by shaded areas (see [Supplementary Table 2](https://github.com/hputnam/Meth_Compare/tree/master/output/supplemental-material) for number of reads in each sample). ([code here]())
- [**Figure 5**](https://github.com/hputnam/Meth_Compare/blob/master/output/figures/Fig_5.jpg): Proportion of CpGs detected by sequencing methods in coding sequences (CDS), introns, 1 Kb flanking regions upstream (Upstream Flank) or downstream of genes (Downstream Flank), and intergenic regions for (**A**) *M. capitata* and (**B**) *P. acuta*. Principal Coordinate Analyses associated with perMANOVA and beta-dispersion tests related to [Supplementary Table 6](https://github.com/hputnam/Meth_Compare/tree/master/output/supplemental-material) that show differences in proportion of CpGs in various genomic locations (CDS, introns, upstream flanks, downstream flanks, and intergenic regions) for (**C**) *M. capitata* and **D**) *P. acuta*. WGBS is represented by green circles, RRBS by purple triangles, and MBDBS by orange diamonds. Percent variation explained by each PCoA axis is included in the axis label. Ellipses depict 95% confidence intervals for each sequencing method. All eigenvectors are significant at the α = 0.05 level. ([code here]())
