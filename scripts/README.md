# scripts

Scripts are organized in the order they should be run.

- [Generating-Genome-Feature-Tracks.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Generating-Genome-Feature-Tracks.ipynb): Create gene, CDS, intron, flanking region, and intergenic feature tracks for *M. capitata* and *P. acuta* genomes
- [Characterizing-CpG-Methylation-5x.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Characterizing-CpG-Methylation-5x.ipynb): Characterize methylation status and genomic locations of CpGs in individual sample data for both species using `bedtools`
- [Characterizing-CpG-Methylation-5x-Union.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Characterizing-CpG-Methylation-5x-Union.ipynb): Characterize methylation status and genomic locations of CpGs from 5x union bedgraphs for both species using `bedtools`
- [Characterizing-CpG-Methylation-5x-Union-Summary-Plots.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Characterizing-CpG-Methylation-5x-Union-Summary-Plots.Rmd): Create summary tables and stacked barplots to understand methylation status and genomic locations of 5x union CpGs in both species
- [Identifying-Genomic-Locations.ipynb](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Identifying-Genomic-Locations.ipynb): Characterize genomic locations of CpGs in upset plots and method-associated DMC using `bedtools`
- [Identifying-Genome-Features-Summary-Plots.Rmd](https://github.com/hputnam/Meth_Compare/blob/master/scripts/Identifying-Genome-Features-Summary-Plots.Rmd): Create summary tables and perform statistical tests to compare CpG locations in various method-associated DMC categories