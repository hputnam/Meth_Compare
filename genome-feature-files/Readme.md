# genome-feature-files

Folder with genome feature tracks for *M. capitata* and *P. acuta* used in downstream analyses.

## *M. capitata*

**Genome information**:

- Genome: `https://gannet.fish.washington.edu/seashell/bu-mox/data/froger/Mcap_Genome/` and `http://gannet.fish.washington.edu/seashell/snaps/Mcap.genome_assembly.fa`
- Genome index: `http://gannet.fish.washington.edu/seashell/snaps/Mcap.genome_assembly.fa.fai` or `https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.genome_assembly.fa.fai`
- Annotated GFF: `https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.GFFannotation.gff`
- All CpG locations: `https://osf.io/wz9j2/`
- CpG location index: `https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap_CpG.gff.idx`

**Genome feature tracks**:

All tracks were derived from the genome in [this script](https://github.com/hputnam/Meth_Compare/blob/master/code/Generating-Genome-Feature-Tracks.ipynb).

- [Genes](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.GFFannotation.gene.gff)
- [Coding Sequences (CDS)](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.GFFannotation.CDS.gff)
- [Introns](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.GFFannotation.intron.gff)
- Flanking regions
	- [All flanking regions](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.GFFannotation.flanks.gff)
	- [Upstream flanking regions](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.GFFannotation.flanks.Upstream.gff)
	- [Downstream flanking regions](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.GFFannotation.flanks.Downstream.gff)
- [Intergenic regions](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.GFFannotation.intergenic.bed)

**Other files**:

- [Mcap.genome_assembly-sequence-lengths.txt](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap.genome_assembly-sequence-lengths.txt): List of chromosome names and lengths used by `bedtools`
- [Mcap-Genome-Feature-Tracks.xml](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Mcap-Genome-Feature-Tracks.xml): IGV session used to verify genome feature tracks.

## *P. acuta*

**Genome information**:

- Genome: `https://gannet.fish.washington.edu/seashell/bu-mox/data/froger/Pact_Genome/` and `https://gannet.fish.washington.edu/seashell/snaps/Pocillopora_acuta_genome_v1.fasta`
- All CpG locations: `https://osf.io/r94xz/`
- CpG location index: `https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact_CpG.gff.idx`

**Genome feature tracks**:

All tracks were derived from the genome in [this script](https://github.com/hputnam/Meth_Compare/blob/master/code/Generating-Genome-Feature-Tracks.ipynb).

- [Genes](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact.GFFannotation.Genes.gff)
- [Coding Sequences (CDS)](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact.GFFannotation.CDS.gff)
- [Introns](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact.GFFannotation.Intron.gff)
- Flanking regions
	- [All flanking regions](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact.GFFannotation.flanks.gff)
	- [Upstream flanking regions](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact.GFFannotation.flanks.Upstream.gff
	- [Downstream flanking regions](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact.GFFannotation.flanks.Downstream.gff)
- [Intergenic regions](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact.GFFannotation.intergenic.bed)

**Other files**:

- [Pact.genome_assembly-sequence-lengths.txt](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact.genome_assembly-sequence-lengths.txt): List of chromosome names and lengths used by `bedtools`
- [Pact-Genome-Feature-Tracks.xml](https://github.com/hputnam/Meth_Compare/blob/master/genome-feature-files/Pact-Genome-Feature-Tracks.xml): IGV session used to verify genome feature tracks.