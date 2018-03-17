# miARma WebPage #

miARma web page that includes guides, installation details and examples has been moved to [http://miarmaseq.com](http://miarmaseq.com)

# miARma #

miARma is a fully customizable pipeline for NGS transcriptome analyses. Including gene/transcripts, miRNAs and circRNAs expression measurements.
Created at Computational Biology and Bioinformatics Group (CbBio)
Institute of Biomedicine of Seville. IBIS (Spain)
Modified and Updated at Bioinformatics Unit at IPBLN-CSIC (Institue for Parasitology and Biomedicine Lopez-Neyra, CSIC).
Granada (Spain). 
Copyright (c) 2017 IBIS & IPBLN. All rights reserved.

### miARma 1.7.2 release (18/Dec/17) ###
Minor bugs fixed.
 
 * New Ensembl BiomaRt URL used
 * Order columns from ReadCount section without checking samples names
 * Fixed a bug in the TargetPrediction
 
### miARma 1.7.1 release (21/Aug/17) ###
Minor bugs fixed eg.
 * No aligned reads in hisat2 paired end analysis added.

Added stuff:
* Unaligned files are compressed
* miRDeeparam added to include parameters to miRDeep execution


### miARma 1.7.0 release (09/Aug/17) ###
 * [Hisat v2.1.0](https://ccb.jhu.edu/software/hisat2/index.shtml) has been included as a mRNA aligner.
 * [STAR v020201](https://github.com/alexdobin/STAR/) has been included as a mRNA aligner.

### 1. Included in miARma ###

#### Quality Software ####
* [Fastqc v0.11.5](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
#### Trimming Softare ####
* [CutAdapt v1.9.1](https://cutadapt.readthedocs.org/en/stable/)
* [Minion v15-065](ftp://ftp.ebi.ac.uk/pub/contrib/enrightlab/kraken/reaper/src/reaper-latest/doc/minion.html)
* [Reaper v15-065](http://www.ebi.ac.uk/~stijn/reaper/reaper.html)
#### Aligners ####
* [Bowtie v1.1.2](http://bowtie-bio.sourceforge.net/index.shtml)
* [Bowtie v2.2.8](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [TopHat v2.1.1](http://ccb.jhu.edu/software/tophat/index.shtml)
* [BWA v0.7.13](http://bio-bwa.sourceforge.net/)
* [Hisat v2.1.0](https://ccb.jhu.edu/software/hisat2/index.shtml)
* [STAR v020201](https://github.com/alexdobin/STAR/)
#### Entity Quantification####
* [feactureCounts v1.5.0-p1](http://bioinf.wehi.edu.au/featureCounts/)
* [Ciri v1.2](http://sourceforge.net/projects/ciri/files/?source=navbar)
* [Ciri v2.0.1](http://sourceforge.net/projects/ciri/files/?source=navbar)
* [miRDeep v2](https://www.mdc-berlin.de/8551903/en/)
#### Others ####
* [RNAfold v2.2.4](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)
* [Samtools v1.3](http://samtools.sourceforge.net/)

### 2. Pre-requisites ###

miARma-Seq is a tool that provides an easy and common interface to various analysis software. It also intends to reduce to the minimum the number of dependencies. Nevertheless, some basic programs listed below must be correctly installed:

* [Perl v5.6.0 or higher.](http://www.cpan.org/src/5.0/perl-5.6.1.tar.gz)
* [Java JDK v.1.6. or higher.](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)
* [R environment v.3.2 or higher.](http://www.r-project.org/)
* [Bioconductor v.1.3 or higher.](https://www.bioconductor.org/install/)

#### Compilers: #####
+ Apple:
    - [Xcode](https://itunes.apple.com/es/app/xcode/id497799835?l=en&mt=12)
+ Linux:
    - [Gcc](https://ftp.gnu.org/gnu/gcc/)
    - [make](https://ftp.gnu.org/gnu/make/)

### How do I get set up? ###

* [miARma can be installed using the following guide](http://idproteins.com/installation)


### Guidelines/How to ###

* [miRNAs guide](http://miarmaseq.idproteins.com/Documentation)
* [mRNAs guide](http://miarmaseq.idproteins.com/Documentation)
* [circRNAs guide](http://miarmaseq.idproteins.com/Documentation)

### Code Documentation ####
* [Perldoc](http://miarmaseq.idproteins.com/PDoc/)

### Who do I talk to? ###

* miARma-devel@idoproteins.com