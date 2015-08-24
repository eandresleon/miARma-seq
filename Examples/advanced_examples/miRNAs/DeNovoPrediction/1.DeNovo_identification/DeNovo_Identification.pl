#!/usr/bin/perl

use strict;

##################################################################################
#1. ReadAlignment variables declaration
##################################################################################
my $read_dir; #Path of the directory with the cutadapt results
my $projectdir; #Directory to save the results
my $threads; #Optional number of threads to perform the analysis
my $bowtie1index; #Indexed genome to align your reads in format .ebwt (Mandatory for analysis. made with bowtie1)
my $verbose; #Option to show the execution data on the screen   
my $logfile; #Run.log file where execution data will be saved
my $statsfile; #Stats.log file where stats data will be saved
my $miARmaPath;#Path to software
my $adapter; #Adapter to trimm at read 3'
my $aligner; #Aligner which will be use in the analysis (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1 and miRDeep)
my $mature_miRNA_file; #a fasta file with all mature sequence from your organism
my $precursor_miRNA_file; #a fasta file with all known pre-miRNa sequence 
my $genome; #fasta file for the cmplete genome of our organism
BEGIN{
	$miARmaPath="../../../";#Path to software. Full path is recommended
	$read_dir="../../reads/"; #Path of the directory with the cutadapt results
	$projectdir="."; #Directory to save the results
	$threads="2"; #Optional number of threads to perform the analysis
	$bowtie1index="../../../Genomes/Indexes/bowtie1/human/hg19"; #Indexed genome to align your reads in format .ebwt (Mandatory for analysis with bowtie1)
	$verbose=""; #Option to show the execution data on the screen   
	$logfile="/run_".$$.".log"; #Run.log file where execution data will be saved
	$adapter="ATCTCGTATGCCGTCTTCTGCTTGAA";#Adapter to trimm at read 3'
	$aligner="miRDeep"; #Aligner which will be use in the analysis (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1 and miRDeep)
	$mature_miRNA_file="hsa_mature_miRBase20.fasta"; #a fasta file with all mature sequence from your organism
	$precursor_miRNA_file="precursors_miRBase20.fasta"; #a fasta file with all known pre-miRNa sequence 
	$statsfile="/stats_".$$.".log"; #Stats.log file where stats data will be saved
	$genome="../../../Genomes/Fasta/hg19.fa"; #fasta file for the cmplete genome of our organism
}

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::Aligner;

##################################################################################
#2. Reading the input directory and executing ReadAlignment function
##################################################################################


# Reading reads directory 
opendir(DIR, $read_dir) || die "$! : $read_dir"; 
my @files= readdir(DIR); 
@files=map("$read_dir$_",@files);

#READALIGNMENT EXECUTION 
# Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	ReadAligment(
		file=>$file,
		aligner=>$aligner,
		threads=>$threads,
		bowtie1index=>$bowtie1index,
		verbose=>$verbose,
		logfile=>$projectdir.$logfile,
		statsfile=>$projectdir.$statsfile,
		projectdir=>$projectdir,
		miARmaPath=>$miARmaPath,
		adapter=>$adapter,
		precursors=>$precursor_miRNA_file,
		mature=>$mature_miRNA_file,
		genome=>$genome,
	);
}
