#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. htseqCount variables declaration
##################################################################################
my $database; #GFF file used to calculate the number of reads in featureCounts analysis
my $logfile; #Run.log file where execution data will be saved
my $verbose; #Option to show the execution data on the screen
my $projectdir; #Directory where htseq_results directory will be created
my $bwa_dir; #Path of the directory with the bowtie2 results
my $threads; #Optional number of threads to perform the analysis
my $miARmaPath;#Path to software
my $Seqtype;
my $bwaindex; #Indexed genome to align your reads (Mandatory for analysis with bwa)

BEGIN{
	$miARmaPath="../../../../";#Path to software. Full path is recommended
	$database="../../../basic_examples/circRNAs/data/Homo_sapiens_GRCh37.74_chr.gtf"; #GFF file used to calculate the number of reads in featureCounts analysis
	$logfile="/run_".$$.".log"; #Run.log file where execution data will be saved
	$verbose=""; #Option to show the execution data on the screen
	$projectdir="."; #Directory where htseq_results directory will be created
	$bwa_dir="../2.Read_alignment/bwa_results/"; #Path of the directory with the bowtie1 results
	$threads=4; #Optional number of threads to perform the analysis
	$Seqtype="Paired";
	$bwaindex="../../../../Genomes/Indexes/BWA/human/hg19.fa"; #Indexed genome to align your reads (Mandatory for analysis with bwa)
}

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::Readcount;

##################################################################################
#2. Reading the input directory and executing htseqCount function
##################################################################################
#Reading CutAdapt results directory, collecting the files and completing with the path
opendir(CUTDIR, $bwa_dir) || die $!; 
my @bwa_files= readdir(CUTDIR);
@bwa_files=map("$bwa_dir$_",@bwa_files);

my @circRNAfiles;

# HTSEQCOUNT EXECUTION
#Reading the array with the path of the files
foreach my $file(@bwa_files){
	#Selecting only the sam files for their processing
	my $result=CIRICount(
	  	file=>$file,
	  	database=>$database,
	  	logfile=>$projectdir.$logfile,
		verbose=>$verbose, 
	  	projectdir=>$projectdir,
		threads=>$threads,
		miARmaPath=>$miARmaPath,
		Seqtype=>$Seqtype,
		bwaindex=>$bwaindex,
	);
	push(@circRNAfiles, $result);
}

CIRIFormat( 
  	input=>\@circRNAfiles, 
  	projectdir=>$projectdir,
  	logfile=>$projectdir.$logfile
);
