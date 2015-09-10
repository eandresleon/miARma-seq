#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. ReadAlignment variables declaration
##################################################################################
my $reads_dir; #Path of the directory with the reaper results
my $projectdir; #Directory to save the results
my $threads; #Optional number of threads to perform the analysis
my $bwaindex; #Indexed genome to align your reads
my $verbose; #Option to show the execution data on the screen   
my $logfile; #Run.log file where execution data will be saved
my $aligner; #Aligner which will be use in the analysis (Allowed values: tophat)
my $miARmaPath;#Path to software
my $statsfile; #Stats.log file where stats data will be saved
my $Seqtype;
BEGIN{
	$miARmaPath="../../../../"; #Path to software: Full path is recommended
	$reads_dir="$miARmaPath/Examples/basic_examples/circRNAs/reads/"; #Path of the directory with the input files
	$projectdir="."; #Directory to save the results
	$threads=4; #Optional number of threads to perform the analysis
	$bwaindex="../../../../Genomes/Indexes/BWA/human/hg19.fa"; #Indexed genome to align your reads (Mandatory for analysis with bwa)
	$verbose=""; #Option to show the execution data on the screen   
	$logfile="/run_".$$.".log"; #Run.log file where execution data will be saved
	$statsfile="/stats_".$$.".log"; #Stats.log file where stats data will be saved
	$aligner="bwa"; #Aligner which will be use in the analysis (Allowed values: tophat)
	$Seqtype="Paired";
}

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::Aligner;

##################################################################################
#2. Reading the input directory and executing ReadAlignment function
##################################################################################

#Reading CutAdapt results directory, collecting the files and completing with the path
opendir(READ_DIR, $reads_dir) || die "$! : Please check the folder where reads are stored\n"; 
my @read_files= readdir(READ_DIR);
@read_files=map("$reads_dir/$_",@read_files);

#READALIGNMENT EXECUTION 
# Reading the array with the names of the files
foreach my $file(@read_files){
	#Selecting only the fastq files for their processing
	ReadAligment(
		file=>$file,
		aligner=>$aligner,
		threads=>$threads,
		bwaindex=>$bwaindex,
		verbose=>$verbose,
		logfile=>$projectdir.$logfile,
		statsfile=>$projectdir.$statsfile,
		projectdir=>$projectdir,
		miARmaPath=>$miARmaPath,
		Seqtype=>$Seqtype,
	);
}
