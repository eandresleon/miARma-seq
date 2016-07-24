#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. ReadAlignment variables declaration
##################################################################################
my $cut_dir; #Path of the directory with the cutadapt results
my $rea_dir; #Path of the directory with the reaper results
my $projectdir; #Directory to save the results
my $threads; #Optional number of threads to perform the analysis
my $bowtie1index; #Indexed genome to align your reads in format .ebwt (Mandatory for analysis with bowtie1)
my $bowtie2index; #Indexed genome to align your reads in format .bt2 (Mandatory for analysis with bowtie2)
my $verbose; #Option to show the execution data on the screen   
my $logfile; #Run.log file where execution data will be saved
my $bowtiemiss; #Max # mismatches in seed alignment in bowtie analysis (0-1)
my $bowtieleng; #Length of seed substrings in bowtie analysis (>3, <32)
my $statsfile; #Stats.log file where stats data will be saved
my $bowtie1parameters; #Other bowtie parameters to perform the analysis using the bowtie1 recommended syntaxis
my $aligner; #Aligner which will be use in the analysis (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1)
my $miARmaPath;#Path to software
my $summary_file;

BEGIN{
	$miARmaPath="../../../../"; #Path to software: Full path is recommended
	$cut_dir="../2.Adapter/cutadapt_results/"; #Path of the directory with the cutadapt results
	$rea_dir="../2.Adapter/Reaper_results/"; #Path of the directory with the reaper results
	$projectdir="."; #Directory to save the results
	$threads="2"; #Optional number of threads to perform the analysis
	$bowtie1index="../../../Genomes/Indexes/bowtie1/human/bw1_homo_sapiens19"; #Indexed genome to align your reads in format .ebwt (Mandatory for analysis with bowtie1)
	$bowtie2index="../../../Genomes/Indexes/bowtie2/human/bw2_homo_sapiens19"; #Indexed genome to align your reads in format .bt2 (Mandatory for analysis with bowtie2)
	$verbose=""; #Option to show the execution data on the screen   
	$logfile="/run_".$$.".log"; #Run.log file where execution data will be saved
	$bowtiemiss="0"; #Max # mismatches in seed alignment in bowtie analysis (0-1)
	$bowtieleng="19"; #Length of seed substrings in bowtie analysis (>3, <32)
	$statsfile="/stats_".$$.".log"; #Stats.log file where stats data will be saved
	$bowtie1parameters=" --best --nomaqround -e 70 -k 1"; #Other bowtie parameters to perform the analysis using the bowtie1 recommended syntaxis
	$aligner="Bowtie1-Bowtie2"; #Aligner which will be use in the analysis (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1)
	$summary_file="Summary_result.xls"; #Variable to save output directory from FastQC function to be use by FastQCStats function
}

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::Aligner;

##################################################################################
#2. Reading the input directory and executing ReadAlignment function
##################################################################################

#Reading CutAdapt results directory, collecting the files and completing with the path
opendir(CUTDIR, $cut_dir) || die $!; 
my @cut_files= readdir(CUTDIR);
@cut_files=map("$cut_dir$_",@cut_files);

#Reading Reaper results directory, collecting the files and completing with the path
opendir(READIR, $rea_dir) || die $!; 
my @rea_files= readdir(READIR);
@rea_files=map("$rea_dir$_",@rea_files);

my @files=(@cut_files, @rea_files);

#READALIGNMENT EXECUTION 
# Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	ReadAligment(
		file=>$file,
		aligner=>$aligner,
		threads=>$threads,
		bowtie2index=>$bowtie2index,
		bowtie1index=>$bowtie1index,
		verbose=>$verbose,
		logfile=>$projectdir.$logfile,
		statsfile=>$projectdir.$statsfile,
		bowtiemiss=> $bowtiemiss,
		bowtieleng=>$bowtieleng,
		projectdir=>$projectdir,
		bowtie1parameters=>$bowtie1parameters,
		miARmaPath=>$miARmaPath,
		summary=>$summary_file,
	);
}