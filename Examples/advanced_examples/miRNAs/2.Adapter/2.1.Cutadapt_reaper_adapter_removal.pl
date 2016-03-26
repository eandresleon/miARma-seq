#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. AdapaterRemoval variables declaration
##################################################################################
my $dir; #Path of the directory with th input files
my $projectdir; #Directory to save the results
my $logfile; #Log file where execution data will be saved
my $statsfile; #Stats file where stats data will be saved
my $verbose; #Optional arguments to show the execution data on screen
my $adapter; #Adapter sequence to be removed in the analysis
my $adaptersoft; #Specific software to remove the adapter from the sequences 
my $reaperparameters; #Parameters to execute reaper.
my $min; #Minimun length of the sequence read to keep with Cutadapt and Reaper Software 
my $max; #Maximun length of the sequence read to use with Cutadapt Software
my $min_quality; #Minimun quality of the sequence read to use with Cutadapt Software
my $miARmaPath;#Path to software
my $summary_file;
BEGIN{
	$miARmaPath="../../../../"; #Path to software: Full path is recommended
	$dir="$miARmaPath/Examples/basic_examples/miRNAs/reads/"; #Path of the directory with the input files
	$projectdir="."; #Directory to save the results
	$logfile="/run_".$$.".log"; #Log file where execution data will be saved
	$statsfile="/stats_".$$.".log"; #Stats file where stats data will be saved
	$verbose=""; #Optional arguments to show the execution data on screen
	$adapter="ATCTCGTATGCCGTCTTCTGCTTGAA"; #Adapter sequence to be removed in the analysis
	$adaptersoft="Reaper-Cutadapt"; #Specific software to remove the adapter from the sequences 
	$reaperparameters="-3p-prefix 12/2/0/0 -dust-suffix-late 20"; #Parameters to execute reaper.
	$min="18"; #Minimun length of the sequence read to keep with Cutadapt and Reaper Software 
	$max="26"; #Maximun length of the sequence read to use with Cutadapt Software
	$min_quality="25"; #Minimun quality of the sequence read to use with Cutadapt Software
	$summary_file="Summary_result.xls"; #Variable to save output directory from FastQC function to be use by FastQCStats function
}

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::Adapt;

##################################################################################
#2. Reading the input directory and executing AdapterRemoval function
##################################################################################

# Reading the directory 
opendir(DIR, $dir) || die "$! : $dir"; 
my @files= readdir(DIR); 

# ADAPTERREMOVAL EXECUTION
AdapterRemoval(
	adaptersoft=>$adaptersoft,
	dir=>$dir,
	files=>\@files,
	adapter=>$adapter,
	logfile=>$projectdir.$logfile,
	statsfile=>$projectdir.$statsfile,
	verbose=>$verbose,
	projectdir=>$projectdir,
	min=>$min,
	max=>$max,
	min_quality=>$min_quality,
	reaperparameters=>$reaperparameters,
	miARmaPath=>$miARmaPath,
	summary=>$summary_file,
);
