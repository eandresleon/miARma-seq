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
my $adaptersoft; #Specific software to remove the adapter from the sequences
my $trimmingnumber; #Number of nucleotides to remove of the sequence read and the quality data
my $readposition; #End of the read to remove the nucleotides (3 or 5). 
my $miARmaPath;#Path to software

BEGIN{
	$miARmaPath="../../../../"; #Path to software: Full path is recommended
	$dir="$miARmaPath/Examples/basic_examples/miRNAs/reads/"; #Path of the directory with the input files
	$projectdir="."; #Directory to save the results
	$logfile="/run_".$$.".log"; #Log file where execution data will be saved
	$statsfile="/stats_".$$.".log"; #Stats file where stats data will be saved
	$verbose="verbose"; #Optional arguments to show the execution data on screen
	$adaptersoft="Adapttrimming"; #Specific software to remove the adapter from the sequences
	$trimmingnumber="12"; #Number of nucleotides to remove of the sequence read and the quality data
	$readposition="5"; #End of the read to remove the nucleotides (3 or 5). 
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
	logfile=>$projectdir.$logfile,
	statsfile=>$projectdir.$statsfile,
	verbose=>$verbose,
	projectdir=>$projectdir,
	trimmingnumber=>$trimmingnumber,
	readposition=>$readposition
);
