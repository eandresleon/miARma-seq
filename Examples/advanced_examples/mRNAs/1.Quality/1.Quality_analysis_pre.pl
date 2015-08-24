#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. FastQC and FastQCStats variables declaration
##################################################################################
my $miARmaPath;
my $dir;
my $projectdir;
my $threads;
my $verbose;
my $prefix;
my $logfile;
my $statsfile;
my $output_dir;

BEGIN { 
	$miARmaPath="../../../../"; #Path to software: Full path is recommended
	$dir="$miARmaPath/Examples/basic_examples/mRNAs/reads/"; #Path of the directory with the input files
	$projectdir="."; #Directory to save the results
	$threads="4"; #Optional number of threads to perform the analysis faster
	$verbose=1; #Optional argument to show the execution data on screen
	$prefix="Pre"; #Label to write in the directory results name
	$logfile="/run_".$$.".log"; #Log file where execution data will be saved
	$statsfile="/stats_".$$.".log"; #Stats file where stats data will be saved
	$output_dir; #Variable to save output directory from FastQC function to be use by FastQCStats function
}
##################################################################################
#2. Reading the input directory and executing FastQC and FastQCStats functions 
##################################################################################

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::Quality;

# Reading the directory collecting the files and completing with the path
opendir(DIR, $dir) || die "$! ($dir)"; 
my @files= readdir(DIR);
@files=map("$dir$_",@files); 

# FASTQC EXECUTION
# Reading the array with the names of the files
foreach my $file(@files){
	#Calling FastQC subroutine of Quality.pm package.
	$output_dir=FastQC(
		miARmaPath=>$miARmaPath,
		file=>$file,
		projectdir=>$projectdir,
		threads=>$threads,
		verbose=>$verbose,
		logfile=>$projectdir.$logfile,
		prefix=>$prefix
	);
}

# FASTQCSTATS EXECUTION
# Calling FastQCStats sobroutine of Quality.pm package. 

FastQCStats(
	dir=>$output_dir, 
	verbose=>$verbose,
	statsfile=>$projectdir.$statsfile,
	logfile=>$projectdir.$logfile
);

