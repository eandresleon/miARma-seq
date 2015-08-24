#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. htseqCount variables declaration
##################################################################################
my $logfile; #Run.log file where execution data will be saved
my $verbose; #Option to show the execution data on the screen
my $projectdir; #Directory where htseq_results directory will be created
my $mirdeep_dir; #Path of the directory with the bowtie2 results
my $miARmaPath;#Path to software

BEGIN{
	$miARmaPath="../../../";#Path to software. Full path is recommended
	$logfile="/run_".$$.".log"; #Run.log file where execution data will be saved
	$verbose=""; #Option to show the execution data on the screen
	$projectdir="."; #Directory where htseq_results directory will be created
	$mirdeep_dir="../1.DeNovo_identification/miRDeep_results/"; #Path of the directory with the bowtie1 results
}

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::Readcount;

##################################################################################
#2. Reading the input directory and executing miRDeepCount function
##################################################################################
#Reading miRDeep results directory, collecting the files and completing with the path
opendir(MIRDIR, $mirdeep_dir) || die $!; 
my @mirdeep_files= readdir(MIRDIR);
@mirdeep_files=map("$mirdeep_dir$_",@mirdeep_files);

my @allRNAfiles;

#Reading the array with the path of the files
foreach my $file(@mirdeep_files){
	#Selecting only the sam files for their processing
	my $result=miRDeepCount(
	  	file=>$file,
	  	logfile=>$projectdir.$logfile,
		verbose=>$verbose, 
	  	projectdir=>$projectdir,
		miARmaPath=>$miARmaPath,
	);
	push(@allRNAfiles, $result) if($result);
}

miRDeepFormat( 
  	input=>\@allRNAfiles, 
  	projectdir=>$projectdir,
  	logfile=>$projectdir.$logfile
);
