#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Adapt;

#Declaring Reaper variables 
my $dir="/Users/apple/Pipeline/examples/Cutadapt/reads";
my $projectdir="/Users/apple/Pipeline/examples/Reaper";
my $verbose;
my $adapter="ATCTCGTATGCCGTCTTCTGCTTGAA";
my $reaperparameters="-3p-prefix 12/2/0/0 -dust-suffix-late 20 -clean-length 18 --nozip";
my $logfile="/run_".$$.".log";

#Reading the directory
opendir(DIR, $dir) || die $!; 
my @files= readdir(DIR);

# REAPER EXECUTION
# Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fq\.gz$/){
		#Printing a message on console about the execution time and the process
		print STDERR "REAPER :: ".date()." Checking $file for Reaper analysis\n";
		Reaper(
			dir=>$dir,
			file=>$file,
			logfile=>$projectdir.$logfile,
			adapter=>$adapter,
			reaperparameters=>$reaperparameters,
			verbose=>$verbose,
			projectdir=>$projectdir	
		);
	}
}
sub date{
	my $dt = DateTime->now;
	return($dt->hms . " [" . $dt->dmy ."]");
}	