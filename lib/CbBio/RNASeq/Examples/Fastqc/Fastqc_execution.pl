#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Quality;

#Declaring FastQC and FastQCStats variables 
my $dir="/Users/apple/Pipeline/examples/Fastq/reads";
my $projectdir="/Users/apple/Pipeline/examples/Fastq";
my $threads="8";
my $verbose;
my $label="Hypoxia_Pre";
my $logfile="/run_".$$.".log";
my $statsfile="/stats_".$$.".log";
my $output_dir;

#Reading the directory
opendir(DIR, $dir) || die $!; 
my @files= readdir(DIR); 


# FASTQC EXECUTION
# Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fq\.gz$/){
		#Printing a message on console about the execution time and the process
		print STDERR "LOG :: ".date()." Checking $file for fastqc analysis\n";
		#Calling fastqc subroutine of Quality.pm package. 
		$output_dir=FastQC(
			file=>$file,
			dir=>$dir,
			projectdir=>$projectdir,
			threads=>$threads,
			verbose=>$verbose,
			logfile=>$projectdir.$logfile,
			label=>$label
		);
	}
}

#FASTQCSTATS EXECUTION
#Calling fastqcStats sobroutine of Quality.pm package. 
FastQCStats(
			dir=>$output_dir, 
			verbose=>$verbose,
			statsfile=>$projectdir.$statsfile,
			logfile=>$projectdir.$logfile
		);

sub date{
	my $dt = DateTime->now;
	return($dt->hms . " [" . $dt->dmy ."]");
}	