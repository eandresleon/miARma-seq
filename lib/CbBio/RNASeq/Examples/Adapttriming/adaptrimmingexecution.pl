#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Adapt;

#Declaring AdaptTriming variables
my $dir="/Users/apple/Pipeline/examples/Cutadapt/reads";
my $projectdir="/Users/apple/Pipeline/examples/Adapttriming";
my $trimmingnumber="12";
my $readposition="3";
my $readstart="\@SRR";
my $logfile="/run_".$$.".log";
my $verbose="verbose";


#Reading the directory
opendir(DIR, $dir) || die $!; 
my @files= readdir(DIR);

# ADAPTTRIMING EXECUTION
# Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fq\.gz$/){
		AdaptTriming(
			verbose=>$verbose,
			dir=>$dir,
			file=>$file,
			projectdir=>$projectdir,
			trimmingnumber=>$trimmingnumber,
			readposition=>$readposition,
			logfile=>$dir.$logfile,
			readstart=>$readstart
		);
	}
}

sub date{
	my $dt = DateTime->now;
	return($dt->hms . " [" . $dt->dmy ."]");
}	