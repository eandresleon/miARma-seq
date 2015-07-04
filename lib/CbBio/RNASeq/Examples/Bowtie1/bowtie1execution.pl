#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Aligner;

#Defining bowtie1 variables
my $dir="/Users/apple/Pipeline/examples/Cutadapt/cutadapt_results";
my $projectdir="/Users/apple/Pipeline/examples/Bowtie1";
my $threads="2";
my $bowtieindex="/Users/apple/Pipeline/examples/Bowtie1/hg19/hg19";
my $verbose="verbose";
my $logfile="/run_".$$.".log"; 
my $bowtiemiss="0";
my $bowtieleng="19";
my $statsfile="/stats_".$$.".log"; 
my $bowtieparameters=" --sam --best --nomaqround -e 70 -k 1";

#Reading the directory 
opendir(DIR, $dir) || die $!; 
my @files= readdir(DIR); 

#BOWTIE1 EXECUTION 
# Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fq\.gz$/){
		bowtie1(
			file=>$file,
			threads=>$threads,
			bowtieindex=>$bowtieindex,
			verbose=>$verbose,
			logfile=>$projectdir.$logfile,
			dir=>$dir,
			statsfile=>$projectdir.$statsfile,
			bowtiemiss=> $bowtiemiss,
			projectdir=>$projectdir,
			bowtieparameters=>$bowtieparameters
		);
	}
}