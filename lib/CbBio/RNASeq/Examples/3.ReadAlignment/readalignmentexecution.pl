#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Aligner;

#Defining bowtie1 variables
my $dir="/Users/apple/Pipeline/examples/Cutadapt/cutadapt_results";
my $projectdir="/Users/apple/Pipeline/examples/ReadAlignment";
my $threads="4";
my $bowtie1index="/Users/apple/Pipeline/examples/Bowtie1/hg19/hg19";
my $bowtie2index="/Users/apple/Pipeline/examples/Bowtie2/hg19/hg19";
my $verbose="verbose";
my $logfile="/run_".$$.".log"; 
my $bowtiemiss="0";
my $bowtieleng="19";
my $statsfile="/stats_".$$.".log"; 
my $bowtie1parameters=" --sam --best --nomaqround -e 70 -k 1";
my $aligner="Bowtie1-Bowtie2";

#Reading the directory 
opendir(DIR, $dir) || die $!; 
my @files= readdir(DIR); 

#BOWTIE1 EXECUTION 
# Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fq\.gz$/){
		ReadAligment(
			file=>$file,
			aligner=>$aligner,
			threads=>$threads,
			bowtie2index=>$bowtie2index,
			bowtie1index=>$bowtie1index,
			verbose=>$verbose,
			logfile=>$projectdir.$logfile,
			dir=>$dir,
			statsfile=>$projectdir.$statsfile,
			bowtiemiss=> $bowtiemiss,
			bowtieleng=>$bowtieleng,
			projectdir=>$projectdir
			bowtie1parameters=>$bowtie1parameters
		);
	}
}