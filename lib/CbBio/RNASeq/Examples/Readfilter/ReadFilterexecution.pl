#!/usr/bin/perl


use strict;
use warnings;
use CbBio::RNASeq::Adapt;

#Declaring ReadFilter variables
my $dir="/Users/apple/Pipeline/examples/Cutadapt/cutadapt_results";
my $projectdir="/Users/apple/Pipeline/examples/Readfilter";
my $logfile="/run_".$$.".log";
my $minreadlength="18";
my $maxreadlength="26";
my $readstart="\@SRR";
my $verbose="verbose";

#This variable collect the name of the output file and will be saved in @readfilterresults
my $results;
my @readfilterresults;

#Opening the directory and reading the files 
opendir(DIR, $dir) || die $!; 
my @files= readdir(DIR); 

# Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	if($file =~ /.*\.fastq$/){
		$results=ReadFilter(
			dir=>$dir,
			projectdir=>$projectdir,
			minreadlength=>$minreadlength,
			maxreadlength=>$maxreadlength,
			readstart=>$readstart,
			logfile=>$projectdir.$logfile,
			file=>$file,
			verbose=>$verbose		
		);
		push (@readfilterresults, $results);
	}
}

sub date{
	my $dt = DateTime->now;
	return($dt->hms . " [" . $dt->dmy ."]");
}	