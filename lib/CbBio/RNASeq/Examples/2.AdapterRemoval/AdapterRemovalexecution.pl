#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Adapt;


#Defining AdapterRemoval variables
my $dir="/Users/apple/Pipeline/examples/Cutadapt/reads";
my $projectdir="/Users/apple/Pipeline/examples/AdapterRemoval";
my $logfile="/run_".$$.".log";
my $statsfile="/stats_".$$.".log";
#my $organism="human";
#my $adaptpredictionnumber="2";
my $verbose;
my $adapter="ATCTCGTATGCCGTCTTCTGCTTGAA";
my $adaptersoft="Reaper-Cutadapt";
my $reaperparameters="-3p-prefix 12/2/0/0 -dust-suffix-late 20 --nozip";
my $min="18";
my $max="26";
my $min_quality="25";

#Declaring variables to keep new files
my @results;

#Reading the directory 
opendir(DIR, $dir) || die "$! : $dir"; 
my @files= readdir(DIR); 

# ADAPTERREMOVAL EXECUTION
#Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fq\.gz$/){
		@results=AdapterRemoval(
			adaptersoft=>$adaptersoft,
			dir=>$dir,
			file=>$file,
			adapter=>$adapter,
			#organism=>$organism,
			logfile=>$projectdir.$logfile,
			statsfile=>$projectdir.$statsfile,
			verbose=>$verbose,
			#adaptpredictionnumber=>$adaptpredictionnumber,
			projectdir=>$projectdir,
			min=>$min,
			max=>$max,
			min_quality=>$min_quality,
			reaperparameters=>$reaperparameters
		);
	}
}

sub date{
	my $dt = DateTime->now;
	return($dt->hms . " [" . $dt->dmy ."]");
}	
