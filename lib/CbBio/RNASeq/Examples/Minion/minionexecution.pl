#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Adapt;


#Defining Minion variables
my $dir="../Cutadapt/cutadapt_results";
my $projectdir=".";
my $logfile="/run_".$$.".log";
my $statsfile="/stats_".$$.".log";
my $organism="human";
my $adaptpredictionnumber="1";
my $verbose;

#This variable collect the name of the predicted adapter
my $adapter;

#Reading the directory 
opendir(DIR, $dir) || die "$! : $dir"; 
my @files= readdir(DIR); 

# MINION EXECUTION
#Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fq\.gz$/){
		$adapter=Minion(
			dir=>$dir,
			file=>$file,
			org=>$organism,
			logfile=>$projectdir.$logfile,
			statsfile=>$projectdir.$statsfile,
			verbose=>$verbose,
			adaptpredictionnumber=>$adaptpredictionnumber
		);
		print "MINION: Predicted adapter for ".$file." : ".$adapter."\n";
	}
}

sub date{
	my $dt = DateTime->now;
	return($dt->hms . " [" . $dt->dmy ."]");
}	


		
