#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Adapt;

#Declaring CutAdapt variables 
my $dir="/Users/apple/Pipeline/examples/Cutadapt/reads";
my $projectdir="/Users/apple/Pipeline/examples/Cutadapt";
my $verbose;
my $logfile="/run_".$$.".log";
my $adapter="ATCTCGTATGCCGTCTTCTGCTTGAA";
my $min="18";
my $max="26";
my $min_quality="25";

#Declaring CutAdaptStats variables 
my $readstart="\@SRR87";
my $statsfile="/stats_".$$.".log";

#This variable collect the name of the output file and will be saved in @results
my $result;
my @results;

#Reading the directory
opendir(DIR, $dir) || die $!; 
my @files= readdir(DIR);

# CUTADAPT EXECUTION
# Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the fastq files for their processing
	if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fq\.gz$/){
		#Printing a message on console about the execution time and the process
		print STDERR "LOG :: ".date()." Checking $file for CutAdapt analysis\n";
		#Calling CutAdapt subroutine of Adapt.pm package. 
		$result=CutAdapt(
			file=>$file,
			dir=>$dir,
			projectdir=>$projectdir,
			adapter=>$adapter,
			verbose=>$verbose,
			logfile=>$projectdir.$logfile,
			min=>$min,
			max=>$max,
			min_quality=>$min_quality
		);
		push (@results, $result);
	}
}

#CUTADADPTSTATS EXECUTION
#Calling fastqcStats sobroutine of Quality.pm package. 
CutAdaptStats(
  	dir=>$dir, 
  	readstart=>$readstart,
  	verbose=>$verbose, 
  	projectdir=>$projectdir,
	statsfile=>$projectdir.$statsfile,
	logfile=>$projectdir.$logfile
);

sub date{
	my $dt = DateTime->now;
	return($dt->hms . " [" . $dt->dmy ."]");
}	

