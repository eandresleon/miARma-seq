#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Quality;

my $file="htseqresults.tab"; #name of the tab file which contains the number of reads from the htseq analysis
my $dir="/Users/apple/Pipeline"; #the path of the directory which contains the files. This directory will be configured as working directory for R
my $projectdir="/Users/apple/Pipeline"; #Path of the directory where will be saved the QC report.
my $targetfile="targetdef.txt"; #path of the target file with the experimental conditions of each sample (or name of the file if it is saved in the provided directory)
my $label="Hypoxia_def"; #string character that will appear in the name the results file
my $filter="yes"; #This value refers to filter processing in the reads (Should be "yes" or "no").
my $ref2="T0h"; #condition in the target file of the reference sample.
my $logfile="/run_".$$.".log"; #Path of run.log file where execution data will be saved
my $verbose="verbose";

QC_EdgeR(
			projectdir=>$projectdir,
			dir=> $dir,
			file=>$file,
			targetfile=>$targetfile,
			label=>$label,
			filter=>$filter,
			ref2=>$ref2,
			logfile=>$dir.$logfile,
			verbose=>$verbose
		);