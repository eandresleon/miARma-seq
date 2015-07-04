#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. DEanalysis variables declaration
##################################################################################
my $dir; #Path of the directory with the tab files
my $projectdir; #Path of the directory where the output files will be saved.
my $targetfile; #Path of the target file.
my $contrastfile; #Path of the contrast file.
my $filter; #This value refers to filter processing in the reads (Should be "yes" or "no").
my $DEsoft; #Specific software to perform the Differential Expression Analysis (Allowed values: edger, noiseq or edger-noiseq)
my $filtermethod; #Method that will be used to filter proccess with Noiseq software. Method must be one of 1,2 or 3 for CPM method, Wilcoxon method and Proportion test, respectively (by default CPM method)
my $cpmvalue; #Cutoff for the counts per million value to be used in filter processing with methods 1 and 3 with Noiseq software and in filter processing with EdgeR (1 cpm by default).
my $noiseq_normethod; #Normalization method to perform the DE analysis with Noiseq. It can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
my $edger_normethod; #Normalization method to perform the DE analysis with EdgeR. It can be one of "TMM" (default), "RLE", "upperquartile" or "none" (no normalization).
my $replicates; #Value to indicate if replicates samples are present in the analysis to perform the DE analysis with EdgeR. It can be "yes" (by default) or "no".
my $logfile; #Run.log file where execution data will be printed
my $verbose; #Optional argument to show the execution data on screen
my $repthreshold; #Number of replicates that have to contains at least a defined number of reads per million to perform the filter with EdgeR software(2 replicates by default)
my $miARmaPath;#Path to software

BEGIN{
	$miARmaPath="../../../";#Path to software. Full path is recommended
	$dir="../2.Read_count/miRNAs_results/"; #Path of the directory with the tab files
	$projectdir="."; #Path of the directory where the output files will be saved.
	$targetfile="targets.txt"; #Path of the target file.
	$contrastfile="contrast.txt"; #Path of the contrast file.
	$filter="yes"; #This value refers to filter processing in the reads (Should be "yes" or "no").
	$DEsoft="EdgeR-Noiseq"; #Specific software to perform the Differential Expression Analysis (Allowed values: edger, noiseq or edger-noiseq)
	$filtermethod=1; #Method that will be used to filter proccess with Noiseq software. Method must be one of 1,2 or 3 for CPM method, Wilcoxon method and Proportion test, respectively (by default CPM method)
	$cpmvalue=2; #Cutoff for the counts per million value to be used in filter processing with methods 1 and 3 with Noiseq software and in filter processing with EdgeR (1 cpm by default).
	$noiseq_normethod="rpkm"; #Normalization method to perform the DE analysis with Noiseq. It can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
	$edger_normethod="upperquartile"; #Normalization method to perform the DE analysis with EdgeR. It can be one of "TMM" (default), "RLE", "upperquartile" or "none" (no normalization).
	$replicates="yes"; #Value to indicate if replicates samples are present in the analysis to perform the DE analysis with EdgeR. It can be "yes" (by default) or "no".
	$logfile="/run_".$$.".log"; #Run.log file where execution data will be printed
	$verbose=""; #Optional argument to show the execution data on screen
	$repthreshold=2; #Number of replicates that have to contains at least a defined number of reads per million to perform the filter with EdgeR software(2 replicates by default)
}
use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::DEAnalysis;

##################################################################################
#2. Reading the input directory and executing DE_Analysis function
##################################################################################
#Reading HTseq results directory, collecting the files and completing with the path
opendir(DIR, $dir) || die $!;
my @files= readdir(DIR);
@files=map("$dir$_",@files);

foreach my $file(@files){
	DE_Analysis(
		projectdir=>$projectdir,
		dir=>$dir,
		file=>$file,
		targetfile=>$targetfile,
		label=>"test",
		filter=>$filter,
		edger_contrastfile=>$contrastfile,
		noiseq_contrastfile=>$contrastfile,
		DEsoft=>$DEsoft,
		filtermethod=>$filtermethod,
		logfile=>$projectdir.$logfile,
		verbose=>$verbose,
		cpmvalue=>$cpmvalue,
		repthreshold=>$repthreshold,
		edger_normethod=>$edger_normethod,
		noiseq_normethod=>$noiseq_normethod,
		replicates=>$replicates,
		miARmaPath=>$miARmaPath
	);
}