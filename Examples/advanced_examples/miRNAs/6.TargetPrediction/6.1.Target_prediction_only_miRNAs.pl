#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. htseqCount variables declaration
##################################################################################
my $organism;#="hattribute to be used as feature ID (default: gene_id) for featureCounts analysis
my $logfile; #Run.log file where execution data will be saved
my $miRNAs; #Format of the input file : miRNA or gene (default:sam)
my $UTRs;
my $verbose; #Option to show the execution data on the screen
my $projectdir; #Directory where htseq_results directory will be created
my $miARmaPath;#Path to software
my $DE_folder_miRNAs;
my $noiseq_cutoff;
my $edger_cutoff;

BEGIN{
	$miARmaPath="../../";#Path to software. Full path is recommended
	$organism="human";
	$logfile="/run_".$$.".log"; #Run.log file where execution data will be saved
	$DE_folder_miRNAs="../7.Differential_expression/";# File from edgeR or Noiseq with differentually expressed miRNAs
	$verbose="verbose"; #Option to show the execution data on the screen
	$projectdir="."; #Directory where htseq_results directory will be created
	$noiseq_cutoff=0.8; #Optional argument to select statistically significant results. 0.8 as default
	$edger_cutoff=0.05; #Optional argument to select statistically significant results. 0.05 as default
}

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::TargetPrediction;

##################################################################################
#2. Reading the input files
##################################################################################

TargetPrediction(
	logfile=>$projectdir.$logfile,
	miRNAs_folder=>$DE_folder_miRNAs,
	verbose=>$verbose, 
	projectdir=>$projectdir,
	organism=>$organism,
	miARmaPath=>$miARmaPath,
	edger_cutoff=>$edger_cutoff,
	noiseq_cutoff=>$noiseq_cutoff,
);
