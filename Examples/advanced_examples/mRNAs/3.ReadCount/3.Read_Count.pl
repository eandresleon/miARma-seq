#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. htseqCount variables declaration
##################################################################################
my $database; #GFF file used to calculate the number of reads in featureCounts analysis
my $seqid; #GFF attribute to be used as feature ID (default: gene_id) for featureCounts analysis
my $parameters; #Other featureCounts parameters to perform the analysis using the featureCounts recommended syntaxis
my $strand; #Whether the data is from a strand-specific assay (yes, no or reverse, yes by default) for featureCounts analysis
my $featuretype; #Feature type (3rd column in GFF file) to be used, all features of other type are ignored (default:exon) for featureCounts analysis
my $logfile; #Run.log file where execution data will be saved
my $verbose; #Option to show the execution data on the screen
my $projectdir; #Directory where htseq_results directory will be created
my $bw1_dir; #Path of the directory with the bowtie1 results
my $bw2_dir; #Path of the directory with the bowtie2 results
my $threads; #Optional number of threads to perform the analysis
my $miARmaPath;#Path to software
my $Seqtype;#Type of sequencing ; could be Paired or Single. [Single by default]

BEGIN{
	$miARmaPath="../../../../";#Path to software. Full path is recommended
	$database="../../../../data/Homo_sapiens_CHR_.GRCh37.74.gtf"; #GFF file used to calculate the number of reads in featureCounts analysis
	$seqid="transcript_id"; #GFF attribute to be used as feature ID (default: gene_id) for featureCounts analysis
	$parameters=" -Q 10"; #Other featureCounts parameters to perform the analysis using the featureCounts recommended syntaxis
	$strand="yes"; #Whether the data is from a strand-specific assay (yes, no or reverse, yes by default) for featureCounts analysis
	$featuretype="exon"; #Feature type (3rd column in GFF file) to be used, all features of other type are ignored (default:exon) for featureCounts analysis
	$logfile="/run_".$$.".log"; #Run.log file where execution data will be saved
	$verbose=""; #Option to show the execution data on the screen
	$projectdir="."; #Directory where htseq_results directory will be created
	$bw1_dir="../2.Aligner/Bowtie1_results/"; #Path of the directory with the bowtie1 results
	$bw2_dir="../2.Aligner/Bowtie2_results/"; #Path of the directory with the bowtie2 results
	$threads=8; #Optional number of threads to perform the analysis
	$Seqtype="Paired";#Type of sequencing ; could be Paired or Single. [Single by default]

}

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::Readcount;

##################################################################################
#2. Reading the input directory and executing htseqCount function
##################################################################################
#Reading CutAdapt results directory, collecting the files and completing with the path
opendir(CUTDIR, $bw1_dir) || die $!; 
my @bw1_files= readdir(CUTDIR);
@bw1_files=map("$bw1_dir$_",@bw1_files);
#Reading Reaper results directory, collecting the files and completing with the path
opendir(READIR, $bw2_dir) || die $!; 
my @bw2_files= readdir(READIR);
@bw2_files=map("$bw2_dir$_",@bw2_files);

my @files=(@bw1_files, @bw2_files);

my @htseqfiles;
my $result; #result file

# HTSEQCOUNT EXECUTION
#Reading the array with the path of the files
foreach my $file(@files){
	#Selecting only the sam files for their processing
	$result=featureCount(
	  	file=>$file,
	  	database=>$database,
	  	seqid=>$seqid,
	  	parameters=>$parameters, 
	  	strand=>$strand, 
	  	featuretype=>$featuretype,
	  	logfile=>$projectdir.$logfile,
		verbose=>$verbose, 
	  	projectdir=>$projectdir,
		threads=>$threads,
		miARmaPath=>$miARmaPath,
		Seqtype=>$Seqtype
	);
	push(@htseqfiles, $result);
}
#HTSEQFORMATEXECUTION
featureFormat( 
  	input=>\@htseqfiles, 
  	projectdir=>$projectdir,
  	logfile=>$projectdir.$logfile
);
