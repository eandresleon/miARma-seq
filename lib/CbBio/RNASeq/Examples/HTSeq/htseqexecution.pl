#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Readcount;

#Declaring htseqCount variables
my $label="Prueba2";
my $htdatabase="/Users/apple/Pipeline/examples/HTSeq/miRBase_Annotation_20_for_hsa_mature_miRNA.gtf";
my $htseqid="transcript_id";
my $htseqparameters=" -a 10"; 
my $htseqstrand="no"; 
my $htseqtype="miRNA",
my $logfile="/run_".$$.".log";
my $htseqformat="bam";
my $verbose="verbose"; 
my $projectdir="/Users/apple/Pipeline/examples/HTSeq";

#Declaring variables needed for the analysis
my $inputdir="/Users/apple/Pipeline/examples/Bowtie2/Bowtie2_results";
my $result;
my @htseqfiles;

#Reading the directory 
opendir(DIR, $inputdir) || die "$! : $inputdir"; 
my @files= readdir(DIR); 

# HTSEQCOUNT EXECUTION
#Reading the array with the names of the files
foreach my $file(@files){
	#Selecting only the sam files for their processing
	if($file =~ /.*\.bam$/ or $file =~ /.*\.sam$/){
		$result=htseqCount(
		  	label=>$label, 
		  	file=>$inputdir."/".$file,
		  	htdatabase=>$htdatabase,
		  	htseqid=>$htseqid,
		  	htseqparameters=>$htseqparameters, 
		  	htseqstrand=>$htseqstrand, 
		  	htseqtype=>$htseqtype,
		  	logfile=>$projectdir.$logfile,
		  	htseqformat=>$htseqformat,
		  	verbose=>$verbose, 
		  	projectdir=>$projectdir
		);
		push(@htseqfiles, $result);
	}
}

#HTSEQFORMATEXECUTION
htseqFormat( 
  		input=>\@htseqfiles, 
  		dir=>$projectdir,
  		label=>$label,
  		logfile=>$projectdir.$logfile
  	);