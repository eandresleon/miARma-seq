#!/usr/bin/perl

use strict;
use warnings;

##################################################################################
#1. ReadAlignment variables declaration
##################################################################################
my $reads_dir; #Path of the directory with the reaper results
my $projectdir; #Directory to save the results
my $threads; #Optional number of threads to perform the analysis
my $bowtie1index; #Indexed genome to align your reads in format .ebwt (Mandatory for analysis with bowtie1)
my $bowtie2index; #Indexed genome to align your reads in format .bt2 (Mandatory for analysis with bowtie2)
my $verbose; #Option to show the execution data on the screen   
my $logfile; #Run.log file where execution data will be saved
my $tophat_multihits; #Max # mismatches in seed alignment in bowtie analysis (0-1)
my $read_mismatches; #Length of seed substrings in bowtie analysis (>3, <32)
my $tophat_seg_mismatches; #Read segments are mapped independently, allowing up to this many mismatches in each segment alignment. The default is 2.
my $tophat_seg_length; #Each read is cut up into segments, each at least this long. These segments are mapped independently. The default is 25.
my $library_type; #The default is unstranded (fr-unstranded). If either fr-firststrand or fr-secondstrand is specified, every read alignment will have an XS attribute tag as explained below. Consider supplying library type options below to select the correct RNA-seq protocol.
my $tophatParameters; #Other tophat parameters to perform the analysis using the recommended syntaxis by tophat
my $aligner; #Aligner which will be use in the analysis (Allowed values: tophat)
my $tophat_aligner;# Tophat uses bowtie2 by default, bowtie1 can be also specified. Especifing bowth, means to repeat the analysis one with each aligner (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1))
my $miARmaPath;#Path to software
my $statsfile; #Stats.log file where stats data will be saved
my $GTF; #is a file used to hold information about gene structure
my $Seqtype;

BEGIN{
	$miARmaPath="../../../../";#Path to software. Full path is recommended
	$reads_dir="../../../basic_examples/mRNAs/reads/"; #Path of the directory with the reaper results
	$projectdir="."; #Directory to save the results
	$threads=14; #Optional number of threads to perform the analysis
	$GTF="../../../../data/Homo_sapiens_CHR_.GRCh37.74.gtf";
	$bowtie2index="../../../../Genomes/Indexes/bowtie2/human/bw2_homo_sapiens19"; #Indexed genome to align your reads in format .bt2 (Mandatory for analysis with bowtie2)
	$verbose=""; #Option to show the execution data on the screen   
	$logfile="/run_".$$.".log"; #Run.log file where execution data will be saved
	$tophat_multihits=5; #Max # mismatches in seed alignment in bowtie analysis (0-1)
	$read_mismatches=2; #Length of seed substrings in bowtie analysis (>3, <32)
	$tophat_seg_mismatches=1; #Read segments are mapped independently, allowing up to this many mismatches in each segment alignment. The default is 2.
	$tophat_seg_length=20; #Each read is cut up into segments, each at least this long. These segments are mapped independently. The default is 25.
	$library_type="fr-firststrand"; #The default is unstranded (fr-unstranded). If either fr-firststrand or fr-secondstrand is specified, every read alignment will have an XS attribute tag. Consider supplying library type options below to select the correct RNA-seq protocol.
	$Seqtype="Paired";
	$statsfile="/stats_".$$.".log"; #Stats.log file where stats data will be saved
	$tophatParameters="--splice-mismatches 0"; #Other tophat parameters to perform the analysis using the recommended syntaxis by tophat
	$aligner="tophat"; #Aligner which will be use in the analysis (Allowed values: tophat)
	$tophat_aligner="Bowtie2" # Tophat uses bowtie2 by default, bowtie1 can be also specified. Especifing bowth, means to repeat the analysis one with each aligner (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1))
}

use lib "$miARmaPath/lib/";
use lib "$miARmaPath/lib/Perl";
use CbBio::RNASeq::Aligner;

##################################################################################
#2. Reading the input directory and executing ReadAlignment function
##################################################################################

#Reading CutAdapt results directory, collecting the files and completing with the path
opendir(READ_DIR, $reads_dir) || die $!; 
my @read_files= readdir(READ_DIR);
@read_files=map("$reads_dir$_",@read_files);

#READALIGNMENT EXECUTION 
# Reading the array with the names of the files
foreach my $file(@read_files){
	#Selecting only the fastq files for their processing
	ReadAligment(
		file=>$file,
		aligner=>$aligner,
		tophat_aligner=>$tophat_aligner,
		threads=>$threads,
		bowtie2index=>$bowtie2index,
		verbose=>$verbose,
		logfile=>$projectdir.$logfile,
		statsfile=>$projectdir.$statsfile,
		tophat_multihits=> $tophat_multihits,
		read_mismatches=>$read_mismatches,
		tophat_seg_mismatches=>$tophat_seg_mismatches,
		tophat_seg_length=>$tophat_seg_length,
		library_type=>$library_type,
		projectdir=>$projectdir,
		tophatParameters=>$tophatParameters,
		miARmaPath=>$miARmaPath,
		GTF=>$GTF,
		Seqtype=>$Seqtype,
	);
}
