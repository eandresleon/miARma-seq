#########################################################################	
#	Check package		 									#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2015 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################

package CbBio::RNASeq::miARma;
#Export package system
require Exporter;
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(run_miARma);

use strict;

=head1 NAME

 check_parameters 

=head1 SYNOPSIS

Adapt package is composed of 9 subroutines: CutAdapt,CutAdaptStats, AdaptTriming, ReadFilter,
Minion, Reaper, check_parameters, AdapterGraph and ReadLengthCount. The aim of this package is to 
process the reads from a fastq file removing a known adapter sequence (CutAdapt, Reaper) or not 
(AdaptTriming), predicting the adapter sequence of the reads (Minion), filtering the reads by length 
(ReadFilter) or printing statistical info after the read processing (CutAdaptStats, AdapterGraph, 
ReadLengthCount). In addition check_parameters acts as main function calling to CutAdapt, CutAdaptStats,  
AdaptTriming, Minion, Reaper, AdapterGraph and ReadLengthCount according to the parameters selected.

=head1 Methods


=head2 check_parameters

  Example    : 
  check_parameters(
  	dir=>"../reads",
  	adaptersoft=>"Cutadapt-Reaper",
  	projectdir=>".",
  	file=>"file.fastq",
  	verbose=>"verbose",
  	logfile=>"run.log",
  	statsfile=>"stats.log",
  	organism=>"human",
  	min=>"18", 
  	max=>"26",
  	min_quality=>"1", 
  	cutadaptparameters=> "-O 4",
  	reaperparameters=>" -3p-prefix 12/2/0/0 -dust-suffix-late 20 --nozip"
  	);
  Description: check_parameters function performs an adapter removal process with cutadapt, reaper software or both, and/or
  the removal of a defined number of nucleotides of the reads (adapttrimming). These analysis generate new fastq files with 
  the conserved reads at Cutadapt_results, Reaper_results or/and AdaptTrimming_results directory in the provided project directory. 
  In addition, if user do not provide the adapter sequence check_parameters can predict it executing Minion software. 
  After processing the reads check_parameters generates statistical data of the process
    Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where the input files are 
  	 [projectdir] Directory where results directory will be created
  	 [file] File which is going to be processing (fasta/fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 [adaptersoft] Specific software to remove the adapter from the sequences (Allowed values: cutadapt, reaper, adapttrimming, 
  	 cutadapt-reaper, cutadapt-adapttrimming, reaper-adapttrimming or cutadapt-reaper-adapttrimming)
  	 
  	 Optional parameters:
  	 [adapter] Adapter sequence to be removed in the analysis (if this sequence is not provided Minion software will be used to predict one)
  	 [organism] To predict an adapter organism is needed to know which genome to blat against (Allowed organism: human, mouse and rat) (Mandatory is the adapter sequence has not been provided)
  	 [adaptpredictionnumber] Number of adapter predictions to show
  	 [minionadaptersequence] Adapter known sequence to compare with the sequence predicted by Minion
  	 [min] Minimun length of the sequence read to keep with Cutadapt and Reaper Software [15 as default to miRNA analysis]
	 [max] Maximun length of the sequence read to use with Cutadapt Software [35 as default to miRNA analysis]
	 [min_quality] Minimun quality of the sequence read to use with Cutadapt Software [0 as default to miRNA analysis]
	 [cutadaptparameters] Other cutadapt parameters to perform the analysis using the cutadapt recommended syntaxis
	 [metafile] Metadata file with the information of the adapter sequences to use with Reaper Software. See reaper instructions to generate this file (Mandatory file for 3p-bc and 5p-bc geometries). 
	 [reaperparameters] Parameters to execute reaper. See reaper instructions to introduce the parameters with the correct syntaxis
  	 [geom] Geometry used in the analysis with Reaper refering to the position of the barcode in the read (No barcode "no-bc" (default value), 3'end "3p-bc" or 5'end "5p-bc").
  	 [tabu] Tabu sequence to remove the read which contain this sequence with Reaper Software (usually this sequence is 5' primer sequence)
  	 [trimmingnumber] Number of nucleotides to remove of the sequence read and the quality data with AdaptTrimming function
  	 [readposition] End of the read to remove the nucleotides with AdaptTrimming function (3 or 5). 
  	 [verbose] Option to show the execution data on the screen

  Returntype : Fastq files and plots with reads distribution at directory Cutadapt_results, Reaper_results and/or AdaptTrimming_results. 
  check_parameters also return the paths of the new generated files
  Requeriments: cutadapt function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Cutadapt v1.2.1 or higher software correctly installed
  	- Reaper software correctly installed
  	- Minion software correctly installed
  	- R v3.1.0 or higher software correctly installed 
  	- Input files on fasta/fastq format (compressed files are accepted)
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub run_miARma{
	#Arguments provided by user are collected by %args. 
	print_header();
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $configuration_file=$args{"config"};
	use lib "$miARmaPath/lib/Perl";
	use Config::IniFiles;
	use DateTime;
	
	my $cfg = Config::IniFiles->new( -file => $configuration_file, -default => "General" );
		
	#Check for general parameters
	if($cfg->SectionExists("General")==1){
		
		#Mandatory parameters: read folder
		if($cfg->exists("General","read_dir") eq "" or ($cfg->val("General","read_dir") eq "")){
			print STDERR "\nERROR " . date() . " read_dir/RNA_read_dir parameter in Section [General] is missing/unfilled. Please check documentation\n";
			help_check_general();
		}
		#Mandatory parameters: label
		elsif($cfg->exists("General","label") eq ""){
			print STDERR "\nERROR " . date() . " label parameter in Section [General] is missing. Please check documentation\n";
			help_check_general();
		}
		#Mandatory parameters: path for binaries
		elsif($cfg->exists("General","miARmaPath") eq "" and $cfg->val("General","miARmaPath") ne ""){
			print STDERR "\nERROR " . date() . " miARmaPath parameter in Section [General] is missing. Please check documentation\n";
			help_check_general();
		}
		#Mandatory parameters: path for results
		elsif($cfg->exists("General","projectdir") eq ""){
			print STDERR "\nERROR " . date() . " projectdir parameter in Section [General] is missing. Please check documentation\n";
			help_check_general();
		}
		#Mandatory parameters: stats folder
		elsif($cfg->exists("General","stats_file") eq ""){
			print STDERR "\nERROR " . date() . " stats_file parameter in Section [General] is missing. Please check documentation\n";
			help_check_general();
		}
		else{
			#First, we are going to check if input files are real fastq files
			check_input_format(-config => $cfg);
			
		}
	}
	else{
		print STDERR "\nERROR " . date() . " Section [General] is missing. Please check documentation\n";
		help_check_general();
	}
	
	#Quality section
	if($cfg->SectionExists("Quality")==1){
		#Mandatory parameters: label
		if($cfg->exists("Quality","prefix") eq ""){
			print STDERR "\nERROR " . date() . " prefix parameter in Section [Quality] is missing. Please check documentation\n";
			help_check_quality();
		}
		else{
			
			#run quality;
			use CbBio::RNASeq::Quality;
			
			# Reading the directory collecting the files and completing with the path
			my $dir=$cfg->val("Quality","read_dir");
			if(opendir(DIR, $dir)){
				print STDERR "miARma :: ".date()." Starting a \"".$cfg->val("Quality","prefix")."\" Quality Analysis\n";
				my @files= readdir(DIR);
				@files=map("$dir/$_",@files);
				my $output_dir;
				
				# # FASTQC EXECUTION
				# # Reading the array with the names of the files
				foreach my $file(@files){
					#Calling FastQC subroutine of Quality.pm package.
					$output_dir=FastQC(
						miARmaPath=>$miARmaPath,
						file=>$file,
						projectdir=>$cfg->val("General","projectdir"),
						threads=>$cfg->val("General","threads"),
						verbose=>$cfg->val("General","verbose"),
						logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile"),
						prefix=>$cfg->val("Quality","prefix"),
					);
				}
				# # FASTQCSTATS EXECUTION
				# # Calling FastQCStats sobroutine of Quality.pm package. 
				FastQCStats(
					dir=>$output_dir, 
					verbose=>$cfg->val("General","verbose"),
					statsfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","stats_file"),
					logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile")
				);
			}
			 else{
			 	print "ERROR :: Please check that your reads are saved in: ($dir)\n";
			 	help_check_quality();
			}
		}
	}
	#Adapter removal
	if($cfg->SectionExists("Adapter")==1){
		#Mandatory parameters: read folder
		if($cfg->exists("Adapter","adaptersoft") eq "" or ($cfg->val("Adapter","adaptersoft") eq "")){
			print STDERR "\nERROR " . date() . " adaptersoft parameter in Section [Adapter] is missing/unfilled. Please check documentation\n";
			help_check_adapter();
		}
		else{
			#run Adapter
			use CbBio::RNASeq::Adapt;
			# Reading the directory collecting the files and completing with the path
					
			my $dir=$cfg->val("General","read_dir");
						
			if(opendir(DIR, $dir)){
				print STDERR "miARma :: ".date()." Starting a adapter removal analysis\n";
				my @files= readdir(DIR);
				AdapterRemoval(
					adaptersoft=>$cfg->val("Adapter","adaptersoft"),
					dir=>$dir,
					files=>\@files,
					adapter=>$cfg->val("Adapter","adapter")|| undef,,
					logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile"),
					statsfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","stats_file"),
					verbose=>$cfg->val("General","verbose")|| undef,,
					projectdir=>$cfg->val("General","projectdir"),
					min=>$cfg->val("Adapter","min")|| undef,,
					max=>$cfg->val("Adapter","max")|| undef,,
					min_quality=>$cfg->val("Adapter","min_quality")|| undef,
					miARmaPath=>$miARmaPath,
					reaperparameters=>$cfg->val("Adapter","reaperparameters") || undef,
					organism=>$cfg->val("Adapter","organism")|| undef,
					trimmingnumber=>$cfg->val("Adapter","trimmingnumber")|| undef,
					readposition=>$cfg->val("Adapter","readposition")|| undef,
					adaptpredictionnumber=>$cfg->val("Adapter","adaptpredictionnumber")|| undef,
					minionadaptersequence=>$cfg->val("Adapter","minionadaptersequence")|| undef,
					cutadaptparameters=>$cfg->val("Adapter","cutadaptparameters")|| undef,
					metafile=>$cfg->val("Adapter","metafile")|| undef,
					reaperparameters=>$cfg->val("Adapter","metafile")|| undef,
					geom=>$cfg->val("Adapter","geom")|| undef,
					tabu=>$cfg->val("Adapter","tabu")|| undef,			
				);
			}	
			else{
				print "ERROR :: Adapt:: Please check that your reads are saved in: ($dir)\n";
			 	help_check_general();
			}
		}
	}
	#Alignment
	if($cfg->SectionExists("Aligner")==1){
		#Mandatory parameters: read folder
		if($cfg->exists("Aligner","aligner") eq "" or ($cfg->val("Aligner","aligner") eq "")){
			print STDERR "\nERROR " . date() . " aligner parameter in Section [Aligner] is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("Aligner","aligner")) eq "bowtie1" and ($cfg->val("Aligner","bowtie1index") eq "")){
			print STDERR "\nERROR " . date() . " Bowtie1 has been selected as aligner but bowtie1index is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("Aligner","aligner")) eq "bowtie2" and ($cfg->val("Aligner","bowtie2index") eq "")){
			print STDERR "\nERROR " . date() . " Bowtie2 has been selected as aligner but bowtie2index is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		else{
			#run Adapter
			use CbBio::RNASeq::Aligner;
			
			my @files;
			if(lc($cfg->val("General","type"))=="mirna"){
			
				#Reading CutAdapt results directory, collecting the files and completing with the path
				my $cut_dir=$cfg->val("General","projectdir")."/cutadapt_results/";
				if($cut_dir){
					opendir(CUTDIR, $cut_dir) || warn "Aligner:: Folder $cut_dir is not found, but cutadapt has been specified as an adaptersoft\n"; 
					my @cut_files= readdir(CUTDIR);
					push(@files,map("$cut_dir$_",@cut_files));
				}
				#Reading Reaper results directory, collecting the files and completing with the path
				my $rea_dir=$cfg->val("General","projectdir")."/reaper_results/";
				if($rea_dir){
					opendir(READIR, $rea_dir) || warn "Aligner:: Folder $rea_dir is not found, but Reaper has been specified as an adaptersoft\n"; 
					my @rea_files= readdir(READIR);
					push(@files,map("$rea_dir$_",@rea_files));
				}
				#Reading Trimming results directory, collecting the files and completing with the path
				my $trim_dir=$cfg->val("General","projectdir")."/AdaptTriming_results/";
				if($trim_dir){
					opendir(READIR, $trim_dir) || warn "Aligner:: Folder $trim_dir is not found, but Adaptrimming has been specified as an adaptersoft\n"; 
					my @trim_files= readdir(READIR);
					push(@files,map("$trim_dir$_",@trim_files));
				}
			}
			else{
				print "ERROR :: You are requesting a miRNA alignment analysis, but no files are found without adapter (needed for alinging)\n";
			 	help_check_adapter();
			}
			if(scalar(@files)>0){
				print STDERR "miARma :: ".date()." Starting a \"".$cfg->val("Aligner","aligner")."\" Alignment Analysis\n";
			
				# Reading the array with the names of the files
				foreach my $file(@files){
					#print STDERR "Reading $file\n";
					ReadAligment(
					    file=>$file,
					    aligner=>$cfg->val("Aligner","aligner"),
						threads=>$cfg->val("General","threads") || 1,
					    bowtie2index=>$cfg->val("Aligner","bowtie2index") || undef,
					    bowtie1index=>$cfg->val("Aligner","bowtie1index") || undef,
						logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile") || undef,
						statsfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","stats_file") || undef,
						verbose=>$cfg->val("General","verbose") || undef,
					    bowtiemiss=>$cfg->val("Aligner","bowtiemiss") || undef,
					    bowtieleng=>$cfg->val("Aligner","bowtieleng") || undef,
						projectdir=>$cfg->val("General","projectdir")|| undef,
					    bowtie1parameters=>$cfg->val("General","bowtie1parameters")|| undef,
					    bowtie2parameters=>$cfg->val("General","bowtie2parameters")|| undef,
					    miARmaPath=>$miARmaPath,
						organism=>$cfg->val("General","organism")|| undef,
					);
				}
			}
		}
	}
	#Count features
	if($cfg->SectionExists("ReadCount")==1){
		#Mandatory parameters: read folder
		if($cfg->exists("ReadCount","database") eq "" or ($cfg->val("ReadCount","database") eq "")){
			print STDERR "\nERROR " . date() . " database parameter in Section [ReadCount] is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		else{
			#Running
			use CbBio::RNASeq::Readcount;
			
			my @files;
			if(lc($cfg->val("General","type"))=="mirna"){
				#Reading CutAdapt results directory, collecting the files and completing with the path
				my $bw1_dir=$cfg->val("General","projectdir")."/Bowtie1_results/";
				if($bw1_dir){
					opendir(BW1DIR, $bw1_dir) || warn "Aligner:: Folder $bw1_dir is not found\n"; 
					my @bw1_files= readdir(BW1DIR);
					push(@files,map("$bw1_dir$_",@bw1_files));
					close BW1DIR;
				}
				#Reading Reaper results directory, collecting the files and completing with the path
				my $bw2_dir=$cfg->val("General","projectdir")."/Bowtie2_results/";
				if($bw2_dir){
					opendir(BW2DIR, $bw2_dir) || warn "Aligner:: Folder $bw2_dir is not found\n"; 
					my @bw2_files= readdir(BW2DIR);
					push(@files,map("$bw2_dir$_",@bw2_files));
					close BW2DIR;
				}
			}
			else{
				print "ERROR :: You are requesting a miRNA readcount analysis, but no aligned files are found (Neither Bowtie1 nor Bowtie2)\n";
			 	help_check_aligner();
			}
			if(scalar(@files)>0){
				print STDERR "miARma :: ".date()." Starting a Readcount Analysis\n";
				my @htseqfiles;
				# Reading the array with the names of the files
				foreach my $file(@files){
					#Selecting only the sam files for their processing
					my $result=featureCount(
						file=>$file,
						database=>$cfg->val("ReadCount","database"),
						seqid=>$cfg->val("ReadCount","seqid") || undef,
						parameters=>$cfg->val("ReadCount","parameters") || undef, 
						strand=>$cfg->val("ReadCount","strand") || undef, 
						featuretype=>$cfg->val("ReadCount","featuretype") || undef,
						logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile") || undef,
						verbose=>$cfg->val("General","verbose") || 0,
						projectdir=>$cfg->val("General","projectdir")|| undef,
						miARmaPath=>$miARmaPath,
						threads=>$cfg->val("General","threads") || 1,
						quality=>$cfg->val("ReadCount","quality") || undef,
					);
					push(@htseqfiles, $result);
				}
				
				#HTSEQFORMATEXECUTION
				featureFormat( 
				  	input=>\@htseqfiles, 
				  	projectdir=>$cfg->val("General","projectdir")|| undef,
				  	logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile") || undef,
				  );
			}
			else{
				print "ERROR :: You are requesting a miRNA readcount analysis, but no aligned files are found (Neither Bowtie1 nor Bowtie2)\n";
				help_check_count();
			}
		}
	}
	
	if($cfg->SectionExists("DEAnalysis")==1){
		#Mandatory parameters: read folder
		if($cfg->exists("DEAnalysis","targetfile") eq "" or ($cfg->val("DEAnalysis","targetfile") eq "")){
			print STDERR "\nERROR " . date() . " targetfile parameter in Section [DEAnalysis] is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("DEAnalysis","contrastfile")) eq "" and ($cfg->val("DEAnalysis","contrastfile") eq "")){
			print STDERR "\nERROR " . date() . " contrastfile parameter in Section [DEAnalysis] is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("DEAnalysis","filter")) eq "" and ($cfg->val("DEAnalysis","filter") eq "")){
			print STDERR "\nERROR " . date() . " filter parameter in Section [DEAnalysis] is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("DEAnalysis","DEsoft")) eq "" and ($cfg->val("DEAnalysis","DEsoft") eq "")){
			print STDERR "\nERROR " . date() . " DEsoft parameter in Section [DEAnalysis] is missing/unfilled. Please check documentation\n";
			help_check_deanalysis();
		}
		else{
			#run DEAnalysis
			use CbBio::RNASeq::DEAnalysis;
			print STDERR "miARma :: ".date()." Differential expression analysis using  ".$cfg->val("DEAnalysis","DEsoft") ." software(s)\n";
			my $dir=$cfg->val("General","projectdir")."/Readcount_results/"; 

			if(opendir(DIR, $dir)){
				opendir(DIR, $dir) || die $!;
				my @files= readdir(DIR);
				@files=map("$dir$_",@files);

				foreach my $file(@files){
					DE_Analysis(
						projectdir=>$cfg->val("General","projectdir")|| undef,
						dir=>$dir,
						file=>$file,
						targetfile=>$cfg->val("DEAnalysis","targetfile"),
						label=>$cfg->val("General","label"),
						filter=>$cfg->val("DEAnalysis","filter")|| undef,
						edger_contrastfile=>$cfg->val("DEAnalysis","contrastfile")|| undef,
						noiseq_contrastfile=>$cfg->val("DEAnalysis","contrastfile")|| undef,
						DEsoft=>$cfg->val("DEAnalysis","DEsoft")|| undef,
						filtermethod=>$cfg->val("DEAnalysis","DEsoft")|| undef,
						logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile") || undef,
						verbose=>$cfg->val("General","verbose") || 0,
						cpmvalue=>$cfg->val("DEAnalysis","cpmvalue")|| undef,
						repthreshold=>$cfg->val("DEAnalysis","repthreshold")|| undef,
						edger_normethod=>$cfg->val("DEAnalysis","edger_normethod")|| undef,
						noiseq_normethod=>$cfg->val("DEAnalysis","noiseq_normethod")|| undef,
						replicates=>$cfg->val("DEAnalysis","replicates") || undef,
						miARmaPath=>$miARmaPath,
					);
				}
			}
			else{
				print "ERROR :: Please check that your tab files are saved in: ($dir)\n";
				help_check_deanalysis();
			}
		}
	}
	if($cfg->SectionExists("TargetPrediction")==1){
		#run DEAnalysis
		my $dir=$cfg->val("General","projectdir");		
		use CbBio::RNASeq::TargetPrediction;
		TargetPrediction(
			miRNAs_folder=>$cfg->val("General","projectdir"),
			logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile") || undef,
			verbose=>$cfg->val("General","verbose") || 0,
			projectdir=>$cfg->val("General","projectdir")|| undef,
			organism=>$cfg->val("General","organism")|| undef,
			miARmaPath=>$miARmaPath,
			edger_cutoff=>$cfg->val("TargetPrediction","edger_cutoff")|| undef,
			noiseq_cutoff=>$cfg->val("TargetPrediction","noiseq_cutoff")|| undef,
		);
		
	}
	if($cfg->SectionExists("DeNovo")==1){
		
		
		if($cfg->exists("DeNovo","bowtie1index") eq "" or ($cfg->val("DeNovo","bowtie1index") eq "")){
			print STDERR "\nERROR " . date() . " bowtie1index parameter in Section [DeNovo] is missing/unfilled. Please check documentation\n";
			help_check_denovo();
		}
		elsif( lc($cfg->val("DeNovo","mature_miRNA_file")) eq "" and ($cfg->val("DeNovo","mature_miRNA_file") eq "")){
			print STDERR "\nERROR " . date() . " mature_miRNA_file parameter in Section [DeNovo] is missing/unfilled. Please check documentation\n";
			help_check_denovo();
		}
		elsif( lc($cfg->val("DeNovo","precursor_miRNA_file")) eq "" and ($cfg->val("DeNovo","precursor_miRNA_file") eq "")){
			print STDERR "\nERROR " . date() . " precursor_miRNA_file parameter in Section [DeNovo] is missing/unfilled. Please check documentation\n";
			help_check_denovo();
		}
		elsif( lc($cfg->val("DeNovo","genome")) eq "" and ($cfg->val("DeNovo","genome") eq "")){
			print STDERR "\nERROR " . date() . " genome parameter in Section [DeNovo] is missing/unfilled. Please check documentation\n";
			help_check_denovo();
		}
		# Read count of new miRNAs
		elsif($cfg->exists("DeNovo","targetfile") eq "" or ($cfg->val("DeNovo","targetfile") eq "")){
			print STDERR "\nERROR " . date() . " targetfile parameter in Section [DeNovo] is missing/unfilled. Please check documentation\n";
			help_check_denovo();
		}
		elsif( lc($cfg->val("DeNovo","contrastfile")) eq "" and ($cfg->val("DeNovo","contrastfile") eq "")){
			print STDERR "\nERROR " . date() . " contrastfile parameter in Section [DeNovo] is missing/unfilled. Please check documentation\n";
			help_check_denovo();
		}
		elsif( lc($cfg->val("DeNovo","filter")) eq "" and ($cfg->val("DeNovo","filter") eq "")){
			print STDERR "\nERROR " . date() . " filter parameter in Section [DeNovo] is missing/unfilled. Please check documentation\n";
			help_check_denovo();
		}
		elsif( lc($cfg->val("DeNovo","DEsoft")) eq "" and ($cfg->val("DeNovo","DEsoft") eq "")){
			print STDERR "\nERROR " . date() . " DEsoft parameter in Section [DeNovo] is missing/unfilled. Please check documentation\n";
			help_check_denovo();
		}
		
		else{
			#run Aligner
			use CbBio::RNASeq::Aligner;
			
			my $read_directory=$cfg->val("General","read_dir") . "/";
			my @fastaq_files;
			if($read_directory){
				opendir(READIR, $read_directory) || warn "DeNovo:: Folder $read_directory is not found\n"; 
				my @files= readdir(READIR);
				push(@fastaq_files,map("$read_directory$_",@files));
			}
			foreach my $file(@fastaq_files){
				ReadAligment(
					file=>$file,
					aligner=>"miRDeep",
					threads=>$cfg->val("General","threads") || 1,
				    bowtie1index=>$cfg->val("DeNovo","bowtie1index") || undef,
					statsfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","stats_file"),
					verbose=>$cfg->val("General","verbose")|| undef,,
					projectdir=>$cfg->val("General","projectdir"),
					logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile"),
					miARmaPath=>$miARmaPath,
					adapter=>$cfg->val("DeNovo","adapter") || undef,
					precursors=>$cfg->val("DeNovo","precursor_miRNA_file"),
					mature=>$cfg->val("DeNovo","mature_miRNA_file"),
					genome=>$cfg->val("DeNovo","genome"),
					organism=>$cfg->val("General","organism"),
				);
			}
			#Once ReadAlignment is finished, is time to count the reads
			my $processed_reads=$cfg->val("General","projectdir") . "/miRDeep_results/";
			my @mirdeep_files;
			if($processed_reads){
				opendir(MIRDE, $processed_reads) || warn "Can't find miRDeep results ($processed_reads)\n"; 
				my @files= readdir(MIRDE);
				push(@mirdeep_files,map("$processed_reads$_",@files));
			}
			else{
				print STDERR "ERROR :: ".date(). " Can't find miRDeep results in \"$processed_reads\"\n"; 
			}
			my @allRNAfiles;
			my $result_file;
			if(scalar(@mirdeep_files)>0){
				foreach my $file(@mirdeep_files){
					#Selecting only the sam files for their processing
					my $result=miRDeepCount(
					  	file=>$file,
						logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile"),
						verbose=>$cfg->val("General","verbose")|| undef,
						projectdir=>$cfg->val("General","projectdir"),
						miARmaPath=>$miARmaPath,
					);
					push(@allRNAfiles, $result) if($result);
				}
			}
			if(scalar(@allRNAfiles)>0){
				$result_file=miRDeepFormat( 
				  	input=>\@allRNAfiles, 
					projectdir=>$cfg->val("General","projectdir"),
					logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile"),
				);
			}
			use CbBio::RNASeq::DEAnalysis;
			use File::Basename;
			
			my($file_name,$path)=fileparse($result_file);

			if($result_file){
				DE_Analysis(
					projectdir=>$cfg->val("General","projectdir")|| undef,
					dir=>$path,
					file=>$result_file,
					targetfile=>$cfg->val("DeNovo","targetfile"),
					label=>"DeNovo_" . $cfg->val("General","label"),
					filter=>$cfg->val("DeNovo","filter")|| undef,
					edger_contrastfile=>$cfg->val("DeNovo","contrastfile")|| undef,
					noiseq_contrastfile=>$cfg->val("DeNovo","contrastfile")|| undef,
					DEsoft=>$cfg->val("DeNovo","DEsoft")|| undef,
					filtermethod=>$cfg->val("DeNovo","DEsoft")|| undef,
					logfile=>$cfg->val("General","projectdir")."/".$cfg->val("General","logfile") || undef,
					verbose=>$cfg->val("General","verbose") || 0,
					cpmvalue=>$cfg->val("DeNovo","cpmvalue")|| undef,
					repthreshold=>$cfg->val("DeNovo","repthreshold")|| undef,
					edger_normethod=>$cfg->val("DeNovo","edger_normethod")|| undef,
					noiseq_normethod=>$cfg->val("DeNovo","noiseq_normethod")|| undef,
					replicates=>$cfg->val("DeNovo","replicates") || undef,
					miARmaPath=>$miARmaPath,
				);
			}
			else{
				print STDERR "Can't find $result_file\n";
			}
		}
	}

sub check_input_format{
	use File::Basename;
	my %args=@_;
	my $cfg=$args{"-config"};#configuration file
	
	my $dir_fastq=$cfg->val("General","read_dir");
	my @fastq_files;
	
	if($dir_fastq){
		opendir(FASTQF, $dir_fastq) || exit "Can't find fastq files on $dir_fastq\n"; 
		my @files= readdir(FASTQF);
		push(@fastq_files,map("$dir_fastq$_",@files));
	}
	my $error=0;
	if(scalar(@fastq_files)>0){
		print STDERR "miARma :: " . date() . " Checking if files in \"$dir_fastq\" are in the correct fastq format\n";
		foreach my $file(@fastq_files){
			
			my ($real_file,$path)=fileparse($file);
			if($real_file !~ /^\./){
				my $compressed_file=0;
				if($real_file =~ /.*\.gz/){
					my $file_tmp="gunzip -d -f -c $file |";
					$file=$file_tmp;
					$compressed_file=1;
				}
				elsif($real_file =~ /.*\.bz2/){
					my $file_tmp="bunzip2 -f -d -c $file |";
					$file=$file_tmp;
					$compressed_file=1;
				}

				open IN,$file or die "Cannot open FASTQ file ($file)\n";
				my $i=0;
				my $mes="Please make sure your file ($real_file) is in accordance with the FASTQ format specifications";
				while(<IN>){
			    	chomp;
			    	$i++;
					if($i == 1){if(/^@\S+/){}else{print STDERR "miARmA ERROR :: ".date() . " First line of FASTQ reads file ($real_file) is not in accordance with the fastq format specifications\n$mes\n";exit;}}                 
					if($i == 2){if(/^\S+$/){}else{print STDERR "miARmA ERROR :: ".date() . " Second line of FASTQ reads file ($real_file) contains whitespace in sequence\n$mes\n";exit;}}     
					if($i == 3){if(/^\+/){}else{print STDERR "miARmA ERROR :: ".date() . " Third line of FASTQ reads file ($real_file) does not start with a '+' character.\n$mes\n";exit;}}     
					if($i == 4){if(/^\S+$/){}else{print STDERR "miARmA ERROR :: ".date() . " Fourth line of FASTQ reads file ($real_file) contains whitespace\n$mes\n";exit;}}     
					
					last if($i == 4);
				}
				close IN;
				#The file seems correct, but cutadapt needs to have a fastq/fq extension
				if($real_file !~ /\.fq(\.bz2)?$/ and $real_file !~ /\.fastq(\.bz2)?$/ and $real_file !~ /\.fq(\.gz)$/ and $real_file !~ /\.fastq(\.gz)$/ and lc($cfg->val("Adapter","adaptersoft")) =~ /cutadapt/){
					print STDERR "ERROR :: Due to cutadapt restrictions, file extension must be .fq or .fastq (or corresponding compressed files .fq.gz/.fq.bz2 or .fastq.gz/fastq.bz2)\nPlease, consider to rename the file $real_file\n";
					$error=1;
					next;
				}
			}
		}
	}
	else{
		print STDERR "Can't find fastq files on $dir_fastq\n";
		exit;
	}
	if($error == 0){
		return();
	}
	else{
		print STDERR "Quitting miARma\n";
		exit;
	}
}	
	
	sub help_check_general{
	    my $usage = qq{
Mandatory parameters:

[General]
read_dir=miRNA_raw_reads_folder/
label=Hypoxia
miARmaPath=.
projectdir=results/
logfile=miARma_analysis_logfile.log
stats_file=stats.log

};

	print STDERR $usage;
	exit();
}
	
	sub help_check_quality{
	    my $usage = qq{
Mandatory parameters:

[Quality]
prefix=Pre

};

	print STDERR $usage;
	return();
}
	
	sub help_check_aligner{
	    my $usage = qq{
Mandatory parameters:

[Aligner]
aligner=Bowtie1
bowtie1index=/genomes/bowtie1/hg19

or

[Aligner]
aligner=Bowtie2
bowtie2index=/genomes/bowtie2/hg19

};

	print STDERR $usage;
	return();
}
	
	sub help_check_adapter{
	    my $usage = qq{
Mandatory parameters:

[Adapter]
adaptersoft=CutAdapt

};

	print STDERR $usage;
	return();
	}
}

	sub help_check_count{
	    my $usage = qq{
Mandatory parameters:

[ReadCount]
database=miRNAs_miRBase20.gft

};

	print STDERR $usage;
	exit();
}

	sub help_check_deanalysis{
		    my $usage = qq{
	Mandatory parameters:

	[DEAnalysis]
	targetfile=targets.txt
	contrastfile=contrast.txt
	filter=yes
	DEsoft=EdgeR-Noiseq

	};

		print STDERR $usage;
		exit();
	}
	
	sub help_check_denovo{
	    my $usage = qq{
Mandatory parameters:

	[DeNovo]
	; Indexed genome to align your reads in format .ebwt (Mandatory for analysis with miRDeep)
	bowtie1index=Genomes/Indexes/bowtie1/human/hg19
	; Adapter to trimm at read 3'
	adapter=ATCTCGTATGCCGTCTTCTGCTTGAA
	; a fasta file with all mature sequence from your organism
	mature_miRNA_file=Examples/miRNAs/DeNovoPrediction/hsa_mature_miRBase20.fasta
	;a fasta file with all known pre-miRNa sequence 
	precursor_miRNA_file=Examples/miRNAs/DeNovoPrediction/precursors_miRBase20.fasta
	;fasta file for the cmplete genome of our organism
	genome=Genomes/Fasta/hg19.fa
	; Complete path of the target file.
	targetfile=Examples/miRNAs/data/targets.txt
	; Path of the contrast file.
	contrastfile=Examples/miRNAs/data/contrast.txt
	#This value refers to filter processing in the reads (Should be "yes" or "no").
	filter=yes
	;Specific software to perform the Differential Expression Analysis (Allowed values: edger, noiseq or edger-noiseq)
	DEsoft=EdgeR-Noiseq
	; providing replicates
	replicates=yes

	};

		print STDERR $usage;
		exit();
	}
	
sub print_header{
	system("clear");
	print "#########################################################################	
#   miARma, miRNA and RNASeq Multiprocess Analysis v1.0 (2015)          #
#                                                                       #
#   Created at Computational Biology and Bioinformatics Group (CbBio)   #
#   Institute of Biomedicine of Seville. IBIS (Spain)                   #
#   Copyright (c) 2015 IBIS. All rights reserved.                       #
#   mail : miARma-devel\@cbbio.es                                        #
#########################################################################\n\n";
}
sub date{
	my $dt = DateTime->now;
	return($dt->hms . " [" . $dt->dmy ."]");
}
1;