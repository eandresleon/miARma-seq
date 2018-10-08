#########################################################################	
#	Check package		 												#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2018 IBIS. All rights reserved.						#
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
	my $check_input=$args{"check"};
	
	use lib "$miARmaPath/lib/Perl";
	use Config::IniFiles;
	use DateTime;
	
	my $cfg = Config::IniFiles->new( -file => $configuration_file, -default => "General" );
	my $log_file=undef;
	my $stat_file=undef;
	my $severe_error=0;
	my $post_qual=0;
	
	my $start_run = time();
	
	#Check for general parameters
	if($cfg->SectionExists("General")==1){
		
		#Mandatory parameters: read folder
		if($cfg->exists("General","read_dir") eq "" or ($cfg->val("General","read_dir") eq "")){
			print STDERR "\nERROR " . date() . " read_dir/RNA_read_dir parameter in Section [General] is missing/unfilled. Please check documentation\n";
			$severe_error=1;			
			help_check_general();
		}
		#Mandatory parameters: label
		elsif($cfg->exists("General","label") eq "" or $cfg->val("General","label") eq ""){
			print STDERR "\nERROR " . date() . " label parameter in Section [General] is missing. Please check documentation\n";
			$severe_error=1;			
			help_check_general();
		}
		#Mandatory parameters: path for binaries
		elsif($cfg->exists("General","miARmaPath") eq "" or $cfg->val("General","miARmaPath") eq ""){
			print STDERR "\nERROR " . date() . " miARmaPath parameter in Section [General] is missing. Please check documentation\n";
			$severe_error=1;			
			help_check_general();
		}
		#Mandatory parameters: path for results
		elsif($cfg->exists("General","output_dir") eq "" or $cfg->val("General","output_dir") eq "" ){
			print STDERR "\nERROR " . date() . " output_dir parameter in Section [General] is missing. Please check documentation\n";
			$severe_error=1;			
			help_check_general();
		}
		#Mandatory parameters: organism
		elsif($cfg->exists("General","organism") eq "" or $cfg->val("General","organism") eq "" ){
			print STDERR "\nERROR " . date() . " organism parameter in Section [General] is missing. Please check documentation\n";
			$severe_error=1;			
			help_check_general();
		}
		#Mandatory parameters: analysis type
		elsif($cfg->exists("General","type") eq "" or $cfg->val("General","type") eq ""){
			print STDERR "\nERROR " . date() . " type parameter in Section [General] is missing . Please check documentation\n";
			$severe_error=1;			
			help_check_general();
		}
		elsif($cfg->val("General","type") ne "circRNA" and $cfg->val("General","type") ne "miRNA" and $cfg->val("General","type") ne "mRNA"){
			print STDERR "\nERROR " . date() . " type parameter in Section [General] is misspelled. Only \"miRNA/mRNA or circRNA\" is accepted and you provide \"".$cfg->val("General","type")."\". Please check documentation\n";
			$severe_error=1;			
			help_check_general();
		}
		
		else{
			#First, we are going to check if input files are real fastq files if no error is found
			#check_input_format(-config => $cfg);
			
			#Time to check that all included files exist
			
			check_input_data(config=>$cfg,check=>$check_input,type=>$cfg->val("General","type"));
			
			system ("mkdir -p " . $cfg->val("General","output_dir"));
			#optional parameters: stats folder
			if($cfg->exists("General","stats_file") ne "" or defined($stat_file)){
				$stat_file=$cfg->val("General","stats_file");
				system("touch $stat_file");
			}
			if($cfg->exists("General","stats_file") eq "" or undef($stat_file)){
				$stat_file=$cfg->val("General","output_dir")."/miARma_stat.$$.log";
				$cfg->newval("General", "stats_file", $stat_file);
				$cfg->RewriteConfig;
				system("touch $stat_file");
			}
			#optional parameters: log_file
			if($cfg->exists("General","logfile") ne "" or defined($log_file)){
				$log_file=$cfg->val("General","logfile");
				system("touch $log_file");
			}
			if($cfg->exists("General","logfile") eq "" or undef($log_file)){
				$log_file=$cfg->val("General","output_dir")."/miARma_logfile.$$.log";
				$cfg->newval("General", "logfile", $log_file);
				$cfg->RewriteConfig;
				system("touch $log_file");
			}
			if($cfg->val("General","read_dir") !~ /\/$/){
				my $real_dir=$cfg->val("General","read_dir"). "/";
				$cfg->newval("General", "read_dir", $real_dir);
				$cfg->RewriteConfig;
			}
			if($cfg->val("General","output_dir") !~ /\/$/){
				my $real_dir=$cfg->val("General","output_dir"). "/";
				$cfg->newval("General", "output_dir", $real_dir);
				$cfg->RewriteConfig;
			}
		}
	}
	else{
		print STDERR "\nERROR " . date() . " Section [General] is missing. Please check documentation\n";
		help_check_general();
	}
	
	#Summary of results
	my $summary_file= $cfg->val("General","output_dir"). "/summary_results_" . $cfg->val("General","label") ."_miARma.xls";
	if(!-e $summary_file){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
		print SUMM date(). " Summary for " . $cfg->val("General","label") . ". ". $cfg->val("General","type") ." analysis using miARma.\n";
		close SUMM;
	}
	
	#Quality section
	if($cfg->SectionExists("Quality")==1){
		#Mandatory parameters: label
		if($cfg->exists("Quality","prefix") eq ""){ 
			#print STDERR "\nERROR " . date() . " prefix parameter in Section [Quality] is missing. Please check documentation\n";
			#help_check_quality();
			$cfg->newval("Quality", "prefix", "pre");
			$cfg->RewriteConfig;
		}
		
		#else{
			
			my $dir=undef;
			if(lc($cfg->val("Quality","prefix")) eq "pre"){
				$dir=$cfg->val("General","read_dir");
			}
			elsif(lc($cfg->val("Quality","prefix")) eq "post"){
				$post_qual=1;
			}
			elsif(lc($cfg->val("Quality","prefix")) eq "both"){
				$dir=$cfg->val("General","read_dir");
				$post_qual=1;
			}
			else{
				print STDERR "\nERROR " . date() . " prefix parameter in Section [Quality] doesn't accept ".$cfg->val("Quality","prefix")." as a value. Only pre, post or both are accepted as correct values\n";
				help_check_quality();
		 	}
			
			if($cfg->SectionExists("DeNovo")==1 and $post_qual==1){
				print date() . "\tWARN:: For DeNovo studies, Post quality analyses is not available\n";
				$post_qual=0;					
			}
			
			my $processed_files;
			#run quality;
			use CbBio::RNASeq::Quality;
			if($dir){
				# Reading the directory collecting the files and completing with the path
				if(opendir(DIR, $dir)){
					my @files= readdir(DIR);
					@files=map("$dir/$_",@files);
					my $output_dir;
				
					open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
					print SUMM "\nQuality [".$dir."/Pre_fastqc_results]\n";
					close SUMM;
					
					# # FASTQC EXECUTION
					# # Reading the array with the names of the files
					print date()." Starting Quality Analysis.\n";
					#$processed_files->{$file}++;
					#Calling FastQC subroutine of Quality.pm package.
					$output_dir=FastQC(
						miARmaPath=>$miARmaPath,
						files=>\@files,
						projectdir=>$cfg->val("General","output_dir"),
						threads=>$cfg->val("General","threads") || 1,
						verbose=>$cfg->val("General","verbose") || 0,
						logfile=>$log_file || $cfg->val("General","logfile"),
						prefix=>"Pre",
					);
					
					# # FASTQCSTATS EXECUTION
					# # Calling FastQCStats sobroutine of Quality.pm package. 
					FastQCStats(
						dir=>$output_dir, 
						verbose=>$cfg->val("General","verbose") || 0,
						statsfile=>$stat_file || $cfg->val("General","stats_file"),
						logfile=>$log_file || $cfg->val("General","logfile"),
						summary=>$summary_file
					);
					print date()." Quality Analysis finished.\n";
				}
				 else{
				 	print "ERROR :: Please check that your reads are saved in: ($dir)\n";
				 	help_check_quality();
				}
			#}
		}
	}
	#Adapter removal
	if($cfg->SectionExists("Adapter")==1){
		#Mandatory parameters: read folder
			if($cfg->exists("Adapter","adaptersoft") eq "" or ($cfg->val("Adapter","adaptersoft") eq "")){
				if(lc($cfg->val("General","type")) eq "mirna"){
					$cfg->newval("Adapter", "adaptersoft", "CutAdapt");
					$cfg->newval("Adapter", "adaptpredictionnumber", "12");
					$cfg->RewriteConfig;
				}
				else{
					print STDERR "\nERROR " . date() . " adaptersoft parameter in Section [Adapter] is missing/unfilled. Please check documentation\n";
					help_check_adapter();
				}
			}
			#run Adapter
			use CbBio::RNASeq::Adapt;
			# Reading the directory collecting the files and completing with the path
					
			my $dir=$cfg->val("General","read_dir");
			my @files;
			if(-e $dir){
				opendir(DIR, $dir) || warn "Adapter:: Folder $dir is not found\n"; 
				my @cut_files= readdir(DIR);
				foreach my $f(sort @cut_files){
					if (!-d "$dir/$f"){
					 	push(@files,"$f");
					 }
				}
			}
			else{
				print "ERROR :: Adapt:: Please check that your reads are saved in: ($dir)\n";
			 	help_check_general();
			}
			check_input_format(-files=>\@files,-dir=>$dir);

			print date()." Starting a Adapter removal analysis\n";
			
			@files=sort(@files);
			my @files_adapter=AdapterRemoval(
				adaptersoft=>$cfg->val("Adapter","adaptersoft"),
				dir=>$dir,
				files=>\@files,
				adapter=>$cfg->val("Adapter","adapter")|| undef,
				adapter_file=>$cfg->val("Adapter","adapter_file")|| undef,
				logfile=>$log_file || $cfg->val("General","logfile"),
				statsfile=>$stat_file || $cfg->val("General","stats_file"),
				verbose=>$cfg->val("General","verbose")|| 0,
				projectdir=>$cfg->val("General","output_dir"),
				min=>$cfg->val("Adapter","min")|| 15,
				max=>$cfg->val("Adapter","max")|| 35,
				min_quality=>$cfg->val("Adapter","min_quality")|| undef,
				miARmaPath=>$miARmaPath,
				reaperparameters=>$cfg->val("Adapter","reaperparameters") || undef,
				organism=>$cfg->val("General","organism")|| undef,
				trimmingnumber=>$cfg->val("Adapter","trimmingnumber")|| undef,
				readposition=>$cfg->val("Adapter","readposition")|| undef,
				adaptpredictionnumber=>$cfg->val("Adapter","adaptpredictionnumber")|| 12,
				minionadaptersequence=>$cfg->val("Adapter","minionadaptersequence")|| undef,
				cutadaptparameters=>$cfg->val("Adapter","cutadaptparameters")|| undef,
				metafile=>$cfg->val("Adapter","metafile")|| undef,
				reaperparameters=>$cfg->val("Adapter","metafile")|| undef,
				geom=>$cfg->val("Adapter","geom")|| undef,
				tabu=>$cfg->val("Adapter","tabu")|| undef,
				summary=>$summary_file,
							
			);
			print date()." Adapter Analysis finished.\n";
				
			#Just in case the user wants to see the quality of the processed reads
			if($post_qual==1){
				print date()." Starting a Post Quality Analysis\n";
				my $output_dir_post;
				open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
				print SUMM "\nQuality [".$dir."/Post_fastqc_results]\n";
				close SUMM;
				#foreach my $processed_files(@files_adapter){
					$output_dir_post=FastQC(
						miARmaPath=>$miARmaPath,
						files=>\@files_adapter,
						projectdir=>$cfg->val("General","output_dir"),
						threads=>$cfg->val("General","threads")|| 1,
						verbose=>$cfg->val("General","verbose")|| 0,
						logfile=>$log_file || $cfg->val("General","logfile"),
						prefix=>"Post",
					);
				#}
				# # FASTQCSTATS EXECUTION
				# # Calling FastQCStats sobroutine of Quality.pm package. 
				FastQCStats(
					dir=>$output_dir_post, 
					verbose=>$cfg->val("General","verbose") || 0,
					statsfile=>$stat_file || $cfg->val("General","stats_file"),
					logfile=>$log_file || $cfg->val("General","logfile"),
					summary=>$summary_file
				);
				print date()." Post Quality Analysis finished.\n";
			#}
		}
	}
	#Alignment
	if($cfg->SectionExists("Aligner")==1){
		#Mandatory parameters: read folder
		if($cfg->exists("Aligner","aligner") eq "" or ($cfg->val("Aligner","aligner") eq "")){
			print STDERR "\nERROR " . date() . " aligner parameter in Section [Aligner] is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("Aligner","aligner")) eq "bowtie1" and ( $cfg->val("Aligner","bowtie1index") eq "" and $cfg->val("Aligner","fasta") eq "")){
			print STDERR "\nERROR " . date() . " Bowtie1 has been selected as aligner but bowtie1index/fasta is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("Aligner","aligner")) eq "bowtie2" and ( $cfg->val("Aligner","bowtie2index") eq "" and $cfg->val("Aligner","fasta") eq "")){
			print STDERR "\nERROR " . date() . " Bowtie2 has been selected as aligner but bowtie2index/fasta is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("Aligner","aligner")) eq "bwa" and ( $cfg->val("Aligner","bwaindex") eq "" and $cfg->val("Aligner","fasta") eq "")){
			print STDERR "\nERROR " . date() . " BWA has been selected as aligner but bwaindex/fasta is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("Aligner","aligner")) eq "star" and ( $cfg->val("Aligner","starindex") eq "" and $cfg->val("Aligner","fasta") eq "")){
			print STDERR "\nERROR " . date() . " STAR has been selected as aligner but starindex/fasta is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("Aligner","aligner")) eq "star" and $cfg->exists("Aligner","startindex") ne "" and ( $cfg->val("Aligner","gft") eq "" )){
			print STDERR "\nERROR " . date() . " STAR has been selected to index a fasta file, but gtf file is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		elsif( lc($cfg->val("Aligner","aligner")) eq "hisat2" and ( $cfg->val("Aligner","hisat2index") eq "" and $cfg->val("Aligner","fasta") eq "")){
			print STDERR "\nERROR " . date() . " Hisat2 has been selected as aligner but hisat2index/fasta is missing/unfilled. Please check documentation\n";
			help_check_aligner();
		}
		else{
			#run Adapter
			use CbBio::RNASeq::Aligner;

			my $do_index=$cfg->val("Aligner","fasta");
			if($do_index and $cfg->val("Aligner","bowtie2index") ne "" and lc($cfg->val("Aligner","aligner")) eq "bowtie2"){
				print STDERR date() . " Inside [Aligner] fasta and bowtie2index are excluyent, discarding fasta parameter\n";
				$do_index=undef;
			}
			if($do_index and $cfg->val("Aligner","bowtie1index") ne "" and lc($cfg->val("Aligner","aligner")) eq "bowtie1"){
				print STDERR date() . " Inside [Aligner] fasta and bowtie1index are excluyent, discarding fasta parameter\n";
				$do_index=undef;
			}
			if($do_index and $cfg->val("Aligner","bwaindex") ne "" and lc($cfg->val("Aligner","aligner")) eq "bwa"){
				print STDERR date() . " Inside [Aligner] fasta and bwaindex are excluyent, discarding fasta parameter\n";
				$do_index=undef;
			}
			if($do_index and $cfg->val("Aligner","hisat2index") ne "" and lc($cfg->val("Aligner","aligner")) eq "hisat2"){
				print STDERR date() . " Inside [Aligner] fasta and hisat2index are excluyent, discarding fasta parameter\n";
				$do_index=undef;
			}
			if($do_index and $cfg->val("Aligner","starindex") ne "" and lc($cfg->val("Aligner","aligner")) eq "star"){
				print STDERR date() . " Inside [Aligner] fasta and bwaiindex are excluyent, discarding fasta parameter\n";
				$do_index=undef;
			}
			
			if($do_index){
				if($cfg->exists("Aligner","indexname") ne ""){
					print STDERR date() . " Indexing ".$cfg->val("Aligner","fasta") ." for a ".$cfg->val("Aligner","aligner")." analysis\n";
					my @index_value=IndexGeneration(
					  	aligner=>$cfg->val("Aligner","aligner"),
					  	fasta=>$cfg->val("Aligner","fasta"),
					  	dir=>$cfg->val("General","output_dir"),
						logfile=>$log_file || $cfg->val("General","logfile"),
					  	indexname=>$cfg->val("Aligner","indexname"),
						miARmaPath=>$miARmaPath,
						threads=>$cfg->val("General","threads") || 1,
						gtf=>$cfg->val("Aligner","gtf") || undef,
						
					 );
					 
 					print STDERR date() . " Index $index_value[0] created\n";

					if(lc($cfg->val("Aligner","aligner")) eq "bowtie1" and $cfg->val("Aligner","bowtie1index") eq ""){
			 			$cfg->newval("Aligner", "bowtie1index", $index_value[0]);
			 			$cfg->RewriteConfig;
					}
					if(lc($cfg->val("Aligner","aligner")) eq "bowtie2" and $cfg->val("Aligner","bowtie2index") eq ""){
			 			$cfg->newval("Aligner", "bowtie2index", $index_value[0]);
			 			$cfg->RewriteConfig;
					}
					if(lc($cfg->val("Aligner","aligner")) eq "bwa" and $cfg->val("Aligner","bwaindex") eq ""){
			 			$cfg->newval("Aligner", "bwaindex", $index_value[0]);
			 			$cfg->RewriteConfig;
					}
					if(lc($cfg->val("Aligner","aligner")) eq "star" and $cfg->val("Aligner","starindex") eq ""){
			 			$cfg->newval("Aligner", "starindex", $index_value[0]);
			 			$cfg->RewriteConfig;
					}
					if(lc($cfg->val("Aligner","aligner")) eq "hisat2" and $cfg->val("Aligner","hisat2index") eq ""){
			 			$cfg->newval("Aligner", "hisat2index", $index_value[0]);
			 			$cfg->RewriteConfig;
					}
				 }
				 else{
		 			print STDERR "\nERROR " . date() . " You are requesting a index but no indexname has been provided. Please check documentation\n";
		 			help_check_index();
				 }
			}
			
			my @files;
			if(lc($cfg->val("General","type")) eq "mirna"){
				#Reading CutAdapt results directory, collecting the files and completing with the path
				my $cut_dir=$cfg->val("General","output_dir")."/cutadapt_results/";
				if(-e $cut_dir){
					opendir(CUTDIR, $cut_dir) || warn "Aligner:: Folder $cut_dir is not found, but cutadapt has been specified as an adaptersoft\n"; 
					my @cut_files= readdir(CUTDIR);
					push(@files,map("$cut_dir$_",@cut_files));
				}
				#Reading Reaper results directory, collecting the files and completing with the path
				my $rea_dir=$cfg->val("General","output_dir")."/Reaper_results/";
				if(-e $rea_dir){
					opendir(READIR, $rea_dir) || warn "Aligner:: Folder $rea_dir is not found, but Reaper has been specified as an adaptersoft\n"; 
					my @rea_files= readdir(READIR);
					push(@files,map("$rea_dir$_",@rea_files));
				}
				#Reading Trimming results directory, collecting the files and completing with the path
				my $trim_dir=$cfg->val("General","output_dir")."/AdaptTriming_results/";
				if(-e $trim_dir){
					opendir(READIR, $trim_dir) || warn "Aligner:: Folder $trim_dir is not found, but Adaptrimming has been specified as an adaptersoft\n"; 
					my @trim_files= readdir(READIR);
					push(@files,map("$trim_dir$_",@trim_files));
				}
				if(!-e $cut_dir and !-e $rea_dir and !-e $trim_dir){
					print STDERR date()." No processed files are found [neither cutadapt, nor reaper nor adaptrimming folders], assuming " .$cfg->val("General","read_dir").  " don't need to be proccessed\n";
					my $dir_reads_all=$cfg->val("General","read_dir");
					if(-e $dir_reads_all){
						opendir(READIR, $dir_reads_all) || die " No reads found to process\n";
						my @dir_reads_all= readdir(READIR);
						push(@files,map("$dir_reads_all$_",@dir_reads_all));
					}
					else{
						print STDERR "Sorry but I couldn't find any read to align in : $cut_dir or in $rea_dir or in $trim_dir or in " . $cfg->val("General","read_dir") ."\n";
						help_check_general();
					}
				}
			}
			elsif(lc($cfg->val("General","type")) eq "mrna"){
				#Mandatory for mRNA analysis
				#if($cfg->exists("Aligner","tophat_aligner") eq "" or ($cfg->val("Aligner","aligner") eq "")){
				if( lc($cfg->val("Aligner","aligner")) eq "tophat" and ( $cfg->val("Aligner","tophat_aligner") eq "")){
					print STDERR "\nERROR " . date() . " aligner parameter in Section [Aligner] is missing/unfilled. Please check documentation\n";
					help_check_aligner();
				}
				elsif( lc($cfg->val("Aligner","tophat_aligner")) eq "bowtie1" and ( $cfg->val("Aligner","bowtie1index") eq "" and $cfg->val("Aligner","fasta") eq "")){
					print STDERR "\nERROR " . date() . " Bowtie1 has been selected as TopHat aligner but bowtie1index/fasta is missing/unfilled. Please check documentation\n";
					help_check_aligner();
				}
				elsif( lc($cfg->val("Aligner","tophat_aligner")) eq "bowtie2" and ( $cfg->val("Aligner","bowtie2index") eq "" and $cfg->val("Aligner","fasta") eq "")){
					print STDERR "\nERROR " . date() . " Bowtie2 has been selected as TopHat aligner but bowtie2index/fasta is missing/unfilled. Please check documentation\n";
					help_check_aligner();
				}
				elsif( lc($cfg->val("Aligner","aligner")) eq "hisat2" and ( $cfg->val("Aligner","hisat2index") eq "" and $cfg->val("Aligner","fasta") eq "")){
					print STDERR "\nERROR " . date() . " hisat2 has been selected as aligner but hisat2index/fasta is missing/unfilled. Please check documentation\n";
					help_check_aligner();
				}
				elsif( lc($cfg->val("Aligner","aligner")) eq "star" and ( $cfg->val("Aligner","starindex") eq "" and $cfg->val("Aligner","fasta") eq "")){
					print STDERR "\nERROR " . date() . " star has been selected as aligner but starindex/fasta is missing/unfilled. Please check documentation\n";
					help_check_aligner();
				}
				elsif( lc($cfg->val("Aligner","aligner")) eq "tophat" and ($cfg->exists("Aligner","gtf") eq "" or $cfg->val("Aligner","gtf") eq "")){
					print STDERR "\nERROR " . date() . " gtf parameter in Section [Aligner] is missing for tophat. Please check documentation\n";
					help_check_aligner();
				}
								
				#Reading read directory, collecting the files and completing with the path
				my $cut_dir=$cfg->val("General","output_dir")."/cutadapt_results/";
				if(-e $cut_dir){
					opendir(CUTDIR, $cut_dir) || warn "Aligner:: Folder $cut_dir is not found, but cutadapt has been specified as an adaptersoft\n"; 
					my @cut_files= readdir(CUTDIR);
					push(@files,map("$cut_dir$_",@cut_files));
				}
				#Reading Reaper results directory, collecting the files and completing with the path
				my $rea_dir=$cfg->val("General","output_dir")."/Reaper_results/";
				if(-e $rea_dir){
					opendir(READIR, $rea_dir) || warn "Aligner:: Folder $rea_dir is not found, but Reaper has been specified as an adaptersoft\n"; 
					my @rea_files= readdir(READIR);
					push(@files,map("$rea_dir$_",@rea_files));
				}
				#Reading Trimming results directory, collecting the files and completing with the path
				my $trim_dir=$cfg->val("General","output_dir")."/AdaptTriming_results/";
				if(-e $trim_dir){
					opendir(READIR, $trim_dir) || warn "Aligner:: Folder $trim_dir is not found, but Adaptrimming has been specified as an adaptersoft\n"; 
					my @trim_files= readdir(READIR);
					push(@files,map("$trim_dir$_",@trim_files));
				}
				if(!-e $cut_dir and !-e $rea_dir and !-e $trim_dir){
					print STDERR date()." No processed files are found [neither cutadapt, nor reaper nor adaptrimming folders], assuming " .$cfg->val("General","read_dir").  " are already processed\n";
					my $dir_reads_all=$cfg->val("General","read_dir");
					if(-e $dir_reads_all){
						opendir(READIR, $dir_reads_all) || die " No reads found to process\n";
						my @dir_reads_all= readdir(READIR);
						push(@files,map("$dir_reads_all$_",@dir_reads_all));
					}
					else{
						print STDERR "Sorry but I couldn't find any read to align in : $cut_dir or in $rea_dir or in $trim_dir or in " . $cfg->val("General","read_dir") ."\n";
						help_check_general();
					}
				}
			}
			elsif(lc($cfg->val("General","type")) eq "circrna"){
				#Mandatory for circRNA analysis
				if($cfg->exists("Aligner","aligner") eq "" or ($cfg->val("Aligner","aligner") eq "")){
					print STDERR "\nERROR " . date() . " aligner parameter in Section [Aligner] is missing/unfilled. Please check documentation\n";
					help_check_aligner();
				}
				elsif( lc($cfg->val("Aligner","aligner")) eq "bwa" and ( $cfg->val("Aligner","bwaindex") eq "" and $cfg->val("Aligner","fasta") eq "")){
					print STDERR "\nERROR " . date() . " Bowtie1 has been selected as aligner but bowtie1index/fasta is missing/unfilled. Please check documentation\n";
					help_check_aligner();
				}
				
				#Reading read directory, collecting the files and completing with the path
				my $cut_dir=$cfg->val("General","output_dir")."/cutadapt_results/";
				if(-e $cut_dir){
					opendir(CUTDIR, $cut_dir) || warn "Aligner:: Folder $cut_dir is not found, but cutadapt has been specified as an adaptersoft\n"; 
					my @cut_files= readdir(CUTDIR);
					push(@files,map("$cut_dir$_",@cut_files));
				}
				#Reading Reaper results directory, collecting the files and completing with the path
				my $rea_dir=$cfg->val("General","output_dir")."/Reaper_results/";
				if(-e $rea_dir){
					opendir(READIR, $rea_dir) || warn "Aligner:: Folder $rea_dir is not found, but Reaper has been specified as an adaptersoft\n"; 
					my @rea_files= readdir(READIR);
					push(@files,map("$rea_dir$_",@rea_files));
				}
				#Reading Trimming results directory, collecting the files and completing with the path
				my $trim_dir=$cfg->val("General","output_dir")."/AdaptTriming_results/";
				if(-e $trim_dir){
					opendir(READIR, $trim_dir) || warn "Aligner:: Folder $trim_dir is not found, but Adaptrimming has been specified as an adaptersoft\n"; 
					my @trim_files= readdir(READIR);
					push(@files,map("$trim_dir$_",@trim_files));
				}
				if(!-e $cut_dir and !-e $rea_dir and !-e $trim_dir){
					print STDERR date()." No processed files are found [neither cutadapt, nor reaper nor adaptrimming folders], assuming " .$cfg->val("General","read_dir").  " are already processed\n";
					my $dir_reads_all=$cfg->val("General","read_dir");
					if(-e $dir_reads_all){
						opendir(READIR, $dir_reads_all) || die " No reads found to process\n";
						my @dir_reads_all= readdir(READIR);
						push(@files,map("$dir_reads_all$_",@dir_reads_all));
					}
					else{
						print STDERR "Sorry but I couldn't find any read to align in : $cut_dir or in $rea_dir or in $trim_dir or in " . $cfg->val("General","read_dir") ."\n";
						help_check_general();
					}
				}
			}
			my $libray_type="fr-firststrand";
			if(lc($cfg->val("General","strand")) eq "no"){
				$libray_type="fr-unstranded";
			}
			if(lc($cfg->val("General","strand")) eq "reverse"){
				$libray_type="fr-secondstrand";
			}
			check_input_format(-files=>\@files,-dir=>"",-log=>$log_file,config=>$cfg);
			
			if(scalar(@files)>0){
				my @alignes;
				print STDERR date()." Starting a \"".$cfg->val("Aligner","aligner")."\" Alignment Analysis\n";
				# Reading the array with the names of the files
				foreach my $file(sort @files){
					my $aligned_file=ReadAligment(
					    file=>$file,
					    aligner=>$cfg->val("Aligner","aligner"),
						threads=>$cfg->val("General","threads") || 1,
					    bowtie2index=>$cfg->val("Aligner","bowtie2index") || undef,
					    bowtie1index=>$cfg->val("Aligner","bowtie1index") || undef,
					    bwaindex=>$cfg->val("Aligner","bwaindex") || undef,
					    hisat2index=>$cfg->val("Aligner","hisat2index") || undef,
					    starindex=>$cfg->val("Aligner","starindex") || undef,
						logfile=>$log_file||$cfg->val("General","logfile"),
						statsfile=>$stat_file|| $cfg->val("General","stats_file"),
						verbose=>$cfg->val("General","verbose") || 0,
					    bowtiemiss=>$cfg->val("Aligner","bowtiemiss") || undef,
					    bowtieleng=>$cfg->val("Aligner","bowtieleng") || undef,
						projectdir=>$cfg->val("General","output_dir")|| undef,
					    bowtie1parameters=>$cfg->val("Aligner","bowtie1parameters")|| undef,
					    bowtie2parameters=>$cfg->val("Aligner","bowtie2parameters")|| undef,
						bwaparameters=>$cfg->val("Aligner","bwaparameters") || undef,
						tophatParameters=>$cfg->val("Aligner","tophatParameters") || undef,
						hisat2parameters=>$cfg->val("Aligner","hisat2parameters") || undef,
						starparameters=>$cfg->val("Aligner","starparameters") || undef,
						miRDeeparameters=>$cfg->val("Aligner","miRDeeparameters") || undef,
					    miARmaPath=>$miARmaPath,
						organism=>$cfg->val("General","organism")|| undef,
						adapter=>$cfg->val("Adapter","adaptersoft")|| undef,
						GTF=>$cfg->val("Aligner","gtf") || undef,
						Seqtype=>$cfg->val("General","seqtype") || "Single",
						tophat_seg_mismatches=>$cfg->val("Aligner","tophat_seg_mismatches") || undef,
						tophat_seg_length=>$cfg->val("Aligner","tophat_seg_length") || undef,
						library_type=>$libray_type || "fr-firststrand",
						tophat_multihits=>$cfg->val("Aligner","tophat_multihits") || undef,
						read_mismatches=>$cfg->val("Aligner","tophat_read_mismatches") || undef,
						tophat_aligner=>$cfg->val("Aligner","tophat_aligner") || undef,
						strand=>$cfg->val("General","strand")|| "yes",
					);
					push(@alignes,$file);
				}
				
				ReadSummary(
				    aligner=>$cfg->val("Aligner","aligner"),
					tophat_aligner=>$cfg->val("Aligner","tophat_aligner") || undef,
					summary=>$summary_file,
					statsfile=>$stat_file|| $cfg->val("General","stats_file"),
					projectdir=>$cfg->val("General","output_dir")|| undef,
				);
				print STDERR date()." \"".$cfg->val("Aligner","aligner")."\" Alignment Analysis finished\n";
				
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
			if(lc($cfg->val("General","type")) ne "circrna"){
				my $bw1_dir=$cfg->val("General","output_dir")."/Bowtie1_results/";
				if(-e $bw1_dir){
					opendir(BW1DIR, $bw1_dir) || warn "3 Aligner:: Folder $bw1_dir is not found\n"; 
					my @bw1_files= readdir(BW1DIR);
					push(@files,map("$bw1_dir$_",@bw1_files));
					close BW1DIR;
				}
				#Reading Reaper results directory, collecting the files and completing with the path
				my $bw2_dir=$cfg->val("General","output_dir")."/Bowtie2_results/";
				if(-e $bw2_dir){
					opendir(BW2DIR, $bw2_dir) || warn "4 Aligner:: Folder $bw2_dir is not found\n"; 
					my @bw2_files= readdir(BW2DIR);
					push(@files,map("$bw2_dir$_",@bw2_files));
					close BW2DIR;
				}
				#Reading results directory, collecting the files and completing with the path
				my $his_dir=$cfg->val("General","output_dir")."/hisat2_results/";
				if(-e $his_dir){
					opendir(HISDIR, $his_dir) || warn "5 Aligner:: Folder $his_dir is not found\n"; 
					my @his_files= readdir(HISDIR);
					push(@files,map("$his_dir$_",@his_files));
					close HISDIR;
				}
				#Reading results directory, collecting the files and completing with the path
				my $star_dir=$cfg->val("General","output_dir")."/star_results/";
				if(-e $star_dir){
					opendir(STARDIR, $star_dir) || warn "5 Aligner:: Folder $star_dir is not found\n"; 
					my @star_files= readdir(STARDIR);
					push(@files,map("$star_dir$_",@star_files));
					close STARDIR;
				}
				if(scalar(@files)>0 and lc($cfg->val("General","type")) ne "circrna"){
					print STDERR date()." Starting a Readcount Analysis\n";
					my @htseqfiles;
					# Reading the array with the names of the files
					foreach my $file(@files){
						#Selecting only the sam files for their processing
						my $result=featureCount(
							file=>$file,
							database=>$cfg->val("ReadCount","database"),
							seqid=>$cfg->val("ReadCount","seqid") || undef,
							parameters=>$cfg->val("ReadCount","parameters") || undef, 
							strand=>$cfg->val("General","strand") || "yes", 
							featuretype=>$cfg->val("ReadCount","featuretype") || undef,
							logfile=>$log_file || $cfg->val("General","logfile"),
							verbose=>$cfg->val("General","verbose") || 0,
							projectdir=>$cfg->val("General","output_dir")|| undef,
							miARmaPath=>$miARmaPath,
							threads=>$cfg->val("General","threads") || 1,
							quality=>$cfg->val("ReadCount","quality") || undef,
							Seqtype=>$cfg->val("General","seqtype") || "Single",
						);
						push(@htseqfiles, $result);
					}
				
					#HTSEQFORMATEXECUTION
					featureFormat( 
					  	input=>\@htseqfiles, 
					  	projectdir=>$cfg->val("General","output_dir")|| undef,
						logfile=>$log_file || $cfg->val("General","logfile"),
						verbose=>$cfg->val("General","verbose") || 0,
					  );
					  featureSummary(
				  		projectdir=>$cfg->val("General","output_dir")|| undef,
						logfile=>$log_file || $cfg->val("General","logfile"),
					  	summary=>$summary_file,
					  	input=>\@htseqfiles, 
					  );
					
  					print STDERR date()." Readcount Analysis finished.\n";
				}
				else{
					print "ERROR :: You are requesting a miRNA readcount analysis, but no aligned files are found (Neither Bowtie1 nor Bowtie2)\n";
					help_check_count();
				}
			}
			elsif(lc($cfg->val("General","type")) eq "circrna"){
				if($cfg->exists("ReadCount","fasta") eq "" or ($cfg->val("ReadCount","fasta") eq "")){
					print STDERR "\nERROR " . date() . " fasta parameter in Section [Readcount] for circRNAs is missing/unfilled. Please check documentation\n";
					$severe_error=1;			
					help_check_count();
				}
				else{
					#Reading CutAdapt results directory, collecting the files and completing with the path
					if($cfg->exists("ReadCount","method") eq ""){
						$cfg->newval("ReadCount", "method", "CIRI1");
						$cfg->RewriteConfig;
					}
					if(lc($cfg->val("ReadCount","method")) ne "ciri1" and lc($cfg->val("ReadCount","method")) ne "ciri2"){
						print STDERR date() . " WARN :: Method for identifying circRNAs is not valid. Using CIRI1.2 as default\n";
						$cfg->newval("ReadCount", "method", "CIRI1");
						$cfg->RewriteConfig;
					}
					my $method=$cfg->val("ReadCount","method");
					
					my $bw2_dir=$cfg->val("General","output_dir")."/bwa_results/";
					if($bw2_dir){
						opendir(BW2DIR, $bw2_dir) || warn "5 Aligner:: Folder $bw2_dir is not found\n"; 
						my @bwa_files= readdir(BW2DIR);
						push(@files,map("$bw2_dir$_",@bwa_files));
						close BW2DIR;
						
						my @circRNAfiles;
						print STDERR date()." Starting a Readcount Analysis\n";
						
						if(lc($method) eq "ciri1"){
							foreach my $file(@files){
								#Selecting only the sam files for their processing
								my $result=CIRICount1(
								  	file=>$file,
									database=>$cfg->val("ReadCount","database")|| undef,
									logfile=>$log_file || $cfg->val("General","logfile"),
									verbose=>$cfg->val("General","verbose") || 0,
									projectdir=>$cfg->val("General","output_dir")|| undef,
									threads=>$cfg->val("General","threads") || 1,
									miARmaPath=>$miARmaPath,
									Seqtype=>$cfg->val("General","seqtype") || "Single",
									fasta=>$cfg->val("ReadCount","fasta") || undef,
									method=>"CIRI1"
								);
								push(@circRNAfiles, $result);
							}
						}
						if(lc($method) eq "ciri2"){
							foreach my $file(@files){
								#Selecting only the sam files for their processing
								my $result=CIRICount2(
								  	file=>$file,
									database=>$cfg->val("ReadCount","database")|| undef,
									logfile=>$log_file || $cfg->val("General","logfile"),
									verbose=>$cfg->val("General","verbose") || 0,
									projectdir=>$cfg->val("General","output_dir")|| undef,
									threads=>$cfg->val("General","threads") || 1,
									miARmaPath=>$miARmaPath,
									Seqtype=>$cfg->val("General","seqtype") || "Single",
									fasta=>$cfg->val("ReadCount","fasta") || undef,
									method=>"CIRI2"
								);
								push(@circRNAfiles, $result);
							}
						}
						
						CIRIFormat( 
						  	input=>\@circRNAfiles, 
							projectdir=>$cfg->val("General","output_dir")|| undef,
							logfile=>$log_file || $cfg->val("General","logfile"),
							verbose=>$cfg->val("General","verbose") || 0,
							summary=>$summary_file,
						);
						@files=undef;
						
	  					print STDERR date()." Readcount Analysis finished.\n";
						
					}
				}
			}
			else{
				print "ERROR :: You are requesting a readcount analysis, but no aligned files are found (Neither Bowtie1 nor Bowtie2 ... )\n";
			 	help_check_aligner();
			}
		}
	}
	#Denovo
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
		
		else{
			#run Aligner
			use CbBio::RNASeq::Aligner;
			use CbBio::RNASeq::Readcount;
			
			my $read_directory=$cfg->val("General","read_dir") . "/";
			my @fastaq_files;
			if($read_directory){
				opendir(READIR, $read_directory) || warn "DeNovo:: Folder $read_directory is not found\n"; 
				my @files= readdir(READIR);
				push(@fastaq_files,map("$read_directory$_",@files));
			}
			print date() . " Starting a De novo identification and quantification of miRNAs\n";
			foreach my $file(@fastaq_files){
				ReadAligment(
					file=>$file,
					aligner=>"miRDeep",
					threads=>$cfg->val("General","threads") || 1,
				    bowtie1index=>$cfg->val("DeNovo","bowtie1index") || undef,
					statsfile=>$stat_file || $cfg->val("General","stats_file"),
					verbose=>$cfg->val("General","verbose")|| 0,
					projectdir=>$cfg->val("General","output_dir"),
					logfile=>$log_file || $cfg->val("General","logfile"),
					miARmaPath=>$miARmaPath,
					adapter=>$cfg->val("DeNovo","adapter") || undef,
					adapter_file=>$cfg->val("DeNovo","adapter_file") || undef,
					precursors=>$cfg->val("DeNovo","precursor_miRNA_file"),
					mature=>$cfg->val("DeNovo","mature_miRNA_file"),
					genome=>$cfg->val("DeNovo","genome"),
					organism=>$cfg->val("General","organism"),
				);
			}
			#Once ReadAlignment is finished, is time to count the reads
			my $processed_reads=$cfg->val("General","output_dir") . "/miRDeep_results/";
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
						logfile=>$log_file || $cfg->val("General","logfile"),
						verbose=>$cfg->val("General","verbose")|| 0,
						projectdir=>$cfg->val("General","output_dir"),
						miARmaPath=>$miARmaPath,
					);
					push(@allRNAfiles, $result) if($result);
				}
			}
			if(scalar(@allRNAfiles)>0){
				$result_file=miRDeepFormat( 
				  	input=>\@allRNAfiles, 
					projectdir=>$cfg->val("General","output_dir"),
					logfile=>$log_file || $cfg->val("General","logfile"),
				);
			}
			if(scalar(@allRNAfiles)>0){
				$result_file=miRDeepSummary( 
				  	input=>\@allRNAfiles, 
					projectdir=>$cfg->val("General","output_dir"),
					summary=>$summary_file,
				);
			}
		}
		print date() . " De novo identification and quantification of miRNAs finished\n";
		
	}	
	#Diferential Expression
	if($cfg->SectionExists("DEAnalysis")==1){
		#Mandatory parameters: read folder
		if($cfg->exists("DEAnalysis","targetfile") eq "" or ($cfg->val("DEAnalysis","targetfile") eq "")){
			print STDERR "\nERROR " . date() . " targetfile parameter in Section [DEAnalysis] is missing/unfilled. Please check documentation\n";
			help_check_deanalysis();
		}
		elsif( lc($cfg->val("DEAnalysis","contrastfile")) eq "" and ($cfg->val("DEAnalysis","contrastfile") eq "")){
			print STDERR "\nERROR " . date() . " contrastfile parameter in Section [DEAnalysis] is missing/unfilled. Please check documentation\n";
			help_check_deanalysis();
		}
		# elsif( lc($cfg->val("DEAnalysis","filter")) eq "" and ($cfg->val("DEAnalysis","filter") eq "")){
		# 	print STDERR "\nERROR " . date() . " filter parameter in Section [DEAnalysis] is missing/unfilled. Please check documentation\n";
		# 	help_check_aligner();
		# }
		# elsif( lc($cfg->val("DEAnalysis","desoft")) eq "" and ($cfg->val("DEAnalysis","desoft") eq "")){
		# 	$cfg->newval("DEAnalysis", "desoft", "edgeR");
		# 	$cfg->RewriteConfig;
		# }
		else{
			
			if( lc($cfg->val("DEAnalysis","desoft")) eq "" and ($cfg->val("DEAnalysis","desoft") eq "")){
				$cfg->newval("DEAnalysis", "desoft", "edgeR");
				$cfg->RewriteConfig;
			}
			#run DEAnalysis
			use CbBio::RNASeq::DEAnalysis;
			print STDERR date()." Starting a differential expression analysis using ".$cfg->val("DEAnalysis","desoft") ." software(s)\n";
			
			my @files;
			#if user did a regular readcount
			my $read_count_dir=$cfg->val("General","output_dir")."/Readcount_results/";
			if(-e $read_count_dir){
				opendir(CUTDIR, $read_count_dir) || warn "DEAnalysis:: Folder $read_count_dir is not found\n"; 
				my @cut_files= readdir(CUTDIR);
				push(@files,map("$read_count_dir/$_",@cut_files));
				
			}
			#if user did a denovo analysis
			my $denovo_dir=$cfg->val("General","output_dir")."/DeNovo_ReadCount/"; 
			if(-e $denovo_dir){
				opendir(READIR, $denovo_dir) ||  warn "DEAnalysis:: Folder $denovo_dir is not found\n"; 
				my @rea_files= readdir(READIR);
				push(@files,map("$denovo_dir/$_",@rea_files));
			}

			#if user did a circRNA analysis
			my $circRNA_dir=$cfg->val("General","output_dir")."/circRNAs_results/"; 
			if(-e $circRNA_dir){
				opendir(READIR, $circRNA_dir) ||  warn "DEAnalysis:: Folder $circRNA_dir is not found\n"; 
				my @rea_files= readdir(READIR);
				push(@files,map("$circRNA_dir/$_",@rea_files));
			}
			
			if(scalar(@files)>0){
				foreach my $file(  sort @files){
					if($file !~ /Size/){
						DE_Analysis(
							projectdir=>$cfg->val("General","output_dir")|| undef,
							file=>$file,
							targetfile=>$cfg->val("DEAnalysis","targetfile"),
							label=>$cfg->val("General","label"),
							filter=>$cfg->val("DEAnalysis","filter")|| "yes",
							edger_contrastfile=>$cfg->val("DEAnalysis","contrastfile")|| undef,
							noiseq_contrastfile=>$cfg->val("DEAnalysis","contrastfile")|| undef,
							DEsoft=>$cfg->val("DEAnalysis","desoft")|| undef,
							filtermethod=>$cfg->val("DEAnalysis","filtermethod")|| undef,
							logfile=>$log_file || $cfg->val("General","logfile"),
							verbose=>$cfg->val("General","verbose") || 0,
							cpmvalue=>$cfg->val("DEAnalysis","cpmvalue")|| undef,
							repthreshold=>$cfg->val("DEAnalysis","repthreshold")|| undef,
							edger_normethod=>$cfg->val("DEAnalysis","edger_normethod")|| undef,
							noiseq_normethod=>$cfg->val("DEAnalysis","noiseq_normethod")|| undef,
							replicates=>$cfg->val("DEAnalysis","replicates") || undef,
							miARmaPath=>$miARmaPath,
							repthreshold=>$cfg->val("DEAnalysis","repthreshold") || undef,
							bcvvalue=>$cfg->val("DEAnalysis","bcvvalue") || undef,
							cpmvalue=>$cfg->val("DEAnalysis","cpmvalue") || undef,
							cutoffvalue=>$cfg->val("DEAnalysis","cutoffvalue") || undef,
							Rdir=>$cfg->val("DEAnalysis","rdir") || undef,
					     	replicatevalue=>$cfg->val("DEAnalysis","replicatevalue") || undef,
							kvalue=>$cfg->val("DEAnalysis","kvalue") || undef,
							lcvalue=>$cfg->val("DEAnalysis","lcvalue") || undef,
							pnrvalue=>$cfg->val("DEAnalysis","pnrvalue") || undef,
							nssvalue=>$cfg->val("DEAnalysis","nssvalue") || undef,
							vvalue=>$cfg->val("DEAnalysis","vvalue") || undef,
							qvalue=>$cfg->val("DEAnalysis","qvalue") || undef,
							rpkm=>$cfg->val("DEAnalysis","rpkm") || undef,
							cpm=>$cfg->val("DEAnalysis","cpm") || undef
						);
					}
				}
				#print STDERR date()." Differential expression Summary.\n";
				DE_AnalysisSummary(
				  	projectdir=>$cfg->val("General","output_dir")|| undef,
					DEsoft=>$cfg->val("DEAnalysis","desoft")|| undef,
					summary=>$summary_file,
					contrastfile=>$cfg->val("DEAnalysis","contrastfile"),
				);
			}
			else{
				print "ERROR :: Please check that your tab files are saved in: ($read_count_dir / $denovo_dir)\n";
				help_check_deanalysis();
			}
			
			print STDERR date()." Differential expression analysis finished.\n";
			
		}
	}
	#Functional Analysis
	if($cfg->SectionExists("FAnalysis")==1){
		 
		use CbBio::RNASeq::FAnalysis;
		 F_Analysis(
		 			projectdir=>$cfg->val("General","output_dir"),
		 			verbose=>$cfg->val("General","verbose") || 0,
		 			logfile=>$log_file || $cfg->val("General","logfile"),
		 			organism=>$cfg->val("General","organism")|| "human",
		 			seqid=>$cfg->val("FAnalysis","seqid") || $cfg->val("Readcount","seqid") ,
		 			edger_cutoff=>$cfg->val("FAnalysis","edger_cutoff") || 0.05,
		 			noiseq_cutoff=>$cfg->val("FAnalysis","noiseq_cutoff") || 0.8,
					miARmaPath=>$miARmaPath,
					fc_threshold=>$cfg->val("FAnalysis","fc_threshold")|| 0,					
		 );
		 
		 F_AnalysisSummary(
			projectdir=>$cfg->val("General","output_dir"),
			summary=>$summary_file,
		 );
		 
	}
	#miRGate prediction
	if($cfg->SectionExists("TargetPrediction")==1){
		#run TargetP
		my $dir=$cfg->val("General","output_dir");		
		use CbBio::RNASeq::TargetPrediction;
		print date()." Starting a target Prediction Analysis using miRGate\n";
		
		if(lc($cfg->val("General","type")) eq "mirna" and $cfg->val("TargetPrediction","genes_folder") eq ""){
			TargetPrediction(
				miRNAs_folder=>$dir,
				logfile=>$log_file || $cfg->val("General","logfile"),
				projectdir=>$cfg->val("General","output_dir")|| undef,
				organism=>$cfg->val("General","organism")|| undef,
				miARmaPath=>$miARmaPath,
				edger_cutoff=>$cfg->val("TargetPrediction","edger_cutoff")|| undef,
				noiseq_cutoff=>$cfg->val("TargetPrediction","noiseq_cutoff")|| undef,
				fc_threshold=>$cfg->val("TargetPrediction","fc_threshold")|| 0,
				verbose=>$cfg->val("General","verbose") || 0,
			);
		}
		elsif(lc($cfg->val("General","type")) eq "mirna" and $cfg->val("TargetPrediction","genes_folder") ne ""){
			TargetPrediction(
				miRNAs_folder=>$dir,
				genes_folder=>$cfg->val("TargetPrediction","genes_folder"),
				logfile=>$log_file || $cfg->val("General","logfile"),
				projectdir=>$cfg->val("General","output_dir")|| undef,
				organism=>$cfg->val("General","organism")|| undef,
				miARmaPath=>$miARmaPath,
				edger_cutoff=>$cfg->val("TargetPrediction","edger_cutoff")|| undef,
				noiseq_cutoff=>$cfg->val("TargetPrediction","noiseq_cutoff")|| undef,
				fc_threshold=>$cfg->val("TargetPrediction","fc_threshold")|| 0,
				verbose=>$cfg->val("General","verbose") || 0,
				method=>$cfg->val("TargetPrediction","method")|| "Pearson"
			);
		}
		elsif(lc($cfg->val("General","type")) eq "mrna" and $cfg->val("TargetPrediction","miRNAs_folder") eq ""){
			TargetPrediction(
				genes_folder=>$dir,
				logfile=>$log_file || $cfg->val("General","logfile"),
				projectdir=>$cfg->val("General","output_dir")|| undef,
				organism=>$cfg->val("General","organism")|| undef,
				miARmaPath=>$miARmaPath,
				edger_cutoff=>$cfg->val("TargetPrediction","edger_cutoff")|| undef,
				noiseq_cutoff=>$cfg->val("TargetPrediction","noiseq_cutoff")|| undef,
				fc_threshold=>$cfg->val("TargetPrediction","fc_threshold")|| 0,
				verbose=>$cfg->val("General","verbose") || 0,		
						
			);
		}
		elsif(lc($cfg->val("General","type")) eq "mrna" and $cfg->val("TargetPrediction","miRNAs_folder") ne "") {
			TargetPrediction(
				genes_folder=>$dir,
				miRNAs_folder=>$cfg->val("TargetPrediction","miRNAs_folder"),
				logfile=>$log_file || $cfg->val("General","logfile"),
				projectdir=>$cfg->val("General","output_dir")|| undef,
				organism=>$cfg->val("General","organism")|| undef,
				miARmaPath=>$miARmaPath,
				edger_cutoff=>$cfg->val("TargetPrediction","edger_cutoff")|| undef,
				noiseq_cutoff=>$cfg->val("TargetPrediction","noiseq_cutoff")|| undef,
				verbose=>$cfg->val("General","verbose") || 0,
				fc_threshold=>$cfg->val("TargetPrediction","fc_threshold")|| 0,
				method=>$cfg->val("TargetPrediction","method")|| "Pearson"				
			);
		}
		
		TargetPrediction_Summary(
		  	projectdir=>$cfg->val("General","output_dir")|| undef,
			summary=>$summary_file
		);
		print date()." Target Prediction Analysis finished.\n";
	}
	
	my $end_run = time();
	my $run_time = $end_run - $start_run;
	print date() . " miARma finished. Job took ". sprintf("%d",$run_time/60) ." minutes\n";
 
sub check_input_data{
	my %args=@_;
	my $cfg=$args{"config"};
	my $check_input=$args{"check"};
	my $type=$args{"type"};
	my $parameter_type=parameter_definition();
	my @sections=$cfg->Sections();
	my @sections_tmp=@sections;
	
	my @print_section=splice(@sections_tmp,1,scalar(@sections_tmp));
	print date()." Starting a miARma analysis for $type\n";
	print date()." Checking provided parameters for: ".join(",",@print_section).". \n";

	foreach my $section (@sections){
		my @parametros=$cfg->Parameters($section);
		foreach my $parametro (@parametros){
			
			if($parameter_type->{$parametro} eq "path"){
				if(-e $cfg->val($section,$parametro)){
				}
				else{
					print date()." Checking $section-$parametro parameter ... error!\n\tERROR Please check that the parameter $parametro inside section [$section] is correct.\n";
					exit;
				}
			}
			if($parameter_type->{$parametro} eq "path_dir" and !$check_input){
				if(-e $cfg->val($section,$parametro)){
					print date()." Checking $section-$parametro parameter ... Exists!\n".date()." The folder specified in ($parametro=".$cfg->val($section,$parametro).") already exists.\n". date(). " Wait 5 seconds to overwrite the folder. Cancel otherwise:\t";
					for (my $var = 5; $var >=0; $var--){
						print "\b" . $var;
						sleep(1);
					}
					print "\n" . date(). " Continue.\n";
				}
			}

			if($parameter_type->{$parametro} eq "file"){
				if(-e $cfg->val($section,$parametro)){
				}
				else{
					print date()." Checking $section-$parametro parameter ... error!\n\tERROR Please check that the parameter $parametro inside section [$section] is correct.\n";
					exit;
				}
			}
			if($parameter_type->{$parametro} eq "bw2"){
				if(-e $cfg->val($section,$parametro) .".1.bt2"){
				}
				else{
					print date()." Checking $section-$parametro parameter ... error!\n\tERROR Please check that the parameter $parametro inside section [$section] is correct.\n";
					exit;
				}
			}
			if($parameter_type->{$parametro} eq "bw1"){
				if(-e $cfg->val($section,$parametro) .".1.ebwt"){
				}
				else{
					print date()." Checking $section-$parametro parameter ... error!\n\tERROR Please check that the parameter $parametro inside section [$section] is correct.\n";
					exit;
				}
			}
			if($parameter_type->{$parametro} eq "bwa"){
				if(-e $cfg->val($section,$parametro) .".bwt"){
				}
				else{
					print date()." Checking $section-$parametro parameter ... error!\n\tERROR Please check that the parameter $parametro inside section [$section] is correct.\n";
					exit;
				}
			}
		}
	}
	if($check_input eq "-check" or $check_input eq "--check" or $check_input eq "check"){
		print date()." All parameters are correct.\n";
		exit;
	}
	else{
		print date()." All parameters are correct.\n";
		return();
	}
}

sub check_input_format{
	use File::Basename;
	my %args=@_;
	
	#In case you provide a config file
	my @fastq_files;
	my $cfg=$args{"-config"};#configuration file
	my $logfile=$args{"-log"} || "/dev/null";#log file
	
	open(LOG,">> ".$logfile) || die $!;

	my @files;
	if($args{"-files"}){
		@files=@{$args{"-files"}};
	}
	my $dir_fastq = $args{"-dir"} || "";
	
	if($cfg and scalar(@files)<0){
		if($dir_fastq){
			if(opendir(FASTQF, $dir_fastq)){
				my @files= readdir(FASTQF);
				push(@fastq_files,map("$dir_fastq$_",@files));
			}
			else{
				print STDERR date() ." ERROR: Can't find fastq files on $dir_fastq\n\n";
				exit;
			}
		}

	}
	#In case yo provide an array of files
	if(scalar(@files)>0){
		if($dir_fastq){
			@fastq_files=map("$dir_fastq/$_",@files);
		}
		else{
			@fastq_files=@files;
		}
	}
	my $exit=0;
	if(scalar(@fastq_files)>0){
		print date() . " Checking if files in folder: \"$dir_fastq\" are in the correct fastq format.\n" if($dir_fastq);
		#print date() . " Checking if files are in the correct fastq format\n" if(!$dir_fastq and $cfg->SectionExists("Adapter")!=1);
		my @wrong_file;
		foreach my $file(@fastq_files){
			if($file){
				my ($real_file,$path)=fileparse($file);
				if($real_file !~ /^\./){
					my $error=0;
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
					#print STDERR "$file <> $real_file\n";
					open IN,$file or die "Cannot open FASTQ file ($file)\n";
					my $i=0;
					my $mes="Please make sure your file ($real_file) is in accordance with the FASTQ format specifications";
					while(<IN>){
				    	chomp;
				    	$i++;
						if($i == 1){if(/^@\S+/){}else{
							$error=1;
							#print STDERR "miARmA ERROR :: ".date() . " First line of FASTQ reads file ($real_file) is not in accordance with the fastq format specifications\n$mes\n";
							push(@wrong_file,$real_file);
						}}                 
						if($i == 2){if(/^\S+$/){}else{
							$error=1;
							push(@wrong_file,$real_file);
							#print STDERR "miARmA ERROR :: ".date() . " Second line of FASTQ reads file ($real_file) contains whitespace in sequence\n$mes\n";
						}}    
						if($i == 3){if(/^\+/){}else{
							$error=1;
							push(@wrong_file,$real_file);
							#print STDERR "miARmA ERROR :: ".date() . " Third line of FASTQ reads file ($real_file) does not start with a '+' character.\n$mes\n";
						}}   
						if($i == 4){if(/^\S+$/){}else{
							$error=1;
							push(@wrong_file,$real_file);
							#print STDERR "miARmA ERROR :: ".date() . " Fourth line of FASTQ reads file ($real_file) contains whitespace\n$mes\n";
						}}  
					
						last if($i == 4);
					}
					close IN;					
				
					if($error == 0){
						#No error, so they are fastq files, but is the extension correct ?
						if($real_file =~ /.*\.fastq$/ and $real_file =~ /.*\.fastq\.lane\.clean$/ and $real_file =~ /.*\.fastq.gz$/ and $real_file =~ /.*\.fq$/ and $real_file =~ /.*\.fq.gz$/ and $real_file =~ /.*\.fq\.bz2$/ and $real_file =~ /.*\.fastq\.bz2$/){
							print STDERR date(). " $real_file is a fastq file, but due to cutadapt and other softwares, this file should include a fastq/fq extension. [Comppresed are also valid as .fq.gz]\n";
							$exit=1;
						}
					}
					else{
						print LOG "Warn ". date() . " $file is not in in accordance with the fastq format\n";
						next;
					}
				}
			}
		}
	}
	else{
		print STDERR "Can't find fastq files on $dir_fastq\n";
		exit;
	}
	#Summary
	if($exit==1){
		print STDERR "Please fix the names of your files\nQuitting\n\n";
		exit;
	}
	close LOG;
}	
	
	sub help_check_general{
	    my $usage = qq{
Mandatory parameters:

[General]
read_dir=miRNA_raw_reads_folder/
label=Hypoxia
miARmaPath=.
projectdir=results/
organism=human
type=miRNA
SetType=Paired
strand=yes
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
	exit();
}
	
	sub help_check_aligner{
	    my $usage = qq{
Mandatory parameters:

For miRNAs:

[Aligner]
aligner=Bowtie1
bowtie1index=/genomes/bowtie1/hg19

or

[Aligner]
aligner=Bowtie2
bowtie2index=/genomes/bowtie2/hg19

for mRNAs:
[Aligner]
aligner=tophat
tophat_aligner=Bowtie2
bowtie2index=Genomes/Indexes/bowtie2/human/hg19
gtf=../../../../data/Homo_sapiens_CHR_.GRCh37.74.gtf
strand=no

};

	print STDERR $usage;
	exit();
}
 sub help_check_index{
	    my $usage = qq{
Mandatory parameters:

[Aligner]
aligner=Bowtie1
fasta=/genomes/bowtie1/hg19.fasta
indexname=hg19
or

[Aligner]
aligner=Bowtie2
fasta=/genomes/bowtie1/hg19.fasta
indexname=hg19

};

	print STDERR $usage;
	exit();
 }	
	sub help_check_adapter{
	    my $usage = qq{
Mandatory parameters:

[Adapter]
adaptersoft=CutAdapt

};

	print STDERR $usage;
	exit();
	}
}

	sub help_check_count{
	    my $usage = qq{
Mandatory parameters:

[ReadCount]
database=miRNAs_miRBase20.gft

for circRNAs:
[ReadCount]
database=HS_nbci_37.2.gft
bwaindex=index/bwa
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
	desoft=EdgeR-Noiseq

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
	desoft=EdgeR-Noiseq
	; providing replicates
	replicates=yes

	};

		print STDERR $usage;
		exit();
	}
	
sub print_header{
	system("clear");
	print "#########################################################################	
#   miARma, miRNA and RNASeq Multiprocess Analysis			#
#                miARma v 1.7.5 (May-2018)                              #
#               		                              		#
#   Created at Computational Biology and Bioinformatics Group (CbBio)   #
#   Institute of Biomedicine of Seville. IBIS (Spain)                   #
#   Modified and Updated at Bioinformatics Unit at IPBLN-CSIC   	#
#   Institue for Parasitology and Biomedicine Lopez-Neyra (IPBLN-CSIC). #
#   Granada (Spain)             				        #
#   Copyright (c) 2018 IBIS & IPBLN. All rights reserved.               #
#   mail : miARma-devel\@idoproteins.com                                 #
#########################################################################\n\n";
}
sub date{
	#my $dt = DateTime->now(time_zone=>'local');
	#return($dt->hms . " [" . $dt->dmy ."]");
	use Time::localtime;
	my $now = ctime();
	return("[$now]");
}

sub parameter_definition{
	
	my $parameter;
	$parameter->{"read_dir"}="path";
	$parameter->{"output_dir"}="path_dir";
	$parameter->{"miARmaPath"}="path";
	$parameter->{"adapter_file"}="file";
	$parameter->{"metafile"}="file";
	$parameter->{"bowtie1index"}="bw1";
	$parameter->{"bowtie2index"}="bw2";
	$parameter->{"bwaindex"}="bwa";
	$parameter->{"fasta"}="file";
	$parameter->{"gtf"}="file";
	$parameter->{"database"}="file";
	$parameter->{"adapter_file"}="file";
	$parameter->{"genome"}="file";
	$parameter->{"mature_miRNA_file"}="file";
	$parameter->{"precursor_miRNA_file"}="file";
	$parameter->{"rdir"}="file";
	$parameter->{"contrastfile"}="file";
	$parameter->{"targetfile"}="file";
	$parameter->{"genes_folder"}="path";
	$parameter->{"miRNAs_folder"}="path";
	
	return($parameter);
}
1;
