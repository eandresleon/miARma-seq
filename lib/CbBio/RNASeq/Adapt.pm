#########################################################################	
#	Adapter processing package		 									#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2015 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################

package CbBio::RNASeq::Adapt;
#Export package system
require Exporter;
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(CutAdapt CutAdaptStats AdaptTriming ReadFilter Minion Reaper AdapterRemoval AdapterGraph ReadLengthCount);

use strict;
use DateTime;
use LWP;
use File::Basename;
use Statistics::R;
use Cwd;
use Cwd 'abs_path';


=head1 NAME

 Adapt 

=head1 SYNOPSIS

Adapt package is composed of 9 subroutines: CutAdapt,CutAdaptStats, AdaptTriming, ReadFilter,
Minion, Reaper, AdapterRemoval, AdapterGraph and ReadLengthCount. The aim of this package is to 
process the reads from a fastq file removing a known adapter sequence (CutAdapt, Reaper) or not 
(AdaptTriming), predicting the adapter sequence of the reads (Minion), filtering the reads by length 
(ReadFilter) or printing statistical info after the read processing (CutAdaptStats, AdapterGraph, 
ReadLengthCount). In addition AdapterRemoval acts as main function calling to CutAdapt, CutAdaptStats,  
AdaptTriming, Minion, Reaper, AdapterGraph and ReadLengthCount according to the parameters selected.

=head1 Methods


=head2 AdapterRemoval

  Example    : 
  AdapterRemoval(
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
  Description: AdapterRemoval function performs an adapter removal process with cutadapt, reaper software or both, and/or
  the removal of a defined number of nucleotides of the reads (adapttrimming). These analysis generate new fastq files with 
  the conserved reads at Cutadapt_results, Reaper_results or/and AdaptTrimming_results directory in the provided project directory. 
  In addition, if user do not provide the adapter sequence AdapterRemoval can predict it executing Minion software. 
  After processing the reads AdapterRemoval generates statistical data of the process
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
  AdapterRemoval also return the paths of the new generated files
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

sub AdapterRemoval{
	#Arguments provided by user are collected by %args. 
	my %args=@_;
	my $adaptersoft=$args{"adaptersoft"}; #Specific software to remove the adapter from the sequences
	my $files=$args{"files"}; #Files which is going to be processing
	my $dir= $args{"dir"}; #Input directory where the files are
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $verbose=$args{"verbose"} || 0; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #Path of stats.log file where stats data will be saved
	my $miARmaPath=$args{"miARmaPath"};
	my $adapter_file=$args{"adapter_file"};
	my $summary_file=$args{"summary"}; #Path to file whre print basic results
	#Declaring the variables to collect the path of the new files and the variables to control the function
	my $cutadapt_result;
	my $reaper_result;
	my $adaptrim_result;
	my $stats_result;
	my @adapter_results;
	my $cutadapt_exec=0;
	my $reaper_exec=0;
	my $adaptrim_exec=0;
	
	my $min=undef;
	if(exists $args{"min"}){
		$min=$args{"min"};
	}
	my $metafile=undef; 
	if(defined($args{"metafile"})){
		$metafile=$args{"metafile"};
	}
	
	#Checking the software provided by the user to perform the adapter removal
	my $adaptererror=0;
	if(lc($adaptersoft) =~ /cutadapt/){
		$cutadapt_exec=1;
		$adaptererror=1;
	}
	if(lc($adaptersoft) =~ /reaper/){
		$reaper_exec=1;
		$adaptererror=1;
	}
	if(lc($adaptersoft) =~ /adapttrimming/){
		$adaptrim_exec=1;
		$adaptererror=1;
	}
	if($adaptererror == 0){
		warn("ADAPTERREMOVAL ERROR :: ".date()." Invalid value for adaptersoft ($adaptersoft). Allowed values are: cutadapt, reaper, adapttrimming, cutadapt-reaper, cutadapt-adapttrimming, reaper-adapttrimming or cutadapt-reaper-adapttrimming");
		help_AdapterRemoval();
	}
	#reading all files
	my $dir_used;
	my $summary_cut;
	my $summary_rea;
	#creating results folder
	system("mkdir -p $projectdir");
	foreach my $file (@$files){
		#Checking if the user has provided the adapter sequence
		my $adapter=undef;
		if(exists $args{"adapter"}){
			$adapter=$args{"adapter"}; #Adapter sequence which is going to be removed from the reads
		}
		
		#Checking mandatory parameters
		if($adaptersoft and $file and $dir and $logfile and $projectdir and $statsfile){
			#Checking hidden files
			if($file !~ /^\./){
				#If the user has not defined the adapter, it will be predicted by Minion function
				if(! defined $adapter and ! defined $adapter_file and ! defined $metafile and lc($adaptersoft) ne "adapttrimming"){

					#Obtaining mandatory parameters for Minion analysis
					my $organism=$args{"organism"}; #Arguments to know which genome to blat against
					#Checking optional parameters
					my $adaptpredictionnumber=undef;
					if(exists $args{"adaptpredictionnumber"}){
						$adaptpredictionnumber=$args{"adaptpredictionnumber"}; 
					}
					my $minionadaptersequence=undef;
					if(exists $args{"minionadaptersequence"}){
						$minionadaptersequence=$args{"minionadaptersequence"}; 
					}
				
					#Calling Minion function to predict the sequence of the adapter
					$adapter=Minion(
						dir=>$dir,
						file=>$file,
						logfile=>$logfile,
						statsfile=>$statsfile,
						adaptpredictionnumber=>$adaptpredictionnumber,
						minionadaptersequence=>$minionadaptersequence,
						org=>$organism,
						verbose=>$verbose,
						miARmaPath=>$miARmaPath,
					);
				}
			
				#CutAdapt analysis will be performed when user provides adaptersoft with cutadapt value
				if($cutadapt_exec == 1){
					#Optional parameters to perform the analysis are predefined as undef variables
					my $max=undef;
					my $min_quality=undef;
					my $cutadaptparameters=undef;
					my $readstart=undef;
					$dir_used="$dir/cutadapt_result";
					#If user provides any optional parameter the value will be collected by the corresponding variable
					if(exists $args{"max"}){
						$max=$args{"max"};
					}
					if(exists $args{"min_quality"}){
						$min_quality=$args{"min_quality"};
					}
					if(defined($args{"cutadaptparameters"})){
						$cutadaptparameters=$args{"cutadaptparameters"};
					}
					if(defined($args{"readstart"})){
						$readstart=$args{"readstart"};
					}

					#Calling CutAdapt function
					$cutadapt_result=CutAdapt(
		  				dir=>$dir,
		  				projectdir=>$projectdir,
		  				file=>$file,
		  				verbose=>$verbose,
		  				logfile=>$logfile,
		  				statsfile=>$statsfile,
		  				adapter=>$adapter,
						adapter_file=>$adapter_file,
		  				min=>$min, 
		  				max=>$max,
		  				min_quality=>$min_quality, 
		  				cutadaptparameters=>$cutadaptparameters,
						miARmaPath=>$miARmaPath,
		  			);
		  			push (@adapter_results, $cutadapt_result);

		  			#Calling ReadLengthCount to count the number of reads of each size
		  			$stats_result=ReadLengthCount(
						file=>$cutadapt_result
		  			);
	  	
		  			#Obtaining the absolute path of the directory and the file to R execution
		  			my $abs_path = abs_path($projectdir);
		  			my $stats_path = abs_path($stats_result);

		  			#Calling AdapterGraph to plot the number of reads with different sizes after adapter removal
		  			AdapterGraph(
						file=>$stats_path,
						dir=>$abs_path."/cutadapt_results/",
						miARmaPath=>$miARmaPath
		  			);

		  			#Calling CutAdaptStats to count the number of the reads before and after the adapter removal
		  			CutAdaptStats(
		  				inputfile=>$dir."/".$file,
		  				outputfile=>$cutadapt_result,  
		  				verbose=>$verbose, 
		  				logfile=> $logfile,
		  				statsfile=>$statsfile,
						summary=>$summary_file,
		  			);
					
					#Reading cutadapt results to summarize processed results
					open(STAT,$statsfile) || die "$! : Cant' read $statsfile\n";
	
					my $file;
					my $processed_reads;
					my $processed_bases;
					my $trimmed_reads;
					my $trimmed_bases;
					my $quality_trimmed="0 (0% of processed reads)";
					my $too_short="0 (0% of processed reads)";
					my $too_long="0 (0% of processed reads)";
	
					while(<STAT>){
						chomp;
						if($_ =~ /^Command line/){
							my @split=split(/\s+/);
							$file=fileparse($split[$#split]);
							my $file;
							my $processed_reads;
							my $processed_bases;
							my $trimmed_reads;
							my $trimmed_bases;
							my $quality_trimmed="0 (0% of processed reads)";
							my $too_short="0 (0% of processed reads)";
							my $too_long="0 (0% of processed reads)";
						}
						if($_ =~ /Processed reads:/){
							$_ =~s/Processed reads:\s+(.*)//g;
							$processed_reads=$1;
						}
						if($_ =~ /Processed bases:/){
							$_ =~s/Processed bases:\s+(.*)//g;
							$processed_bases=$1;
						}
						if($_ =~ /Trimmed reads:/){
							$_ =~s/Trimmed reads:\s+(.*)//g;
							$trimmed_reads=$1;
						}
						if($_ =~ /Quality-trimmed:/){
							$_ =~s/Quality-trimmed:\s+(.*)//g;
							$quality_trimmed=$1;
						}
						if($_ =~ /Trimmed bases:/){
							$_ =~s/Trimmed bases:\s+(.*)//g;
							$trimmed_bases=$1;
						}
						if($_ =~ /Too short reads:/){
							$_ =~s/Too short reads:\s+(.*)//g;
							$too_short=$1;
						}
						if($_ =~ /Too long reads:/){
							$_ =~s/Too long reads:\s+(.*)//g;
							$too_long=$1;		
						}
						if($file){
							$summary_cut->{$file}="$processed_reads\t$processed_bases\t$trimmed_reads\t$trimmed_bases\t$quality_trimmed\t$too_short\t$too_long";
						}
					}
					close STAT;
				}
				if($reaper_exec == 1){
					#Optional parameters to perform the analysis are predefined as undef variables
					my $reaperparameters=undef;
					my $geom=undef;
					my $tabu=undef; 
					$dir_used="$dir/Reaper_results";
					
					#If user provides any optional parameter the value will be collected by the corresponding variable
					if(defined($args{"reaperparameters"})){
						$reaperparameters=$args{"reaperparameters"};
					}
					if(defined($args{"geom"})){
						$geom=$args{"geom"};
					}
					#If user provides a metafile with adapter sequences this option will be executed and adapter sequence in the adapter variable will not be used
					if(defined($args{"metafile"})){
						$metafile=$args{"metafile"};
						$adapter=undef;
					}
					if(defined($args{"tabu"})){
						$tabu=$args{"tabu"};
					}

					#Calling Reaper function
					$reaper_result=Reaper(
						dir=>$dir,
						file=>$file,
						logfile=>$logfile,
						adapter=>$adapter,
						geom=>$geom,
						metafile=>$metafile,
						tabu=>$tabu,
						min=>$min,
						reaperparameters=>$reaperparameters,
						verbose=>$verbose,
						projectdir=>$projectdir,
						miARmaPath=>$miARmaPath
					);
					push (@adapter_results, $reaper_result);

					#Calling ReadLengthCount to count the number of reads of each size
		  			$stats_result=ReadLengthCount(
						file=>$reaper_result,
		  			);

		  			#Obtaining the absolute path of the directory and the file to R execution
		  			my $abs_path = abs_path($projectdir);
		  			my $stats_path = abs_path($stats_result);


		  			#Calling AdapterGraph to plot the number of reads with different sizes after adapter removal
		  			AdapterGraph(
						file=>$stats_path,
						dir=>$abs_path."/Reaper_results/",
						miARmaPath=>$miARmaPath
		  			);
					print STDERR date()." Please check $dir for results.\nA summary can be consulted in $statsfile\n" if($verbose);
					
					my $stat_file=$reaper_result;
					$stat_file=~s/\.lane\.clean/\.sumstat/g;
					
					my $processed_reads=0;
					my $accepted=0;
					my $discarded=0;
					my $trimmed_bases=0;
					my $quality_trimmed=0;
					my $too_short=0;
					open(STAT,$stat_file) || die "$! : Cant $stat_file\n";
					while(<STAT>){
						chomp;
						if($_ =~ /^total_input=/){
							$_=~s/total_input=(.*)//g;
							$processed_reads=$1;
						}
						if($_ =~ /^total_accepted=/){
							$_=~s/total_accepted=(.*)//g;
							$accepted=$1;
						}
						if($_ =~ /^total_discarded=/){
							$_=~s/total_discarded=(.*)//g;
							$discarded=$1;
						}
						if($_ =~ /^discarded_length_cutoff=/){
							$_ =~s/discarded_length_cutoff=(.*)//g;
							$too_short=$1;
						}
						if($_ =~ /^bases_removed_by_quality=/){
							$_ =~s/bases_removed_by_quality=(.*)//g;
							$quality_trimmed=$1;		
						}
						if($_ =~ /^discarded_low_complexity=/){
							$_ =~s/discarded_low_complexity=(.*)//g;
							$trimmed_bases=$1;		
						}
						if($file){
							$summary_rea->{$file}="$processed_reads\t$accepted\t$discarded\t$too_short\t$quality_trimmed\t$trimmed_bases";
						}
					}
					close STAT;
					
				}
				if($adaptrim_exec == 1){
					#Optional parameters to perform the analysis are predefined as undef variables
					my $trimmingnumber=undef;
					my $readposition=undef;
					$dir_used="$dir/AdaptTriming_results";
					
					#If user provides any optional parameter the value will be collected by the corresponding variable
					if(defined($args{"trimmingnumber"})){
						$trimmingnumber=$args{"trimmingnumber"};
					}
					if(defined($args{"readposition"})){
						$readposition=$args{"readposition"};
					}

					my $adapt_file=AdaptTriming(
						dir=>$dir,
						file=>$file,
						trimmingnumber=>$trimmingnumber,
						readposition=>$readposition,
						verbose=>$verbose,
						logfile=>$logfile,
						projectdir=>$projectdir
					);
					push (@adapter_results, $adapt_file);
					print STDERR date()." Please check $dir for results.\nA summary can be consulted in $statsfile\n" if($verbose);
					
				}
		  	}else{
				if($file !~ /^\./){
			  		die ("ADAPTERREMOVAL ERROR :: ".date()."File($file) has an invalid format. AdapterRemoval only accepts .fastq, .fq, fastq.gz or .fq.gz files ");
				}
		  	}
		}else{
			#Registering error
	   		open(LOG,">> ".$logfile) || die "ADAPTERREMOVAL ERROR :: ".date()."Can't open '$logfile': $!";
	    	print LOG "ADAPTERREMOVAL ERROR :: ".date()." Directory($dir), projectdir($projectdir), file($file), statsfile($statsfile) and/or logfile($logfile) have not been provided";
	    	close(LOG);

			#If mandatory parameters have not been provided program will die and show error message
			warn("ADAPTERREMOVAL ERROR :: ".date()." Directory($dir), projectdir($projectdir), file($file), statsfile($statsfile) and/or logfile($logfile) have not been provided");
			help_AdapterRemoval();
		}
	}

	#Printing summary results
	if(scalar(keys %$summary_cut)>0 and $cutadapt_exec==1){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
		print SUMM "\nAdapter Removal [".$projectdir."/cutadapt_result]\n";
		print SUMM "Filename\tProcessed Reads\tProcessed Bases\tTrimmed reads\tTrimmed bases\tQuality-Discarded\tToo short reads\tToo long reads\n";
		foreach my $processed_file (sort keys %$summary_cut){
			print SUMM $processed_file ."\t". $summary_cut->{$processed_file}."\n";
		}
		close SUMM;
	}
	
	if(scalar(keys %$summary_rea)>0 and $reaper_exec==1){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
		print SUMM "\nAdapter Removal [".$projectdir."/Reaper_results]\n";
		print SUMM "Filename\tProcessed Reads\tAccepted Reads\tDiscarded Reads\tSize-Discarded\tQuality-Discarded\tLow Complexity-Discarded\n";
		foreach my $processed_file (sort keys %$summary_rea){
			print SUMM $processed_file ."\t". $summary_rea->{$processed_file}."\n";
		}
		close SUMM;
	}
	#Returning the names of the new files
	
	return(@adapter_results);
	
	sub help_AdapterRemoval{
	    my $usage = qq{
		  	$0 

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
  	 		               
			Examples:
			AdapterRemoval(dir=>"../reads", adaptersoft=>"Cutadapt-Reaper", projectdir=>".", file=>"file.fastq", verbose=>"verbose", logfile=>"run.log", statsfile=>"stats.log", organism=>"human", min=>"18", max=>"26", min_quality=>"1", cutadaptparameters=> "-O 4", reaperparameters=>" -3p-prefix 12/2/0/0 -dust-suffix-late 20 --nozip");
	};

	print STDERR $usage;
	exit(); 
	}
}

=head2 AdapterGraph

  Example    : 
  	AdapterGraph(
		file=>"file.fastq.lane.report.clean.len",
		dir=>"."
  	);
 
  Description: AdapterGraph generates graphs with the reads distribution after adapter removal
  Input parameters: 
	Mandatory parameters:
  	 [dir] Directory where the results plot will be generated 
  	 [file] File with the stats of the adapter removal from reaper or ReadLengthCount function. This file contains two columns: the size of the read
  	 and the number of the reads with that size. 
  Returntype : PDF file with the plot in the provided directory
  Requeriments: AdapterGraph function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- R v3.1.0 or higher software correctly installed
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut


sub AdapterGraph{
	#Arguments provided by user are collected by %args. File, dir, adapter, projectdir and logfile
	# are mandatory arguments while verbose is optional. 
	my %args=@_;
	my $file=$args{"file"}; #Path of the file with the stats of the adapter removal from reaper
	my $dir= $args{"dir"}; #Directory where the results plot will be generated
	my $miARmaPath=$args{"miARmaPath"};

	if($dir and $file){
		if($file =~ /.lane.report.clean.len$/){

			#Obatining the name of the file
			my $name=fileparse($file, qr{\.f.*});

			#Calling R from perl
			my $R;
			
			#If user has defined the directory where R is installed the bridge will be created since that directory
			if(defined($args{"Rdir"})){
				my $Rdir=$args{"Rdir"}; #path where R software is installed
				$R = Statistics::R->new($Rdir);
			}
			else
			{
				$R = Statistics::R->new() ;
			}

			#Starting R 
			$R->startR;

			#Declaring R instructions for the quality control analysis. QC_EdgeR R function is needed 
			my $cmds = <<EOF;
			source("$miARmaPath/lib/CbBio/RNASeq/R-Scripts/AdapterGraph.R") 
			AdapterGraph("$dir", "$file", "$name")
EOF

			#R commands execution
			my $out2 = $R->run($cmds);

		}else{
			die("ADAPTERGRAPH ERROR :: Provided file($file) is not the correct stats file. The correct file has \".lane.report.clean.len\" as extension");
		}
	}else{
		warn("ADAPTERGRAPH :: Directory($dir) and/or file($file) have not been provided");
		help_AdapterGraph();
	}
	sub help_AdapterGraph{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
  	 		[dir] Directory where the results plot will be generated 
  	 		[file] File with the stats of the adapter removal from reaper
  	 		               
			Examples:
			AdapterGraph(file=>"file.fastq.lane.report.clean.len", dir=".");
			

		};
		print STDERR $usage;
		exit();  
	}
}

=head2 ReadLengthCount

  Example    : 
  	ReadLengthCount(
		file=>"file.fastq"
  	);
 
  Description: ReadLengthCount counts the number of reads with the different sizes and prints the stats in a new file 
  with the same name in the same directory but with the ".lane.report.clean.len" extension
  Input parameters: 
	Mandatory parameters:
  	 [file] Fastq file with .fastq or .fq extension
  Returntype : File with the statistical data of the number of reads and the path of the new file. 
  Requeriments: ReadLengthCount function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub ReadLengthCount{
	#Arguments provided by user are collected by %args. File, dir, adapter, projectdir and logfile
	# are mandatory arguments while verbose is optional. 
	my %args=@_;
	my $file=$args{"file"}; #Path of the fastq file after adapter removal

	#Checking the extension of the file
	if($file =~ /.*\.fastq$/ or $file=~ /.*\.fq$/ or $file=~ /.*\lane.clean$/){
		open(FASTQ, $file) || die "READLENGTHCOUNT ERROR :: ".date()."Can't open '$file': $!";
		#Defining the variables needed to control the process
		my $linebeftseq=0;
		my $chars;
		my $read_hash;
		my $firstline=0;
		my $readstart;
		#Reading each line of the file
		while(<FASTQ>){
			#Deleting the last character (\n)
			chomp;
			if($_ =~ /^(@.{3,4})/ and $firstline == 0){
				$readstart= $1;
				$firstline=1; 
			}
			#Checking the first line of the reads 
			if($_ =~ /^$readstart/){
				$linebeftseq=1;
			}
			#Else the count adds one number
			else{
				$linebeftseq++;
			}
			#The line 2 contain the sequence of the read 
			if($linebeftseq == 2){
				#Keeping the length of the read
				$chars = length($_);
				#Keeping the length ina a hash
				$read_hash->{$chars}++;
			}
		}
		close FASTQ;
		#Opening the new stats file
		open(NEW,"> ".$file.".lane.report.clean.len") || die "READLENGTHCOUNT ERROR :: ".date()."Can't open '$file': $!";
		print NEW "length\tcount\n";
		#Printing the stats
 		foreach my $read (sort {$a<=>$b} keys %$read_hash){
    		print NEW "$read\t$read_hash->{$read}\n";
		}
		close NEW;
		return($file.".lane.report.clean.len");
	}else{
		warn("READLENGTHCOUNT ERROR :: Provided file($file) has an invalid format. Allowed extensions: .fastq or .fq");
		help_ReadLengthCount();
	}	
	sub help_ReadLengthCount{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
  	 		[file] Fastq file with .fastq or .fq extension
  	 		[readstart] String of characters which defines the begining of the reads
  	 		               
			Examples:
			ReadLengthCount(file=>"file.fastq");
			

		};
		print STDERR $usage;
		exit();  
	}
}



=head2 CutAdapt

  Example    : 
  CutAdapt(
  	dir=>"./reads",
  	projectdir=>".",
  	file=>"file.fastq",
  	verbose=>"verbose",
  	logfile=>"run.log",
  	adapter=>"ATCTCGTATGCCGTCTTCTGCTTGAA", 
  	min=>"18", 
  	max=>"26",
  	min_quality=>"1", 
  	cutadaptparameters=> "-O 4"
  	statsfile=>"stats.log"
  	);
  Description: CutAdapt function performs a cutadapt analysis consisting on the adapter 
  removal of the reads contained in the fastq file provided. This analysis generates a new fastq 
  file with the conserved reads which will be saved on a new cutadapt_results directory on the 
  project directory.
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where the input files are 
  	 [projectdir] Directory where results directory will be created
  	 [file] File which is going to be processing (fasta/fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [adapter] Adapter sequence to be removed in the analysis
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 Optional parameters:
  	 [min] Minimun length of the sequence read [15 as default to miRNA analysis]
	 [max] Maximun length of the sequence read [35 as default to miRNA analysis]
	 [min_quality] Minimun quality of the sequence read [0 as default to miRNA analysis]
	 [cutadaptparameters] Other cutadapt parameters to perform the analysis using the cutadapt recommended syntaxis
  	 [verbose] Option to show the execution data on the screen   
  Returntype : Fastq file at directory cutadapt_results. cutadapt returns also the path of 
  the results file.
  Requeriments: cutadapt function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Input files on fasta/fastq format (compressed files are accepted)
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub CutAdapt{

	#Arguments provided by user are collected by %args. File, dir, adapter, projectdir and logfile
	# are mandatory arguments while verbose is optional. 
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/cutadapt/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/cutadapt/";
	}

	#First, check that cutadapt is in path:
	my @cutadapt_bin=`which cut_adapt`;
	#Executing the command
	if(scalar(@cutadapt_bin)<1){
		die "CUTADAPT ERROR :: system args failed: $? : Is cutadapt installed and exported to \$PATH ?";
	}

	my $file=$args{"file"}; #File which is going to be processing
	my $adapter=$args{"adapter"}; #Adapter sequence which is going to be removed from the reads
	my $adapter_file=$args{"adapter_file"}; #Adapter sequence which is going to be removed from the reads
	my $dir= $args{"dir"} ."/"; #Input directory where the files are
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}."/"; #Input directory where results directory will be created
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #Path of the statsfile to write the stats data

	#Defining the results directory
	my $output_dir="/cutadapt_results/";
	
	#The minimun and maximum length has been established by default as 15 and 35 for miRNAseq
	#analysis. However these values will be changed if the user provide any other.
	my $min=15;
	my $max=35;
	my $min_quality=0;
	#If user has provided the minimun or maximum length will be collected by args and
	#the min and max variables will be overwritten
	if($args{"min"}){
		$min=$args{"min"};
	}
	if(defined $args{"max"}){
		$max=$args{"max"};
	}
	if(defined $args{"min_quality"}){
		$min_quality=$args{"min_quality"};
	}
	#Collecting cutadapt parameters provided by the user
	my $cutadaptpardef= " -m $min -M $max -q $min_quality ";
	#The number of missmatches can be provided by the user
	if(defined($args{"cutadaptparameters"})){
		my $cutadaptparameters=$args{"cutadaptparameters"};
		$cutadaptpardef.=" $cutadaptparameters ";
	}

	#Checking the mandatory parameters
	if ($file and $dir and $adapter and $projectdir and $statsfile and ! $adapter_file){
		#Checking hidden files
		if($file !~ /^\./){

			#Printing process data on screen
			print  STDOUT "CUTADAPT :: ".date()." Checking $file for CutAdapt analysis\n" if($verbose);
			#Cutadapt execution command 
			my $command="cut_adapt -b ".$adapter." ".$cutadaptpardef.$dir."/".$file;
			#Extracting the name of the file
			my $name=fileparse($file, qr{\.f.*});

			#commandef is the command that will be executed by system composed of the results
			#directory creation and the cutadapt execution. The output data 
			#will be saved on a .fastq file in the output directory. We use this extension to avoid 
			#problems reading gz files in the subsequent analysis. The error data will be 
			#redirected to the stats.log file
			my $commanddef="mkdir -p ".$projectdir.$output_dir." ; ".$command. " > ".$projectdir.$output_dir.$name."_cut.fastq 2>> ".$statsfile;
			
			#Execution data will be register on the run.log file. Opening the run.log file 
			open (LOG,">> ".$logfile) || die "CUTADAPT ERROR :: ".date()."Can't open '$logfile': $!";
			#Printing the date and command execution on the run.log file
			print LOG "CUTADAPT :: ".date()." Executing $commanddef\n";
			close(LOG);
			#If verbose option has been provided by user print the data on screen
			if($verbose){
				print STDOUT "CUTADAPT :: ".date()." Executing $commanddef\n";
			}
			#Executing the command or if system can't be executed die showing the error.
			system($commanddef) == 0
			or die "CUTADAPT ERROR :: system args failed: $? ($commanddef)";
			
			#The name of the output file is returned to main program
			return($projectdir.$output_dir.$name."_cut.fastq");
		}
		else{
			die ("CUTADAPT ERROR :: ".date()."File($file) has an invalid format. Cutadapt only accepts .fastq, .fq, fastq.gz or .fq.gz files ");
		}
	}
	elsif ($file and $dir and ! $adapter and $projectdir and $statsfile and $adapter_file){
		
		#Read the file where each adapter is specified per read file
		my $current_adapter=undef;
		open(ADAPT_FILE,$adapter_file) || die "Can't file adapt_file at $adapter_file ($!)\n";
		while(<ADAPT_FILE>){
			chomp;
			my ($current_file,$adapter)=split(/\t/);
			if($current_file eq $file){
				$current_adapter=$adapter;
			}
			if($current_file eq $file.".bz2"){
				$current_adapter=$adapter;
			}
			if($current_file eq $file.".gz"){
				$current_adapter=$adapter;
			}
		}
		close ADAPT_FILE;
		#checking that we have found the adapter for the specified file;
		
		if(defined $current_adapter){
			print STDOUT "CUTADAPT :: ".date()." Checking $file for CutAdapt analysis\n" if($verbose);
			my $command="cut_adapt -b ".$current_adapter." ".$cutadaptpardef.$dir."/".$file;
			#Extracting the name of the file
			my $name=fileparse($file, qr{\.f.*});

			#commandef is the command that will be executed by system composed of the results
			#directory creation and the cutadapt execution. The output data 
			#will be saved on a .fastq file in the output directory. We use this extension to avoid 
			#problems reading gz files in the subsequent analysis. The error data will be 
			#redirected to the stats.log file
			my $commanddef="mkdir -p ".$projectdir.$output_dir." ; ".$command. " > ".$projectdir.$output_dir.$name."_cut.fastq 2>> ".$statsfile;
			
			#Execution data will be register on the run.log file. Opening the run.log file 
			open (LOG,">> ".$logfile) || die "CUTADAPT ERROR :: ".date()."Can't open '$logfile': $!";
			#Printing the date and command execution on the run.log file
			print LOG "CUTADAPT :: ".date()." Executing $commanddef\n";
			close(LOG);
			#If verbose option has been provided by user print the data on screen
			if($verbose){
				print STDOUT "CUTADAPT :: ".date()." Executing $commanddef\n";
			}
			#Executing the command or if system can't be executed die showing the error.
			system($commanddef) == 0
			or die "CUTADAPT ERROR :: system args failed: $? ($commanddef)";
			
			#The name of the output file is returned to main program
			return($projectdir.$output_dir.$name."_cut.fastq");
		}
		else{
			warn("CUTADAPT ERROR :: ".date()." For the file $file, we have not found any correct adapter in $adapter_file. Please check read_file names\n");
			next;
		}
	}
	
	else 
	{
		#Registering error
   		open(LOG,">> ".$logfile) || die "CUTADAPT ERROR :: ".date()."Can't open '$logfile': $!";
    	print LOG "CUTADAPT ERROR :: ".date()." Directory($dir), projectdir($projectdir), file($file), adapter_file ($adapter_file), adapter ($adapter), statsfile($statsfile) and/or logfile($logfile) have not been provided";
    	close(LOG);


		#If mandatory parameters have not been provided program will die and show error message
		warn("CUTADAPT ERROR :: ".date()." Directory($dir), projectdir($projectdir), file($file), adapter_file ($adapter_file), adapter ($adapter), statsfile($statsfile) and/or logfile($logfile) have not been provided");
		help_CutAdapt();

	}

	sub help_CutAdapt{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[dir] Input directory where the input files are 
  	 		[projectdir] Directory where results directory will be created
  	 		[file] File in fastq format which is going to be processing (Extensions allowed: .fastq, .fq, .fq.gz, .fastq.gz)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[adapter] Adapter sequence to be removed in the analysis
  	 		[statsfile] Path of stats.log file where stats data will be saved
						 
			Optional parameters:
  	 		[min] Minimun length of the sequence read [15 as default to miRNA analysis]
	 		[max] Maximun length of the sequence read [35 as default to miRNA analysis]
	 		[min_quality] Minimun quality of the sequence read [0 as default to miRNA analysis]
	 		[cutadaptparameters] Other cutadapt parameters to perform the analysis using the cutadapt recommended syntaxis
  	 		[verbose] Option to show the execution data on the screen  
  	 		               
			Examples:
			CutAdapt(dir=>"./reads", projectdir=>".", file=>"file.fatsq", verbose=>"verbose", adapter=>"ATCTCGTATGCCGTCTTCTGCTTGAA", min=>"18", max=>"26", min_quality=>"1", cutadaptparameters=> "-O 4", logfile=> "run.log", statsfile=>"stats.log"); 

	};

	print STDERR $usage;
	exit(); 
	}


}

=head2 CutAdaptStats

  Example    : 
  CutAdaptStats(
  	inputfile=>"./reads/file.fastq",
  	outputfile=>"./cutadapt_results/file.fastq",  
  	verbose=>"verbose", 
  	logfile=> "run.log",
  	statsfile=>"stats.log",
  	);
  Description: This function collects the mandatory parameters inputfile, outputfile, statsfile and logfile 
  to count the number of reads before and after cutadapt analysis. It also 
  calculates the percent of conserved reads. These data will be printed on statsfile or on 
  the screen if verbose option was selected.  
  Input parameters: 
	Mandatory parameters:
  	 [inputfile] Path of the fastq file before adapter removal (Allowed extensions: .fastq, .fq, .fastq.gz or .fq.gz)
  	 [outputfile] Path of the fastq file after adapter removal (Allowed extensions: .fastq)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 Optional parameters:
  	 [verbose] Option to show the execution data on the screen   
  Returntype : Print stats data on stats.log file and on the screen.
  Requeriments: CutAdaptStats function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub CutAdaptStats{

	use File::Basename;
	#Arguments provided by user are collected by %args. Dir, readstart,statsfile and logfile
	#are mandatory arguments while verbose is optional.
	my %args=@_;
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #Path of the statsfile to write the stats data
	my $logfile=$args{"logfile"}; # Path of the logfile to write the execution data
	my $inputfile=$args{"inputfile"}; #Fastq file before adapter removal
	my $outputfile=$args{"outputfile"}; #Fastq file after adapter removal
	my $summary_file=$args{"summary"}; #Path to file whre print basic results
	#Variable declaration
	my $readstart;
	my $firstline=0;
	my $percent;
	my %reads;
	my $pre;
	my $post;
	
	#print STDERR "Recibo in $inputfile y out $outputfile\n";
	#Checking mandatory parameters
	if($inputfile and $outputfile and $logfile and $statsfile){
		
	}
    else{
    	#Registering error
   		open(LOG,">> ".$logfile) || die "CUTADAPTSTATS ERROR :: ".date()."Can't open '$logfile': $!";
    	print LOG "CUTADAPTSTATS ERROR :: ".date()."Inputfile($inputfile), outputfile($outputfile), logfile($logfile) and/or statsfile($statsfile) have not been provided";
    	close LOG;

    	#If mandatory parameters have not been provided program will die and show error message
    	warn ("CUTADAPTSTATS ERROR :: ".date()."Inputfile($inputfile), outputfile($outputfile), logfile($logfile) and/or statsfile($statsfile) have not been provided");
    	help_CutAdaptStats();
		
    }
    sub help_CutAdaptStats{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[inputfile] Path of the fastq file before adapter removal (Allowed extensions: .fastq, .fq, .fastq.gz or .fq.gz)
  	 		[outputfile] Path of the fastq file after adapter removal (Allowed extensions: .fastq)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[statsfile] Path of stats.log file where stats data will be saved
						 
			Optional parameters:
  	 		[verbose] Option to show the execution data on the screen  
  	 		               
			Examples:
			CutAdaptStats(inputfile=>"./reads/file.fastq", outputfile=>"./cutadapt_results/file.fastq", verbose=>"verbose", logfile=> "run.log", statsfile=>"stats.log");
	};

	print STDERR $usage;
	exit(); 
	}

}

=head2 AdaptTriming

  Example    : 
  	AdaptTriming(
		dir=>"./reads",
		file=>"file.fastq",
		trimmingnumber=>"12",
		position=>"3",
		verbose=>"verbose",
		logfile=>"run.log",
		projectdir=>"/."
				
	);
  Description: AdaptTriming is a function to remove a defined number of nucleotides provided 
  by trimmingnumber variable of the 5 or 3 end of each read of the provided fastq file. The 
  quality data of the removed nucleotides will be removed too. The output data will be saved 
  in a new fastqfile with the same name in the subdirectory AdaptTriming_results on the project directory.
  This function also return the name of the results file. 
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory with the fastqfiles to be processed
  	 [file] File which is going to be processing (fasta/fastq format)
  	 [trimmingnumber] Number of nucleotides to remove of the sequence read and the quality data
  	 [readposition] End of the read to remove the nucleotides (3 or 5).
  	 [logfile] Path of run.log file where execution data will be saved
  	 [projectdir] Directory where results directory will be created
  	Optional parameters:
  	 [verbose] Option to show the execution data on the screen   
  	 
  Returntype : Fastq file with the same name in the subdirectory AdaptTriming_results on the project 
  directory. This function also return the name of the results file. 
  Requeriments: AdaptTrimming function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Input fastq files.
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub AdaptTriming{

	#Arguments provided by user are collected by %args. Dir, file, readstart, trimmingnumber, 
	#readposition, projectdir and logfile are mandatory arguments while verbose is optional.
	my %args=@_;
	my $dir=$args{"dir"}; #directory to read the fastq files
	my $file=$args{"file"}; # File which is going to be processing
	my $trimmingnumber =$args{"trimmingnumber"}; #number of nucleotides to trimm from the reads.
	my $readposition=$args{"readposition"}; #position of the sequence to be trimmed.
	my $logfile=$args{"logfile"}; # path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	
	
	#Checking the mandatory arguments
	if($file and $dir and $trimmingnumber and $readposition and $projectdir){
		#Checking hidden files
		if($file !~ /^\./){
			#}
		##Checking the file extension
		#elsif($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq.gz$/ or $file =~ /.*\.fq$/ or $file =~ /.*\.fq.gz$/ or $file =~ /.*\.fq.bz2$/ or $file =~ /.*\.fastq.bz2$/){
			print STDERR date()." Trimming $trimmingnumber nt of the $readposition end of $file\n" if($verbose);
    		
			#Printing process information
			open(LOG,">> ".$logfile) || die "ADAPTTRIMING ERROR :: ".date()."Can't open '$logfile': $!";
    		print LOG "ADAPTTRIMING :: ".date()." Trimming $trimmingnumber nt of the $readposition end of $file\n";
    		#If verbose option has been provided program will print the data on screen too.
    		
			close LOG;
			#Defining the results directory
			my $output_dir= $projectdir."/AdaptTriming_results/";
			#Creating results directory
			system("mkdir -p ".$output_dir);
			
			my $fake_file;
			if($file =~ /\.gz$/){
				#Opening the fastq file
				open(FILE, "gunzip -c $dir/$file |") || die "ADAPTTRIMING ERROR :: ".date()."Can't open '$dir/$file': $!";
				$fake_file=$file;
				$fake_file=~s/\.gz$//g;
			}
			elsif($file =~ /\.bz2$/){
				#Opening the fastq file
				open(FILE, "bunzip2 -c $dir/$file |") || die "ADAPTTRIMING ERROR :: ".date()."Can't open '$dir/$file': $!";
				$fake_file=$file;
				$fake_file=~s/\.bz2$//g;
			}
			else{
				#Opening the fastq file
				open(FILE, $dir."/".$file) || die "ADAPTTRIMING ERROR :: ".date()."Can't open '$dir/$file': $!";
				$fake_file=$file;
			}
			#Extracting the name of the file
			#my $name=fileparse($fake_file, qr{\.f.*});
			my $name=fileparse($fake_file);
			$name=~s/\.\w+$//g;
			my $outputfile = $output_dir.$name."_at.fastq";
			
			#Opening the results file where the reads will be printed
			open(NEW,">> ".$outputfile) || die "ADAPTTRIMING ERROR :: ".date()."Can't open '$outputfile': $!";
			#Setting the counter
			my $linebeftseq=0;
			my $readstart;
			my $firstline=0;
			#Reading each line of the file
			while(<FILE>){
				#Obtaining the begining of the reads
				if($_ =~ /^(@.{3,4})/ and $firstline == 0){
					$readstart= $1;
					$firstline=1; 
				}
				#Checking the first line of the reads and printing it in the new file 
				if($_ =~ /^$readstart/){
					$linebeftseq=1;
					print NEW $_;
				}
				#Else the count adds one number
				else{
					$linebeftseq++;
				}
				#The lines 2 and 4 contains the sequence and the quality data of the bases
				#These lines will be trimmed 
				if($linebeftseq == 2 or $linebeftseq == 4){
					#For position 5 the sequence of the 5' end (start of the line) will not be printed in the newfile
					if($readposition == 5){
						$_ =~s/^.{$trimmingnumber}//;
					}
					#For position 3 the sequence of the 3' end (end of the line) will not be printed in the new file
					elsif($readposition == 3){
						$_ =~s/.{$trimmingnumber}$//;
					}
					else{
						die ("ADAPTTRIMING ERROR :: Invalid value for readposition ($readposition). Allowed values are 3 or 5");
					}
					print NEW $_;
					#After the last line of the read the counter will be reseted
					if($linebeftseq == 4){
						$linebeftseq=0;
					}
				}
				if($linebeftseq == 3 ){
					print NEW $_;
				}
			}
			close FILE;
			close NEW;	
			return($outputfile);
		}
		else 
		{
			die ("ADAPTTRIMING ERROR :: ".date()."File($file) has an invalid format. AdaptTriming only accepts .fastq or .fq files ");
		}
	}
	else
	{
		#Registering error
   		open(LOG,">> ".$logfile) || die "ADAPTTRIMING ERROR :: ".date()."Can't open '$logfile': $!";
    	print LOG "ADAPTTRIMING :: ".date()."Directory($dir), projectdir($projectdir), file($file), trimmingnumber($trimmingnumber), readposition($readposition) and/or logfile($logfile) have not been provided";
    	close LOG;

    	#If mandatory parameters have not been provided program will die and show error message
    	warn ("ADAPTTRIMING :: ".date()."Directory($dir), projectdir($projectdir), file($file), trimmingnumber($trimmingnumber), readposition($readposition) and/or logfile($logfile) have not been provided");
    	help_AdaptTriming();
		
	}
	sub help_AdaptTriming{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[dir] Input directory with the fastqfiles to be processed
  	 		[file] File which is going to be processing (fasta/fastq format)
  	 		[trimmingnumber] Number of nucleotides to remove of the sequence read and the qualitydata
  	 		[readposition] End of the read to remove the nucleotides (3 or 5).
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[projectdir] Directory where results directory will be created
						 
			Optional parameters:
  	 		[verbose] Option to show the execution data on the screen  
  	 		               
			Examples:
			AdaptTriming(dir=>"./reads", file=>"file.fastq", trimmingnumber=>"12", position=>"3", verbose=>"verbose", logfile=> "run.log", projectdir=>".")
	};

	print STDERR $usage;
	exit(); 
	}


}

=head2 ReadFilter

  Example    : 
  	ReadFilter(
			dir=>"./reads",
			minreadlength=>"18",
			maxreadlength=>"26",
			logfile=>"run.log",
			file=>"file.fastq",
			verbose=>"verbose",
			projectdir=>"."	
	);
  Description: ReadFilter is a function to keep the reads with lengths between the defined 
  minimun read length (minreadlength) and the maximun (maxreadlength). ReadFilter takes a
  fastq file and keep only the reads with the correct size. These reads will be printed in 
  a new fastq file with the same name in the subdirectory Readfilter_results on the 
  project directory. This function also return the name of the new file. 
  Input parameters: 
	Mandatory parameters:
  	 [file] File which is going to be processing (fastq format)
  	 [dir] Input directory with the fastq files to be processed
  	 [minreadlength] Number of the minimun read length to keep, the reads less than minreadlength will be discarded
  	 [maxreadlength] Number of the maximun read length to keep, the reads greater than maxreadlength will be discarded
  	 [logfile] Path of run.log file where execution data will be saved
  	 [projectdir] Directory where results directory will be created
  	Optional parameters:
  	 [verbose] Option to show the execution data on the screen   
  	  
  Returntype : Fastq file with the same name in the subdirectory Readfilter_results on the project 
  directory and the name of the new file fastq file.
  Requeriments: adapttrimming function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub ReadFilter{

	#Arguments provided by user are collected by %args. Dir, minreadlength, maxreadlength,
	#readstart, logfile and file are mandatory arguments while verbose is optional.
	my %args=@_;
	my $dir=$args{"dir"}; #directory to read the fastq files 
	my $minreadlength =$args{"minreadlength"}; #Number of the minimun read length to keep, the reads less than minreadlength will be discarded
	my $maxreadlength=$args{"maxreadlength"}; #Number of the maximun read length to keep, the reads greater than maxreadlength will be discarded
	my $logfile=$args{"logfile"}; # Path of run.log file to write the execution data
	my $file=$args{"file"}; #File which is going to be processing (fastq format)
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	
	#Defining the results directory
	my $outputfile = $projectdir."/Readfilter_results/".$file;
	
	#Checking the mandatory arguments
	if($dir and $projectdir and $logfile and $file and $minreadlength and $maxreadlength){
		#Checking the file extension
		#if($file =~ /.*\.fastq$/ or $file=~ /.*\.fq$/){
		if($file !~ /^\./){
			#Printing the process information on the log file
			open(LOG,">> ".$logfile) || die "READFILTER ERROR :: ".date()."Can't open '$logfile': $!";
			print LOG "READFILTER :: ".date()." Keeping the reads between $minreadlength and $maxreadlength of $file \n";
			close LOG;
			#If verbose option has been provided program will print the data on screen too.
	    	if($verbose){
				print STDOUT "READFILTER :: ".date()." Keeping the reads between $minreadlength and $maxreadlength of $file \n";
			}
			#Opening the file
			open(FILE, $dir."/".$file) || die "READFILTER ERROR :: ".date()."Can't open '$dir/$file': $!";
			#Making the output directory
			system("mkdir -p ".$projectdir."/Readfilter_results/");
			#Opening the new fastq file where the reads with the defined length will be printed
			open(NEW,">> ".$outputfile) || die "READFILTER ERROR :: ".date()."Can't open '$outputfile': $!";
			#Defining the variables needed to control the process
			my $linebeftseq=0;
			my $firstline;
			my $readfilter=0;
			my $readstart;
			my $readbegin=0;
			#Reading each line of the file
			while(<FILE>){
				#Deleting the last character (\n)
				chomp;
				#Obtaining the begining of the reads
				if($_ =~ /^(@.{3,4})/ and $readbegin == 0){
					$readstart= $1;
					$readbegin=1; 
				}
				#Checking the first line of the reads and keeping it in the firstline variable 
				if($_ =~ /^$readstart/){
					$linebeftseq=1;
					$firstline=$_;
				}
				#Else the count adds one number
				else{
					$linebeftseq++;
				}
				#The line 2 contain the sequence of the read 
				if($linebeftseq == 2){
					#Checking if the size of the read have the defined size.
					if(length($_)>=$minreadlength and length($_)<=$maxreadlength){
						#Printing the first and second line of the reads with a length between the defined limits. 
						print NEW $firstline."\n";
						print NEW $_."\n";
						#Setting a variable to select the next lines of the selected reads
						$readfilter=1;
					}
				}
				#Printing the third line of the selected reads
				if($linebeftseq == 3 and $readfilter == 1){
					print NEW $_."\n";
				}
				if($linebeftseq == 4 ){
					#Re-setting the counter of lines and the firstline varaible
					$linebeftseq=0;
					$firstline=undef;
					#Printing the forth line of the selected reads
					if($readfilter == 1){
						print NEW $_."\n";
						#Re-setting the selected read variable
						$readfilter=0;
					}
				}
			}
			#Closing the files
			close FILE;
			close NEW;
			
			#Returning the name of the file	
			return($file);
		}
		else
		{
			if($file !~ /^\./){
				die ("READFILTER ERROR :: ".date()."File($file) has an invalid format. ReadFilter only accepts .fastq or .fq files ");
			}
		}
	}
	else
	{
		#Registering error
   		open(LOG,">> ".$logfile) || die "READFILTER ERROR :: ".date()."Can't open '$logfile': $!";
    	print LOG "READFILTER ERROR :: ".date()."Directory($dir), projectdir($projectdir), file($file), minreadlength($minreadlength), maxreadlength($maxreadlength) and/or logfile($logfile) have not been provided";
    	close LOG;

    	#If mandatory parameters have not been provided program will die and show error message
    	warn ("READFILTER ERROR :: ".date()."Directory($dir), projectdir($projectdir), file($file), minreadlength($minreadlength), maxreadlength($maxreadlength) and/or logfile($logfile) have not been provided");
    	help_ReadFilter();
	}

	sub help_ReadFilter{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[file] File which is going to be processing (fastq format)
  	 		[dir] Input directory with the fastq files to be processed
  	 		[minreadlength] Number of the minimun read length to keep, the reads less than minreadlength will be discarded
  			[maxreadlength] Number of the maximun read length to keep, the reads greater than maxreadlength will be discarded
  	 		[logfile] Path of run.log file where execution data will be saved
  			[projectdir] Directory where results directory will be created
						 
			Optional parameters:
  	 		[verbose] Option to show the execution data on the screen  
  	 		               
			Examples:
			ReadFilter(dir=>"./reads", minreadlength=>"18", maxreadlength=>"26", logfile=>"run.log",file=>"file.fastq", verbose=>"verbose", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}
}

=head2 Minion

  Example    : 
  	Minion(
			dir=>"./reads",
			file=>"file.fastq",
			logfile=>"run.log",
			statsfile=>"stats.log",
			adaptpredictionnumber=>"4",
			minionadaptersequence=>"ACGCGCCATGTGATGAGCAGTGTCACGTAC",
			organism=>"human",
			verbose=>"verbose",
	);
  Description: Minion is a function to predict adapter sequence from fastq files. For this purpouse Minion 
  function executes Minion software and perform a Blat analysis. Minion is used to predict the 
  3' adapter of the provided fastq files. After execution Minion prints 2 alternative predictions 
  of the adapter sequences (by default) in the stats file on the provided directory. The number 
  of predictions can be defined with the optional adaptpredictionnumber variable. These sequences are evaluated 
  by Blat analysis to check if any of them is a bilogical sequence. If not Minion function returns the adapter sequence.
  Also Minion allows the comparison of a defined adapter sequence with all candidate sequences. Minion ouptputs 
  the candidate sequence that best matches the query sequence and shows the alignment between 
  the two sequences. 
  Input parameters: 
	Mandatory parameters:
  	 [file] Name of the fastq file to make the adapter prediction
  	 [dir] Input directory with the fastq files to be processed
  	 [organism] Arguments to know which genome to blat against (Allowed organism: human, mouse and rat)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	Optional parameters:
  	 [adaptpredictionnumber] Number of predictions to show
  	 [minionadaptersequence] Adapter known sequence to compare with Minion
  	 [verbose] Option to show the execution data on the screen   
  Returntype : Adapter prediction in the stats file on the provided directory
  Requeriments: adapttrimming function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Minion software correctly installed
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub Minion{

	#Arguments provided by user are collected by %args. Dir, file, organism, statsfile and logfile
	#are mandatory arguments while adaptpredictionnumber, minionadaptersequence and verbose
	#are optional.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/kraken/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/kraken/";
	}
	
	#First, check that minion is in path:
	my @minion_bin=`which minion`;
	#Executing the command
	if(scalar(@minion_bin)<1){
		die "MINION ERROR :: system args failed: $? : Is minion installed and exported to \$PATH ?";
	}
	my $dir=$args{"dir"}; #Directory to read the fastq files 
	my $file=$args{"file"}; #Name of the file which is going to be processing (fastq format)
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $statsfile=$args{"statsfile"}; #Path of the stats file to write the execution data
	my $adaptpredictionnumber =$args{"adaptpredictionnumber"}; #Number of predictions to show
	my $minionadaptersequence=$args{"minionadaptersequence"}; #Adapter known sequence to compare with Minion
	my $organism= lc ($args{"org"}); #Arguments to know which genome to blat against
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	
	#Defining the variables needed for the process 
	my $minionparameters="";
	
	#Matching the organism with the correspondant genome
	my $orgDB->{"human"}="hg19";
	$orgDB->{"mouse"}="mm10";
	$orgDB->{"rat"}="rn5";

	
	#Checking the mandatory arguments
	if($file and $logfile and $statsfile and $dir and $organism){

		##Printing the process information
		print STDOUT date()." Predicting an adapter sequence for $file \n" if($verbose);

		if(!exists($orgDB->{$organism})){
			print STDERR "MINION :: The provided organism ($organism) is not compatible with our software (". join(",",(keys %$orgDB)) . ").\n Please write an email to miARma-devel\@cbbio.es to include it\n\n";
			exit;
		}
		#Checking the optional parameters
		#To show a defined number of predictions
		if(defined($args{"adaptpredictionnumber"})){
			my $adaptpredictionnumber=$args{"adaptpredictionnumber"};
			$minionparameters.=" -show $adaptpredictionnumber";
		}
		#To compare a defined sequence to all candidate sequences. 
		if(defined($args{"minionadaptersequence"})){
			my $minionadaptersequence=$args{"minionadaptersequence"};
			$minionparameters.=" -adapter $minionadaptersequence";
		}
		#Printing the file info in the stats file
		open(STATS,">> ".$statsfile) || die "MINION ERROR :: ".date()."Can't open '$statsfile': $!";
		print STATS "MINION :: ".date()." Adapter prediction for $file sample\n";
		#Defining the execution command
		my $bzfiles;
		if($file =~ /\.bz2$/){
			$bzfiles="bunzip2 -c  $dir/$file | minion search-adapter $minionparameters";
		}
		else{
			$bzfiles="minion search-adapter -i ".$dir."/".$file.$minionparameters ;
		}

		my $command="$bzfiles  1>/tmp/minion.sq 2>> ".$statsfile;
		#Executing the command
		system($command) == 0
		or die "MINION ERROR :: system args failed: $? ($command)\nAre you sure this fastq file contains and adapter sequence?";
		#Printing the process information on the log file
		open(LOG,">> ".$logfile) || die "MINION ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "MINION :: ".date()." Executing $command \n";
		close LOG;
		close STATS;
		#If verbose option has been provided program will print the data on screen too.
    	if($verbose){
			print STDOUT "MINION :: ".date()." Executing $command \n";
		}

		#Obtaining Minion Predictions to perform Blat analysis
		my $prediction=0;
		my $results;
		my $criterion;
		my $seq;
		open(MINION,"/tmp/minion.sq") || die "MINION ERROR :: Can't find Minion predictions: $!\n";
		while(<MINION>){
			chomp;
			if(/criterion=/){
				$prediction++;
				$criterion=$_;
				$criterion=~s/criterion=(.+)/$1/g;
			}
			elsif(/sequence=/){
				$seq=$_;
				$seq=~s/sequence=(.+)/$1/g;
			}
			if($criterion and $seq){
				$results->{$prediction}->{$criterion}=$seq;
				$criterion="";
				$seq="";
			}
		}
		close MINION;

		my $number_possible_results=scalar(keys %$results);
		if($number_possible_results==0){
			print STDERR date(). " MINION ERROR :: Minion predictions for $file didn't comply our requeriments\n";
			next;
		}
		#Performing Blat analysis with the predicted adapters
		my $adapter=BlatResults(
			results=>$results,
			org=>$organism,
			orgDB=>$orgDB->{$organism},
			verbose=>$verbose
		);
		if($adapter){
			return($adapter);
		}
		else{
			print STDERR "MINION WARN :: None of the minion predictions (total : $number_possible_results) were valid. Asumming already trimmed (use adapter=no in ini file for trimmed input)\n";
			return();
		}
		
	}
	else{
		#Registering error
   		open(LOG,">> ".$logfile) || die "MINION ERROR :: ".date()." Can't open '$logfile': $!";
    	print LOG "MINION ERROR :: ".date()." Directory($dir), file($file), organism($organism), statsfile($statsfile) and/or logfile($logfile) have not been provided";
    	close LOG;

    	#If mandatory parameters have not been provided program will die and show error message
    	warn ("MINION ERROR :: ".date()." Directory($dir), file($file), organism($organism), statsfile($statsfile) and/or logfile($logfile) have not been provided");
    	help_Minion();
	}
	sub help_Minion{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[file] Name of the fastq file to make the adapter prediction
  	 		[dir] Input directory with the fastq files to be processed
  	 		[organism] Name of the organism to know which genome to blat against (Allowed organism: human, mouse and rat)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[statsfile] Path of stats.log file where stats data will be saved
						 
			Optional parameters:
  	 		[adaptpredictionnumber] Number of predictions to show
  	 		[minionadaptersequence] Adapter known sequence to compare with Minion
  	 		[verbose] Option to show the execution data on the screen    
  	 		               
			Examples:
			Minion(dir=>"./reads", file=>"file.fastq", logfile=>"run.log", statsfile=>"stats.log", adaptpredictionnumber=>"4", minionadaptersequence=>"ACGCGCCATGTGATGAGCAGTGTCACGTAC", organism=>"human", verbose=>"verbose");
	};

	print STDERR $usage;
	exit(); 
	}
}

=head2 BlatResults

  Example    : 
  	BlatResults(
			result=>$result,
			org=>"human",
			orgDB=>"hg19"		
	);
  Description: BlatResults is a function to execute to blat the adapter predicted by Minion against the genome of interest. 
  A good adapter must not align to the reference genome.
  Input parameters: 
	Mandatory parameters:
  	 [results] Adapter prediction from Minion. This is a hash with 3 dimensions: Number of predictions, Minion criterion and predicted adapter sequence
  	 [org] Name of the organism to know which genome to blat against (Allowed organism: human, mouse and rat)
  	 [orgDB] Genomic build of the selected organism (Allowed builds: hg19, mm10 and rn5)
  Returntype : Valid sequence
  Requeriments: BlatResults function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- LWP
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut
sub BlatResults{
	
	#Arguments provided by user are collected by %args. Results, org and dB are mandatory arguments.
	my %args=@_;
	my $results=$args{"results"}; #Minion predictions 
	my $org=$args{"org"}; #Organism to blat with
	my $dB=$args{"orgDB"}; #Genomic build to be used
	my $ok=0; #Checking results
	my $verbose=$args{"verbose"};
	#Establishing conditions to conect Blat server
	my $browser = LWP::UserAgent->new;

	#Checking the mandatory arguments
	if($results and $org and $dB){
		#Accesing to each prediction for each sample
		foreach my $pred(sort {$a <=> $b} keys %$results){
			foreach my $criteria (keys %{$results->{$pred}}){
				my $sequence=$results->{$pred}->{$criteria};
				if($criteria eq "sequence-density" and length($sequence)>=20){
					print STDERR "MINION :: Checking $sequence\n" if($verbose);
					#Printing process information
					print STDOUT "MINION :: " . date() . " Checking $pred prediction against ".$org.":".$dB." genome based on $criteria criteria\n" if($verbose);
					#Establishing the Blat web address 
					my $URL = "http://genome.ucsc.edu/cgi-bin/hgBlat";

					#Connecting to Blat web to match the Minion adapter prediction to the corresponding genomic build 
					my $response = $browser->post( $URL,
					    [ 'org' => $org, 
						  'db' => $dB,
					      'userSeq' => $sequence,
					      'type' => "BLAT's guess",
						  'sort' => "query,score",
					      'output' => "hyperlink"
					    ]
					);
					#Obtaining the blat results and printing in a temp file
					if ($response->is_success) {
		     			open(BLAT,">/tmp/blat_result") || die "MINION ERROR :: ".date()."Can't find Minion predictions: $!";
						print BLAT $response->decoded_content; 
						close BLAT;

						#Checking Blat results for genomic matches
						open(FILE,"/tmp/blat_result") || die "MINION ERROR :: ".date()."Can't find Minion predictions: $!";
						while(<FILE>){
							chomp;
							if(/YourSeq/){
								my $data=$_;
								$data =~s/.* (YourSeq.*)/$1/g;
								my (undef,$score,undef,undef,undef,$identity,$chr,$strand,$start,$end)=split(/\s+/,$data);
								print STDOUT "Your sequence is mapped to the following genomic Coordinates : [$chr:$start-$end;$strand], identity:$identity, score:$score\n" if($verbose);
								#print STDOUT "So, returning :" . $results->{$pred +1}->{"fanout-score"} ."\n";								
								#return($results->{$pred +1}->{"fanout-score"});
								next;
							}elsif(/Sorry/){
								$ok=1;
							}
						}
						close FILE;
						if($ok==1){
							return($sequence);
						}
		 			}
		 			else {
		     			print STDERR $response->status_line;
		     			exit;
					}
				}
			}
		}
		if($ok==0){
			return("");
		}
	}
	else{
		warn ("MINION ERROR :: Results, org($org) and/or orgdB($dB) variables have not been provided to the BlatResults subroutine");
		help_BlatResults();
		
	}
	sub help_BlatResults{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[results] Adapter prediction from Minion. This is a hash with 3 dimensions: Number of predictions, Minion criterion and predicted adapter sequence
  	 		[org] Arguments to know which genome to blat against (Allowed organism: human, mouse and rat)
  	 		[orgDB] Genomic build of the selected organism (Allowed builds: ?????)
  	 	            
			Examples:
			BlatResults(result=>\$result, org=>"human", orgDB=>"hg19");
	};

	print STDERR $usage;
	exit(); 
	}
}

=head2 Reaper

  Example    : 
  	Reaper(
			dir=>"./reads",
			file=>"file.fastq",
			logfile=>"run.log",
			adapter=>"ATCTCGTATGCCGTCTTCTGCTTGAA",
			reaperparameters=>" -3p-prefix 12/2/0/0 -dust-suffix-late 20 --nozip",
			min=>"18",
			verbose=>"verbose",
			projectdir=>"."		
	);
  Description: Reaper is a function to execute Reaper software. Reaper software is a program for 
  demultiplexing, trimming and filtering short read sequencing data. Reaper takes the sequence 
  of the adapter and barcode at 5' or 3' and removes them from the reads. In addtion you can filter 
  by length or discard low complexity regions. Output files will be saved in the Reaper_results on the project
  directory. See the instructions of Reaper software for more information (ftp://ftp.ebi.ac.uk/pub/contrib/enrightlab/kraken/reaper/src/reaper-latest/doc/reaper.html)
  Input parameters: 
	Mandatory parameters:
  	 [file] Fastq file to make the adapter prediction
  	 [dir] Input directory with the fastq files to be processed
  	 [logfile] Path of run.log file where execution data will be saved
  	 [projectdir] Directory where results directory will be created
  	 [metafile] Metadata file with the information of the adapter sequences. See reaper instructions to generate this file (Mandatory file for 3p-bc and 5p-bc geometries). 
  	 [adapter] Adapter sequence to remove from the reads with no barcode geometry (Mandatory if metafile have not been provided).
  	Optional parameters:
  	 [reaperparameters] Parameters to execute reaper. See reaper instructions to introduce the parameters with the correct syntaxis
  	 [geom] Geometry used in the analysis refering to the position of the barcode in the read (No barcode "no-bc" (default value), 3'end "3p-bc" or 5'end "5p-bc").
  	 [tabu] Tabu sequence to remove the read which contain this sequence (usually this sequence is 5' primer sequence) 
  	 [min] Minimun length of the sequence read [15 as default to miRNA analysis]
  Returntype : Fastq file in the Reaper_results directory. Reaper also return the path of the new file
  Requeriments: Reaper function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Reaper software correctly installed
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub Reaper{

	#Arguments provided by user are collected by %args. Dir, file, logfile and reaperparameters
	#are mandatory arguments while verbose is optional.
	my %args=@_;
	
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/kraken/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/kraken/";
	}
	
	#First, check that reaper is in path:
	my @reaper_bin=`which reaper`;
	#Executing the command
	if(scalar(@reaper_bin)<1){
		die "REAPER ERROR ::system args failed: $? : Is reaper installed and exported to \$PATH ?";
	}


	my $dir=$args{"dir"}; #directory to read the fastq files 
	my $file=$args{"file"}; #File which is going to be processing (fastq format)
	my $logfile=$args{"logfile"}; # path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	
	#Defining default parameters for the analysis
	my $geom="no-bc"; 
	my $metafile=undef; 
	my $tabu="-"; 
	my $adapter=undef;
	my $reaperpardef=" --nozip";
	my $min="15";

	#If any additional parameter has been provided will be collected in the corresponding variable
	if(defined($args{"geom"})){
		$geom=$args{"geom"}; #Geometry used in the analysis refering to the position of the barcode in the read (No barcode "no-bc", 3'end "3p-bc" or 5'end "5p-bc") 
	}
	if(defined($args{"metafile"})){
		$metafile=$args{"metafile"}; #metadata file with the information of the adapter sequences. See reaper instructions to generate this file. 
	}
	if(defined($args{"adapter"})){
		$adapter=$args{"adapter"}; #Adapter sequence to remove from the reads with no barcode geometry. 
	}
	if(defined($args{"tabu"})){
		$tabu=$args{"tabu"}; #Tabu sequence to remove the read which contain this sequence.
	}
	if(defined($args{"reaperparameters"})){
		my $reaperparameters=$args{"reaperparameters"}; #Parameters to execute reaper. See reaper instructions to introduce the parameters with the correct syntaxis
		$reaperpardef.=" $reaperparameters ";
	}
	if(defined($args{"min"})){
		my $min=$args{"min"}; #Minimun length of the sequence read [15 as default to miRNA analysis]
		$reaperpardef.=" -clean-length $min ";
	}
	
	#Defining the output directory
	my $resultsdir= $projectdir."/Reaper_results";
	my $command;
	
	#Checking the mandatory arguments
	if($file and $logfile and $geom and $dir and $projectdir){
		#Checking hidden files
		if($file !~ /^\./){
			#}
		#elsif($file =~ /.*\.fastq$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq\.gz$/){
			#Printing process data on screen
			print STDERR date()." Checking $file for Reaper analysis\n" if($verbose);
			#Extracting the name of the file
			my $name=fileparse($file, qr{\.f.*});
			
			my $bzfiles;
			if($file =~ /\.bz2$/){
				$bzfiles="bunzip2 -c $dir/$file | reaper ";
			}
			else{
				$bzfiles="reaper -i $dir/$file ";
			}
			#Defining the execution command
			if($metafile){
				$command="mkdir -p ".$resultsdir." ; $bzfiles -basename ".$resultsdir."/".$name."_rea.fastq -geom ".$geom." -meta ".$metafile. $reaperpardef." 2>> ".$logfile;
			}elsif($geom eq "no-bc" and $adapter){
				$command="mkdir -p ".$resultsdir." ; $bzfiles -basename ".$resultsdir."/".$name."_rea.fastq -geom ".$geom." -3pa ".$adapter." -tabu ".$tabu.$reaperpardef." 2>> ".$logfile;
			}else{
				die("REAPER ERROR :: ".date()."Reaper needs a metafile ($metafile) with the adapter information or the adapter($adapter) in the case of no-bc geometry. No one of these options has been provided");
			}

			#Executing the command
			system($command) == 0
			or die "REAPER ERROR :: system args failed: $? ($command)";
			
			#Printing the process information on the log file
			open(LOG,">> ".$logfile) || die $!;
			print LOG "REAPER :: ".date()." Executing $command \n";
			#If verbose option has been provided program will print the data on screen too.
	    	if($verbose){
				print STDERR "REAPER :: ".date()." Executing $command \n";
			}

			close LOG;
			return($resultsdir."/".$name."_rea.fastq.lane.clean")
		}else{
			if($file !~ /^\./){
				die ("REAPER ERROR :: ".date()."File($file) has an invalid format. Reaper only accepts .fastq, .fq, fastq.gz or .fq.gz files ");
			}
		}
	}
	else 
	{
		#Registering error
   		open(LOG,">> ".$logfile) || die $!;
    	print LOG "REAPER ERROR :: ".date()." Directory($dir), projectdir($projectdir), file($file), geometry($geom) and/or logfile($logfile) have not been provided";
    	close LOG;

    	#If mandatory parameters have not been provided program will die and show error message
    	warn ("REAPER ERROR :: ".date()." Directory($dir), projectdir($projectdir), file($file), geometry($geom) and/or logfile($logfile) have not been provided");
    	help_Reaper();
	}

	sub help_Reaper{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
  	 		[file] Fastq file to make the adapter prediction
  	 		[dir] Input directory with the fastq files to be processed
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[projectdir] Directory where results directory will be created
  	 		[metafile] Metadata file with the information of the adapter sequences. See reaper instructions to generate this file (Mandatory file for 3p-bc and 5p-bc geometries). 
  			[adapter] Adapter sequence to remove from the reads with no barcode geometry (Mandatory if metafile have not been provided).
  			Optional parameters:
  	 		[reaperparameters] Parameters to execute reaper. See reaper instructions to introduce the parameters with the correct syntaxis
  	 		[geom] Geometry used in the analysis refering to the position of the barcode in the read (No barcode "no-bc" (default value), 3'end "3p-bc" or 5'end "5p-bc").
  			[tabu] Tabu sequence to remove the read which contain this sequence (usually this sequence is 5' primer sequence) 
  	 	            
			Examples:
			Reaper(dir=>"./reads", file=>"file.fastq", logfile=>"run.log", adapter=>"ATCTCGTATGCCGTCTTCTGCTTGAA", reaperparameters=>" -3p-prefix 12/2/0/0 -dust-suffix-late 20 -clean-length 18 --nozip", verbose=>"verbose", projectdir=>".")
	};

	print STDERR $usage;
	exit(); 
	}
}
sub date{
	#my $dt = DateTime->now(time_zone=>'local');
	#return($dt->hms . " [" . $dt->dmy ."]");
	use Time::localtime;
	my $now = ctime();
	return("[$now]");
}
1;