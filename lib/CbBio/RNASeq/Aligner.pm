#########################################################################	
#	Bowtie processing package		 									#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2015 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################

package CbBio::RNASeq::Aligner;
require Exporter;
#Export package system
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(bowtie1_index bowtie1 bowtie2_index bowtie2 IndexGeneration ReadAligment bwa tophat ReadSummary hisat2 hisat2_index star star_index);

use strict;
use DateTime;
use File::Basename;
use Cwd 'abs_path';
=head1 NAME

 Aligner

=head1 SYNOPSIS

Aligner is composed of 6 subroutines:bowtie1_index, bowtie1, bowtie2_index and bowtie2, IndexGeneration and
ReadAlignment. The numbers 1/2 refers to the version of bowtie software 1.0.0. or 2.2.0. The index subroutines 
build the respective index from the genome sequence in fasta format. The bowtie subroutines 
perform an alignment of the reads to a genome index with the respective bowtie programme. IndexGeneration is 
a fuction to create the corresponding index (bowtie1/2 or both) for the alignment of the reads and ReadAligment performs 
the alignment with bowtie1/2 or both.

=head1 Methods

=head2 IndexGeneration

  Example    : 
  index_generation(
  	aligner=>"Bowtie1",
  	fasta=>"genome.fasta",
  	dir=>".",
  	logfile=>"run.log",
  	indexname=>"hg19"
  );
  Description: IndexGeneration is a fuction to create the corresponding index for the alignment 
  of the reads. This function allow the execution of bowtie1_index and/or bowtie2_index. The corresponding 
  index genomes will be saved on the provided by the user.
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where new index will be saved
  	 [fasta] Path of the genomic fasta sequence to build the index
  	 [logfile] Path of run.log file where execution data will be saved
  	 [indexname] Name to write in the index files
  	 [aligner] Aligner which will be use in the analysis to generate the corresponding index (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1)
  Returntype : Directory with index genome and the path of the new index 
  Requeriments: bowtie1_index function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v1.0.0 or higher software correctly installed
  	- Input genomic sequence on fasta format  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut


sub IndexGeneration{

	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname, aligner and logfile
	# are mandatory arguments.
	my %args=@_;
	my $fasta=$args{"fasta"}; #the path of genome sequence in fasta format
	my $dir=$args{"dir"}; #directory to create the genomeindex directory where new index will be saved
	my $logfile=$args{"logfile"}; #path of the logfile to write the execution data
	my $indexname=$args{"indexname"}; #name to write in the index files
	my $aligner=$args{"aligner"}; #Aligner which will be use in the analysis to generate the corresponding index
	my $miARmaPath=$args{"miARmaPath"}; #path to miARma stuff
	my $threads=$args{"threads"}; # number of threads to use (for new aligners such us hisat2 and star)
	my $gtf=$args{"gtf"}; # GTF file to create a trasncriptme in star
	
	
	#Declaring the variables to collect the path of the new index
	my $bowtie1index;
	my $bowtie2index;
	my $bwaindex;
	my @index;

	if ($fasta and $dir and $indexname and $logfile and $aligner){
		#Checking the aligner which will be used to generate the correspondig index
		if (lc($aligner) eq "bowtie1"){
			#Calling bowtie1_index function
			$bowtie1index= bowtie1_index(
	  			fasta=>$fasta,
	  			dir=>$dir,
	  			logfile=>$logfile,
	  			indexname=>$indexname,
				miARmaPath=>$miARmaPath
				
	  		);
	  		#Saving the path of the new bowtie1 index
	  		push(@index, $bowtie1index);
		}
		elsif (lc($aligner) eq "bowtie2"){
			#Calling bowtie2_index function
			$bowtie2index=bowtie2_index(
	  			fasta=>$fasta,
	  			dir=>$dir,
	  			logfile=>$logfile,
	  			indexname=>$indexname,
				miARmaPath=>$miARmaPath,
				threads=>$threads
				
	  		);
	  		#Saving the path of the new bowtie2 index
	  		push(@index, $bowtie2index);
		}
		elsif (lc($aligner) eq "bowtie1-bowtie2" or lc($aligner) eq "bowtie2-bowtie1"){
			#Calling bowtie1_index function
			$bowtie1index=bowtie1_index(
	  			fasta=>$fasta,
	  			dir=>$dir,
	  			logfile=>$logfile,
	  			indexname=>$indexname,
				miARmaPath=>$miARmaPath
				
	  		);
	  		#Saving the path of the Bowtie1 index in an array
	  		push(@index, $bowtie1index);

			#Calling bowtie2_index function
			$bowtie2index=bowtie2_index(
	  			fasta=>$fasta,
	  			dir=>$dir,
	  			logfile=>$logfile,
	  			indexname=>$indexname,
				miARmaPath=>$miARmaPath,
				threads=>$threads
				
	  		);
	  		#Saving the path of the Bowtie2 index in an array
	  		push(@index, $bowtie2index);	
		}
		elsif(lc($aligner) eq "bwa"){
			#Calling bowtie2_index function
			$bwaindex=bwa_index(
	  			fasta=>$fasta,
	  			dir=>$dir,
	  			logfile=>$logfile,
	  			indexname=>$indexname,
				miARmaPath=>$miARmaPath
				
	  		);
	  		#Saving the path of the new bowtie2 index
	  		push(@index, $bwaindex);
		}
		elsif(lc($aligner) eq "hisat2"){
			#Calling bowtie2_index function
			$bwaindex=hisat2_index(
	  			fasta=>$fasta,
	  			dir=>$dir,
	  			logfile=>$logfile,
	  			indexname=>$indexname,
				miARmaPath=>$miARmaPath,
				threads=>$threads
				
	  		);
	  		#Saving the path of the new bowtie2 index
	  		push(@index, $bwaindex);
		}
		elsif(lc($aligner) eq "star"){
			#Calling bowtie2_index function
			$bwaindex=star_index(
	  			fasta=>$fasta,
	  			dir=>$dir,
	  			logfile=>$logfile,
	  			indexname=>$indexname,
				miARmaPath=>$miARmaPath,
				threads=>$threads,
				gtf=>$gtf
				
	  		);
	  		#Saving the path of the new bowtie2 index
	  		push(@index, $bwaindex);
		}
		else{
			die("INDEX_GENERATION ERROR:: ".date()."Provided aligner argument ($aligner) has an invalid value. Allowed values for this parameters: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1");
		}
		#Returning the path of the new index in an array
	  	return(@index);
	}
	else{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "INDEX_GENERATION ERROR:: ".date()."Can't open '$logfile': $!";
    	print LOG "INDEX_GENERATION ERROR:: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile), aligner($aligner) and/or fasta file($fasta) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program will die and show error message
		warn ("INDEX_GENERATION ERROR:: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile), aligner($aligner) and/or fasta file($fasta) have not been provided");
		help_IndexGeneration();
	}
	sub help_IndexGeneration{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[dir] Input directory where new index will be saved
  	 		[fasta] Path of the genomic fasta sequence to build the index
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[indexname] Name to write in the index files
  	 		[aligner] Aligner which will be use in the analysis to generate the corresponding index (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1)
						             
			Examples:
			IndexGeneration(fasta=>"genome.fasta", dir=>".", logfile=>"run.log", indexname=>"hg19", aligner=>"Bowtie1");

	};

	print STDERR $usage;
	exit(); 
	}

}

=head2 ReadAligment

  Example    : 
  ReadAligment(
		file=>"file.fastq",
		aligner=>"Bowtie1-Bowtie2",
		threads=>"4",
		bowtie1index=>"./genomeindex1/hg19",
		bowtie1index=>"./genomeindex2/hg19",
		verbose=>"verbose", 
		logfile=>"run.log", 
		bowtiemiss=>"0", 
		bowtieleng=>"19",
		statsfile=>"stats.log", 
		bowtie1parameters=>" --sam --best --nomaqround -e 70 -k 1",
		projectdir=>"."
  	);
  Description: ReadAlignment is a fuction to align the reads with the corresponding index. This function
  executes bowtie1 and/or bowtie2 functions. For this reason ReadAligment collects the arguments needed to execute 
  these functions as well as the parameter aligner which defines the aligner which will be employed in the analysis.  
  Execution and stats data will be saved on run.log and stats.log file and will show on the screen 
  if verbose option is selected. 
  Input parameters: 
	Mandatory parameters:
  	 [file] Path of the file which is going to be aligned (fasta/fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 [projectdir] Directory where results directory will be created
  	 [aligner] Aligner which will be use in the analysis (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1)
  	 Optional parameters:
  	 [bowtie1index] Indexed genome to align your reads in format .ebwt (Mandatory for analysis with bowtie1)
  	 [bowtie2index] Indexed genome to align your reads in format .bt2 (Mandatory for analysis with bowtie2)
  	 [threads] Optional number of threads to perform the analysis  
  	 [bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 [bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 [bowtie1parameters] Other bowtie parameters to perform the analysis using the bowtie1 recommended syntaxis
  	 [bowtie2parameters] Other bowtie parameters to perform the analysis using the bowtie2 recommended syntaxis
  	 [verbose] Option to show the execution data on the screen   
  Returntype : File at directory bowtie1_results and/or bowtie2_results according to the aligner selected for the analysis.
  Also returns the path of the output file 
  Requeriments: bowtie1 function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v1.0.0 or higher software correctly installed
  	- Input files on fastq format on the provided directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut


sub ReadAligment{

	#Arguments provided by user are collected by %args. Dir, file, aligner, statsfile, projectdir 
	#and logfile are mandatory arguments while verbose and threads are optional.
	my %args=@_;
	my $file=$args{"file"}; #Name of the file which is going to be aligned
	my $aligner=$args{"aligner"}; #Aligner which will be use in the analysis 
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #Path of the statsfile to write the stats data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $miARmaPath=$args{"miARmaPath"};
	my $organism=$args{"organism"}; #Organism to align
	my $adapter=$args{"adapter"}; #are this reads processed by miArma 
	my $adapter_file=$args{"adapter_file"}; #are this reads processed by miArma 
	#Declaring the variables to collect the path of the new files
	my $output_file1;
	my $output_file2;
	my @results;
	if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.lane\.clean$/ or $file =~ /.*\.fastq.gz$/ or $file =~ /.*\.fq$/ or $file =~ /.*\.fq.gz$/ or $file =~ /.*\.fq\.bz2$/ or $file =~ /.*\.fastq\.bz2$/){
		if($file and ($aligner or $adapter_file) and $logfile and $statsfile and $projectdir){
			#Checking the selected aligner for the analysis
			if(lc($aligner) eq "bowtie1"){
				#Collecting specific bowtie1 index
				my $bowtie1index=$args{"bowtie1index"}; #Genome index to align your reads in .ebwt format
				#Optional parameters are predefined as undef
				my $bowtiemiss=undef;
				my $bowtieleng=undef;
				my $threads=undef;
				my $bowtie1parameters=undef;
				my $Seqtype;
				#If any optional parameter is provided by the user will be collected
				if(defined($args{"bowtiemiss"})){
					$bowtiemiss=$args{"bowtiemiss"};
				}
				if(defined($args{"bowtieleng"})){
					$bowtieleng=$args{"bowtieleng"};
				}
				if(defined($args{"threads"})){
					$threads=$args{"threads"}; 
				}
				if(defined($args{"bowtie1parameters"})){
					$bowtie1parameters=$args{"bowtie1parameters"};
				}
				if(defined($args{"Seqtype"})){
					$Seqtype=$args{"Seqtype"};
				}else{
					$Seqtype="SingleEnd";
				}
				
				if($bowtie1index){
					#Calling Bowtie1 function
					$output_file1=bowtie1( 
						file=>$file,
						threads=>$threads,
						bowtieindex=>$bowtie1index,
						verbose=>$verbose, 
						logfile=>$logfile, 
						bowtiemiss=>$bowtiemiss, 
						bowtieleng=>$bowtieleng,
						statsfile=>$statsfile, 
						bowtieparameters=>$bowtie1parameters,
						projectdir=>$projectdir,
						miARmaPath=>$miARmaPath,
						adapter=>$adapter,	
						Seqtype=>$Seqtype
			  		);
			  		return($output_file1);
			  	}
			  	else{
			  		die("READALIGNMENT ERROR:: ".date()."Index argument ($bowtie1index) has not been provided");
			  	}
			}
			elsif(lc($aligner) eq "bowtie2"){
				#Collecting specific bowtie2 index
				my $bowtie2index=$args{"bowtie2index"}; #Genome index to align your reads in .bt2 format
				#Optional parameters are predefined as undef
				my $bowtiemiss=undef;
				my $bowtieleng=undef;
				my $threads=undef;
				my $bowtie2parameters=undef;
				my $Seqtype;
				#If any optional parameter is provided by the user will be collected
				if(defined($args{"bowtiemiss"})){
					$bowtiemiss=$args{"bowtiemiss"};
				}
				if(defined($args{"bowtieleng"})){
					$bowtieleng=$args{"bowtieleng"};
				}
				if(defined($args{"threads"})){
					$threads=$args{"threads"}; 
				}
				if(defined($args{"bowtie2parameters"})){
					$bowtie2parameters=$args{"bowtie2parameters"};
				}
				if(defined($args{"Seqtype"})){
					$Seqtype=$args{"Seqtype"};
				}else{
					$Seqtype="SingleEnd";
				}
				#Calling Bowtie2 function
				$output_file2=bowtie2( 
					file=>$file,
					threads=>$threads,
					bowtieindex=>$bowtie2index,
					verbose=>$verbose, 
					logfile=>$logfile, 
					bowtiemiss=>$bowtiemiss, 
					bowtieleng=>$bowtieleng,
					statsfile=>$statsfile, 
					bowtieparameters=>$bowtie2parameters,
					projectdir=>$projectdir,
					miARmaPath=>$miARmaPath,
					adapter=>$adapter,	
					Seqtype=>$Seqtype
								
		  		);
		  		return($output_file2);

			}
			elsif(lc($aligner) eq "bowtie1-bowtie2" or lc($aligner) eq "bowtie2-bowtie1"){
				#Collecting specific bowtie index
				my $bowtie1index=$args{"bowtie1index"}; #Genome index to align your reads in .ebwt format
				my $bowtie2index=$args{"bowtie2index"}; #Genome index to align your reads in .bt2 format

				#Optional parameters are predefined as undef
				my $bowtiemiss=undef;
				my $bowtieleng=undef;
				my $threads=undef;
				my $bowtie1parameters=undef;
				my $bowtie2parameters=undef;
				my $Seqtype;
				#If any optional parameter is provided by the user will be collected
				if(defined($args{"bowtiemiss"})){
					$bowtiemiss=$args{"bowtiemiss"};
				}
				if(defined($args{"bowtieleng"})){
					$bowtieleng=$args{"bowtieleng"};
				}
				if(defined($args{"threads"})){
					$threads=$args{"threads"}; 
				}
				if(defined($args{"bowtie1parameters"})){
					$bowtie1parameters=$args{"bowtie1parameters"};
				}
				if(defined($args{"bowtie2parameters"})){
					$bowtie2parameters=$args{"bowtie2parameters"};
				}
				if(defined($args{"Seqtype"})){
					$Seqtype=$args{"Seqtype"};
				}else{
					$Seqtype="SingleEnd";
				}
				
				#Calling Bowtie1 function
				$output_file1=bowtie1(
					file=>$file,
					threads=>$threads,
					bowtieindex=>$bowtie1index,
					verbose=>$verbose, 
					logfile=>$logfile, 
					bowtiemiss=>$bowtiemiss, 
					bowtieleng=>$bowtieleng,
					statsfile=>$statsfile, 
					bowtieparameters=>$bowtie1parameters,
					projectdir=>$projectdir,
					miARmaPath=>$miARmaPath,
					adapter=>$adapter,
					Seqtype=>$Seqtype,
				
			  	);
			  	push(@results, $output_file1);

			  	#Calling Bowtie2 function
				$output_file2=bowtie2( 
					file=>$file,
					threads=>$threads,
					bowtieindex=>$bowtie2index,
					verbose=>$verbose, 
					logfile=>$logfile, 
					bowtiemiss=>$bowtiemiss, 
					bowtieleng=>$bowtieleng,
					statsfile=>$statsfile, 
					bowtieparameters=>$bowtie2parameters,
					projectdir=>$projectdir,
					miARmaPath=>$miARmaPath,
					adapter=>$adapter,
					Seqtype=>$Seqtype,
				
		  		);
		  		push(@results, $output_file2);
		
			}
			elsif(lc($aligner) eq "tophat"){
				my $library_type;
				my $tophat_seg_mismatches;
				my $tophat_seg_length;
				my $threads;
				my $tophatParameters;
				my $tophat_multihits;
				my $read_mismatches;
				my $GTF=$args{"GTF"};
				my $tophat_aligner=$args{"tophat_aligner"};
				my $Seqtype;
				if($GTF){
					if(!$tophat_aligner){
						$tophat_aligner="bowtie";
					}
					if(defined($args{"Seqtype"})){
						$Seqtype=$args{"Seqtype"};
					}else{
						$Seqtype="SingleEnd";
					}
					print STDOUT "\t". date()." Starting a $Seqtype analysis with $tophat_aligner\n" if($verbose);
					
					#If any optional parameter is provided by the user will be collected
					if(defined($args{"library_type"})){
						$library_type=$args{"library_type"};
					}
					if(defined($args{"tophat_seg_mismatches"})){
						$tophat_seg_mismatches=$args{"tophat_seg_mismatches"};
					}
					if(defined($args{"threads"})){
						$threads=$args{"threads"}; 
					}
					if(defined($args{"tophat_seg_length"})){
						$tophat_seg_length=$args{"tophat_seg_length"};
					}
					if(defined($args{"tophatParameters"})){
						$tophatParameters=$args{"tophatParameters"};
					}
					if(defined($args{"tophat_multihits"})){
						$tophat_multihits=$args{"tophat_multihits"};
					}
					if(defined($args{"read_mismatches"})){
						$read_mismatches=$args{"read_mismatches"};
					}
				
					my $bowtieindex;
					my $output_file2;
				
					
					if(lc($tophat_aligner) eq "bowtie2"){
						#Collecting specific bowtie2 index
						$bowtieindex=$args{"bowtie2index"}; #Genome index to align your reads in .bt2 format
						#Calling Bowtie2 function
						$output_file2=TopHat( 
							file=>$file,
							threads=>$threads,
							bowtieindex=>$bowtieindex,
							verbose=>$verbose, 
							logfile=>$logfile, 
							library_type=>$library_type, 
							tophat_seg_length=>$tophat_seg_length,
							tophat_seg_mismatches=>$tophat_seg_mismatches,
							tophat_multihits=>$tophat_multihits,
							read_mismatches=>$read_mismatches,
							statsfile=>$statsfile, 
							tophatParameters=>$tophatParameters,
							projectdir=>$projectdir,
							miARmaPath=>$miARmaPath,
							tophat_aligner=>$tophat_aligner,
							GTF=>$GTF,
							Seqtype=>$Seqtype,
				  		);
					}
					elsif(lc($tophat_aligner) eq "bowtie" or lc($tophat_aligner) eq "bowtie1"){
						#Collecting specific bowtie2 index
						$bowtieindex=$args{"bowtie1index"}; #Genome index to align your reads in .ebwt format
						#Calling Bowtie2 function
						$output_file2=TopHat( 
							file=>$file,
							threads=>$threads,
							bowtieindex=>$bowtieindex,
							verbose=>$verbose, 
							logfile=>$logfile, 
							library_type=>$library_type, 
							tophat_seg_length=>$tophat_seg_length,
							tophat_seg_mismatches=>$tophat_seg_mismatches,
							tophat_multihits=>$tophat_multihits,
							read_mismatches=>$read_mismatches,
							statsfile=>$statsfile, 
							tophatParameters=>$tophatParameters,
							projectdir=>$projectdir,
							miARmaPath=>$miARmaPath,
							tophat_aligner=>$tophat_aligner,
							GTF=>$GTF,
							Seqtype=>$Seqtype,
							
				  		);
					}
					elsif(lc($tophat_aligner) eq "bowtie1-bowtie2" or lc($tophat_aligner) eq "bowtie2-bowtie1"){
						
						#Collecting specific bowtie2 index
						$bowtieindex=$args{"bowtie2index"}; #Genome index to align your reads in .bt2 format
						$tophat_aligner="bowtie";
						#Calling Bowtie2 function
						$output_file2=TopHat( 
							file=>$file,
							threads=>$threads,
							bowtieindex=>$bowtieindex,
							verbose=>$verbose, 
							logfile=>$logfile, 
							library_type=>$library_type, 
							tophat_seg_length=>$tophat_seg_length,
							tophat_seg_mismatches=>$tophat_seg_mismatches,
							tophat_multihits=>$tophat_multihits,
							read_mismatches=>$read_mismatches,
							statsfile=>$statsfile, 
							tophatParameters=>$tophatParameters,
							projectdir=>$projectdir,
							miARmaPath=>$miARmaPath,
							tophat_aligner=>$tophat_aligner,
							GTF=>$GTF,
							Seqtype=>$Seqtype,
							
				  		);
						#Calling Bowtie1 function
						$bowtieindex=$args{"bowtie1index"}; #Genome index to align your reads in .ebwt format
						$tophat_aligner="bowtie1";
						my $output_file=TopHat( 
							file=>$file,
							threads=>$threads,
							bowtieindex=>$bowtieindex,
							verbose=>$verbose, 
							logfile=>$logfile, 
							library_type=>$library_type, 
							tophat_seg_length=>$tophat_seg_length,
							tophat_seg_mismatches=>$tophat_seg_mismatches,
							tophat_multihits=>$tophat_multihits,
							read_mismatches=>$read_mismatches,
							statsfile=>$statsfile, 
							tophatParameters=>$tophatParameters,
							projectdir=>$projectdir,
							miARmaPath=>$miARmaPath,
							tophat_aligner=>$tophat_aligner,
							GTF=>$GTF,
							Seqtype=>$Seqtype,
							
				  		);
						push(@results,$output_file2,$output_file)
					}
				
					if(!defined $bowtieindex){
						print STDERR "ERROR : You are requiesting to use ($aligner) $tophat_aligner, but no Indexed files are found on $bowtieindex. Please fill the bowtie2index/bowtie1index properly\n\n";
						exit;
					}
			  		return($output_file2);
				}
				else{
					die("READALIGNMENT ERROR :: ".date()." Provided GTF argument ($GTF) has an invalid value");
				}
			}
			elsif(lc($aligner) eq "bwa"){
				#Collecting specific bowtie1 index
				my $bwaindex=$args{"bwaindex"}; #Genome index to align your reads in .ebwt format
				#Optional parameters are predefined as undef
				my $threads=undef;
				my $Seqtype=$args{"Seqtype"};
				if(defined($args{"threads"})){
					$threads=$args{"threads"}; 
				}
				if(defined($args{"Seqtype"})){
					$Seqtype=$args{"Seqtype"};
				}else{
					$Seqtype="SingleEnd";
				}

				if($bwaindex){
					#Calling Bowtie1 function
					$output_file1=bwa( 
						file=>$file,
						threads=>$threads,
						bwaindex=>$bwaindex,
						verbose=>$verbose, 
						logfile=>$logfile, 
						statsfile=>$statsfile, 
						projectdir=>$projectdir,
						miARmaPath=>$miARmaPath,
						Seqtype=>$Seqtype,
						bwaparameters=>$args{"bwaparameters"} || undef
			  		);
			 		return($output_file1);
			  	}
			  	else{
			  		die("READALIGNMENT ERROR:: ".date()."Index argument ($bwaindex) has not been provided");
			  	}
			}
			elsif(lc($aligner) eq "mirdeep"){
				
				#Collecting specific bowtie1 index
				my $bowtie1index=$args{"bowtie1index"}; #Genome index to align your reads in .ebwt format
				#Optional parameters are predefined as undef
				my $threads=undef;
				my $adapter=undef;
				my $adapter_file=$args{"adapter_file"};
				
				my $mature;
				my $precursors;
				my $genome;
				if(defined($args{"threads"})){
					$threads=$args{"threads"}; 
				}
				if(defined($args{"adapter"})){
					$adapter=$args{"adapter"}; 
				}
				if(defined($args{"precursors"})){
					$precursors=$args{"precursors"}; 
				}
				if(defined($args{"mature"})){
					$mature=$args{"mature"}; 
				}
				if(defined($args{"genome"})){
					$genome=$args{"genome"}; 
				}
				#if we have defined a different adapter for each read-file
				my $file_name=fileparse($file);
				
				if($adapter_file){
					open(ADAPT_FILE,$adapter_file) || die "miRDeeP ERROR :: Can't file adapt_file at $adapter_file ($!)\n";
					while(<ADAPT_FILE>){
						chomp;
						my ($current_file,$current_adapter)=split(/\t/);
						if($current_file eq $file_name){
							$adapter=$current_adapter;
						}
						if($current_file eq $file_name.".bz2"){
							$adapter=$current_adapter;
						}
						if($current_file eq $file_name.".gz"){
							$adapter=$current_adapter;
						}
					}
					close ADAPT_FILE;
				}
				
				if($bowtie1index and $mature and $precursors and $genome and $organism){
					#Calling Bowtie1 function
					$output_file1=miRDeep( 
						file=>$file,
						threads=>$threads,
						bowtie1index=>$bowtie1index,
						verbose=>$verbose, 
						logfile=>$logfile, 
						statsfile=>$statsfile, 
						projectdir=>$projectdir,
						miARmaPath=>$miARmaPath,
						adapter=>$adapter || undef,
						mature=>$mature,
						precursors=>$precursors,
						genome=>$genome,
						organism=>$organism,
						miRDeeparameters=>$args{"miRDeeparameters"}
			  		);
			 		return($output_file1);
			  	}
			  	else{
			  		die("READALIGNMENT ERROR:: ".date()."Index argument ($bowtie1index) , adapter ($adapter), mature fasta file ($mature), organism($organism), precursor fasta file ($precursors) or fasta genome ($genome) has not been provided");
			  	}
			}
			elsif(lc($aligner) eq "hisat2"){
				#Collecting specific hisat2 index
				my $hisat2index=$args{"hisat2index"}; #Genome index to align your reads in .bt2 format
				#Optional parameters are predefined as undef
				my $threads=undef;
				my $hisat2parameters=undef;
				my $Seqtype;
				#If any optional parameter is provided by the user will be collected
				if(defined($args{"threads"})){
					$threads=$args{"threads"}; 
				}
				if(defined($args{"hisat2parameters"})){
					$hisat2parameters=$args{"hisat2parameters"};
				}
				if(defined($args{"Seqtype"})){
					$Seqtype=$args{"Seqtype"};
				}else{
					$Seqtype="SingleEnd";
				}
				#Calling hisat2 function
				$output_file2=hisat2( 
					file=>$file,
					threads=>$threads,
					hisat2index=>$hisat2index,
					verbose=>$verbose, 
					logfile=>$logfile, 
					statsfile=>$statsfile, 
					hisat2parameters=>$hisat2parameters,
					projectdir=>$projectdir,
					miARmaPath=>$miARmaPath,
					adapter=>$adapter,	
					Seqtype=>$Seqtype,
					strand=>$args{"strand"}
		  		);
		  		return($output_file2);
			}
			elsif(lc($aligner) eq "star"){
				#Collecting specific star2 index
				my $starindex=$args{"starindex"}; #Genome index to align your reads in .bt2 format
				#Optional parameters are predefined as undef
				my $threads=undef;
				my $starparameters=undef;
				my $Seqtype;
				#If any optional parameter is provided by the user will be collected
				if(defined($args{"threads"})){
					$threads=$args{"threads"}; 
				}
				if(defined($args{"starparameters"})){
					$starparameters=$args{"starparameters"};
				}
				if(defined($args{"Seqtype"})){
					$Seqtype=$args{"Seqtype"};
				}else{
					$Seqtype="SingleEnd";
				}
				#Calling star function
				$output_file2=star( 
					file=>$file,
					threads=>$threads,
					starindex=>$starindex,
					verbose=>$verbose, 
					logfile=>$logfile, 
					statsfile=>$statsfile, 
					starparameters=>$starparameters,
					projectdir=>$projectdir,
					miARmaPath=>$miARmaPath,
					adapter=>$adapter,	
					Seqtype=>$Seqtype
								
		  		);
		  		return($output_file2);
			}
			else{
				die("READALIGNMENT ERROR :: ".date()." Provided aligner argument ($aligner) has an invalid value. Allowed values for this parameters: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1, tophat and bwa");
			}
		}
		else{
			#Registering the error
	   		open(LOG,">> ".$logfile) || die "READALIGNMENT ERROR :: ".date()." Can't open logfile '$logfile': $!";
	    	print LOG "READALIGNMENT ERROR :: ".date()." File ($file), logfile($logfile), statsfile ($statsfile), aligner($aligner) and/or projectdir($projectdir) have not been provided";
	    	close LOG;

			#If mandatory parameters have not been provided program will die and show error message
			warn ("READALIGNMENT ERROR :: ".date()." File ($file), logfile($logfile), statsfile ($statsfile), aligner($aligner) and/or projectdir($projectdir) have not been provided");
			help_ReadAligment();
		}
	}
	else{
		next;
	}
	return(@results);
	
	sub help_ReadAligment{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[file] Path of the file which is going to be aligned (fasta/fastq format)
			[logfile] Path of run.log file where execution data will be saved
			[statsfile] Path of stats.log file where stats data will be saved
			[projectdir] Directory where results directory will be created
			[aligner] Aligner which will be use in the analysis (Allowed values: Bowtie1, Bowtie2 and Bowtie1-Bowtie2/Bowtie2-Bowtie1)
			Optional parameters:
			[bowtie1index] Indexed genome to align your reads in format .ebwt (Mandatory for analysis with bowtie1)
			[bowtie2index] Indexed genome to align your reads in format .bt2 (Mandatory for analysis with bowtie2)
			[bwaindex] Indexed genome to align your reads in format .bwt (Mandatory for analysis with bwa)
			[threads] Optional number of threads to perform the analysis  
			[bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
			[bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
			[bowtie1parameters] Other bowtie parameters to perform the analysis using the bowtie1 recommended syntaxis
			[bowtie2parameters] Other bowtie parameters to perform the analysis using the bowtie2 recommended syntaxis
			
			[verbose] Option to show the execution data on the screen   
						             
			Examples:
			ReadAligment(file=>"file.fastq", aligner=>"Bowtie1-Bowtie2", threads=>"4", bowtie1index=>"./genomeindex1/hg19", bowtie2index=>"./genomeindex2/hg19",verbose=>"verbose", logfile=>"run.log", bowtiemiss=>"0", bowtieleng=>"19",statsfile=>"stats.log", bowtie1parameters=>" --sam --best --nomaqround -e 70 -k 1", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}
	
}



=head2 bowtie1_index

  Example    : 
  bowtie1_index(
  	fasta=>"genome.fasta",
  	dir=>".",
  	logfile=>"run.log",
  	indexname=>"hg19"
  );
  Description: This function collects the directory and the path of genome sequence in fasta
  format to build a new index. This index can be used to perform a bowtie1 analysis. The 
  execution data will be printed on run.log file. This fuction returns the path of the new 
  bowtie1 index and saves the new index on the directory named Bowtie1_index at directory 
  provided by the user.
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where new index will be saved
  	 [fasta] Path of the genomic fasta sequence to build the index
  	 [logfile] Path of run.log file where execution data will be saved
  	 [indexname] Name to write in the index files
  Returntype : Directory named Bowtie1_index with index genome. bowtie1_index also returns the path of the new bowtie1 index 
  Requeriments: bowtie1_index function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v1.0.0 or higher software correctly installed
  	- Input genomic sequence on fasta format  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub bowtie1_index{

	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/bowtie1/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/bowtie1/";
	}
	#First, check that bowtie1 is in path:
	my @bowtie1_bin=`which bowtie`;
	#Executing the command
	if(scalar(@bowtie1_bin)<1){
		die "BOWTIE1_INDEX ERROR :: system args failed: $? : Is bowtie1 installed and exported to \$PATH ?";
	}
	my $fasta=$args{"fasta"}; #the path of genome sequence in fasta format
	my $dir=$args{"dir"}; #directory to create the genomeindex directory where new index will be saved
	my $logfile=$args{"logfile"}; #path of the logfile to write the execution data
	my $indexname=$args{"indexname"}; #name to write in the index files
	
	#Variable declaration
	my $index_output;
	my $command;
	my $commanddef;
	
	#Checking the mandatory arguments
	if ($fasta and $dir and $indexname and $logfile){
		print STDERR date()." Generating the index genome $indexname from $fasta. This process could take some hours.\n";
		#bowtie-build execution command from a fasta file. The output index will be saved
		#in the genomeindex1 directory with the name index 
		my $command="bowtie-build -f ".$fasta." ".$dir."/Bowtie1_index/".$indexname;
		#commandef is the command that will be executed by system composed of the index
		#directory creation, the module loading, and the bowtie-build execution. The error 
		#data will be redirected to the run.log file
		$commanddef="mkdir -p ".$dir."/Bowtie1_index ; ".$command ." >> ".$logfile." 2>&1";
		#Printing the date and command execution on the run.log file
		open (LOG,">> ".$logfile) || die "BOWTIE1_INDEX ERROR :: Can't open $logfile: $!";
		print LOG "BOWTIE1_INDEX :: ".date()." Executing $commanddef\n";
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		or die "BOWTIE1_INDEX ERROR :: system args failed: $? ($commanddef)";
		close LOG;
		#Returning the path of the new bowtie1 index 
		return($dir."/Bowtie1_index/".$indexname);
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "BOWTIE1_INDEX ERROR :: Can't open $logfile: $!";
    	print LOG "BOWTIE1_INDEX ERROR :: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile) and/or fasta file($fasta) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program will die and show error message
		warn ("BOWTIE1_INDEX ERROR:: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile) and/or fasta file($fasta) have not been provided");
		help_bowtie1_index();
	}	

	sub help_bowtie1_index{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[dir] Input directory where new index will be saved
  	 		[fasta] Path of the genomic fasta sequence to build the index
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[indexname] Name to write in the index files
						             
			Examples:
			bowtie1_index(fasta=>"genome.fasta", dir=>".", logfile=>"run.log", indexname=>"hg19");

	};

	print STDERR $usage;
	exit(); 
	}
}

=head2 bowtie1

  Example    : 
  bowtie1( 
		file=>"./file.fastq",
		threads=>"4",
		bowtieindex=>"./hg19",
		verbose=>"verbose", 
		logfile=>"run.log", 
		bowtiemiss=>"0", 
		bowtieleng=>"19",
		statsfile=>"stats.log", 
		bowtieparameters=>"--best --nomaqround -e 70 -k 1",
		projectdir=>"."
  	);
  Description: Bowtie1 takes the reads of the provided fastq file and align them with the bowtie 
  index to generate a sam file which will be saved on Bowtie1_results on the project directory. 
  Execution and stats data will be saved on run.log and stats.log file and will show on the screen 
  if verbose option is selected. In this case the bowtieindex has to be built with bowtie1_index function.  
  Input parameters: 
	Mandatory parameters:
  	 [file] Path of the file which is going to be aligned (fasta/fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 [bowtieindex]  Indexed genome to align your files in format .ebwt
  	 [projectdir] Directory where results directory will be created
  	 Optional parameters:
  	 [threads] Optional number of threads to perform the analysis  
  	 [bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 [bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 [bowtieparameters] Other bowtie parameters to perform the analysis using the bowtie recommended syntaxis
  	 [verbose] Option to show the execution data on the screen
	 [Seqtype] Sequencing method. SingleEnd by default. Acepted values : [Paired-End|Single-End]
  Returntype : File at directory Bowtie1_results.Also returns the path of the output file 
  Requeriments: bowtie1 function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v1.0.0 or higher software correctly installed
  	- Input files on fastq format on the provided directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub bowtie1{

	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/bowtie1/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/bowtie1/";
	}
	
	#First, check that fastqc is in path:
	my @bowtie1_bin=`which bowtie`;
	#Executing the command
	if(scalar(@bowtie1_bin)<1){
		die "BOWTIE1_INDEX ::system args failed: $? : Is bowtie1 installed and exported to \$PATH ?";
	}

	my $file=$args{"file"}; #Name of the file which is going to be aligned
	my $bowtieindex=$args{"bowtieindex"}; #Genome index in format .ebwt 
	my $threads=$args{"threads"}; #Optional number of threads to perform the analysis  
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #Path of the statsfile to write the stats data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $Seqtype=$args{"Seqtype"}; #Sequencing method. SingleEnd by default. Acepted values : [Paired-End|Single-End]
	my $adapter=$args{"adapter"}; #are this reads processed by miArma 
	
	# Variable declaration and describing results directory 
	my $commanddef;
	my $output_dir="/Bowtie1_results/";
	my $mate_file=$file;
	
	#Variable to collect the optional parameters 
	my $bowtiepardef= "--sam ";
	#The number of missmatches can be provided by the user
	if(defined($args{"bowtiemiss"})){
		my $bowtiemiss=$args{"bowtiemiss"};
		$bowtiepardef.=" -n $bowtiemiss";
	}
	#Seed length can be provided by the user
	if(defined($args{"bowtieleng"})){
		my $bowtieleng=$args{"bowtieleng"};
		$bowtiepardef.=" -l $bowtieleng";
	}
	#Number of threads can be provided by the user
	if($threads>0){
		$bowtiepardef.= " -p $threads";
	}
	#Any other bowtie parameter can be provided by the user using the correct sintaxis
	if(defined($args{"bowtieparameters"})){
		my $bowtieparameters=$args{"bowtieparameters"};
		$bowtiepardef.=" $bowtieparameters";
	}
	#Checking the mandatory parameters
	if ($file and $projectdir and $bowtieindex and $logfile and $statsfile){ 
		#Extracting the name of the file
		my $name=fileparse($file, qr{\.f.*});
		my $output_file_final=$name;
		$output_file_final=~s/_1//g;
		#Bowtie execution command
		my $command;
		my $compressed_file=0;
		
		my $output_file_bw;
		if($adapter){
			$output_file_bw=$projectdir.$output_dir.$output_file_final."_bw1.bam";
		}
		else{
			$output_file_bw=$projectdir.$output_dir.$output_file_final."_nat_bw1.bam";
		}
		
		if(lc($Seqtype) eq "pairedend" or lc($Seqtype) eq "paired" or lc($Seqtype) eq "paired-end"){
			#Check if the file is a paired-end file
			if($file =~ /.*_R?1.*/){
				#it contains the _1 label
				my $mate_file=$file;
				if($file =~  /.*_R1.*/){
					$mate_file=~s/_R1(.+)/_R2$1/g;
				}
				else{
					$mate_file=~s/_1(.+)/_2$1/g;
				}
				if(-e $mate_file){
					if($file ne $mate_file){
						print STDERR date()." Checking $file for bowtie1 (Paired-End) analysis\n" if($verbose);
						if($file =~ /\.gz$/){
							print STDERR date()." Uncompressing $file\n" if($verbose);
							#In case gzip
							system("gunzip -f -k $file");
							system("gunzip -f -k $mate_file");
							#Changing new extension
							$file=~s/\.gz$//g;
							$mate_file=~s/\.gz$//g;
							$command="bowtie ".$bowtiepardef." ".$bowtieindex." -1 $file -2 $mate_file ". $output_file_bw;
							$compressed_file=1;
						}
						elsif($file =~ /\.bz2$/){
							#in case bzip2
							print STDERR date()." Uncompressing $file\n" if($verbose);
							system("bunzip2 -f -d -k $file");
							system("bunzip2 -f -d -k $mate_file");
							#New extension
							$file=~s/\.bz2$//g;
							$mate_file=~s/\.bz2$//g;
							
							#$command="bowtie ".$bowtiepardef." ".$bowtieindex." -1 - ".`< bunzip2 -d -c -k $file`." -2 - ". `< bunzip2 -d -c -k $mate_file` ." " . $output_file_bw;
							$command="bowtie ".$bowtiepardef." ".$bowtieindex." -1 ".$file." -2 ". $mate_file ." " . $output_file_bw;
							$compressed_file=1;
						}
						else{
							$command="bowtie ".$bowtiepardef." ".$bowtieindex." -1 ".$file." -2 ". $mate_file ." " . $output_file_bw;
							$compressed_file=0;
						}
					}
				}
				else{
					print STDERR "ERROR:: You have requested a Paired-End analysis, so for the file $file a $mate_file file is needed\n";
					last;
				}
			}
			elsif($file =~ /.*_R?2.*/){
				#The mate pair
			}else{
				print STDERR "ERROR :: Your input file ($file) does not seem to be a paired end. Check bowtie naming files for paired end analysis\n";
				return();
			}
		}
		else{
			print STDOUT "\tBOWTIE 1 :: ".date()." Checking $file for bowtie1 (single-end) analysis\n" if($verbose);			
			if($file =~ /\.gz$/){
				print STDOUT "\tBOWTIE 1 :: ".date()." Uncompressing $file\n" if($verbose);
				#In case gzip
				$command="gunzip -f -d -k -c $file | bowtie ".$bowtiepardef." ".$bowtieindex." - ". $output_file_bw;
				$file=~s/\.gz$//g;
				$compressed_file=1;
			}
			elsif($file =~ /\.bz2$/){
				#in case bzip2
				#New extension
				$command="bunzip2 -f -d -k -c $file | bowtie ".$bowtiepardef." ".$bowtieindex." - ". $output_file_bw ;
				$compressed_file=1;
				$file=~s/\.bz2$//g;
			}
			else{
				$command="bowtie ".$bowtiepardef." ".$bowtieindex." ".$file." ".$output_file_bw;
				$compressed_file=0;
			}			
		}
		if($verbose){
			#commandef is the command will be executed by system composed of the results directory creation 
			#and the bowtie execution. The stats data will be printed on the screen
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command;
			#Printing on the screen the date and the execution data
			print STDOUT "\tBOWTIE 1 :: ".date()." Executing $commanddef\n";
		}   	
		else
		{
			#Bowtie1 execution without verbose option, stats data will be printed on statsfile
			#and execution data on run.log file.
			#Opening stats file
			open (STATS,">> ".$statsfile) || die "BOWTIE1 ERROR :: Can't open $statsfile: $!";
			print STATS "BOWTIE1 :: File:".$file."\n";
			#commandef is the command will be executed by system composed of the results directory creation 
			#and the bowtie execution. The stats data will be redirected to the stats.log file
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command." 2>> ".$statsfile;
			
		}

		#Opening the run.log and printing the execution data
		open (LOG,">> ".$logfile) || die "BOWTIE1 ERROR :: Can't open $logfile: $!";
		print LOG "BOWTIE 1 :: ".date()." Executing $commanddef\n File:".$file."\n";
		close LOG;
			
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		or die "BOWTIE1 ERROR :: system args failed: $? ($commanddef)";
		close STATS;
		#The path of the output file is returned to the main program
		
		return($output_file_bw);
	}
	else
	{
		#Registering the error on run.log file
   		open(LOG,">> ".$logfile) || die "BOWTIE1 ERROR :: Can't open $logfile: $!";
    	print LOG "BOWTIE1 ERROR :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($bowtieindex) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program will die and show error message
		warn ("BOWTIE1 ERROR:: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($bowtieindex) have not been provided");
		help_bowtie1();
		

	}	

	sub help_bowtie1{
	    my $usage = qq{
		  	$0 

			Needed parameters:
  	 		[file] Path of the file which is going to be aligned (fasta/fastq format)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[statsfile] Path of stats.log file where stats data will be saved
  	 		[bowtieindex]  Indexed genome to align your files in format .ebwt
  	 		[projectdir] Directory where results directory will be created

  	 		Optional parameters:
  	 		[threads] Optional number of threads to perform the analysis  
  	 		[bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 		[bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 		[bowtieparameters] Other bowtie parameters to perform the analysis using the bowtie recommended syntaxis
  	 		[verbose] Option to show the execution data on the screen   
						             
			Examples:
			bowtie1(file=>"./file.fastq", threads=>"4", bowtieindex=>"./hg19", verbose=>"verbose", logfile=>"run.log", bowtiemiss=>"0", bowtieleng=>"19",statsfile=>"stats.log", bowtieparameters=>" --sam --best --nomaqround -e 70 -k 1", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}	
}

=head2 bowtie2_index

  Example    : 
  bowtie2_index(
  	fasta=>"genome.fasta", 
  	dir=>".",
  	logfile=>"run.log",
  	indexname=>"hg19"
  )
  Description: This function collects the directory and the path of genome sequence in fasta
  format to build a new index. This index can be used to perform a bowtie2 analysis. The 
  execution data will be printed on run.log file. This fuction returns the path of the new 
  bowtie2 index and saves the new index on the directory named Bowtie2_index at directory 
  provided by the user.
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where new index will be saved
  	 [fasta] Genomic fasta sequence to build the index
  	 [logfile] Path of run.log file where execution data will be saved
  	 [indexname] Name to write in the index files 
  Returntype : Directory named Bowtie2_index with index genome. Also returns the path of the new bowtie2 index 
  Requeriments: cutadapt function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v2.2.0 or higher software correctly installed
  	- Input genomic sequence on fasta format 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub bowtie2_index{
	
	##Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $threads=$args{"threads"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/bowtie2/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/bowtie2/";
	}
	
	#First, check that bowtie2 is in path:
	my @bowtie2_bin=`which bowtie2`;
	#Executing the command
	if(scalar(@bowtie2_bin)<1){
		die "BOWTIE2 ERROR ::system args failed: $? : Is bowtie2 installed and exported to \$PATH (".$ENV{PATH}.") ?";
	}
	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname 
	#and logfile are mandatory arguments.
	my $fasta=$args{"fasta"}; #the path of genome sequence in fasta format
	my $dir=$args{"dir"}; #directory to create the genomeindex directory where new index will be saved
	my $logfile=$args{"logfile"}; #path of the logfile to write the execution data
	my $indexname=$args{"indexname"}; #path of the logfile to write the execution data
	
	#Variable declaration
	my $command;
	my $commanddef;

	#Checking the mandatory arguments
	if ($fasta and $dir and $logfile and $indexname){
		print STDERR date()." Generating the index genome $indexname from $fasta. This process could take some hours\n";
		#bowtie2-build execution command from a fasta file. The output index will be saved
		#in the genomeindex2 directory with the name index 
		$command= "bowtie2-build --threads $threads -f ".$fasta." ".$dir."/Bowtie2_index/".$indexname;
		#commandef is the command that will be executed by system composed of the index
		#directory creation, the module loading, and the bowtie2-build execution. The error 
		#data will be redirected to the run.log file
		$commanddef="mkdir -p ".$dir."/Bowtie2_index ; ".$command ." >> ".$logfile." 2>&1";
		#Printing the date and command execution on the run.log file
		open (LOG,">> ".$logfile) || die "BOWTIE2_INDEX ERROR :: Can't open $logfile: $!";
		print LOG "BOWTIE2_INDEX :: ".date()." Executing $commanddef\n";
		
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		or die "BOWTIE2_INDEX ERROR :: system args failed: $? ($command)";	
		close LOG;
		#Returning the path of the new bowtie2 index 
		return($dir."/Bowtie2_index/".$indexname);
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "BOWTIE2_INDEX ERROR :: Can't open $logfile: $!";
    	print LOG "BOWTIE2_INDEX ERROR :: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile) and/or fasta file($fasta) have not been provided";

		#If mandatory parameters have not been provided program dies and shows error message
		warn ("BOWTIE2_INDEX ERROR :: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile) and/or fasta file($fasta) have not been provided");
		help_bowtie2_index();
	}

	sub help_bowtie2_index{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[dir] Input directory where new index will be saved
  	 		[fasta] Path of the genomic fasta sequence to build the index
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[indexname] Name to write in the index files
			[threads] Number of threads to use			             
			Examples:
			bowtie2_index(fasta=>"genome.fasta", dir=>".", logfile=>"run.log", indexname=>"hg19", threads=>4);

	};

	print STDERR $usage;
	exit(); 
	}
}


=head2 bowtie2

  Example    : 
  bowtie2( 
	file=>"./file.fastq",
	threads=>"4",
	bowtieindex=>"./hg19",
	verbose=>"verbose", 
	logfile=>"run.log", 
	bowtiemiss=>"0", 
	bowtieleng=>"19",
	statsfile=>"stats.log", 
	bowtieparameters=>" -I 50 -X 200",
	projectdir=>"."
  )
  Description: Bowtie2 takes the reads and align them with the bowtie index to generate a sam file which 
  will be saved on Bowtie2_results on the project directory. Execution and stats data will be saved on directory stats.log 
  file and will show on screen if verbose option is selected. In this case the bowtieindex 
  has to be built with bowtie2_index 
  function. 
  Input parameters: 
	Mandatory parameters:
  	 [file] Name of the file which is going to be align (fasta/fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 [bowtieindex]  Indexed genome to align your files in format .bt2
  	 [projectdir] Directory where bowtie2_results directory will be created
  	 Optional parameters:
  	 [threads] Optional number of threads to perform the analysis
  	 [bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 [bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 [bowtieparameters] Other bowtie parameters to perform the analysis using the bowtie recommended syntaxis
  	 [verbose] Option to show the execution data on the screen   
  Returntype : File at directory Bowtie2_results. Also returns the path of the output file
  Requeriments: bowtie2 function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v2.2.0 or higher software correctly installed
  	- Input files on fastq format on the provided directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub bowtie2{
	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/bowtie2/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/bowtie2/";
	}
	
	#First, check that bowtie2 is in path:
	my @bowtie2_bin=`which bowtie2`;
	#Executing the command
	if(scalar(@bowtie2_bin)<1){
		die "BOWTIE2 ERROR ::system args failed: $? : Is bowtie2 installed and exported to \$PATH (".$ENV{PATH}.") ?";
	}

	my $file=$args{"file"}; #File which is going to be aligned
	my $bowtieindex=$args{"bowtieindex"}; #Genome index in format .bt2
	my $threads=$args{"threads"}; #Optional number of threads to perform the analysis faster
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #path of the statsfile to write the stats data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $Seqtype=$args{"Seqtype"}; #Sequencing method. SingleEnd by default. Acepted values : [Paired-End|Single-End]
	my $adapter=$args{"adapter"}; #are this reads processed by miArma 
	
	#Variable declaration and describing results directory 
	my $commanddef;
	my $output_dir="/Bowtie2_results/";

	#Collecting bowtie2 parameters provided by the user
	my $bowtiepardef= "";
	#The number of missmatches can be provided by the user
	if(defined($args{"bowtiemiss"})){
		my $bowtiemiss=$args{"bowtiemiss"};
		$bowtiepardef.=" -N $bowtiemiss";
	}
	#Seed length can be provided by the user
	if(defined($args{"bowtieleng"})){
		my $bowtieleng=$args{"bowtieleng"};
		$bowtiepardef.=" -L $bowtieleng";
	}
	#Number of threads can be provided by the user
	if($threads>0){
		$bowtiepardef.= " -p $threads";
	}
	#Any other bowtie parameter can be provided by the user using the correct sintaxis
	if(defined($args{"bowtieparameters"})){
		my $bowtieparameters=$args{"bowtieparameters"};
		$bowtiepardef.=" $bowtieparameters";
	}
	
	#Checking the mandatory parameters
	if ($file and $projectdir and $bowtieindex and $logfile and $statsfile){ 
		#Extracting the name of the file
		my $name=fileparse($file, qr{\.f.*});
		my $output_file_bw2;
		if($adapter){
			$output_file_bw2=$projectdir.$output_dir.$name."_bw2.sam";
		}
		else{
			$output_file_bw2=$projectdir.$output_dir.$name."_nat_bw2.bam";
		}
		#Bowtie2 execution command
		my $command;
		if(lc($Seqtype) eq "pairedend" or lc($Seqtype) eq "paired" or lc($Seqtype) eq "paired-end"){
			print STDOUT "\tBOWTIE 2 :: ".date()." Checking $file for bowtie2 (paired-end) analysis\n" if($verbose);
			
			#Check if the file is a paired-end file
			if($file =~ /.*_R?1.*/){
				#it contains the _1 label
				my $mate_file=$file;
				if($file =~  /.*_R1.*/){
					$mate_file=~s/_R1(.+)/_R2$1/g;
				}
				else{
					$mate_file=~s/_1(.+)/_2$1/g;
				}
				if(-e $mate_file){
					if($file ne $mate_file){
						$command="bowtie2 ".$bowtiepardef." -x ".$bowtieindex." -1 ".$file." -2 ". $mate_file ." --met-file ".$projectdir.$output_dir.$name.".metrics -S --un ".$projectdir.$output_dir.$name."_no_aligned.fastq -S ". $output_file_bw2;
					}
				}
				else{
					print STDERR "ERROR:: You have requested a Paired-End analysis, so for the file $file a $mate_file file is needed\n";
					last;
				}
			}
			else{
				return();
			}
		}
		else{
			print STDOUT "\tBOWTIE 2 :: ".date()." Checking $file for bowtie2 (single-end) analysis\n" if($verbose);
			$command="bowtie2 ".$bowtiepardef." -x ".$bowtieindex." -U ".$file." --met-file ".$projectdir.$output_dir.$name.".metrics --un ".$projectdir.$output_dir.$name."_no_aligned.fastq -S ". $output_file_bw2;
		}
		
		#Bowtie execution with verbose option
		if($verbose){
			#commandef is the command will be executed by system composed of the results directory 
			#creation and the bowtie2 execution. The stats data will be printed on the screen
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command;
			#Printing on the screen the date and the execution data
			print STDOUT "\tBOWTIE2 :: ".date()." Executing $commanddef\n";
		}   
		#Bowtie execution printing stats on stats file	
		else{
			#Opening and printing stats on statsfile
			open (STATS,">> ".$statsfile) || die "BOWTIE2 ERROR :: Can't open $statsfile: $!";
			print STATS "BOWTIE2 :: File:".$file."\n";
			#commandef is the command will be executed by system composed of the results directory creation
			#and the bowtie execution. The stats data will be redirected to the stats.log file
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command." 2>> ".$statsfile;
		}
		#Opening the run.log and printing the execution data
		open (LOG,">> ".$logfile) || die "BOWTIE2 ERROR :: Can't open $logfile: $!";
		print LOG "BOWTIE2 :: ".date()." Executing $commanddef\n";
		close LOG;

		#print STDERR "$command\n";
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		 or die "BOWTIE2 ERROR :: system args failed: $? ($commanddef)";
		close STATS;
		#The path of the output file is returned to the main program
		return ($output_file_bw2);		
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "BOWTIE2 ERROR :: Can't open $logfile: $!";
    	print LOG "BOWTIE2 :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($bowtieindex) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn ("BOWTIE2 :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($bowtieindex) have not been provided");
		help_bowtie2();
	}
	sub help_bowtie2{
	    my $usage = qq{
		  	$0 

			Needed parameters:
  	 		[file] Name of the file which is going to be align (fasta/fastq format)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[statsfile] Path of stats.log file where stats data will be saved
  	 		[bowtieindex]  Indexed genome to align your files in format .bt2
  	 		[projectdir] Directory where bowtie2_results directory will be created

  	 		Optional parameters:
  	 		[threads] Optional number of threads to perform the analysis
  	 		[bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 		[bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 		[bowtieparameters] Other bowtie parameters to perform the analysis using the bowtie recommended syntaxis
  	 		[verbose] Option to show the execution data on the screen   
						             
			Examples:
			bowtie2(file=>"file.fastq", threads=>"4", bowtieindex=>"./hg19", verbose=>"verbose", logfile=>"run.log", bowtiemiss=>"0", bowtieleng=>"19",statsfile=>"stats.log", bowtieparameters=>" -I 50 -X 200", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}	
}

=head2 TopHat

  Example    : 
  TopHat( 
	file=>"./file.fastq",
	threads=>"4",
	bowtieindex=>"./hg19",
	verbose=>"verbose", 
	logfile=>"run.log", 
	projectdir=>".".
  )
  Description: TopHat takes the reads and align them with the bowtie index to generate a sam file which 
  will be saved on Bowtie2_results (or Bowtie1_results) on the project directory. Execution and stats data will be saved on directory stats.log 
  file and will show on screen if verbose option is selected. In this case the bowtieindex 
  has to be built with bowtie2_index 
  function. 
  Input parameters: 
	Mandatory parameters:
  	 [file] Name of the file which is going to be align (fasta/fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 [bowtieindex]  Indexed genome to align your files in format .bt2
  	 [projectdir] Directory where bowtie2_results directory will be created
  	 Optional parameters:
  	 [threads] Optional number of threads to perform the analysis
  	 [verbose] Option to show the execution data on the screen   
  	 [tophat_aligner] Aligner to use within tophat [bowtie,bowtie2 or bowtie1-bowtie2| default=bowtie2]   
  	 [tophat_multihits] Instructs TopHat to allow up to this many alignments to the reference for a given read, and choose the alignments based on their alignment scores if there are more than this number. The default is 20 for read mapping    
  	 [tophat_seg_mismatches] Read segments are mapped independently, allowing up to this many mismatches in each segment alignment. The default is 2.  
  	 [tophat_seg_length] Each read is cut up into segments, each at least this long. These segments are mapped independently. The default is 25.
	 [tophat_read_mismatches] Final read alignments having more than these many mismatches are discarded. The default is 2.
  	 [tophatParameters] Other parameter as explained in http://ccb.jhu.edu/software/tophat/manual.shtml   
  Returntype : File at directory Bowtie2_results. Also returns the path of the output file
  Requeriments: TopHat function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Input files on fastq format on the provided directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub TopHat{
	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $tophat_aligner=$args{"tophat_aligner"};
	my $output_dir;
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/tophat/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/tophat/";
	}

	if(lc($tophat_aligner) eq "bowtie1"){
		$tophat_aligner="bowtie1";
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/bowtie1/";
		$output_dir="Bowtie1_results/";
	}
	else{
		$tophat_aligner="bowtie";
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/bowtie2/";
		$output_dir="Bowtie2_results/";
	}
	#First, check that bowtie2 is in path:
	my @tophat_bin=`which tophat`;
	#Executing the command
	if(scalar(@tophat_bin)<1){
		die "TOPHAT ERROR ::system args failed: $? : Is tophat installed and exported to \$PATH ?";
	}

	my $file=$args{"file"}; #File which is going to be aligned
	my $bowtieindex=$args{"bowtieindex"}; #Genome index in format .bt2
	my $threads=$args{"threads"}; #Optional number of threads to perform the analysis faster
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #path of the statsfile to write the stats data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $GTF=$args{"GTF"};
	my $Seqtype=$args{"Seqtype"}; #Sequencing method. SingleEnd by default. Acepted values : [Paired-End|Singe-End]
	#Variable declaration and describing results directory 
	my $commanddef;

	#Collecting bowtie2 parameters provided by the user
	my $bowtiepardef= "";
	#The number of missmatches can be provided by the user
	if(defined($args{"library_type"})){
		my $library_type=$args{"library_type"};
		$bowtiepardef.=" --library-type $library_type";
	}
	#Seed length can be provided by the user
	if(defined($args{"tophat_seg_length"})){
		my $tophat_seg_length=$args{"tophat_seg_length"};
		$bowtiepardef.=" --segment-length $tophat_seg_length";
	}
	#Seed length can be provided by the user
	if(defined($args{"tophat_seg_mismatches"})){
		my $tophat_seg_mismatches=$args{"tophat_seg_mismatches"};
		$bowtiepardef.=" --segment-mismatches $tophat_seg_mismatches";
	}
	#Seed length can be provided by the user
	if(defined($args{"read_mismatches"})){
		my $read_mismatches=$args{"read_mismatches"};
		$bowtiepardef.=" --read-mismatches $read_mismatches";
	}
	#Seed length can be provided by the user
	if(defined($args{"tophat_multihits"})){
		my $tophat_multihits=$args{"tophat_multihits"};
		$bowtiepardef.=" --max-multihits $tophat_multihits";
	}
	#Number of threads can be provided by the user
	if($threads>0){
		$bowtiepardef.= " -p $threads";
	}
	#Any other bowtie parameter can be provided by the user using the correct sintaxis
	if(defined($args{"tophatParameters"})){
		my $tophatParameters=$args{"tophatParameters"};
		$bowtiepardef.=" $tophatParameters";
	}
	
	if(lc($args{"tophat_aligner"}) eq "bowtie1"){
		my $tophat_aligner=lc($args{"tophat_aligner"});
		$bowtiepardef.=" --$tophat_aligner";
	}
	if(lc($args{"tophat_aligner"}) eq "bowtie2"){
		$bowtiepardef.=" --b2-very-fast";
	}

	#Checking the mandatory parameters
	if ($file and $projectdir and $bowtieindex and $logfile and $statsfile and $GTF and $bowtieindex){ 
		#Extracting the name of the file
		my $name=fileparse($file, qr{\.f.*});
		my $command;
		if(lc($Seqtype) eq "pairedend" or lc($Seqtype) eq "paired" or lc($Seqtype) eq "paired-end"){
			#Check if the file is a paired-end file
			if($file =~ /.*_R?1.*/){
				#it contains the _1 label
				my $mate_file=$file;
				if($file =~  /.*_R1.*/){
					$mate_file=~s/_R1(.+)/_R2$1/g;
				}
				else{
					$mate_file=~s/_1(.+)/_2$1/g;
				}
				if(-e $mate_file){
					if($file ne $mate_file){
						print STDOUT "\tTOPHAT :: ".date()." Checking $file for TopHat analysis using $tophat_aligner\n" if($verbose);
						$command="tophat --no-coverage-search ".$bowtiepardef." --GTF ". abs_path($GTF) . " --transcriptome-index=$projectdir/$output_dir/transcriptome/tr -o $projectdir/$output_dir" ." ".$bowtieindex." ".$file ." ". $mate_file;
					}
				}
				else{
					print STDERR "ERROR:: You have requested a Paired-End analysis, so for the file $file a $mate_file file is needed\n";
					last;
				}
			}
			else{
				return();
			}
		}
		else{
			print STDOUT "\tTOPHAT :: ".date()." Checking $file for TopHat analysis using $tophat_aligner\n" if($verbose);
			$command="tophat --no-coverage-search ".$bowtiepardef." --GTF ". abs_path($GTF) . " --transcriptome-index=$projectdir/$output_dir/transcriptome/tr -o $projectdir/$output_dir" ." ".$bowtieindex." ".$file;
		}
		#tophat execution command
		#tophat execution with verbose option
		if($verbose){
			#commandef is the command will be executed by system composed of the results directory 
			#creation and the bowtie2 execution. The stats data will be printed on the screen
			$commanddef= "mkdir -p ".$projectdir."/".$output_dir." ;".$command;
			#Printing on the screen the date and the execution data
			print STDOUT "\tTOPHAT :: ".date()." Executing $commanddef\n" ;
		}   
		#tophat execution printing stats on stats file	
		else{
			#Opening and printing stats on statsfile
			open (STATS,">> ".$statsfile) || die "TOPHAT ERROR :: Can't open $statsfile: $!";
			print STATS "TOPHAT :: File:".$file."\n";
			#commandef is the command will be executed by system composed of the results directory creation
			#and the tophat execution. The stats data will be redirected to the stats.log file
			$commanddef= "mkdir -p ".$projectdir."/".$output_dir." ;".$command." 2>> ".$statsfile;
		}
		#Opening the run.log and printing the execution data
		open (LOG,">> ".$logfile) || die "TOPHAT ERROR :: Can't open $logfile: $!";
		print LOG "TOPHAT :: ".date()." Executing $commanddef\n";
		close LOG;

		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		or die "TOPHAT ERROR :: system args failed: $? ($commanddef)";
		close STATS;
		
		#Renaming
		my $name=fileparse($file, qr{\..*$});
		$name=~s/_R?1//;
		$name=~s/_R?2//;
		if(-e "$projectdir/$output_dir/accepted_hits.bam"){
			my $output_file_mapped;
			my $output_file_unmapped;
			my $output_file_stats;
			my $output_file_ins;
			my $output_file_del;
			my $output_file_junct;
			
			#Renaming results
			if(lc($tophat_aligner) eq "bowtie1"){
				$output_file_mapped="$projectdir/$output_dir/$name\_top_bw1.bam";
				$output_file_unmapped="$projectdir/$output_dir/$name\_top_no_aligned.bam";
				$output_file_stats="$projectdir/$output_dir/$name\_top_bw1_align_summary.txt";
				$output_file_del="$projectdir/$output_dir/$name\_top_bw1_deletions.bed";
				$output_file_ins="$projectdir/$output_dir/$name\_top_bw1_insertions.bed";
				$output_file_junct="$projectdir/$output_dir/$name\_top_bw1_junctions.bed";
			}
			else{
				$output_file_mapped="$projectdir/$output_dir/$name\_top_bw2.bam";
				$output_file_unmapped="$projectdir/$output_dir/$name\_top_no_aligned.bam";
				$output_file_stats="$projectdir/$output_dir/$name\_top_bw2_align_summary.txt";
				$output_file_del="$projectdir/$output_dir/$name\_top_bw2_deletions.bed";
				$output_file_ins="$projectdir/$output_dir/$name\_top_bw2_insertions.bed";
				$output_file_junct="$projectdir/$output_dir/$name\_top_bw2_junctions.bed";
			}
			
			system("mv $projectdir/$output_dir/accepted_hits.bam $output_file_mapped");
			system("mv $projectdir/$output_dir/align_summary.txt $output_file_stats");
			system("mv $projectdir/$output_dir/deletions.bed $output_file_del");
			system("mv $projectdir/$output_dir/insertions.bed $output_file_ins");
			system("mv $projectdir/$output_dir/junctions.bed $output_file_junct");
			
			if($output_file_unmapped){
				system("mv $projectdir/$output_dir/unmapped.bam $output_file_unmapped");
			}
			#The path of the output file is returned to the main program
			return ($output_file_mapped);
		}
		else{
			print STDERR "ERROR :: Missing file : $projectdir/$output_dir/accepted_hits.bam\n";
			return();
		}
		return ("$projectdir/$output_dir/$name.bam");		
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "TOPHAT ERROR :: Can't open $logfile: $!";
    	print LOG "TOPHAT :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($bowtieindex) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn ("TOPHAT :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile),GTF ($GTF) and/or index($bowtieindex) have not been provided");
		help_tophat();
	}
	sub help_tophat{
	    my $usage = qq{
		  	$0 

			Needed parameters:
  	 		[file] Name of the file which is going to be align (fasta/fastq format)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[statsfile] Path of stats.log file where stats data will be saved
  	 		[bowtieindex]  Indexed genome to align your files in format .bt2
  	 		[projectdir] Directory where bowtie2_results directory will be created

  	 		Optional parameters:
  	 		[threads] Optional number of threads to perform the analysis
  	 		[bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 		[bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 		[bowtieparameters] Other bowtie parameters to perform the analysis using the bowtie recommended syntaxis
  	 		[verbose] Option to show the execution data on the screen   
						             
			Examples:
			bowtie2(file=>"file.fastq", threads=>"4", bowtieindex=>"./hg19", verbose=>"verbose", logfile=>"run.log", bowtiemiss=>"0", bowtieleng=>"19",statsfile=>"stats.log", bowtieparameters=>" -I 50 -X 200", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}	
}

=head2 bwa_index

  Example    : 
  bwa_index(
  	fasta=>"genome.fasta",
  	dir=>".",
  	logfile=>"run.log",
  	indexname=>"hg19"
  );
  Description: This function collects the directory and the path of genome sequence in fasta
  format to build a new index. This index can be used to perform a bwa analysis. The 
  execution data will be printed on run.log file. This fuction returns the path of the new 
  bwa index and saves the new index on the directory named Bowtie1_index at directory 
  provided by the user.
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where new index will be saved
  	 [fasta] Path of the genomic fasta sequence to build the index
  	 [logfile] Path of run.log file where execution data will be saved
  	 [indexname] Name to write in the index files
  Returntype : Directory named Bowtie1_index with index genome. bwa_index also returns the path of the new bwa index 
  Requeriments: bwa_index function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v1.0.0 or higher software correctly installed
  	- Input genomic sequence on fasta format  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub bwa_index{

	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/bwa/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/bwa/";
	}
	#First, check that bwa is in path:
	my @bwa_bin=`which bwa`;
	#Executing the command
	if(scalar(@bwa_bin)<1){
		die "BWA_INDEX ERROR :: system args failed: $? : Is bwa installed and exported to \$PATH ?";
	}
	my $fasta=$args{"fasta"}; #the path of genome sequence in fasta format
	my $dir=$args{"dir"}; #directory to create the genomeindex directory where new index will be saved
	my $logfile=$args{"logfile"}; #path of the logfile to write the execution data
	my $indexname=$args{"indexname"}; #name to write in the index files
	
	#Variable declaration
	my $index_output;
	my $command;
	my $commanddef;
	
	#Checking the mandatory arguments
	if ($fasta and $dir and $indexname and $logfile){
		print STDERR date()." Generating the index genome $indexname from $fasta. This process could take some hours\n";
		#bowtie-build execution command from a fasta file. The output index will be saved
		#in the genomeindex1 directory with the name index 
		my $command="bwa index -a bwtsw -p ".$dir."/BWA_index/$indexname ".$fasta;
		#commandef is the command that will be executed by system composed of the index
		#directory creation, the module loading, and the bowtie-build execution. The error 
		#data will be redirected to the run.log file
		$commanddef="mkdir -p ".$dir."/BWA_index ; ".$command ." >> ".$logfile." 2>&1";
		#Printing the date and command execution on the run.log file
		open (LOG,">> ".$logfile) || die "BWA_INDEX ERROR :: Can't open $logfile: $!";
		print LOG "BWA_INDEX :: ".date()." Executing $commanddef\n";
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		or die "BWA_INDEX ERROR :: system args failed: $? ($commanddef)";
		close LOG;
		#Returning the path of the new bwa index 
		return($dir."/BWA_index/".$indexname);
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "BWA_INDEX ERROR :: Can't open $logfile: $!";
    	print LOG "BWA_INDEX ERROR :: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile) and/or fasta file($fasta) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program will die and show error message
		warn ("BWA_INDEX ERROR:: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile) and/or fasta file($fasta) have not been provided");
		help_bwa_index();
	}	

	sub help_bwa_index{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[dir] Input directory where new index will be saved
  	 		[fasta] Path of the genomic fasta sequence to build the index
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[indexname] Name to write in the index files
						             
			Examples:
			bwa_index(fasta=>"genome.fasta", dir=>".", logfile=>"run.log", indexname=>"hg19");

	};

	print STDERR $usage;
	exit(); 
	}
}

sub bwa{

	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/bwa/".":$miARmaPath/bin/".$arch."/samtools/:";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/bwa/".":$miARmaPath/bin/".$arch."/samtools/:";
	}
	#First, check that bwa is in path:
	my @bwa_bin=`which bwa`;
	#Executing the command
	if(scalar(@bwa_bin)<1){
		die "BWA ::system args failed: $? : Is bwa installed and exported to \$PATH ?";
	}
	#First, check that bwa is in path:
	my @sam_bin=`which samtools`;
	#Executing the command
	if(scalar(@sam_bin)<1){
		die "BWA ::system args failed: $? : Is samtools installed and exported to \$PATH ?";
	}
	
	my $file=$args{"file"}; #Name of the file which is going to be aligned
	my $bwaindex=$args{"bwaindex"}; #Genome index in format .bwt 
	my $threads=$args{"threads"}; #Optional number of threads to perform the analysis  
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #Path of the statsfile to write the stats data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $Seqtype=$args{"Seqtype"}; #Sequencing method. SingleEnd by default. Acepted values : [Paired-End|Singe-End]
	my $bwapardef=$args{"bwaparameters"}; #additional parameter for bwa
	
	# Variable declaration and describing results directory 
	my $output_dir="/bwa_results/";

	#Variable to collect the optional parameters 
	$bwapardef.= " -T 19 ";
	#The number of missmatches can be provided by the user

	#Number of threads can be provided by the user
	if($threads>0){
		$bwapardef.= " -t $threads";
	}
	#Any other bowtie parameter can be provided by the user using the correct sintaxis
	#Checking the mandatory parameters
	open (LOG,">> ".$logfile) || die "BWA ERROR :: Can't open $logfile: $!";
	
	if ($file and $projectdir and $bwaindex and $logfile and $statsfile){ 
		#Extracting the name of the file
		my $name=fileparse($file, qr{\.f.*});
		my $command;
		my $stat_command;
		if(lc($Seqtype) eq "pairedend" or lc($Seqtype) eq "paired" or lc($Seqtype) eq "paired-end"){
			print STDOUT "\tBWA :: ".date()." Checking $file for bwa analysis (Paired End)\n" if($verbose);
			print LOG "BWA :: ".date()." Checking $file for bwa analysis (Paired End)\n";
			
			#Check if the file is a paired-end file
			if($file =~ /.*_R?1.*/){
				#it contains the _1 label
				my $mate_file=$file;
				if($file =~  /.*_R1.*/){
					$mate_file=~s/_R1(.+)/_R2$1/g;
				}
				else{
					$mate_file=~s/_1(.+)/_2$1/g;
				}
				if(-e $mate_file){
					if($file ne $mate_file){						
						my $real_name=$projectdir.$output_dir.$name."_bwa.sam";
						$real_name=~s/_\d.+bwa/_bwa/g;
						my $stat_file_bam=$real_name;
						$stat_file_bam=~s/\.sam/_align_summary.txt/g;
						if($file =~ /\.bz2$/){							
							$command="bwa mem ".$bwapardef." ".$bwaindex."  <(bunzip2 -c $file) <(bunzip2 -c $mate_file)>".$real_name;
							$stat_command="samtools flagstat $real_name > $stat_file_bam";
						}
						else{
							$command="bwa mem ".$bwapardef." ".$bwaindex." ".$file." $mate_file >".$real_name;
							$stat_command="samtools flagstat $real_name > $stat_file_bam";
							
						}
					}
				}
				else{
					print STDERR "ERROR:: You have requested a Paired-End analysis, so for the file $file a $mate_file file is needed\n";
					last;
				}
			}
			else{
				return();
			}
		}
		else{
			print STDOUT "\tBWA :: ".date()." Checking $file for bwa analysis (Single End)\n" if($verbose);
			print LOG "BWA :: ".date()." Checking $file for bwa analysis (Single End)\n";
			#bwa execution command
			$command="bwa mem ".$bwapardef." ".$bwaindex." ".$file." >".$projectdir.$output_dir.$name."_bwa.sam";
			my $stat_file_bam=$projectdir.$output_dir.$name."_bwa_align_summary.txt";
			$stat_command="samtools flagstat ".$projectdir.$output_dir.$name."_bwa.sam"." > $stat_file_bam";
			
		}

		my $commanddef;
		#bwa execution with verbose option
		if($verbose){
			#commandef is the command will be executed by system composed of the results directory creation 
			#and the bowtie execution. The stats data will be printed on the screen
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command;
			#Printing on the screen the date and the execution data
			print STDOUT "\tBWA :: ".date()." Executing $commanddef\n";
		}   	
		else
		{
			#bwa execution without verbose option, stats data will be printed on statsfile
			#and execution data on run.log file.
			#Opening stats file
			open (STATS,">> ".$statsfile) || die "BWA ERROR :: Can't open $statsfile: $!";
			print STATS "\tBWA :: ".date()." Executing $commanddef\n";
			#commandef is the command will be executed by system composed of the results directory creation 
			#and the bowtie execution. The stats data will be redirected to the stats.log file
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command. " 2>> ".$statsfile . "; " .$stat_command." 2>> ".$statsfile;
			
		}

		#Opening the run.log and printing the execution data
		print LOG "BWA :: ".date()." Executing $commanddef\n File:".$file."\n";
		close LOG;
			
		#Executing the command or if system can't be executed die showing the error.
		system_bash($commanddef) == 0
		or die "BWA ERROR :: system args failed: $? ($commanddef)";
		close STATS;
		#The path of the output file is returned to the main program
		return($projectdir.$output_dir.$name."_bwa.sam");
	}
	else
	{
		#Registering the error on run.log file
   		open(LOG,">> ".$logfile) || die "BWA ERROR :: Can't open $logfile: $!";
    	print LOG "BWA ERROR :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($bwaindex) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program will die and show error message
		warn ("BWA ERROR:: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($bwaindex) have not been provided");
		help_bwa();
	}	

	sub help_bwa{
	    my $usage = qq{
		  	$0 

			Needed parameters:
  	 		[file] Path of the file which is going to be aligned (fasta/fastq format)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[statsfile] Path of stats.log file where stats data will be saved
  	 		[bwaindex]  Indexed genome to align your files in format .bwt
  	 		[projectdir] Directory where results directory will be created

  	 		Optional parameters:
  	 		[threads] Optional number of threads to perform the analysis  
   	 		[verbose] Option to show the execution data on the screen   
						             
			Examples:
			bwa(file=>"./file.fastq", threads=>"4", bwaindex=>"./hg19.fa", verbose=>"verbose", logfile=>"run.log", statsfile=>"stats.log", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}	
}

=head2 miRDeep

  Example    : 
  miRDeep( 
	file=>"./file.fastq",
	threads=>"4",
	bowtieindex=>"./hg19",
	verbose=>"verbose", 
	logfile=>"run.log", 
	bowtiemiss=>"0", 
	bowtieleng=>"19",
	statsfile=>"stats.log", 
	bowtieparameters=>" -I 50 -X 200",
	projectdir=>"."
  )
  Description: miRDeep takes the reads and align them with the bowtie index to generate a sam file which 
  will be saved on miRDeep_results on the project directory. Execution and stats data will be saved on directory stats.log 
  file and will show on screen if verbose option is selected. In this case the bowtieindex 
  has to be built with miRDeep_index 
  function. 
  Input parameters: 
	Mandatory parameters:
  	 [file] Name of the file which is going to be align (fasta/fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 [bowtieindex]  Indexed genome to align your files in format .bt2
  	 [projectdir] Directory where  miRDeep_results directory will be created
  	 Optional parameters:
  	 [threads] Optional number of threads to perform the analysis
  	 [bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 [bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 [bowtieparameters] Other bowtie parameters to perform the analysis using the bowtie recommended syntaxis
  	 [verbose] Option to show the execution data on the screen   
  Returntype : File at directory miRDeep_results. Also returns the path of the output file
  Requeriments: miRDeep function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v2.2.0 or higher software correctly installed
  	- Input files on fastq format on the provided directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub miRDeep{
	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/common/mirdeep/:.:$miARmaPath/bin/".$arch."/bowtie1/:$miARmaPath/bin/".$arch."/Vienna/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/common/mirdeep/:.:$miARmaPath/bin/".$arch."/bowtie1/:$miARmaPath/bin/".$arch."/Vienna/";
	}
	
	#First, check that miRDeep is in path:
	my @miRDeep_bin=`which miRDeep2.pl`;
	#Executing the command
	if(scalar(@miRDeep_bin)<1){
		die "miRDeep ERROR ::system args failed: $? : Is miRDeep installed and exported to \$PATH ?";
	}

	my $file=$args{"file"}; #File which is going to be aligned
	my $bowtieindex=$args{"bowtie1index"}; #Genome index in format .bt2
	my $threads=$args{"threads"}; #Optional number of threads to perform the analysis faster
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #path of the statsfile to write the stats data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $adapter=$args{"adapter"}; #Sequencing method. SingleEnd by default. Acepted values : [Paired-End|Single-End]
	my $adapter_file=$args{"adapter_file"}; #Sequencing method. SingleEnd by default. Acepted values : [Paired-End|Single-End]
	my $mature_miRNA_file=$args{"mature"}; #a fasta file with all mature sequence from your organism
	my $precursor_miRNA_file=$args{"precursors"}; #a fasta file with all known pre-miRNa sequence 
	my $genome=$args{"genome"}; # fasta file for the complete genome of our organism
	my $organism=$args{"organism"}; # fasta file for the complete genome of our organism
	
	#Variable declaration and describing results directory 
	my $commanddef_mapper;
	my $commanddef_novo;
	my $output_dir="/miRDeep_results/";
	
	system("mkdir -p $projectdir");
	#Collecting miRDeep parameters provided by the user
	my $bowtiepardef= "";

	#Any other bowtie parameter can be provided by the user using the correct sintaxis
	if(defined($args{"miRDeeparameters"})){
		my $miRDeeparameters=$args{"miRDeeparameters"};
		$bowtiepardef.=" $miRDeeparameters";
	}
	
	#Checking the mandatory parameters
	if ($file and $projectdir and $bowtieindex and $logfile and $statsfile and $genome and $mature_miRNA_file and $precursor_miRNA_file){
		if($file =~ /\.gz$/){
			print STDOUT "\tmiRDeep :: ".date()." Uncompressing $file\n" if($verbose);
			#In case gzip
			system("gunzip -f $file");
			#Changing new extension
			$file=~s/\.gz$//g;
		}
		elsif($file =~ /\.bz2$/){
			#in case bzip2
			print STDOUT "\tmiRDeep :: ".date()." Uncompressing $file\n" if($verbose);
			system("bunzip2 -f -d $file");
			#New extension
			$file=~s/\.bz2$//g;			
		}
		 
		print STDOUT "\tmiRDeep :: ".date()." Checking $file for miRDeep analysis\n" if($verbose);
		#Extracting the name of the file
		my $name=fileparse($file, qr{\.f.*});
		
		#if no adapter is known, we rely on minion
		if(!defined $adapter){
			use CbBio::RNASeq::Adapt;
			$adapter=Minion(
				dir=>".",
				file=>$file,
				logfile=>$logfile,
				statsfile=>$statsfile,
				adaptpredictionnumber=>4,
				org=>$organism,
				verbose=>$verbose,
				miARmaPath=>$miARmaPath,
			);				
		}
		#miRDeep execution command
		
		my $command_mapper;
		if(defined($adapter) and lc($adapter) ne "no"){
			$command_mapper="export PERL5LIB=$miARmaPath/lib/Perl/; mapper.pl ".$file." -e -h -i -j -n -m -k ".$adapter ." -o ". $threads ." -p " . $bowtieindex." -s ". $projectdir.$output_dir.$name.".fa -t ".$projectdir.$output_dir.$name."_vs_genome.arf";
		}
		#if no adaptder is provided or found by minion, try without adapter
		elsif(defined($adapter) and lc($adapter) eq "no"){
			$command_mapper="export PERL5LIB=$miARmaPath/lib/Perl/; mapper.pl ".$file." -e -h -i -j -n -m -o ". $threads ." -p " . $bowtieindex." -s ". $projectdir.$output_dir.$name.".fa -t ".$projectdir.$output_dir.$name."_vs_genome.arf";
		}
		else{
			$command_mapper="export PERL5LIB=$miARmaPath/lib/Perl/; mapper.pl ".$file." -e -h -i -j -n -m -o ". $threads ." -p " . $bowtieindex." -s ". $projectdir.$output_dir.$name.".fa -t ".$projectdir.$output_dir.$name."_vs_genome.arf";
		}
		my $command_novo=" export PERL5LIB=$miARmaPath/lib/Perl/;miRDeep2.pl ".$projectdir.$output_dir.$name.".fa ".$genome. " ". $projectdir.$output_dir.$name."_vs_genome.arf $mature_miRNA_file none $precursor_miRNA_file -r $name -P -d -c -v $bowtiepardef ";
		#Bowtie execution with verbose option
		if($verbose){
			#commandef is the command will be executed by system composed of the results directory 
			#creation and the miRDeep execution. The stats data will be printed on the screen
			$commanddef_mapper= "mkdir -p ".$projectdir.$output_dir." ;".$command_mapper;
			$commanddef_novo= "mkdir -p ".$projectdir.$output_dir." ;".$command_novo;
			#Printing on the screen the date and the execution data
			print STDOUT "\tmiRDeep :: ".date()." Executing $command_mapper\n$command_novo";
		}   
		#Bowtie execution printing stats on stats file	
		else{
			#Opening and printing stats on statsfile
			system("touch $statsfile");
			open (STATS,">> ".$statsfile) || die "miRDeep ERROR :: Can't open $statsfile: $!";
			print STATS "miRDeep :: File:".$file."\n";
			#commandef is the command will be executed by system composed of the results directory creation
			#and the bowtie execution. The stats data will be redirected to the stats.log file
			$commanddef_mapper= "mkdir -p ".$projectdir.$output_dir." ;".$command_mapper." >> ".$statsfile ." 2>&1";
			$commanddef_novo= "mkdir -p ".$projectdir.$output_dir." ;".$command_novo." >> ".$statsfile ." 2>&1";
			
		}
		#Opening the run.log and printing the execution data
		open (LOG,">> ".$logfile) || die "miRDeep ERROR :: Can't open $logfile: $!";
		print LOG "miRDeep :: ".date()." Executing $commanddef_mapper\n";
		print LOG "miRDeep :: ".date()." Executing $commanddef_novo\n";

		print STDOUT "\tmiRDeep :: ".date()." Mapping $file with miRDeep mapper module\n" if($verbose);
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef_mapper) == 0
		 or die "miRDeep ERROR :: system args failed: $? ($commanddef_mapper)";

		print STDOUT "\tmiRDeep :: ".date()." Identifying miRNAs (including novel) in $file with miRDeep module\n" if($verbose);
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef_novo) == 0
		 or die "miRDeep ERROR :: system args failed: $? ($commanddef_novo)";

		#Clean all temporary files by miRDeep
		my $clean_cmd= "mv result_*.csv $projectdir/$output_dir/$name\.xls && rm -rf $miARmaPath/dir_* $miARmaPath/*.bed $miARmaPath/*.html $miARmaPath/mirna_results_* $miARmaPath/expression_analyses $miARmaPath/mapper.log $miARmaPath/bowtie.log $miARmaPath/error_* $miARmaPath/mirdeep_runs";
		print LOG "miRDeep :: ".date()." $clean_cmd\n";
		system($clean_cmd) == 0
		 or die "miRDeep ERROR :: system args failed: $? ($clean_cmd)";
		
		close STATS;
		close LOG;
		
		return();
		
		#The path of the output file is returned to the main program
		#my $miRDeep_counts=parse_miRDeep(
		#	file=>$name,
		#	path=>$projectdir.$output_dir,
		#);
		#return ($miRDeep_counts);		
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "miRDeep ERROR :: Can't open $logfile: $!";
    	print LOG "miRDeep :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($bowtieindex) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn ("miRDeep :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($bowtieindex) have not been provided");
		help_miRDeep();
	}
	sub help_miRDeep{
	    my $usage = qq{
		  	$0 

			Needed parameters:
  	 		[file] Name of the file which is going to be align (fasta/fastq format)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[statsfile] Path of stats.log file where stats data will be saved
  	 		[bowtieindex]  Indexed genome to align your files in format .bt2
  	 		[projectdir] Directory where miRDeep_results directory will be created

  	 		Optional parameters:
  	 		[threads] Optional number of threads to perform the analysis
  	 		[bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 		[bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 		[bowtieparameters] Other bowtie parameters to perform the analysis using the bowtie recommended syntaxis
  	 		[verbose] Option to show the execution data on the screen   
						             
			Examples:
			miRDeep(file=>"file.fastq", threads=>"4", bowtieindex=>"./hg19", verbose=>"verbose", logfile=>"run.log", bowtiemiss=>"0", bowtieleng=>"19",statsfile=>"stats.log", bowtieparameters=>" -I 50 -X 200", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}	
}

=head2 hisat2_index

  Example    : 
  hisat2_index(
  	fasta=>"genome.fasta",
  	dir=>".",
  	logfile=>"run.log",
  	indexname=>"hg19"
  );
  Description: This function collects the directory and the path of genome sequence in fasta
  format to build a new index. This index can be used to perform a hisat2 analysis. The 
  execution data will be printed on run.log file. This fuction returns the path of the new 
  hisat2 index and saves the new index on the directory named Bowtie1_index at directory 
  provided by the user.
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where new index will be saved
  	 [fasta] Path of the genomic fasta sequence to build the index
  	 [logfile] Path of run.log file where execution data will be saved
  	 [indexname] Name to write in the index files
  Returntype : Directory named Bowtie1_index with index genome. hisat2_index also returns the path of the new hisat2 index 
  Requeriments: hisat2_index function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v1.0.0 or higher software correctly installed
  	- Input genomic sequence on fasta format  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub hisat2_index{

	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $threads=$args{"threads"};
	
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/hisat2/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/hisat2/";
	}
	#First, check that hisat2 is in path:
	my @hisat2_bin=`which hisat2`;
	#Executing the command
	if(scalar(@hisat2_bin)<1){
		die "HISAT2_INDEX ERROR :: system args failed: $? : Is hisat2 installed and exported to \$PATH ?";
	}
	my $fasta=$args{"fasta"}; #the path of genome sequence in fasta format
	my $dir=$args{"dir"}; #directory to create the genomeindex directory where new index will be saved
	my $logfile=$args{"logfile"}; #path of the logfile to write the execution data
	my $indexname=$args{"indexname"}; #name to write in the index files
	
	#Variable declaration
	my $index_output;
	my $command;
	my $commanddef;
	
	#Checking the mandatory arguments
	if ($fasta and $dir and $indexname and $logfile){
		print STDERR date()." Generating the index genome $indexname from $fasta. This process could take some hours\n";
		#bowtie-build execution command from a fasta file. The output index will be saved
		#in the genomeindex1 directory with the name index 
		my $command="hisat2-build -p $threads -f ".$fasta." ".$dir."/hisat2_index/".$indexname;
		#commandef is the command that will be executed by system composed of the index
		#directory creation, the module loading, and the bowtie-build execution. The error 
		#data will be redirected to the run.log file
		$commanddef="mkdir -p ".$dir."/hisat2_index ; ".$command ." >> ".$logfile." 2>&1";
		#Printing the date and command execution on the run.log file
		open (LOG,">> ".$logfile) || die "HISAT2_INDEX ERROR :: Can't open $logfile: $!";
		print LOG "HISAT2_INDEX :: ".date()." Executing $commanddef\n";
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		or die "HISAT2_INDEX ERROR :: system args failed: $? ($commanddef)";
		close LOG;
		#Returning the path of the new hisat2 index 
		return($dir."/hisat2_index/".$indexname);
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "HISAT2_INDEX ERROR :: Can't open $logfile: $!";
    	print LOG "HISAT2_INDEX ERROR :: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile) and/or fasta file($fasta) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program will die and show error message
		warn ("HISAT2_INDEX ERROR:: ".date()." Directory ($dir), indexname ($indexname), logfile($logfile) and/or fasta file($fasta) have not been provided");
		help_hisat2_index();
	}	

	sub help_hisat2_index{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[dir] Input directory where new index will be saved
  	 		[fasta] Path of the genomic fasta sequence to build the index
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[indexname] Name to write in the index files
			[threads] Number of threads
						             
			Examples:
			hisat2_index(fasta=>"genome.fasta", dir=>".", logfile=>"run.log", indexname=>"hg19");

	};

	print STDERR $usage;
	exit(); 
	}
}


=head2 hisat2

  Example    : 
  hisat2( 
	file=>"./file.fastq",
	threads=>"4",
	bowtieindex=>"./hg19",
	verbose=>"verbose", 
	logfile=>"run.log", 
	bowtiemiss=>"0", 
	bowtieleng=>"19",
	statsfile=>"stats.log", 
	bowtieparameters=>" -I 50 -X 200",
	projectdir=>"."
  )
  Description: Bowtie2 takes the reads and align them with the bowtie index to generate a sam file which 
  will be saved on Bowtie2_results on the project directory. Execution and stats data will be saved on directory stats.log 
  file and will show on screen if verbose option is selected. In this case the bowtieindex 
  has to be built with hisat2_index 
  function. 
  Input parameters: 
	Mandatory parameters:
  	 [file] Name of the file which is going to be align (fasta/fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 [bowtieindex]  Indexed genome to align your files in format .bt2
  	 [projectdir] Directory where hisat2_results directory will be created
  	 Optional parameters:
  	 [threads] Optional number of threads to perform the analysis
  	 [bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 [bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 [bowtieparameters] Other bowtie parameters to perform the analysis using the bowtie recommended syntaxis
  	 [verbose] Option to show the execution data on the screen   
  Returntype : File at directory Bowtie2_results. Also returns the path of the output file
  Requeriments: hisat2 function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v2.2.0 or higher software correctly installed
  	- Input files on fastq format on the provided directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub hisat2{
	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/hisat2/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/hisat2/";
	}
	
	#First, check that hisat2 is in path:
	my @hisat2_bin=`which hisat2`;
	#Executing the command
	if(scalar(@hisat2_bin)<1){
		die "HISAT2 ERROR ::system args failed: $? : Is hisat2 installed and exported to \$PATH (".$ENV{PATH}.") ?";
	}

	my $file=$args{"file"}; #File which is going to be aligned
	my $hisat2index=$args{"hisat2index"}; #Genome index in format .bt2
	my $threads=$args{"threads"}; #Optional number of threads to perform the analysis faster
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #path of the statsfile to write the stats data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $Seqtype=$args{"Seqtype"}; #Sequencing method. SingleEnd by default. Acepted values : [Paired-End|Single-End]
	my $adapter=$args{"adapter"}; #are this reads processed by miArma 
	my $strand=$args{"strand"}; #are this reads processed by miArma 
	
	#Variable declaration and describing results directory 
	my $commanddef;
	my $output_dir="/hisat2_results/";
	my $hisatpardef;
	#Collecting hisat2 parameters provided by the user
	if(defined($args{"hisat2parameters"})){
		my $hisatparameters=$args{"hisat2parameters"};
		$hisatpardef.=" hisat2parameters";
	}
	#Number of threads can be provided by the user
	if($threads>0){
		$hisatpardef.= " -p $threads";
	}

	#Checking the mandatory parameters
	if ($file and $projectdir and $hisat2index and $logfile and $statsfile){ 
		#Extracting the name of the file
		my $name=fileparse($file, qr{\.f.*});
		my $output_file_bw2;
		if($adapter){
			$output_file_bw2=$projectdir.$output_dir.$name."_his";
			$output_file_bw2=~s/_R1(.+)/$1/g;
			$output_file_bw2=~s/_1(.+)/$1/g;
		}
		else{
			$output_file_bw2=$projectdir.$output_dir.$name."_nat_his";
			$output_file_bw2=~s/_R1(.+)/$1/g;
			$output_file_bw2=~s/_1(.+)/$1/g;
		}
		#Bowtie2 execution command
		my $command;
		if(lc($Seqtype) eq "pairedend" or lc($Seqtype) eq "paired" or lc($Seqtype) eq "paired-end"){
			print STDOUT "\tHISAT 2 :: ".date()." Checking $file for hisat2 (paired-end) analysis\n" if($verbose);
			
			if(lc($strand) eq "reverse"){
				$hisatpardef.= " --rna-strandness RF";
			}
			elsif(lc($strand) eq "yes"){
				$hisatpardef.= " --rna-strandness FR";
			}

			#Check if the file is a paired-end file
			if($file =~ /.*_R?1.*/){
				#it contains the _1 label
				my $mate_file=$file;
				if($file =~  /.*_R1.*/){
					$mate_file=~s/_R1(.+)/_R2$1/g;
				}
				else{
					$mate_file=~s/_1(.+)/_2$1/g;
				}
				if(-e $mate_file){
					if($file ne $mate_file){
						$command="hisat2 -q ".$hisatpardef." -x ".$hisat2index." -1 ".$file." -2 ". $mate_file ." --met-file ".$output_file_bw2.".metrics --un-conc-gz ".$output_file_bw2."_no_aligned.fastq.gz -S ". $output_file_bw2 . ".sam";
					}
				}
				else{
					print STDERR "ERROR:: You have requested a Paired-End analysis, so for the file $file a $mate_file file is needed\n";
					last;
				}
			}
			else{
				return();
			}
		}
		else{
			if(lc($strand) eq "reverse"){
				$hisatpardef.= " --rna-strandness R";
			}
			elsif(lc($strand) eq "yes"){
				$hisatpardef.= " --rna-strandness F";
			}
			
			print STDOUT "\tHISAT 2 :: ".date()." Checking $file for hisat2 (single-end) analysis\n" if($verbose);
			$command="hisat2 -q ".$hisatpardef." -x ".$hisat2index." -U ".$file." --met-file ". $output_file_bw2.".metrics --un-gz ".$output_file_bw2."_no_aligned.fastq.gz -S ". $output_file_bw2.".sam";
		}
		
		#Bowtie execution with verbose option
		if($verbose){
			#commandef is the command will be executed by system composed of the results directory 
			#creation and the hisat2 execution. The stats data will be printed on the screen
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command;
			#Printing on the screen the date and the execution data
			print STDOUT "\tHISAT2 :: ".date()." Executing $commanddef\n";
		}   
		#Bowtie execution printing stats on stats file	
		else{
			#Opening and printing stats on statsfile
			open (STATS,">> ".$statsfile) || die "HISAT2 ERROR :: Can't open $statsfile: $!";
			print STATS "HISAT2 :: File:".$file."\n";
			#commandef is the command will be executed by system composed of the results directory creation
			#and the hisat execution. The stats data will be redirected to the stats.log file
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command." 2>> ".$statsfile;
		}
		#Opening the run.log and printing the execution data
		open (LOG,">> ".$logfile) || die "HISAT2 ERROR :: Can't open $logfile: $!";
		print LOG "HISAT2 :: ".date()." Executing $commanddef\n";
		close LOG;

		#print STDERR "$command\n";
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		 or die "HISAT2 ERROR :: system args failed: $? ($commanddef)";
		close STATS;
		#The path of the output file is returned to the main program
		return ($output_file_bw2);		
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "HISAT2 ERROR :: Can't open $logfile: $!";
    	print LOG "HISAT2 :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($hisat2index) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn ("HISAT2 :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($hisat2index) have not been provided");
		help_hisat2();
	}
	sub help_hisat2{
	    my $usage = qq{
		  	$0 

			Needed parameters:
  	 		[file] Name of the file which is going to be align (fasta/fastq format)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[statsfile] Path of stats.log file where stats data will be saved
  	 		[hisatindex]  Indexed genome to align your files in format .bt2
  	 		[projectdir] Directory where hisat2_results directory will be created

  	 		Optional parameters:
  	 		[threads] Optional number of threads to perform the analysis
  	 		[hisatmiss] Max # mismatches in seed alignment in hisat analysis (0-1)
  	 		[hisatlength] Length of seed substrings in hisat analysis (>3, <32)
  	 		[hisatparameters] Other hisat parameters to perform the analysis using the hisat recommended syntaxis
  	 		[verbose] Option to show the execution data on the screen   
						             
			Examples:
			hisat2(file=>"file.fastq", threads=>"4", hisatindex=>"./hg19", verbose=>"verbose", logfile=>"run.log", hisatmiss=>"0", hisatleng=>"19",statsfile=>"stats.log", hisatparameters=>" -I 50 -X 200", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}	
}

=head2 star_index

  Example    : 
  star_index(
  	fasta=>"genome.fasta",
  	dir=>".",
  	logfile=>"run.log",
  	indexname=>"hg19"
  );
  Description: This function collects the directory and the path of genome sequence in fasta
  format to build a new index. This index can be used to perform a star analysis. The 
  execution data will be printed on run.log file. This fuction returns the path of the new 
  star index and saves the new index on the directory named Bowtie1_index at directory 
  provided by the user.
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where new index will be saved
  	 [fasta] Path of the genomic fasta sequence to build the index
  	 [logfile] Path of run.log file where execution data will be saved
  	 [indexname] Name to write in the index files
  Returntype : Directory named Bowtie1_index with index genome. star_index also returns the path of the new star index 
  Requeriments: star_index function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v1.0.0 or higher software correctly installed
  	- Input genomic sequence on fasta format  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub star_index{

	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $threads=$args{"threads"};
	
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/STAR/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/STAR/";
	}
	#First, check that star is in path:
	my @star_bin=`which STAR`;
	#Executing the command
	if(scalar(@star_bin)<1){
		die "STAR_INDEX ERROR :: system args failed: $? : Is star installed and exported to \$PATH ?";
	}
	my $fasta=$args{"fasta"}; #the path of genome sequence in fasta format
	my $dir=$args{"dir"}; #directory to create the genomeindex directory where new index will be saved
	my $logfile=$args{"logfile"}; #path of the logfile to write the execution data
	my $gtf=$args{"gtf"}; #name to write in the index files
	
	#Variable declaration
	my $index_output;
	my $command;
	my $commanddef;
	
	#Checking the mandatory arguments
	if ($fasta and $dir and $gtf and $logfile){
		print STDERR date()." Generating the index genome from $fasta. This process could take some hours and it could need >=32 GB of RAM memory\n";
		#bowtie-build execution command from a fasta file. The output index will be saved
		#in the genomeindex1 directory with the name index 
		my $command="STAR --runThreadN $threads --runMode genomeGenerate --genomeDir $dir/star_index/ --genomeFastaFiles $fasta --sjdbGTFfile $gtf";
		#commandef is the command that will be executed by system composed of the index
		#directory creation, the module loading, and the bowtie-build execution. The error 
		#data will be redirected to the run.log file
		$commanddef="mkdir -p ".$dir."/star_index ; ".$command ." >> ".$logfile." 2>&1";
		#Printing the date and command execution on the run.log file
		open (LOG,">> ".$logfile) || die "STAR_INDEX ERROR :: Can't open $logfile: $!";
		print LOG "STAR_INDEX :: ".date()." Executing $commanddef\n";
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		or die "STAR_INDEX ERROR :: system args failed: $? ($commanddef)";
		close LOG;
		#Returning the path of the new star index 
		return($dir."/star_index/");
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "STAR_INDEX ERROR :: Can't open $logfile: $!";
    	print LOG "STAR_INDEX ERROR :: ".date()." Directory ($dir), GTF ($gtf), logfile($logfile) and/or fasta file($fasta) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program will die and show error message
		warn ("STAR_INDEX ERROR:: ".date()." Directory ($dir), GTF ($gtf), logfile($logfile) and/or fasta file($fasta) have not been provided");
		help_star_index();
	}	

	sub help_star_index{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[dir] Input directory where new index will be saved
  	 		[fasta] Path of the genomic fasta sequence to build the index
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[indexname] Name to write in the index files
			[threads] Number of threads
						             
			Examples:
			star_index(fasta=>"genome.fasta", dir=>".", logfile=>"run.log", indexname=>"hg19");

	};

	print STDERR $usage;
	exit(); 
	}
}


=head2 star

  Example    : 
  star( 
	file=>"./file.fastq",
	threads=>"4",
	bowtieindex=>"./hg19",
	verbose=>"verbose", 
	logfile=>"run.log", 
	bowtiemiss=>"0", 
	bowtieleng=>"19",
	statsfile=>"stats.log", 
	bowtieparameters=>" -I 50 -X 200",
	projectdir=>"."
  )
  Description: Bowtie2 takes the reads and align them with the bowtie index to generate a sam file which 
  will be saved on Bowtie2_results on the project directory. Execution and stats data will be saved on directory stats.log 
  file and will show on screen if verbose option is selected. In this case the bowtieindex 
  has to be built with star_index 
  function. 
  Input parameters: 
	Mandatory parameters:
  	 [file] Name of the file which is going to be align (fasta/fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [statsfile] Path of stats.log file where stats data will be saved
  	 [bowtieindex]  Indexed genome to align your files in format .bt2
  	 [projectdir] Directory where star_results directory will be created
  	 Optional parameters:
  	 [threads] Optional number of threads to perform the analysis
  	 [bowtiemiss] Max # mismatches in seed alignment in bowtie analysis (0-1)
  	 [bowtielength] Length of seed substrings in bowtie analysis (>3, <32)
  	 [bowtieparameters] Other bowtie parameters to perform the analysis using the bowtie recommended syntaxis
  	 [verbose] Option to show the execution data on the screen   
  Returntype : File at directory Bowtie2_results. Also returns the path of the output file
  Requeriments: star function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Bowtie v2.2.0 or higher software correctly installed
  	- Input files on fastq format on the provided directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub star{
	#Arguments provided by user are collected by %args. Dir, path of fasta file, indexname and logfile
	# are mandatory arguments.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/STAR/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/STAR/";
	}
	
	#First, check that star is in path:
	my @star_bin=`which STAR`;
	#Executing the command
	if(scalar(@star_bin)<1){
		die "STAR ERROR ::system args failed: $? : Is star installed and exported to \$PATH (".$ENV{PATH}.") ?";
	}

	my $file=$args{"file"}; #File which is going to be aligned
	my $starindex=$args{"starindex"}; #Genome index in format .bt2
	my $threads=$args{"threads"}; #Optional number of threads to perform the analysis faster
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #path of the statsfile to write the stats data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $Seqtype=$args{"Seqtype"}; #Sequencing method. SingleEnd by default. Acepted values : [Paired-End|Single-End]
	my $adapter=$args{"adapter"}; #are this reads processed by miArma 
	
	#Variable declaration and describing results directory 
	my $commanddef;
	my $output_dir="/star_results/";
	my $starpardef;
	#Collecting star parameters provided by the user
	if(defined($args{"starparameters"})){
		my $hisatparameters=$args{"starparameters"};
		$starpardef.=" starparameters";
	}
	#Number of threads can be provided by the user
	if($threads>0){
		$starpardef.= " --runThreadN $threads";
	}
	
	#Checking the mandatory parameters
	if ($file and $projectdir and $starindex and $logfile and $statsfile){ 
		#Extracting the name of the file
		my $name=fileparse($file, qr{\.f.*});
		my $output_file_bw2;
		if($adapter){
			$output_file_bw2=$projectdir.$output_dir.$name."_str_";
			$output_file_bw2=~s/_R1(.+)/$1/g;
			$output_file_bw2=~s/_1(.+)/$1/g;
		}
		else{
			$output_file_bw2=$projectdir.$output_dir.$name."_nat_str_";
			$output_file_bw2=~s/_R1(.+)/$1/g;
			$output_file_bw2=~s/_1(.+)/$1/g;
		}
		#Bowtie2 execution command
		my $command;
		if(lc($Seqtype) eq "pairedend" or lc($Seqtype) eq "paired" or lc($Seqtype) eq "paired-end"){
			print STDOUT "\tSTAR :: ".date()." Checking $file for star (paired-end) analysis\n" if($verbose);
			
			#Check if the file is a paired-end file
			if($file =~ /.*_R?1.*/){
				#it contains the _1 label
				my $mate_file=$file;
				if($file =~  /.*_R1.*/){
					$mate_file=~s/_R1(.+)/_R2$1/g;
				}
				else{
					$mate_file=~s/_1(.+)/_2$1/g;
				}
				if(-e $mate_file){
					if($file ne $mate_file){
						if($file =~ /\.gz$/){
							$starpardef.= " --readFilesCommand gunzip -c ";
						}
						elsif($file =~ /\.bz2$/){
							$starpardef.= " --readFilesCommand bunzip2 -c ";
						}
						else{
							#none
						}
						#$command="star -q ".$starpardef." -x ".$starindex." -1 ".$file." -2 ". $mate_file ." --met-file ".$projectdir.$output_dir.$name.".metrics --un ".$projectdir.$output_dir.$name."_no_aligned.fastq -S ". $output_file_bw2;
						$command="STAR --runMode alignReads --genomeDir $starindex --readFilesIn $file $mate_file  --outFileNamePrefix $output_file_bw2 $starpardef --outSAMtype BAM SortedByCoordinate"
					}
				}
				else{
					print STDERR "ERROR:: You have requested a Paired-End analysis, so for the file $file a $mate_file file is needed\n";
					last;
				}
			}
			else{
				return();
			}
		}
		else{
			print STDOUT "\tSTAR :: ".date()." Checking $file for star (single-end) analysis\n" if($verbose);
			if($file =~ /\.gz$/){
				$starpardef.= " --readFilesCommand gunzip -c ";
			}
			elsif($file =~ /\.bz2$/){
				$starpardef.= " --readFilesCommand bunzip2 -c ";
			}
			else{
				#none
			}
			#$command="star -q ".$starpardef." -x ".$starindex." -U ".$file." --met-file ".$projectdir.$output_dir.$name.".metrics --un ".$projectdir.$output_dir.$name."_no_aligned.fastq -S ". $output_file_bw2;
			$command="STAR --runMode alignReads --genomeDir $starindex --readFilesIn $file  --outFileNamePrefix $output_file_bw2 $starpardef --outSAMtype BAM SortedByCoordinate"
			
		}
		
		#Bowtie execution with verbose option
		if($verbose){
			#commandef is the command will be executed by system composed of the results directory 
			#creation and the star execution. The stats data will be printed on the screen
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command;
			#Printing on the screen the date and the execution data
			print STDOUT "\tSTAR :: ".date()." Executing $commanddef\n";
		}   
		#Bowtie execution printing stats on stats file	
		else{
			#Opening and printing stats on statsfile
			open (STATS,">> ".$statsfile) || die "STAR ERROR :: Can't open $statsfile: $!";
			print STATS "STAR :: File:".$file."\n";
			#commandef is the command will be executed by system composed of the results directory creation
			#and the hisat execution. The stats data will be redirected to the stats.log file
			$commanddef= "mkdir -p ".$projectdir.$output_dir." ;".$command." >> ".$logfile." 2>&1";
		}
		#Opening the run.log and printing the execution data
		open (LOG,">> ".$logfile) || die "STAR ERROR :: Can't open $logfile: $!";
		print LOG "STAR :: ".date()." Executing $commanddef\n";
		close LOG;

		#print STDERR "$command\n";
		#Executing the command or if system can't be executed die showing the error.
		system($commanddef) == 0
		 or die "STAR ERROR :: system args failed: $? ($commanddef)";
		close STATS;
		#The path of the output file is returned to the main program
		return ($output_file_bw2);		
	}
	else
	{
		#Registering the error
   		open(LOG,">> ".$logfile) || die "STAR ERROR :: Can't open $logfile: $!";
    	print LOG "STAR :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($starindex) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn ("STAR :: ".date()." File($file), logfile($logfile), projectdir ($projectdir), statsfile($statsfile) and/or index($starindex) have not been provided");
		help_star();
	}
	sub help_star{
	    my $usage = qq{
		  	$0 

			Needed parameters:
  	 		[file] Name of the file which is going to be align (fasta/fastq format)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[statsfile] Path of stats.log file where stats data will be saved
  	 		[hisatindex]  Indexed genome to align your files in format .bt2
  	 		[projectdir] Directory where star_results directory will be created

  	 		Optional parameters:
  	 		[threads] Optional number of threads to perform the analysis
  	 		[hisatmiss] Max # mismatches in seed alignment in hisat analysis (0-1)
  	 		[hisatlength] Length of seed substrings in hisat analysis (>3, <32)
  	 		[hisatparameters] Other hisat parameters to perform the analysis using the hisat recommended syntaxis
  	 		[verbose] Option to show the execution data on the screen   
						             
			Examples:
			star(file=>"file.fastq", threads=>"4", hisatindex=>"./hg19", verbose=>"verbose", logfile=>"run.log", hisatmiss=>"0", hisatleng=>"19",statsfile=>"stats.log", hisatparameters=>" -I 50 -X 200", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}	
}


sub ReadSummary{
	use File::Basename;
	#Arguments provided by user are collected by %args. Dir, file, aligner, statsfile, projectdir 
	#and logfile are mandatory arguments while verbose and threads are optional.
	my %args=@_;
	my $input=$args{"files"}; #Input is an referenced array with the paths of the files which have to be joined
	my $summary_file=$args{"summary"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}; #Optional arguments to show the execution data on screen
	my $statsfile=$args{"statsfile"}; #Path of the statsfile to write the stats data
	my $aligner=$args{"aligner"}; #Input directory where results directory will be created
	my $tophat_aligner=$args{"tophat_aligner"};
	
	my $summary;
	my $summary_path;
	my $summary_bw2;
	
	if(lc($aligner) eq "bowtie1" ){
		$summary_path->{$aligner}=$projectdir ."/Bowtie1_results/";
		open(STAT,$statsfile);
		my $real_file;
		my $processed;
		my $aligned;
		my $failed;
		while(<STAT>){
			chomp;
			if($_ =~/^BOWTIE1/){
				my $file=$_;
				$file=~s/BOWTIE1 :: File://g;
				$real_file=fileparse($file);
				$processed=0;
				$aligned=0;
				$failed=0;
			}
			if($_ =~/^# reads processed:/){
				$processed=$_;
				$processed=~s/.*: //g;
			}
			if($_ =~/^# reads with at/){
				$aligned=$_;
				$aligned=~s/.*: //g;
			}
			if($_ =~/^# reads that failed/){
				$failed=$_;
				$failed=~s/.*: //g;
			}
			if($real_file and $failed){
				$summary->{$real_file}="$processed\t$aligned\t$failed";
			}
		}
	}
	if(lc($aligner) eq "bowtie2" ){
		$summary_path->{$aligner}=$projectdir ."/Bowtie2_results/";
		open(STAT,$statsfile);
		my $real_file;
		my $processed;
		my $aligned;
		my $failed;
		my $multimapping;
		my $overall;
		while(<STAT>){
			chomp;
			if($_ =~/^BOWTIE2/){
				my $file=$_;
				$file=~s/BOWTIE2 :: File://g;
				$real_file=fileparse($file);
				$processed=0;
				$aligned=0;
				$failed=0;
				$multimapping=0;
				$overall=0;
			}
			if($_ =~/reads; of these:/){
				$processed=$_;
				$processed=~s/^(\d+) reads; of these:/$1/g;
			}
			if($_ =~/ aligned 0 times/){
				$failed=$_;
				$failed=~s/\s+(\d+) .*/$1/g;
			}
			if($_ =~/aligned exactly 1 time/){
				$aligned=$_;
				$aligned=~s/\s+(\d+) .*/$1/g;
			}
			if($_ =~/aligned >1 times/){
				$multimapping=$_;
				$multimapping=~s/\s+(\d+) .*/$1/g;
			}
			if($_ =~/overall alignment rate/){
				$overall=$_;
				$overall=~s/^(\d+\.\d+)% .*/$1%/g;
			}
			if($real_file and $failed){
				$summary_bw2->{$real_file}="$processed\t$aligned\t$multimapping\t$overall\t$failed";
			}
		}
	}
	if(lc($aligner) eq "bowtie2-bowtie1" or lc($aligner) eq "bowtie1-bowtie2"){
		ReadSummary(
		    aligner=>"bowtie1",
			summary=>$summary_file,
			statsfile=>$statsfile,
			projectdir=>$projectdir,
		);
		ReadSummary(
		    aligner=>"bowtie2",
			summary=>$summary_file,
			statsfile=>$statsfile,
			projectdir=>$projectdir,
		);
	}
	if(lc($aligner) eq "tophat"){
		my $dirname;
		if(lc($tophat_aligner) eq "bowtie2"){
			#test if tophat was used with bw2
			$dirname = $projectdir ."/Bowtie2_results/";
		}
		if(lc($tophat_aligner) eq "bowtie1"){
			#test if tophat was used with bw2
			$dirname = $projectdir ."/Bowtie1_results/";
		}
		if(lc($tophat_aligner) eq "bowtie2-bowtie1" or lc($tophat_aligner) eq "bowtie1-bowtie2"){
			ReadSummary(
			    aligner=>"tophat",
				tophat_aligner=>"bowtie1",
				summary=>$summary_file,
				statsfile=>$statsfile,
				projectdir=>$projectdir,
			);
			ReadSummary(
		    	aligner=>"tophat",
				tophat_aligner=>"bowtie2",
				summary=>$summary_file,
				statsfile=>$statsfile,
				projectdir=>$projectdir,
			);
		}
		
		if($dirname){
			$summary_path->{$aligner}=$dirname;
		
			opendir ( DIR, $dirname ) || die "Error in opening dir $dirname\n";
			my @files= readdir DIR;
			foreach my $filename (sort @files){		     
				if($filename =~ /_align_summary.txt$/){
					open(RESULTS, $dirname ."/". $filename) || warn $!;	
					my $real_file=fileparse($filename);
					$real_file=~s/_align_summary.txt$//g;
					my $processed;
					my $aligned;
					my $failed;
					my $multimapping;
					my $overall;
					while(<RESULTS>){
						chomp;
						if($_ =~/Input/){
							$processed=$_;
							$processed=~s/\s+Input\s+:\s+(\d+)/$1/g;
						}
						if($_ =~/Mapped/){
							$aligned=$_;
							$aligned=~s/\s+Mapped\s+:\s+(\d+) .+/$1/g;
						}
						if($_ =~/of these/){
							$multimapping=$_;
							$multimapping=~s/\s+of these:\s+(.+)/$1/g;
							$multimapping=~s/have multiple alignments.+//g;
							$multimapping=~s/have multiple alignments//g;
						}
						if($_ =~/overall read mapping rate/){
							$overall=$_;
							$overall=~s/^(\d+\.\d+)% .*/$1%/g;
						}
					
						$failed=$processed - $aligned;
					
						if($real_file and $failed){
							$summary_bw2->{$real_file}="$processed\t$aligned\t$multimapping\t$overall\t$failed";
						}
					}
					close STAT;
				}
			}
		}
	}
	if(lc($aligner) eq "bwa" ){
		my $dirname = $projectdir ."/bwa_results/";
		$summary_path->{$aligner}=$dirname;
		opendir ( DIR, $dirname ) || die "Error in opening dir $dirname\n";
		my @files= readdir DIR;
	
		foreach my $filename (sort @files){		     
			if($filename =~ /_align_summary.txt$/){
				open(RESULTS, $dirname ."/". $filename) || warn $!;	
				my $real_file=fileparse($filename);
				$real_file=~s/_align_summary.txt//g;
				my $processed;
				my $aligned;
				my $failed;
				my $multimapping;
				my $overall;
				my $line=1;
				while(<RESULTS>){
					chomp;
					if($line ==1){
						$processed=$_;
						$processed=~s/^(\d+) .+/$1/g;
					}
					if($line == 5){
						$aligned=$_;
						$aligned=~s/^(\d+) .+/$1/g;
					}
					$line++;
					my $failed=$processed-$aligned;
					my $perc= sprintf("%.2f",($aligned*100)/$processed);
					if($processed and $aligned){
						$summary->{$real_file}="$processed\t$aligned ($perc%)\t$failed";
						next;
					}
				}
				
			}
		}
	}
	if(lc($aligner) eq "hisat2" ){
		$summary_path->{$aligner}=$projectdir ."/hisat2_results/";
		open(STAT,$statsfile);
		my $real_file;
		my $processed;
		my $aligned;
		my $failed;
		my $multimapping;
		my $overall;
		while(<STAT>){
			chomp;
			if($_ =~/^HISAT2/){
				my $file=$_;
				$file=~s/HISAT2 :: File://g;
				$real_file=fileparse($file);
				$processed=0;
				$aligned=0;
				$failed=0;
				$multimapping=0;
				$overall=0;
			}
			if($_ =~/reads; of these:/){
				$processed=$_;
				$processed=~s/^(\d+) reads; of these:/$1/g;
			}
			if($_ =~/ aligned 0 times/){
				$failed=$_;
				$failed=~s/\s+(\d+) .*/$1/g;
			}
			if($_ =~/aligned exactly 1 time/){
				$aligned=$_;
				$aligned=~s/\s+(\d+) .*/$1/g;
			}
			if($_ =~/aligned >1 times/){
				$multimapping=$_;
				$multimapping=~s/\s+(\d+) .*/$1/g;
			}
			if($_ =~/overall alignment rate/){
				$overall=$_;
				$overall=~s/^(\d+\.\d+)% .*/$1%/g;
			}
			if($real_file and $failed){
				$summary_bw2->{$real_file}="$processed\t$aligned\t$multimapping\t$overall\t$failed";
			}
		}
	}
	if(lc($aligner) eq "star" ){
		$summary_path->{$aligner}=$projectdir ."/star_results/";

		opendir ( DIR, $summary_path->{$aligner} ) || die "Error in opening dir ".$summary_path->{$aligner} ."\n";
		while( (my $real_file = readdir(DIR))){
		    if($real_file =~/_Log\.final\.out$/){
			   open(STAT,$summary_path->{$aligner} ."/". $real_file) || die "$! Can't find ".$summary_path->{$aligner}."/$real_file\n";
			   my $processed;
			   my $aligned;
			   my $failed;
			   my $multimapping;
			   my $overall;
			   while(<STAT>){
				   chomp;
				   if($_ =~/Number of input reads/){
			   			$processed=$_;
			   			$processed=~s/\s+Number of input reads\s+\|\s+(\d+)/$1/g;
			   	   }
				   if($_ =~/Uniquely mapped reads %/){
			   			$overall=$_;
			   			$overall=~s/\s+Uniquely mapped reads %\s+\|\s+(\d+)/$1/g;
			   	   }
				   if($_ =~/Uniquely mapped reads number/){
						$aligned=$_;
				   		$aligned=~s/\s+Uniquely mapped reads number\s+\|\s+(\d+)/$1/g;
				   }
				   if($_ =~/Number of reads mapped to multiple loci/){
						$multimapping=$_;
						$multimapping=~s/\s+Number of reads mapped to multiple loci\s+\|\s+(\d+)/$1/g;
				    }
		 		   if($_ =~/% of reads unmapped: other/){
		 				$failed=$_;
		 				$failed=~s/\s+% of reads unmapped: other\s+\|\s+(\d+)/$1/g;
		 		    }
					if($real_file and $failed){
						my $used_file=$real_file;
						$used_file=~s/_Log.*//g;
						$summary_bw2->{$used_file}="$processed\t$aligned\t$multimapping\t$overall\t$failed";
					}
				}
			   
		   }
		}
		close STAT;
		closedir(DIR);
	}
	
	#saving stats
	if(scalar(keys %$summary)>0){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
		print SUMM "\nAlignment [".$summary_path->{$aligner}."]\n";
		print SUMM "Filename\tProcessed Reads\tAligned reads\tFailed to align\n";
		foreach my $processed_file (sort keys %$summary){
			print SUMM $processed_file ."\t". $summary->{$processed_file}."\n";
		}
		close SUMM;
	}
	if(scalar(keys %$summary_bw2)>0){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
		print SUMM "\nAlignment [".$summary_path->{$aligner}."]\n";
		print SUMM "Filename\tProcessed Reads\tAligned reads\tMultimapping reads\tOverall alignment\tFail to align\n";
		foreach my $processed_file (sort keys %$summary_bw2){
			print SUMM $processed_file ."\t". $summary_bw2->{$processed_file}."\n";
		}
		close SUMM;
	}
}
sub system_bash {
  my @args = ( "bash", "-c", shift );
  system(@args);
}

sub date{
	#my $dt = DateTime->now(time_zone=>'local');
	#return($dt->hms . " [" . $dt->dmy ."]");
	use Time::localtime;
	my $now = ctime();
	return("[$now]");
}
1;

