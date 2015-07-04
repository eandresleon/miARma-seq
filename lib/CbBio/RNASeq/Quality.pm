#########################################################################	
#	Quality package					 									#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2015 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################

package CbBio::RNASeq::Quality;
#Export package system
require Exporter;
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(FastQC FastQCStats QC_EdgeR);

use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use DateTime;
use Statistics::R;

=head1 NAME

 Quality

=head1 SYNOPSIS

Quality package is composed by three subroutines: FastQC, FastQCStats and QC_EdgeR. The aim of this 
package is perform a quality analysis of the fastq files by FastQC function or the tabulated files 
with the count of the reads (from htseq-count software) by QC_EdgeR. This package also contains the 
FastQCStats subroutine to extract the main stats of FastQC analysis. 


=head1 Methods

=head2 FastQC

  Example    : 
  FastQC( 
  	projectdir=>".", 
  	file=>"../reads/file.fatsq", 
  	verbose=>"verbose", 
  	threads=>"4", 
  	logfile=> "run.log"
  	)
  Description: FastQC function performs a fastqc analysis on the provided file to check the 
  quality of the raw data. The output files will be saved on a new directory called 
  fastqc_results on the project directory. Execution data will be written on a logfile 
  in the provided path. Optionally, number of threads can be provided to run the process 
  faster and verbose to show execution data on the screen. 
  Input parameters: 
	Mandatory parameters: 
  	 [file] Name of the file which is going to be processing (fastq format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [projectdir] Directory where results directory will be created
  	 [label] Character string to put in the name of the results directory
  	 Optional parameters:
  	 [threads] Number of threads to perform the analysis faster (1 thread by default)
  	 [verbose] Option to show the execution data on the screen
  Returntype : Directory named as input file at directory fastqc_results. FastQC function also
  returs the path of the results directory.
  Requeriments: FastQC function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Fastqc v0.10.1 or higher software correctly installed
  	- Input files on fastq format (compressed files are accepted)
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub FastQC{
	
	#Arguments provided by user are collected by %args. File, dir, label, projectdir and logfile are mandatory
	#arguments while threads and verbose are optional.
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};

	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/common/fastqc/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/common/fastqc/";
	}
	
	my @fastqc_bin=`which fastqc`;
	#Executing the command
	if(scalar(@fastqc_bin)<1){
		die "FASTQC ERROR :: system args failed: $? : Is fastqc installed and exported to \$PATH ?";
	}
	my @java_bin=`which java`;
	#Executing the command
	if(scalar(@java_bin)<1){
		die "FASTQC ERROR :: system args failed: $? : Is java installed and exported to \$PATH ?";
	}
	
	#my $fastqc_bin="$miARmaPath/bin/common/fastqc/fastqc";
	my $file=$args{"file"}; #Name of the file to perform the analysis
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be saved
	my $threads=$args{"threads"}; #Optional number of threads to perform the analysis faster
	my $verbose=$args{"verbose"}; #Optional argument to show the execution data on screen
	my $label=$args{"prefix"}; #Label to write in the directory results name
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	#Describing results directory 
	my $output_dir= "/".$label."_fastqc_results/";
	
	#Checking the mandatory arguments 	
	if($file and $logfile and $label and $projectdir){
		#Obtaining the absolute path of the file
		if($file =~ /.*\.fastq$/ or $file =~ /.*\.fastq\.gz$/ or $file=~ /.*\.fq$/ or $file=~ /.*\.fq\.gz$/ or $file=~ /.*\.fastq\.lane\.clean$/ or $file=~ /.*\.fq\.bz2$/ or $file=~ /.*\.fastq\.bz2$/){
			$file = abs_path($file);
		  	system("mkdir -p $projectdir");
			#Printing process data on screen and in log file
			print STDERR "FASTQC :: ".date()." Checking $file for FastQC analysis\n";
	
			#Variable declarations
			my $command;
			my $commanddef;

			#Checking if user has defined the number of threads to perform the analysis.
			if($threads>0){
				#Fastqc execution command with a defined number of threads and the output directory
				$command="fastqc -f fastq -t ".$threads." ".$file." -o ".$projectdir.$output_dir ;
				#commandef is the command will be executed by system composed of the results directory creation 
				#and the fastqc execution. The error data will be redirected to the run.log file
				$commanddef="mkdir -p ".$projectdir.$output_dir." ; ".$command." >> ".$logfile." 2>&1";
			}
			else{
				#Fastqc execution command with the default number of threads (1) and the output directory 
				$command="fastqc -f fastq ".$file." -o ".$projectdir.$output_dir;
				#commandef is the command that will be executed by system composed of the results directory creation
				#and the fastqc execution. The error data will be redirected to the run.log file 
				$commanddef="mkdir -p ".$projectdir.$output_dir." ;".$command ." >> ".$logfile." 2>&1";
			}
			
			#Execution data will be register on the run.log file. Opening the run.log file 
			open(LOG,">> ".$logfile) || die "FASTQC ERROR :: Can't open [$logfile]: $!";
			# Printing the date and command execution on the run.log file
			print LOG "FASTQC :: ".date()." Checking $file for FastQC analysis\n";
			print LOG "FASTQC :: ".date()." Executing $commanddef\n";
			#Executing the command or if system can't be executed die showing the error.
			system($commanddef) == 0
	   	    or die "FASTQC ERROR :: system args failed: $? ($commanddef)";
	   	    #If verbose option has been provided by user print the data on screen
			if($verbose){
				print STDOUT "FASTQC :: ".date()." Executing $commanddef\n" if($verbose);
			}
			#The path of output results is returned to main program
		}
		return($projectdir.$output_dir);
		
		#else{
		#	warn("FASTQC ERROR :: ".date()." File($file) has an invalid format.");
		#	help_FastQC();
		#}
	}
	else
	{
		#Registering the error in run.log file
   		open(LOG,">> ".$logfile) || die "FASTQC ERROR :: Can't open $logfile: $!";
    	print LOG "FASTQC ERROR :: ".date()." Projectdir($projectdir), file($file), label($label) and/or logfile($logfile) have not been provided";
    	close LOG;

		#If mandatory parameters have not been provided program will die and show error message
		warn ("FASTQC ERROR :: ".date()." Projectdir($projectdir), file($file), label($label) and/or logfile($logfile) have not been provided.");
		help_FastQC();
	}

	#Help function shows the usage of FastQC function
	sub help_FastQC{
	    my $usage = qq{
		  	$0 

			Needed parameters: 
  	 		[file] Name of the file in fastq format which is going to be processing (Extensions allowed: .fastq, .fq, .fq.gz, .fastq.gz)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[projectdir] Directory where results directory will be created
  	 		[label] Character string to put in the name of the results directory
						 
			Optional parameters:
			[threads] Number of threads to perform the analysis faster (1 thread by default)
  	 		[verbose] Option to show the execution data on the screen.
  	 		               
			Examples:
			FastQC(projectdir=>".", file=>"../reads/file.fastq", verbose=>"verbose", threads=>"4", logfile=> "run.log"); 

	};

	print STDERR $usage;
	exit(); 
	}
}

=head2 FastQCStats

  Example    : 
  FastQCStats(
  	dir=>"./Pre_fastqc_results",
  	verbose=>"verbose",
  	statsfile=> "stats.log",
  	logfile=> "run.log");
  Description: FastQCStats function shows the main stats from the fastqc analysis. This 
  function collects the directory to access to the previously created by fastqc function 
  fastqc_results directory and open the directory corresponding with each file analyzed. 
  Then it opens the fastqc_report.html and looks for the number of total sequences, 
  the sequence lenght, the enconding type and the GC content. These data will be printed 
  on statsfile or on the screen if verbose option was selected. 
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where fastqc results are placed
  	 [statsfile] Path of stats.log file where stats data will be printed
  	 [logfile] Path of run.log file where execution data will be printed
  	 Optional parameters:
  	 [verbose] Option to show the stats data on the screen
  Requeriments: FastQCStats function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Previous fastqc results on the fastqc_results directory on the provided directory
  Returntype : print data on stats.log or on the screen
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub FastQCStats{
	#Arguments provided by user are collected by %args. Dir, statsfile and logfile are mandatory
	#arguments while verbose is optional. 
	my %args=@_;
	my $dir=$args{"dir"}; #Input directory where fastqc results are placed
	my $verbose=$args{"verbose"}; #Optional parameter to show the stats data on the screen
	my $statsfile=$args{"statsfile"}; #Path of stats.log file where stats data will be printed
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be printed
	
	#Variable declaration
	my @filesres; 
	
	#Checking the dir argument 
	if($dir and $statsfile and $logfile){
		#Opening the fastqc_results directory and reading the directories.	
		opendir(DIRRES, $dir) || die "FASTQCSTATS ERROR :: Can't open $dir : $!";
		my @filesres=readdir(DIRRES);
		#Accessing to each directory
		foreach my $fileres(@filesres){
			#Declaring variables to keep the stats 
			my $seqnumber=0;
			my $seqlength=0;
			my $encoding;
			my $gc;
			#Selecting the results directory
			if ($fileres !~ /\.zip$/ and $fileres !~ /^\./){
				#Opening the report file of each directory
				open(HTML, $dir.$fileres."/fastqc_report.html") || die "FASTQCSTATS ERROR :: Can't open $dir/$fileres/fastqc_report.html : $!";;
				#Declaring the line counters. These counters are used because the 
				#data appears in a different lines than the regular expression. When the
				#regular expression is found the counter changes and in the next line 
				#evaluate the regular expresion to obtain the value. 
				my $linebeftseq=0;
				my $linebefseqlen=0;
				my $linebefenc=0;
				my $linebefgc=0;
				#Reading the lines of the report
				while(<HTML>){
					chomp;
					#If the line contains Total Sequences the counter increases and enter
					#in the next clause obtaining the number of total sequences 
					if($_ =~ /Total Sequences/){
						$linebeftseq=1;
					}
					elsif($linebeftseq == 1)
					{
						$_ =~s/<td>(\d+)<\/td>//;
						$seqnumber=$1;
						$linebeftseq=0;
					}
					#If the line contains Sequence length the counter increases and enter
					#in the next clause obtaining the lenght of the Sequence
					elsif($_ =~ /Sequence length</){
						$linebefseqlen=1;
					}
					elsif($linebefseqlen == 1)
					{
						$_ =~s/<td>(.+)<\/td>//;
						$seqlength=$1;
						$linebefseqlen=0;
					}
					#If the line contains Encoding the counter increases and enter
					#in the next clause obtaining the encoding system
					elsif($_ =~ /Encoding/){
						$linebefenc=1;
					}
					elsif($linebefenc == 1)
					{
						$_ =~s/<td>(.+)<\/td>//;
						$encoding=$1;
						$linebefenc=0;
					}
					#If the line contains %GC the counter increases and enter
					#in the next clause obtaining the GC content
					elsif($_ =~ /%GC/){
						$linebefgc=1;
						
					}
					elsif($linebefgc == 1)
					{
						$_ =~s/<td>(\d+)<\/td>//;
						$gc=$1;
						$linebefgc=0;
					}
					#If the line don�t contain any of the regular expressions the counters
					# will keep the 0 value.
					else
					{
						$linebeftseq=0;
						$linebefseqlen=0;
						$linebefenc=0;
						$linebefgc=0;
					}
				}
				#The statistical data is printed on logfile provided by user. Opening the log file.
				open(STATS, ">> ".$statsfile) || die "FASTQCSTATS ERROR :: Can't open $statsfile : $!";
				#Printing the previous variables in logfile
				print STATS "FASTQCSTATS :: ".date()."\n Name\t $fileres\n 
				Total Sequences:\t $seqnumber\t Sequence length:\t $seqlength\n 
				Encoding:\t $encoding\t GCcontent:\t $gc%\n";
				close STATS;
				#If verbose option has been provided program will print the previous variables
				# on screen too. 
				if($verbose){
					print STDOUT "FASTQCSTATS :: ".date()."\n Name\t $fileres\n 
					Total Sequences:\t $seqnumber\t Sequence length:\t $seqlength\n 
					Encoding:\t $encoding\t GCcontent:\t $gc%\n";
				}
			#Closing the directory and the file
			close HTML;
			close DIRRES;
			}
		}
		print STDERR "\n\nFASTQC :: ".date()." Please check $dir/Pre_fastqc_results for detailed quality information about each sample.\nA summary can be consulted in $statsfile\n\n";
	}
	else
	{
		#Registering error in logfile
   		open(LOG,">> ".$logfile) || die "FASTQCSTATS ERROR :: Can't open $logfile : $!";
    	print LOG "FASTQCSTATS ERROR :: ".date()." Directory ($dir), logfile($logfile) and/or statsfile($statsfile) have not been provided";
    	close LOG;


		#If mandatory parameters have not been provided program will die and show error mesagge
		warn("FASTQCSTATS ERROR :: ".date()." Directory ($dir), logfile($logfile) and/or statsfile($statsfile) have not been provided to the FastQCStats subroutine");
		help_FastQCStats();
	}

	#Help function shows the usage of FastQCStats function
	sub help_FastQCStats{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[dir] Input directory where fastqc results are placed
  	 		[statsfile] Path of stats.log file where stats data will be printed
  	 		[logfile] Path of run.log file where execution data will be printed
						 
			Optional parameters:
  	 		[verbose] Option to show the execution data on the screen.
  	 		               
			Examples:
			FastQCStats(dir=>"./Pre_fastqc_results", logfile=> "run.log", statsfile=> "stats.log",verbose=>verbose); 

	};

	print STDERR $usage;
	exit(); 
	}
	
}

=head2 QC_EdgeR

  Example    : 
  QC_EdgeR(
			projectdir=>".",
			dir=>"/Hypoxia_htseq_results",
			file=>"htseqresults.tab",
			targetfile=>"targets.txt",
			label=>"Hypoxia",
			filter=>"yes",
			cpmvalue=>2,
			repthreshold=>2,
			ref1=>2,
			ref2=>"T0,
			logfile=>"run.log",
			Rdir=>"/Users/Applications/R",
			verbose=>"verbose"
		);
  
  Description: QC_EdgeR function calls a R function called QC_EdgeR which takes a tab file with the 
  number of the reads from htseq-count analysis and a target file with the experimental conditions 
  of the samples to analyze the quality of the samples. For this purpose QC_EdgeR prints in a pdf file 
  differents plots: boxplots and density plots with the distribution of the reads in each sample with 
  and without normalization, MDS and PCA plots of the samples and heatmaps of the samples with the 
  more expressed elements(genes, miRNAs...) and between the samples. 
  Input parameters: 
	Mandatory parameters:
	[projectdir] Path of the directory where will be saved the QC report.
	[dir] the path of the directory which contains the files. This directory will be configured as working directory for R
	[file] name of the tab file which contains the number of reads from the htseq analysis
	[targetfile] Path of the tabulated file which contains the experimental condiction of each sample. First column must 
	contain the names to be used to the plots, and the next columns the condition of each factor. QC_EdgeR works with 1 or 2 factors.
	[label] Name of the experiment to appear in the title of the plots and in the name of the pdf file results.
	[filter] This value refers to filter processing in the reads (Should be "yes" or "no").
	[ref2] Reference condition to use in the analysis consisting on a character string that defines the reference (ie: if drug 
	conditions are with and without, ref2 should be without).
	[logfile] Path of run.log file where execution data will be printed
  	 Optional parameters:
  	 [cpmvalue] Cutoff for the counts per million value to be used in filter processing (1 cpm by default).
  	 [repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default).
  	 [ref1] Optional parameter only needed when 2 factors are analysed. Reference condition to use in the analysis consisting on 
	the number of the column of the target file which contains the reference condition (ie: if target contains drug and time as second 
	and third column and we want to use drug as reference ref1 should be 2). By default this argument is 2 for one factor experiments
	[Rdir] Path where R software is installed
	[verbose] Optional argument to show the execution data on screen
  Requeriments: QC_EdgeR function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- R v3.1.0 or higher software correctly installed
  	- Output file from htseq-count analysis with the number of the reads of each sample
  	
  Returntype : pdf file on the given directory with the differents plots.
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub QC_EdgeR{
	#First, check that R is in path:
	my @R_bin=`which R`;
	#Executing the command
	if(scalar(@R_bin)<1){
		die "QC_EDGER ERROR :: system args failed: $? : Is R installed and exported to \$PATH ?";
	}

	#Arguments provided by user are collected by %args. File, dir, projectdir, targetfile, label, filter
	#logfile is mandatory arguments while cpmvalue, repthreshold, verbose and Rdir are optional.
	my %args=@_;
	my $file=$args{"file"}; #name of the tab file which contains the number of reads from the htseq analysis
	my $dir=$args{"dir"}; #the path of the directory which contains the files. This directory will be configured as working directory for R
	my $projectdir=$args{"projectdir"}; #Path of the directory where will be saved the QC report.
	my $targetfile=$args{"targetfile"}; #path of the target file with the experimental conditions of each sample
	my $label=$args{"label"}; #string character that will appear in the name the results file
	my $filter=$args{"filter"}; #This value refers to filter processing in the reads (Should be "yes" or "no").
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be saved

	if($file and $dir and $projectdir and $targetfile and $label and $filter and $logfile){
	
		#Optional parameters
		my $verbose=$args{"verbose"}; #Optional argument to show the execution data on screen
		#These parameters have defaults values but can be modified if the user provides other ones. 
		my $cpmvalue=1; #Cutoff for the counts per million value to be used in filter processing.
		my $repthreshold=2; #Number of replicates that have to contains at least a defined number of reads per million to perform the filter.


		#If user has provided the cpmvalue, repthreshold and/or ref1 will be collected by args and these variables will be overwritten
		if(exists $args{"cpmvalue"}){
			$cpmvalue=$args{"cpmvalue"};
		}
		if(exists $args{"repthreshold"}){
			$repthreshold=$args{"repthreshold"};
		}


		# Printing the date and command execution on screen
		print STDERR "QC_EDGER :: ".date()." Starting Quality Control Analysis of $file \n"; 

		#Calling R from perl
		my $R;
		
		#If user has defined the directory where R is installed the bridge will be created since that directory
		if(defined($args{"Rdir"})){
			my $Rdir=$args{"Rdir"}; #path where R software is installed
			$R = Statistics::R->new($Rdir) ;
		}
		else
		{
			$R = Statistics::R->new() ;
		}

		#Starting R 
		$R->startR;

		#Declaring R instructions for the quality control analysis. QC_EdgeR R function is needed 
		my $cmds = <<EOF;
		setwd("$dir")
		source("http://valkyrie.us.es/CbBio/RNASeq/R-Scripts/QC_EdgeR.R") 
		filenames<-QC_EdgeR("$projectdir","$dir", "$file", "$targetfile", "$label", "$filter", $cpmvalue, $repthreshold)
EOF
		
		#Printing the execution data on log file and on the screen if verbose parameter is defined 
		open (LOG,">> ".$logfile) || die "QC_EDGER ERROR :: Can't open $logfile : $!";
		print LOG "QC_EDGER :: ".date()." Executing $cmds\n";
		close LOG;

		if($verbose){
				print STDOUT "QC_EDGER :: ".date()." Executing $cmds\n";
		}

		#R commands execution
		my $out2 = $R->run($cmds);

		#Obtaining the names of the generated files
		my $media = $R->get('filenames');
		print "QC_EDGER :: ".date()." The file ". $media." has been generated \n";
	}
	else
	{
		warn("QC_EDGER ERROR :: File($file), dir($dir), projectdir($projectdir), targetfile($targetfile), label($label) filter($filter), and/or logfile($logfile) have not been provided");
		help_QC_EdgeR();
		
	}

	#Help function shows the usage of QC_EdgeR function
	sub help_QC_EdgeR{
	    my $usage = qq{
		  	$0 

			Needed parameters:
			[projectdir] Path of the directory where will be saved the QC report.
			[dir] the path of the directory which contains the files. This directory will be configured as working directory for R
			[file] name of the tab file which contains the number of reads from the htseq analysis
			[targetfile] Path of the tabulated file which contains the experimental condiction of each sample. First column must 
			contain the names to be used to the plots, and the next columns the condition of each factor. QC_EdgeR works with 1 or 2 factors.
			[label] Name of the experiment to appear in the title of the plots and in the name of the pdf file results.
			[filter] This value refers to filter processing in the reads (Should be "yes" or "no").
			[logfile] Path of run.log file where execution data will be printed
						 
			Optional parameters:
			[cpmvalue] Cutoff for the counts per million value to be used in filter processing (1 cpm by default).
  	 		[repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default).
 			[Rdir] Path where R software is installed
			[verbose] Optional argument to show the execution data on screen
  	 		               
			Examples:
			QC_EdgeR(projectdir=>".", dir=>"/Hypoxia_htseq_results", file=>"htseqresults.tab", targetfile=>"targets.txt", label=>"Hypoxia", filter=>"yes", cpmvalue=>2, repthreshold=>2,
			ref1=>2,ref2=>"T0h",logfile=>"run.log", Rdir=>"/Users/Applications/R", verbose=>"verbose");
	};

	print STDERR $usage;
	exit(); 
	}

}

sub date{
	my $dt = DateTime->now;
	return($dt->hms . " [" . $dt->dmy ."]");
}	
1;