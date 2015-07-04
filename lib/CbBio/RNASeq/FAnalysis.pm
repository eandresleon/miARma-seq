#########################################################################	
#	Functinal analysis package							#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2015 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################

package CbBio::RNASeq::FAnalysis;
#Export package system
require Exporter;
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(F_Analysis);

use strict;
use DateTime;
use Statistics::R;
use Cwd;
use Cwd 'abs_path';

=head1 NAME

 FAnalysis 

=head1 SYNOPSIS

FAnalysis package is composed by three subroutines: DE_noiseq, F_Analysis and F_Analysis. The aim of this package is 
perform the functional analysis (DE) analysis from the tabulated files with the count of the reads 
(from htseq-count software) with Noiseq software (DE_Noiseq) or EdgeR software (F_Analysis). F_Analysis is a common function
to execute one of these functions or both. 


=head1 Methods

=head2 DE_noiseq

  Example    : 
  F_Analysis(
	projectdir=>".",
	noiseq_dir=>"../noiseq_result/",
	edger_dir=>"../edgeR_result/",
	logfile=> "logfile.log",
	verbose=>1,
	organism=>"mouse",
 );
  
  Description: F_Analysis performs a functional analysis (DE) Analysis from the tabulated files with the count of the reads 
  (from htseq-count software) with Noiseq and/or EdgeR software. These analysis generate tabulated files with the results of each condition evaluated
  as well as pdf files with descriptive plots.   
  Input parameters: 
	
	Mandatory parameters:
	[projectdir] Path of the directory where the output files will be saved.
	[noiseq_dir] Path of the directory which contains the file with the reads. This directory will be configured as working directory for R
	[edger_dir] Name of the tab file which contains the number of reads from the htseq analysis
	[organism] Path of the target file. This is a tabulated file which contains the experimental condiction of each sample. 
	[logfile] Path of run.log file where execution data will be printed
  	 
  	 Optional parameters:
  	 [verbose] Optional argument to show the execution data on screen
  	 [edger_cutoff] Optional argument to show the execution data on screen
  	 [noiseq_cutoff] Optional argument to show the execution data on screen
  
  Requeriments: F_Analysis function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- R v3.1.0 or higher software correctly installed
  	- EdgeR package correctly installed
  	- Noiseq package correctly installed
  	- Output file from htseq-count analysis with the number of the reads of each sample
  	
  Returntype : pdf file with descriptive plots of the analysis and tab files (xls extension) for condition evaluated in Noiseq_results and or EdgeR-results directory.  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub F_Analysis{
	
	#Arguments provided by user are collected by %args. 
	my %args=@_;
	my $projectdir=$args{"projectdir"}; #Path of the directory where will be saved the results.
	my $noiseq_dir=$args{"noiseq_dir"}; #Character string that will appear in the name the results file
	my $edger_dir=$args{"edger_dir"}; #This value refers to filter processing in the reads (Should be "yes" or "no").
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be saved
	my $organism=$args{"organism"}; #Specific software to perform the functional analysis Analysis
	my $verbose=$args{"verbose"};
	my $seq_id=$args{"seq_id"}; #Type of entities used in differential expresion, genes (gene_id) or transcripts (transcript_id)
	
	my $dataset;
	if(defined($args{"dataset"})){
		$dataset=$args{"dataset"};
	}
	else{
		$dataset=undef;
	}
	open(LOG,">$logfile") || warn "Cant create logfile ($logfile)\n";
	
	#Checking mandatory parameters
	if($projectdir and $logfile and $organism and ($noiseq_dir or $edger_dir) and $seq_id){
		#Obtaining the absolute path of the directory and the file to R execution
	  	my $projectdir = abs_path($projectdir);
		#for edgeR
		if($edger_dir){
			print STDERR "LOG :: " . date() . " Starting a Functional Analysis based on edgeR data\n";
			print LOG "LOG :: " . date() . " Starting a Functional Analysis based on edgeR data\n";
			
			my $edger_cutoff;
			if(defined $args{"edger_cutoff"}){
				$edger_cutoff=$args{"edger_cutoff"};
			}
			else{
				$edger_cutoff=0.05;
			}
		  	opendir(DIR, $edger_dir) || die $!;
			my @files= readdir(DIR);
			#@files=map("$edger_dir/$_",@files);
			my $universe_edger_data;
			my $up_edger_data;
			my $down_edger_data;
			
			my $up_number=0;
			my $down_number=0;
			my $universe_number=0;
			foreach my $edger_files(@files){
				if($edger_files =~ /\.xls$/){
					print STDERR "LOG :: " . date() . " Reading $edger_files and filtering by FDR <=$edger_cutoff\n";
					print LOG "LOG :: " . date() . " Reading $edger_files\n" if($verbose);
					#Read edgeR data
					open(EDGER,$edger_dir."/" . $edger_files) || warn "cant find $edger_dir/$edger_files\n";
					while(<EDGER>){
						chomp;
						my($feature,$fc,$cpm,$pvalue,$fdr)=split(/\t/);
						$universe_edger_data->{$feature}++;
						#filtering over-expressed
						if($fdr<=$edger_cutoff and $fc>0){
							$up_edger_data->{$feature}++;
							$up_number++;
							#last if scalar(keys %$up_edger_data)==5;
						}
						#filtering down-expressed
						if($fdr<=$edger_cutoff and $fc<0){
							$down_edger_data->{$feature}++;
							$down_number++;
						}
						$universe_number++;
					}
					close EDGER;
					print LOG "Number of genes at universe: " . scalar(keys %$universe_edger_data) ."\n" if($verbose);
					print LOG "Number of Up-regulated genes (edgeR cutoff=$edger_cutoff): " . scalar(keys %$up_edger_data)  ."\n" if($verbose);
					print LOG "Number of Down-regulated genes (edgeR cutoff=$edger_cutoff): " . scalar(keys %$down_edger_data)  ."\n" if($verbose);
					
					if($universe_number>0 and ($up_number>0 or $down_number>0)){
						goseq(
							universe=>$universe_edger_data,
							up=>$up_edger_data,
							down=>$down_edger_data,
							org=>$organism,
							method=>"edgeR",
							verbose=>$verbose,
							logfile=>$logfile,
							projectdir=>$projectdir,
							organism=>$organism,
							cut_off=>$edger_cutoff,
							seq_id=>$seq_id,
							dataset=>$dataset,
						);
					}
					else{
						warn "WARN :: ".date()." The number of selected proteins is 0. No functional analysis can be done\nWARN :: ".date()." Upregulated genes n=$up_number, Downregulated genes n=$down_number\nPlease check $logfile\n\n";
						return();
					}
				}
			}
		}		
		if($noiseq_dir){
			print STDERR "\nLOG :: " . date() . " Starting a Functional Analysis based on NoiSeq data\n";
			print LOG "LOG :: " . date() . " Starting a Functional Analysis based on NoiSeq data\n";
			my $noiseq_cutoff;
			if(defined $args{"noiseq_cutoff"}){
				$noiseq_cutoff=$args{"noiseq_cutoff"};
			}
			else{
				$noiseq_cutoff=0.8;
			}
		  	opendir(DIR, $noiseq_dir) || die $! . " Noiseq dir: [$noiseq_dir]\n";
			my @files= readdir(DIR);

			my $universe_noiseq_data;
			my $up_noiseq_data;
			my $down_noiseq_data;
			
			my $up_number=0;
			my $down_number=0;
			my $universe_number=0;
			foreach my $noiseq_files(@files){
				if($noiseq_files =~ /\.xls$/){
					print STDERR "LOG :: " . date() . " Reading $noiseq_files\n";
					print LOG "LOG :: " . date() . " Reading $noiseq_files\n" if($verbose);
					#Read edgeR data
					open(NOISEQ,$noiseq_dir."/" . $noiseq_files) || warn "cant find $noiseq_dir/$noiseq_files\n";
					while(<NOISEQ>){
						chomp;
						my($feature,undef,undef,$fc,undef,$fdr)=split(/\t/);
						#print LOG "$feature\t$fc\t$fdr\n";
						$universe_noiseq_data->{$feature}++;
						#filtering over-expressed
						if($fdr>=$noiseq_cutoff and $fc>0){
							$up_noiseq_data->{$feature}++;
							$up_number++;
							#last if scalar(keys %$up_edger_data)==5;
						}
						#filtering down-expressed
						if($fdr>=$noiseq_cutoff and $fc<0){
							$down_noiseq_data->{$feature}++;
							$down_number++;
						}
						$universe_number++;
					}
					close NOISEQ;
					print LOG "Number of genes at universe: " . scalar(keys %$universe_noiseq_data) ."\n" if($verbose);
					print LOG "Number of Up-regulated genes (NoiSeq cutoff=$noiseq_cutoff): " . scalar(keys %$up_noiseq_data)  ."\n" if($verbose);
					print LOG "Number of Down-regulated genes (NoiSeq cutoff=$noiseq_cutoff): " . scalar(keys %$down_noiseq_data)  ."\n" if($verbose);
				}
			}

			#Once all fila are read
			goseq(
				universe=>$universe_noiseq_data,
				up=>$up_noiseq_data,
				down=>$down_noiseq_data,
				org=>$organism,
				method=>"noiseq",
				verbose=>$verbose,
				logfile=>$logfile,
				projectdir=>$projectdir,
				organism=>$organism,
				cut_off=>$noiseq_cutoff,
				seq_id=>$seq_id,
				dataset=>$dataset,
			);
		}
		else{
			print STDERR "WARN :: No files provided for functional anaysis, please check your edgeR folder ($edger_dir) and/or your NoiSeq folder ($noiseq_dir)\n";
			return();
		}
	}
	else{
		#If mandatory parameters have not been provided program will die and show error message
		warn("F_Analysis ERROR :: ".date()." projectdir($projectdir), organism($organism), logfile ($logfile), seq_id ($seq_id) have not been provided");
		help_F_Analysis();
	}
	close LOG;
	sub help_F_Analysis{
	    my $usage = qq{
		  	$0 

			Mandatory arguments:
			  [projectdir] Path of the directory where will be saved the QC report.
			  [up] Path of the tab file which contains the number of reads from the htseq analysis
			  [down] Path of the tabulated file which contains the experimental condiction of each sample. 
			  [organism] Path of the tabulated file which contains the experimental condiction of each sample. 
			  [method] Name of the method: two posibilities availables: edgeR or noiseq
			  [seq_id] gene_id or transcript_id
			    Example:

			 goseq(projectdir=\"./Project\",up=\"up\",down=\"down\",universe=\"universe\",organism=\"human\",method=\"edgeR\",seq_id=\"gene_id\")"  	 
	};

	print STDERR $usage;
	exit(); 
	}
}


sub goseq{
	#Arguments provided by user are collected by %args. 
	my %args=@_;
	my $projectdir=$args{"projectdir"}; #Path of the directory where will be saved the results.
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be saved
	my $organism=$args{"organism"}; #Specific software to perform the functional analysis Analysis
	my $universe=$args{"universe"}; #all detected genes in my sample
	my $up=$args{"up"}; #Genes upregulated in my experiment
	my $down=$args{"down"}; #Genes downregulated in my experiment 
	my $method=$args{"method"}; # Method used for DE calculation, can be either edgeR or NoiSeq
	my $verbose=$args{"verbose"};
	my $cut_off=$args{"cut_off"};
	my $seq_id;
	if(defined($args{"seq_id"})){
		$seq_id=$args{"seq_id"};
	}
	else{
		$seq_id="gene_id";
	}
	my @R_bin=`which R`;
	#Executing the command
	if(scalar(@R_bin)<1){
		die "Functional Analysis ERROR :: system args failed: $? : Is R installed and exported to \$PATH ?";
	}
	open(LOG,">>$logfile") || warn "Cant create logfile ($logfile)\n";
	
	#Checking mandatory parameters
	if($projectdir and $logfile and $organism and $up and $universe and $down and $method and $seq_id){
		my $output_dir=$projectdir."/Functional_Analysis_results";
		my $output_file="$output_dir/upRegulated_$method\_results.xls";

		if (scalar(keys %$up)<1 or scalar(keys %$down)<1 or scalar(keys %$universe)<1){
			warn "WARN :: ".date()." The number of selected proteins is 0. No functional analysis can be done\nWARN :: ".date()." Upregulated genes n=".scalar(keys %$up).", Downregulated genes n=".scalar(keys %$down)."\nPlease check $logfile\n\n";
			return();
		}
		if(scalar(keys %$up)<75 or scalar(keys %$down)<75){
			print STDERR "\nWARN :: The number of differentially expressed genes are too small. Some errors could appear. Try to use a less restringet cut_off (>$cut_off)\n\n";
		}
			
		my $Rcommand="projectdir=\"$output_dir\",up=c(".join(",", keys (%$up))."),down=c(".join(",", keys (%$down))."),universe=c(".join(",", keys (%$universe))."),organism=\"$organism\",method=\"$method\",seq_id=\"$seq_id\"";
		# Printing the date and command execution on screen
		print STDERR "GOSeqR :: ".date()." Starting Functional Analysis\n"; 

		#Creating results directory
		my $command="mkdir -p ".$output_dir;
		system($command) == 0
		   or die "F_Analysis ERROR :: ".date()." System args failed: $? ($command)";

		#Calling R from perl
		my $R;
		#If user has defined the directory where R is installed the bridge will be created since that directory
		if(defined($args{"Rdir"})){
			my $Rdir=$args{"Rdir"};#path where R software is installed
			$R = Statistics::R->new($Rdir) ;
		}
		else
		{
			$R = Statistics::R->new() ;
		}

		#Starting R 
		$R->startR;
		#Declaring R instructions for the functional analysis analysis. DE_noiseq R function is needed 
		#source("http://valkyrie.us.es/CbBio/RNASeq/R-Scripts/F_Analysis.R")
		
		my $cmds = <<EOF;
		setwd("$output_dir")
		source("../../../lib/CbBio/RNASeq/R-Scripts/F_Analysis.R")
		resultsfiles<-NA
		resultsfiles<-F_Analysis($Rcommand)
EOF
		
		#Printing the execution data on log file and on the screen if verbose parameter is defined 
		open (LOG,">> ".$logfile) || die "F_Analysis ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "F_Analysis :: ".date()." Executing $cmds\n";
		if($verbose){
			#print STDOUT "F_Analysis :: ".date()." Executing $cmds\n";
		}
		
		#R commands execution
		my $out2 = $R->run($cmds);
		#Obtaining the files generated in the analysis
		my @media = $R->get('resultsfiles');
				
		foreach my $result (@media){
			if($result eq "NA"){
				print STDERR "ERROR :: Nothing significative was found\n";
			}
			else{
				foreach my $res_file (@{$result}){
					print "F_Analysis :: ".date()." The file ".$res_file." has been generated \n";
				}
			}
		}
		close LOG;
		
	}
	else{
		#Registering error
	   	open(LOG,">> ".$logfile) || die "F_Analysis ERROR :: ".date()."Can't open '$logfile': $!";
	    #print LOG "F_Analysis ERROR :: ".date()." projectdir($projectdir), up($up), down($down), universe($universe), organism($organism), logfile ($logfile) have not been provided");
	    close(LOG);
		#If mandatory parameters have not been provided program will die and show error message
		warn("F_Analysis ERROR :: ".date()." projectdir($projectdir), up($up), down($down), universe($universe), organism($organism), logfile ($logfile) have not been provided");
		help_goSeq();
	}
	close LOG;
	sub help_goSeq{
	    my $usage = qq{
		  	$0 

			Mandatory arguments:
			  [projectdir] Path of the directory where will be saved the QC report.
			  [up] Path of the tab file which contains the number of reads from the htseq analysis
			  [down] Path of the tabulated file which contains the experimental condiction of each sample. 
			  [organism] Path of the tabulated file which contains the experimental condiction of each sample. 
			  [method] Name of the method: two posibilities availables: edgeR or noiseq
			    Example:

			 goseq(projectdir=\"./Project\",up=\"up\",down=\"down\",universe=\"universe\",organism=\"human\")"  	 

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

