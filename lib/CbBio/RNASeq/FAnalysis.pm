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
@EXPORT=qw(F_Analysis F_AnalysisSummary);

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
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be saved
	my $organism=$args{"organism"}; #Specific software to perform the functional analysis Analysis
	my $verbose=$args{"verbose"};
	my $seq_id=$args{"seqid"}; #Type of entities used in differential expresion, genes (gene_id) or transcripts (transcript_id)
	my $miARmaPath=$args{"miARmaPath"};
	my $fc_threshold=$args{"fc_threshold"};
	
	my $dataset;
	if(defined($args{"dataset"})){
		$dataset=$args{"dataset"};
	}
	else{
		$dataset=undef;
	}
	open(LOG,">>$logfile") || warn "Cant create logfile ($logfile)\n";
	
	my $noiseq_dir=undef;
	my $edger_dir=undef;
	if(-e $projectdir ."/" . "Noiseq_results/"){
		$noiseq_dir=$projectdir ."/" . "Noiseq_results/";
	}
	if(-e $projectdir ."/" . "EdgeR_results/"){
		$edger_dir=$projectdir ."/" . "EdgeR_results/";
	}
	#Checking mandatory parameters
	my $check=0;
	if($projectdir and $logfile and $organism and ($noiseq_dir or $edger_dir) and $seq_id){
		print date() . " Starting a Functional Analysis.\n";
		
		#Obtaining the absolute path of the directory and the file to R execution
	  	my $projectdir = abs_path($projectdir);
		#for edgeR
		if($edger_dir){
			print LOG date() . " Starting a Functional Analysis based on edgeR data\n";
			print STDOUT date() . " Starting a Functional Analysis based on edgeR data\n" if($verbose);
			$check=1;
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
					my $label=$edger_files;
					$label=~s/(.*)EdgeR_results_(.*)\.xls/$1$2/g;
					print STDOUT date() . " Reading $edger_files and filtering by FDR <=$edger_cutoff\n" if($verbose);
					print LOG date() . " Reading $edger_files and filtering by FDR <=$edger_cutoff\n";
					#Read edgeR data
					open(EDGER,$edger_dir."/" . $edger_files) || warn "cant find $edger_dir/$edger_files\n";
					while(<EDGER>){
						chomp;
						my($feature,$length,$fc,$cpm,$pvalue,$fdr)=split(/\t/);
						$universe_edger_data->{$feature}++;
						#filtering over-expressed
						if($fdr<=$edger_cutoff and $fc>=$fc_threshold){
							$up_edger_data->{$feature}++;
							$up_number++;
							#last if scalar(keys %$up_edger_data)==5;
						}
						#filtering down-expressed
						if($fdr<=$edger_cutoff and $fc <($fc_threshold * -1)){
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
							label=>$label,
							miARmaPath=>$miARmaPath
						);
					}
					else{
						warn date()." WARN :: The number of selected proteins is 0. No functional analysis can be done\n". date()." Upregulated genes n=$up_number, Downregulated genes n=$down_number\n".date()." Please check $logfile\n\n";
						return();
					}
				}
			}
			#return();
		}		
		if($noiseq_dir){
			$check=1;
			print STDOUT date() . " Starting a Functional Analysis based on NoiSeq data\n" if($verbose);
			print LOG date() . " Starting a Functional Analysis based on NoiSeq data\n";
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
					my $label=$noiseq_files;
					$label=~s/(.*)NOISeq_results_(.*)\.xls/$1$2/g;
					
					print STDOUT date() . " Reading $noiseq_files\n" if($verbose);
					print LOG date() . " Reading $noiseq_files\n" if($verbose);
					#Read edgeR data
					open(NOISEQ,$noiseq_dir."/" . $noiseq_files) || warn "cant find $noiseq_dir/$noiseq_files\n";
					while(<NOISEQ>){
						chomp;
						my($feature,undef,undef,$fc,undef,$fdr)=split(/\t/);
						#$feature=~s/\"//g;
						#print STDERR "$feature\t$fc\t$fdr\n";
						
						$universe_noiseq_data->{$feature}++;
						#filtering over-expressed
						if($fdr>=$noiseq_cutoff and $fc>=$fc_threshold){
							$up_noiseq_data->{$feature}++;
							$up_number++;
							#last if scalar(keys %$up_edger_data)==5;
						}
						#filtering down-expressed
						if($fdr>=$noiseq_cutoff and $fc <($fc_threshold * -1)){
							$down_noiseq_data->{$feature}++;
							$down_number++;
						}
						$universe_number++;
					}
					close NOISEQ;
					print STDERR "Number of genes at universe: " . scalar(keys %$universe_noiseq_data) ."\n" if($verbose);
					print STDERR "Number of Up-regulated genes (NoiSeq cutoff=$noiseq_cutoff): " . scalar(keys %$up_noiseq_data)  ."\n" if($verbose);
					print STDERR "Number of Down-regulated genes (NoiSeq cutoff=$noiseq_cutoff): " . scalar(keys %$down_noiseq_data)  ."\n" if($verbose);
					
					#Once all fila are read
					if($universe_number>0 and ($up_number>0 or $down_number>0)){
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
							label=>$label,
							miARmaPath=>$miARmaPath
						);
					}
					else{
						warn date()." WARN :: The number of selected proteins is 0. No functional analysis can be done\n". date()." Upregulated genes n=$up_number, Downregulated genes n=$down_number\n".date()." Please check $logfile\n\n";
						return();
					}
				}
			}
			#return();
		}
		else{
			if($check ==0){
				print STDERR date() ." WARN :: No files provided for functional anaysis, please check your edgeR folder ($edger_dir) and/or your NoiSeq folder ($noiseq_dir)\n";
				return();
			}
		}
		print date() . " Functional Analysis finished.\n";
		return();
	}
	else{
		#If mandatory parameters have not been provided program will die and show error message
		warn(date()." F_Analysis ERROR :: projectdir($projectdir), organism($organism), logfile ($logfile), seq_id ($seq_id), noiseq_dir ($noiseq_dir) or edgeR_dir ($edger_dir) have not been provided");
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
	my $label=$args{"label"}; #label for saving files
	my $miARmaPath=$args{"miARmaPath"};
	
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

		#Creating results directory
		my $command="mkdir -p ".$output_dir;
		system($command) == 0
		   or die "F_Analysis ERROR :: ".date()." System args failed: $? ($command)";

		
		#save entities to read with R		
		my $up_file="$output_dir/.up_entities_$method.txt";
		open (UP,">$up_file") || die "$! $up_file";
		print UP join("\n", keys (%$up));
		close UP;
		
		my $down_file="$output_dir/.down_entities_$method.txt";;
		(open DW,">$down_file") || die "$! $down_file";
		print DW join("\n", keys (%$down));
		close DW;
		
		my $universe_file="$output_dir/.universe_entities_$method.txt";
		open (UN,">$universe_file")  || die "$! $universe_file";
		print UN join("\n", keys (%$universe));
		close UN;
		
		my $Rcommand="projectdir=\"$output_dir\",up=\"$up_file\",down=\"$down_file\",universe=\"$universe_file\",organism=\"$organism\",method=\"$method\",seq_id=\"$seq_id\",label=\"$label\"";
		# Printing the date and command execution on screen
		print STDOUT "GOSeqR :: ".date()." Starting Functional Analysis with $Rcommand\n" if($verbose);

		
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
		#source("/Users/eandres/Proyectos/EduardoAndres/miARma/New/lib/CbBio/RNASeq/R-Scripts/F_Analysis.R")

		my $cmds = <<EOF;
		source("$miARmaPath/lib/CbBio/RNASeq/R-Scripts/F_Analysis.R")
		setwd("$output_dir")
		resultsfiles<-NA
		resultsfiles<-F_Analysis($Rcommand)
EOF
		# for testing: 	source("/Users/eandres/Proyectos/miARma/lib/CbBio/RNASeq/R-Scripts/F_Analysis.R")	
		# 		source("http://valkyrie.us.es/CbBio/RNASeq/R-Scripts/F_Analysis.R")

		#Printing the execution data on log file and on the screen if verbose parameter is defined 
		open (LOG,">> ".$logfile) || die "F_Analysis ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG date()." Executing $cmds\n";
		#R commands execution
		my $out2 = $R->run($cmds);
		#Obtaining the files generated in the analysis
		my @media = $R->get('resultsfiles');
				
		foreach my $result (@media){
			if($result eq "NA"){
				if (scalar(keys %$up)<1 or scalar(keys %$down)<1 or scalar(keys %$universe)<1){
					warn "WARN :: ".date()." The number of selected proteins is 0. No functional analysis can be done\nWARN :: ".date()." Upregulated genes n=".scalar(keys %$up).", Downregulated genes n=".scalar(keys %$down)."\nPlease check $logfile\n";
					return();
				}
				if(scalar(keys %$up)<51 or scalar(keys %$down)<51){
					print STDERR date(). " WARN :: The number of differentially expressed genes is too small (<50). Some errors could appear. Try to use a less restringet cut_off for $method(>$cut_off)\n";
				}
				print STDERR date(). " ERROR :: Nothing significative was found using $method data. No Functional analisis result files were generated\n\n";
			}
			else{
				foreach my $res_file (@{$result}){
					print "F_Analysis :: ".date()." The file ".$res_file." has been generated \n" if($verbose);
				}
			}
		}
		close LOG;
		return();
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

sub F_AnalysisSummary{
	use File::Basename;
	#Arguments provided by user are collected by %args. Dir, file, aligner, statsfile, projectdir 
	#and logfile are mandatory arguments while verbose and threads are optional.
	my %args=@_;
	my $summary_file=$args{"summary"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}; #Optional arguments to show the execution data on screen

	my $data;
	my $category->{"CC"}="GO Cellular Component";
	$category->{"BP"}="GO Biological Process";
	$category->{"MF"}="GO Molecular Function";
	$category->{"KEGG"}="KEGG Pathways";
	
	opendir(GOSEQ, $projectdir ."/Functional_Analysis_results/") || warn "F_Analysis:: Folder \"$projectdir/Functional_Analysis_results/\" is not found\n"; 
	my @goseq_files= readdir(GOSEQ);
	foreach my $goseq_file (sort @goseq_files){
		if($goseq_file =~ /\.xls$/){	
			open(GOSEQFILE,$projectdir ."/Functional_Analysis_results/$goseq_file") || die "$! $projectdir/Functional_Analysis_results/$goseq_file";
			my $pval_threshold=0.05;
			while(<GOSEQFILE>){
				chomp;
				$_=~s/\"//g;
				if($_ !~ /category/){
					
					
					my($goterm,$pval_over,$pval_under,undef,undef,$term,$ontology)=split(/\t/);
					if($pval_under<=$pval_threshold){
						$data->{$goseq_file}->{$ontology}->{DOWN}->{$goterm}=$term;
					}
					elsif($pval_over<=$pval_threshold){
						$data->{$goseq_file}->{$ontology}->{UP}->{$goterm}=$term;
					}
					else{
						$data->{$goseq_file}->{"Z"}->{DOWN}->{""}=1;
						$data->{$goseq_file}->{"Z"}->{UP}->{""}=1;
					}
				}
			}
			close GOSEQFILE;
		}
	}
	close GOSEQ;
	
	if(scalar(keys %$data)>0){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
		print SUMM "\nFunctional Analysis by GoSeq [".$projectdir ."/Functional_Analysis_results/]\n";
		print SUMM "File\tNumber of Over Represented Terms (Pval<0.05)\tNumber of Under Represented Terms (Pval<0.05)\tOntology\n";
		foreach my $file (sort keys %$data){
			my $printed=0;
			foreach my $ontology (sort keys %{$data->{$file}}){
				if($ontology ne "Z"){
					print SUMM "$file\t". scalar(keys %{$data->{$file}->{$ontology}->{UP}}) ."\t" . scalar(keys %{$data->{$file}->{$ontology}->{DOWN}}) ."\t". $category->{$ontology}."\n";
					#print SUMM join("\t",keys %{$data->{$file}->{$ontology}->{UP}}) ."\n";
					$printed=1;
				}
				else{
					print SUMM "$file\t0\t0\t-\n" if($printed==0);
				}
			}
		}
		close SUMM;
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

