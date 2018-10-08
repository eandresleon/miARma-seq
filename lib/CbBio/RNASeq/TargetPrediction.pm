
#########################################################################	
#	TargetPrediction processing package		 									#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2015 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################

package CbBio::RNASeq::TargetPrediction;
#Export package system
require Exporter;
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(TargetPrediction TargetPrediction_Summary correlation);

use strict;
use LWP;
use File::Basename;
use Cwd;
use Cwd 'abs_path';

=head1 NAME

 TargetPrediction 

=head1 SYNOPSIS

TargetPrediction package is composed of 9 subroutines: CutTargetPrediction,CutTargetPredictionStats, TargetPredictionTriming, ReadFilter,
Minion, Reaper, TargetPredictionerRemoval, TargetPredictionerGraph and ReadLengthCount. The aim of this package is to 
process the reads from a fastq file removing a known TargetPredictioner sequence (CutTargetPrediction, Reaper) or not 
(TargetPredictionTriming), predicting the TargetPredictioner sequence of the reads (Minion), filtering the reads by length 
(ReadFilter) or printing statistical info after the read processing (CutTargetPredictionStats, TargetPredictionerGraph, 
ReadLengthCount). In addition TargetPredictionerRemoval acts as main function calling to CutTargetPrediction, CutTargetPredictionStats,  
TargetPredictionTriming, Minion, Reaper, TargetPredictionerGraph and ReadLengthCount according to the parameters selected.

=head1 Methods


=head2 TargetPredictionerRemoval

  Example    : 
  TargetPredictionerRemoval(
  	dir=>"../reads",
  	TargetPredictionersoft=>"CutTargetPrediction-Reaper",
  	projectdir=>".",
  	file=>"file.fastq",
  	verbose=>"verbose",
  	logfile=>"run.log",
  	statsfile=>"stats.log",
  	organism=>"human",
  	min=>"18", 
  	max=>"26",
  	min_quality=>"1", 
  	cutTargetPredictionparameters=> "-O 4",
  	reaperparameters=>" -3p-prefix 12/2/0/0 -dust-suffix-late 20 --nozip"
  	);
  Description: TargetPredictionerRemoval function performs an TargetPredictioner removal process with cutTargetPrediction, reaper software or both, and/or
  the removal of a defined number of nucleotides of the reads (TargetPredictiontrimming). These analysis generate new fastq files with 
  the conserved reads at CutTargetPrediction_results, Reaper_results or/and TargetPredictionTrimming_results directory in the provided project directory. 
  In addition, if user do not provide the TargetPredictioner sequence TargetPredictionerRemoval can predict it executing Minion software. 
  After processing the reads TargetPredictionerRemoval generates statistical data of the process
    Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where the input files are 
   	 
  	 Optional parameters:
  	 [TargetPredictioner] TargetPredictioner sequence to be removed in the analysis (if this sequence is not provided Minion software will be used to predict one)

  Returntype : Fastq files and plots with reads distribution at directory CutTargetPrediction_results, Reaper_results and/or TargetPredictionTrimming_results. 
  TargetPredictionerRemoval also return the paths of the new generated files
  Requeriments: cutTargetPrediction function requires for a correct analysis:	my @files= readdir DIR;

  	- Perl v5.10.0 or higher software correctly installed
  	- CutTargetPrediction v1.2.1 or higher software correctly installed
  	- Reaper software correctly installed
  	- Minion software correctly installed
  	- R v3.1.0 or higher software correctly installed 
  	- Input files on fasta/fastq format (compressed files are accepted)
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub TargetPrediction{
	#Arguments provided by user are collected by %args. 
	my %args=@_;
	my $dir= $args{"dir"}; #Input directory where the files are
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $organism=$args{"organism"};
	my $fc_threshold=$args{"fc_threshold"};
	my $method=$args{"method"};
	
	my $miARmaPath=$args{"miARmaPath"};
	use lib "$miARmaPath/lib/Perl";
	use DateTime;
	
	my $edgeR_cutoff;
	if(defined($args{"edger_cutoff"})){
		$edgeR_cutoff=$args{"edger_cutoff"};
	}
	else{
		$edgeR_cutoff=0.05;
	}

	my $NoiSeq_cutoff;
	if(defined($args{"noiseq_cutoff"})){
		$NoiSeq_cutoff=$args{"noiseq_cutoff"};
	}
	else{
		$NoiSeq_cutoff=0.8;
	}
	
	#miRNA target prediction
	if(defined($args{"miRNAs_folder"}) and !defined($args{"genes_folder"})){
		print STDOUT "LOG :: " . date() . " Starting a miRNA target prediction from ".$args{"miRNAs_folder"}."\n" if($verbose);
		opendir(DIR, $args{"miRNAs_folder"}) || die $!;
		my @files= readdir(DIR);
		system("mkdir -p $projectdir/miRGate_results/");
		foreach my $file (@files){
			if(-d $args{"miRNAs_folder"} ."/$file"){
				#check if data comes from edgeR or Noiseq
				if(lc($file) =~ /edger/){
					opendir(DIR, $args{"miRNAs_folder"} ."/$file") || die "$! : ".$args{"miRNAs_folder"} ."/$file";
					my @files= readdir(DIR);
					@files=map($args{"miRNAs_folder"} ."/$file/$_",@files);
					foreach my $edger_files(@files){
						if($edger_files =~ /\.xls$/){							
							my $output_file=fileparse($edger_files);
							$output_file=~s/\.xls/_miRNAs_miRGate.xls/g;							
							my $data=get_edgeR_data(
							    file=>$edger_files,
							    edger_cutoff=>$edgeR_cutoff,
							    verbose=>$verbose,
								fc_threshold=>$fc_threshold,
								
							);							
							print STDOUT "LOG :: " . date() . " Searching targets from edgeR file [$edger_files] for " . scalar(keys %$data) ." miRNAs using miRGate\n" if($verbose);
							 miRGate(
							 	miRNAs=>$data,
							 	organism=>$organism,
							 	verbose=>$verbose,
								output=>$projectdir ."/miRGate_results/". $output_file,
								miARmaPath=>$miARmaPath,
							);
						}
					}
				}
				if(lc($file) =~ /noiseq/){
					opendir(DIR, $args{"miRNAs_folder"} ."/$file") || die "$! : ".$args{"miRNAs_folder"} ."/$file";
					my @files= readdir(DIR);
					@files=map($args{"miRNAs_folder"} ."/$file/$_",@files);
					foreach my $noise_files(@files){
						if($noise_files =~ /\.xls$/){
							my $output_file=fileparse($noise_files);
							$output_file=~s/\.xls/_miRNAs_miRGate.xls/g;
							my $data=get_NoiSeq_data(
					            file=>$noise_files,
					            NoiSeq_cutoff=>$NoiSeq_cutoff,
					            verbose=>$verbose,
								fc_threshold=>$fc_threshold,
							);
						print STDOUT "LOG :: " . date() . " Searching targets from NoiSeq file [$noise_files] for " . scalar(keys %$data) ." miRNAs using miRGate\n" if($verbose);
							 miRGate(
								miRNAs=>$data,
								organism=>$organism,
								verbose=>$verbose,
								output=>$projectdir . "/miRGate_results/". $output_file,
								miARmaPath=>$miARmaPath,
							);
						}
					}
				}
			}
		}
	}
	#Gene target prediction
	if(defined($args{"genes_folder"}) and !defined($args{"miRNAs_folder"}) ){
		print STDOUT "LOG :: " . date() . " Starting a gene target prediction from ".$args{"genes_folder"}."\n" if($verbose);
		opendir(DIR, $args{"genes_folder"}) || die $!;
		my @files= readdir(DIR);
		system("mkdir -p $projectdir/miRGate_results/");
		foreach my $file (@files){
			if(-d $args{"genes_folder"} ."/$file"){
				#check if data comes from edgeR or Noiseq
				if(lc($file) =~ /edger/){
					opendir(DIR, $args{"genes_folder"} ."/$file") || die "$! : ".$args{"genes_folder"} ."/$file";
					my @files= readdir(DIR);
					@files=map($args{"genes_folder"} ."/$file/$_",@files);
					foreach my $edger_files(@files){
						if($edger_files =~ /\.xls$/){
							my $output_file=fileparse($edger_files);
							$output_file=~s/\.xls/_genes_miRGate.xls/g;
							my $data=get_edgeR_data(
							    file=>$edger_files,
							    edger_cutoff=>$edgeR_cutoff,
							    verbose=>$verbose,
								fc_threshold=>$fc_threshold,
								
							);
							print STDOUT "LOG :: " . date() . " Searching targets from edgeR file [$edger_files] for " . scalar(keys %$data) ." genes using miRGate\n" if($verbose);
							 miRGate(
							 	UTRs=>$data,
							 	organism=>$organism,
							 	verbose=>$verbose,
								output=>$projectdir. "/miRGate_results/". $output_file
							);
						}
					}
				}
				if(lc($file) =~ /noiseq/){
					opendir(DIR, $args{"genes_folder"} ."/$file") || die "$! : ".$args{"genes_folder"} ."/$file";
					my @files= readdir(DIR);
					@files=map($args{"genes_folder"} ."/$file/$_",@files);
					foreach my $noise_files(@files){
						if($noise_files =~ /\.xls$/){
							my $output_file=fileparse($noise_files);
							$output_file=~s/\.xls/_genes_miRGate.xls/g;
							my $data=get_NoiSeq_data(
					            file=>$noise_files,
					            NoiSeq_cutoff=>$NoiSeq_cutoff,
					            verbose=>$verbose,
								fc_threshold=>$fc_threshold,
								
							);
						print STDOUT "LOG :: " . date() . " Searching targets from NoiSeq file [$noise_files] for " . scalar(keys %$data) ." genes using miRGate\n" if($verbose);
							 miRGate(
								UTRs=>$data,
								organism=>$organism,
								verbose=>$verbose,
								output=>$projectdir ."/miRGate_results/". $output_file
							);
						}
					}
				}
			}
		}
	}
	#miRNA-gene target prediction
	if(defined($args{"miRNAs_folder"}) and defined($args{"genes_folder"})){
		
		my $miRNA_edgeR_results;
		my $miRNA_noiseq_results;
		my $genes_edgeR_results;
		my $genes_noiseq_results;
		
		my $miRNA_edgeR_RPKMs;
		my $miRNA_noiseq_RPKMs;
		my $genes_edgeR_RPKMs;
		my $genes_noiseq_RPKMs;
		
		####### miRNAs ########		
		print STDOUT "LOG :: " . date() . " Starting a miRNA-gene target prediction\n" if($verbose);
		opendir(DIR, $args{"miRNAs_folder"}) || die $!;
		my @files= readdir(DIR);
		system("mkdir -p $projectdir/miRGate_results/");

		foreach my $file (@files){
			if(-d $args{"miRNAs_folder"} ."/$file"){
				#check if data comes from edgeR or Noiseq
				if(lc($file) =~ /edger/){
					opendir(DIR, $args{"miRNAs_folder"} ."/$file") || die "$! : ".$args{"miRNAs_folder"} ."/$file";
					my @files= readdir(DIR);
					@files=map($args{"miRNAs_folder"} ."/$file/$_",@files);
					
					foreach my $edger_files(@files){
						if($edger_files =~ /\.xls$/){
							my $output_file=fileparse($edger_files);
							$output_file=~s/\.xls/_miRNAs_miRGate.xls/g;
							$miRNA_edgeR_results->{$edger_files}=get_edgeR_data(
							    file=>$edger_files,
							    edger_cutoff=>$edgeR_cutoff,
							    verbose=>$verbose,
								fc_threshold=>$fc_threshold,
							);
						}
					}
					foreach my $edger_files(@files){
						if($edger_files =~ /RPKM\.tsv$/){
							$miRNA_edgeR_RPKMs=get_RPKM(
								file=>$edger_files
							);
						}
					}
				}
				if(lc($file) =~ /noiseq/){
					opendir(DIR, $args{"miRNAs_folder"} ."/$file") || die "$! : ".$args{"miRNAs_folder"} ."/$file";
					my @files= readdir(DIR);
					@files=map($args{"miRNAs_folder"} ."/$file/$_",@files);
					
					foreach my $noise_files(@files){
						if($noise_files =~ /\.xls$/){
							$miRNA_noiseq_results->{$noise_files}=get_NoiSeq_data(
								file=>$noise_files,
								NoiSeq_cutoff=>$NoiSeq_cutoff,
								verbose=>$verbose,
								fc_threshold=>$fc_threshold,
								
							);
						}
					}
					foreach my $noise_files(@files){
						if($noise_files =~ /RPKM\.tsv$/){
							$miRNA_noiseq_RPKMs=get_RPKM(
								file=>$noise_files
							);
						}
					}
				}
			}
		}
		
		####### genes ######
		
		opendir(DIR, $args{"genes_folder"}) || die $!;
		my @files= readdir(DIR);
		foreach my $file (@files){
			if(-d $args{"genes_folder"} ."/$file"){
				#check if data comes from edgeR or Noiseq
				if(lc($file) =~ /edger/){
					opendir(DIR, $args{"genes_folder"} ."/$file") || die "$! : ".$args{"genes_folder"} ."/$file";
					my @files= readdir(DIR);
					@files=map($args{"genes_folder"} ."/$file/$_",@files);
					foreach my $edger_files(@files){
						if($edger_files =~ /\.xls$/){
							$genes_edgeR_results->{$edger_files}=get_edgeR_data(
							    file=>$edger_files,
							    edger_cutoff=>$edgeR_cutoff,
							    verbose=>$verbose,
								fc_threshold=>$fc_threshold,
							);
						}
					}
					foreach my $edger_files(@files){
						if($edger_files =~ /RPKM\.tsv$/){
							$genes_edgeR_RPKMs=get_RPKM(
								file=>$edger_files
							)	
						}
					}
				}
				if(lc($file) =~ /noiseq/){
					opendir(DIR, $args{"genes_folder"} ."/$file") || die "$! : ".$args{"genes_folder"} ."/$file";
					my @files= readdir(DIR);
					@files=map($args{"genes_folder"} ."/$file/$_",@files);
					foreach my $noise_files(@files){
						if($noise_files =~ /\.xls$/){
							$genes_noiseq_results->{$noise_files}=get_NoiSeq_data(
					            file=>$noise_files,
					            NoiSeq_cutoff=>$NoiSeq_cutoff,
					            verbose=>$verbose,
								fc_threshold=>$fc_threshold,
							);
						}
					}
					foreach my $noise_files(@files){
						if($noise_files =~ /RPKM\.tsv$/){
							$genes_noiseq_RPKMs=get_RPKM(
								file=>$noise_files
							);
						}
					}
				}
			}
		}

		###### Both in miRGate ########
		
		#edgeR
		my $all_possible_pairs;
		
		foreach my $file2 (keys %$miRNA_edgeR_results){
			foreach my $m (keys %{$miRNA_edgeR_results->{$file2}}){
				foreach my $file1 (keys %$genes_edgeR_results){
					foreach my $g (keys %{$genes_edgeR_results->{$file1}}){
						$all_possible_pairs->{$m}->{$g}++;
					}
				}
			}
		}

		foreach my $gene_edger_data (sort keys %$genes_edgeR_results){
			my $g_o_edger=fileparse($gene_edger_data);
			$g_o_edger=~s/\.xls$/_genes/;
			foreach my $miRNA_files (sort keys %$miRNA_edgeR_results){
				my $m_o_edger=fileparse($miRNA_files);
				$m_o_edger=~s/\.xls$/_miRNAs/;
				print STDOUT "LOG :: ".date()." Gathering Prediction from ".scalar(keys %{$miRNA_edgeR_results->{$miRNA_files}})." miRNAs in " . fileparse($miRNA_files) ." and ".scalar(keys %{$genes_edgeR_results->{$gene_edger_data}})." genes from " . fileparse($gene_edger_data) . " using miRGate\n" if($verbose);
				
				 my $pares=miRGate(
				 	UTRs=>$genes_edgeR_results->{$gene_edger_data},
				 	miRNAs=>$miRNA_edgeR_results->{$miRNA_files},
				 	organism=>$organism,
				 	verbose=>$verbose,
				 	output=>$projectdir ."/miRGate_results/". $g_o_edger ."_" . $m_o_edger .".xls",
				 );

				 my $stat_corr=correlation(
			 		pares=>$pares,
					all_pares=>$all_possible_pairs,
					UTRs_rpkm=>$genes_edgeR_RPKMs,
					miRNAs_rpkm=>$miRNA_edgeR_RPKMs,
				 	organism=>$organism,
				 	verbose=>$verbose,
				 	output=>$projectdir ."/miRGate_results/". $g_o_edger,
					miARmaPath=>$miARmaPath,
					method=>$method,
					DE_type=>"edgeR",
					logfile=>$logfile
				 );
				 Integration(
			 		pares=>$pares,
				 	corr_pares=>$stat_corr,
				 	verbose=>$verbose,
				 	output=>$projectdir ."/miRGate_results/". $g_o_edger,
					logfile=>$logfile					
				 );
			}
		}
		#Noiseq
		
		my $all_possible_pairs;
		
		foreach my $file2 (keys %{$miRNA_noiseq_results}){
			foreach my $m (keys %{$miRNA_noiseq_results->{$file2}}){
				foreach my $file1 (keys %{$genes_noiseq_results}){
					foreach my $g (keys %{$genes_noiseq_results->{$file1}}){
						$all_possible_pairs->{$m}->{$g}++;
					}
				}
			}
		}
		
		foreach my $gene_noiseq_data (sort keys %$genes_noiseq_results){
			my $g_o_noiseq=fileparse($gene_noiseq_data);
			$g_o_noiseq=~s/\.xls$/_genes/;
			foreach my $miRNA_files (sort keys %$miRNA_noiseq_results){
				my $m_o_noiseq=fileparse($miRNA_files);
				$m_o_noiseq=~s/\.xls$/_miRNAs/;
				print STDOUT "LOG :: ".date()." Gathering Prediction from ".scalar(keys %{$miRNA_noiseq_results->{$miRNA_files}})." miRNAs in " . fileparse($miRNA_files) ." and ".scalar(keys %{$genes_noiseq_results->{$gene_noiseq_data}})." genes from " . fileparse($gene_noiseq_data) . " using miRGate\n" if($verbose);
				 my $pares=miRGate(
				 	UTRs=>$genes_noiseq_results->{$gene_noiseq_data},
				 	miRNAs=>$miRNA_noiseq_results->{$miRNA_files},
				 	organism=>$organism,
				 	verbose=>$verbose,
				 	output=>$projectdir ."/miRGate_results/". $g_o_noiseq ."_" . $m_o_noiseq . ".xls",
				 );
				 my $stat_corr=correlation(
				 	pares=>$pares,
					all_pares=>$all_possible_pairs,
					UTRs_rpkm=>$genes_noiseq_RPKMs,
					miRNAs_rpkm=>$miRNA_noiseq_RPKMs,
				 	organism=>$organism,
				 	verbose=>$verbose,
				 	output=>$projectdir ."/miRGate_results/". $g_o_noiseq,
					miARmaPath=>$miARmaPath,
					method=>$method,
					DE_type=>"NOISeq",
					logfile=>$logfile
				 );
				 Integration(
			 		pares=>$pares,
				 	corr_pares=>$stat_corr,
				 	verbose=>$verbose,
				 	output=>$projectdir ."/miRGate_results/". $g_o_noiseq,
					logfile=>$logfile					
				 );
			}
		}
	}
	
	sub help_TargetPredictionerRemoval{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
  	 		[dir] Input directory where the input files are 
  	 
  	 		Optional parameters:
  	 		[TargetPredictioner] TargetPredictioner sequence to be removed in the analysis (if this sequence is not provided Minion software will be used to predict one)
  	 		               
			Examples:
			TargetPredictionerRemoval(dir=>"../reads", TargetPredictionersoft=>"CutTargetPrediction-Reaper", projectdir=>".", file=>"file.fastq", verbose=>"verbose", logfile=>"run.log", statsfile=>"stats.log", organism=>"human", min=>"18", max=>"26", min_quality=>"1", cutTargetPredictionparameters=> "-O 4", reaperparameters=>" -3p-prefix 12/2/0/0 -dust-suffix-late 20 --nozip");
	};

	print STDERR $usage;
	exit(); 
	}
	return();
}
sub get_edgeR_data{
	my %args=@_;
	my $file= $args{"file"}; #Input directory where the files are
	my $edger_cutoff=$args{"edger_cutoff"};
	my $verbose=$args{"verbose"};
	my $organism=$args{"organism"};
	my $fc_threshold=$args{"fc_threshold"};

	my $data;
	print STDERR "LOG :: " . date() . " Reading $file (cut_off $edger_cutoff and fc_threshold $fc_threshold)\n" if($verbose);
	
	#Read edgeR data
	open(EDGER,$file) || warn "cant find $file\n";
	while(<EDGER>){
		chomp;
		if($_ !~/logCPM/){
			$_=~s/\"//g;
			my @data=split(/\t/);
			my $feature=$data[0];
			my $fc=$data[($#data)-3];
			my $fdr=$data[$#data];
			#filtering over-expressed
			if($fdr<=$edger_cutoff and $fc>=$fc_threshold){
				$data->{$feature}="$fc\t$fdr";
			}
			#filtering down-expressed
			if($fdr<=$edger_cutoff and $fc<($fc_threshold * -1)){
				$data->{$feature}="$fc\t$fdr";
			}
		}
	}
	close EDGER;
	return($data);
}	
sub get_NoiSeq_data{
	my %args=@_;
	my $file= $args{"file"}; #Input directory where the files are
	my $noiseq_cutoff=$args{"NoiSeq_cutoff"};
	my $verbose=$args{"verbose"};
	my $fc_threshold=$args{"fc_threshold"};
	my $data;
	print STDOUT "LOG :: " . date() . " Reading $file (cut_off $noiseq_cutoff and fc_threshold $fc_threshold)\n" if($verbose);
	#Read edgeR data
	open(NOISEQ,$file) || warn "cant find $file\n";
	while(<NOISEQ>){
		chomp;
		if($_ !~/ranking/){
			$_=~s/\"//g;
			my($feature,undef,undef,$fc,undef,$fdr)=split(/\t/);
			#filtering over-expressed
			if($fdr>=$noiseq_cutoff and $fc>=$fc_threshold){
				$data->{$feature}="$fc\t$fdr";
			}
			#filtering down-expressed
			if($fdr>=$noiseq_cutoff and $fc<($fc_threshold * -1 )){
				$data->{$feature}="$fc\t$fdr";
			}
		}
	}
	close NOISEQ;
	return($data);
}
sub miRGate{

	my %args=@_;
	my $miRNAs= $args{"miRNAs"};
	my $UTRs= $args{"UTRs"};
	my $organism=$args{"organism"};
	my $output_file=$args{"output"};
	my $verbose=$args{"verbose"};

	my $pares;
	
	my $miARmaPath=$args{"miARmaPath"};
	use lib "$miARmaPath/lib/Perl";
	use LWP::Simple;
	my $parser;
	#print STDERR "REcibidos un total de " . scalar(keys %$miRNAs) ." miRNAs y " . scalar(keys %$UTRs) ." genes y los guard en $output_file\n";
	open (OUT,">$output_file") || warn "Cant create output ($output_file)\n";
	print OUT "#miRNA\tmiRNA FC\tmiRNA FDR\tEnsembl Gene\tGene Symbol\tEnsembl Transcript\tGene FC\tGene FDR\tMethod\tTarget Site\tScore\tEnergy\n";
	my $correct=0;
	if($miRNAs and !$UTRs){
		$correct=1;
		my $source="miRNA_predictions";
		my $cont=0;
		foreach my $miRNA (sort keys %$miRNAs){
			$cont++;
			my $result;
			my $SourceListingURL = "http://mirgate.bioinfo.cnio.es/ResT/miARma/$organism/$source/$miRNA";
			sleep 5;
			my $xmldoc = get($SourceListingURL);
			if($xmldoc){
				my @lines=split(/\n/,$xmldoc);
				foreach my $line (@lines) {
					chomp;
					my($id,$found_miRNA,$ensembl,$hgnc,$utr,$method,$target_site,$score,$energy)=split(/\t/,$line);
					my $fc;
					if(exists $miRNAs->{$found_miRNA}){
						$fc=$miRNAs->{$found_miRNA};
					}
					else{
						$fc=$miRNAs->{$found_miRNA};
					}
					if(lc($miRNA) eq lc($found_miRNA)){
						print OUT $miRNA ."\t" . $fc ."\t". $ensembl ."\t" . $hgnc ."\t" . $utr ."\tNA\tNA";
						print OUT "\t" . $method ."\t" . $target_site ."\t" . $score ."\t" . $energy ."\n";
					}
				}
			}
		}
		
	}
	if($UTRs and !$miRNAs){
		$correct=1;
		my $source="gene_predictions";
		my $cont=0;
		foreach my $gene (sort keys %$UTRs){
			$cont++;
			my $result;
			my $SourceListingURL = "http://mirgate.bioinfo.cnio.es/ResT/miARma/$organism/$source/$gene";
			sleep 5;
			my $xmldoc = get($SourceListingURL);
			if($xmldoc){
				my @lines=split(/\n/,$xmldoc);
				foreach my $line (@lines) {
					chomp;
					my($id,$found_miRNA,$ensembl,$hgnc,$utr,$method,$target_site,$score,$energy)=split(/\t/,$line);
					my $fc;
					if(exists $UTRs->{$ensembl}){
						$fc=$UTRs->{$ensembl};
					}
					else{
						$fc=$UTRs->{$ensembl};
					}
					
					print OUT $found_miRNA ."\tNA\tNA\t" . $ensembl ."\t" . $hgnc ."\t" . $utr ."\t". $fc;
					print OUT "\t" . $method ."\t" . $target_site ."\t" . $score ."\t" . $energy ."\n";
				}
			}
			else{
				sleep 5;
				next;
			}
		}
	}
	if($miRNAs and $UTRs){
		$correct=1;
		my $source="miRNA_predictions";
		foreach my $miRNA (keys %$miRNAs){
			my $result;
			my $SourceListingURL = "http://mirgate.bioinfo.cnio.es/ResT/miARma/$organism/$source/$miRNA";
			sleep 5;
			my $xmldoc = get($SourceListingURL);
			if($xmldoc){
				my @lines=split(/\n/,$xmldoc);
				foreach my $line (@lines) {
					chomp;
					my($id,$found_miRNA,$ensembl,$hgnc,$utr,$method,$target_site,$score,$energy)=split(/\t/,$line);
					my $fc;
					if(exists $UTRs->{$ensembl}){
						$fc=$UTRs->{$ensembl};
					}
					elsif(exists $UTRs->{$hgnc}){
						$fc=$UTRs->{$hgnc};
					}
					elsif(exists $UTRs->{$hgnc}){
						$fc=$UTRs->{$hgnc};
					}
					else{
						$fc=$UTRs->{$ensembl};
					}
					if(lc($miRNA) eq lc($found_miRNA)){
						$result->{$id}->{miRNA}=$found_miRNA;
						$result->{$id}->{EnsEMBL}=$ensembl;
						$result->{$id}->{HGNC}=$hgnc;
						$result->{$id}->{utr}=$utr;
						$result->{$id}->{method}=$method;
						$result->{$id}->{target_site}=$target_site;
						$result->{$id}->{score}=$score;
						$result->{$id}->{energy}=$energy;
					}
				}
			}
			else{
				sleep 2;
				next;
			}
			foreach my $id (keys %$result){
				foreach my $gene (keys %$UTRs){
					if(($gene eq $result->{$id}->{EnsEMBL}) or ($gene eq $result->{$id}->{HGNC}) or ($gene eq $result->{$id}->{utr})){
						if(exists ($UTRs->{$gene}) and $miRNAs->{$miRNA}){
							my ($miRFC,$miRFDR)=split(/\t/,$miRNAs->{$miRNA});
							my ($gnFC,$gnFDR)=split(/\t/,$UTRs->{$gene});
							
							
							$pares->{$miRNA}->{$result->{$id}->{HGNC}}->{"$gnFC\t$miRFC"}->{$result->{$id}->{method}}=$result->{$id}->{target_site};
							print OUT $miRNA ."\t" . $miRNAs->{$miRNA} ."\t" . $result->{$id}->{EnsEMBL} ."\t" . $result->{$id}->{HGNC} ."\t" . $result->{$id}->{utr} ."\t" . $UTRs->{$gene};
							print OUT "\t" . $result->{$id}->{method} ."\t" . $result->{$id}->{target_site} ."\t" . $result->{$id}->{score} ."\t".$result->{$id}->{energy} ."\n";
						}
					}
				}
			}
		}
	}
	
	close OUT;
	if($correct == 0){
		print STDOUT "\t". date(). " miRGate:: No data for miRGate found in [ $miRNAs $UTRs ]\n" if($verbose);
		return();
	}
	
	return($pares);
}
sub TargetPrediction_Summary{
	my %args=@_;
	my $summary_file=$args{"summary"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	
	my $dirname = $projectdir ."/miRGate_results/";
	opendir ( DIR, $dirname ) || die "Error in opening dir $dirname\n";
	my @files= readdir DIR;
	
	my $miRNAs;
	my $genes;
		
	foreach my $filename (sort @files){		     
		if($filename =~ /\.xls$/){
			#Filtering for xls (miRgate output format)

			my $connections=0;
			#Parding filename to gather the tumor name
			open(RESULTS, $dirname ."/". $filename) || die $!;	
				
			while(<RESULTS>){
				chomp;
				if($_ !~ /FDR/){
					my @data=split("\t");
			
					#Variable asignment
					my $miRNA=$data[0];
					my $miRNA_FC=$data[1];
					my $miRNA_FDR=$data[2];
					my $symbol=$data[4];
					my $transcript=$data[5];
					my $transcript_FC=$data[6];
					my $transcript_FDR=$data[7];
					my $methods=$data[8];

					$miRNAs->{$filename}->{$miRNA}++;
					$genes->{$filename}->{$symbol}++;
				}
			}
			close RESULTS;
		}
	}
	my $cont=0;
	
	if(scalar(keys %$miRNAs)>0 and scalar(keys %$genes)>0){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n\n";
		print SUMM "\nmiRNA-mRNA Target Predictions by miRGate [$dirname]\n";
		#print SUMM "File\tNumber of DE elements (Pval <=0.05)\tNumber of DE elements (FDR <=0.05)\n";
	
		print SUMM "\nmiRNAs with more associations\n";
		print SUMM "File\tmiRNA\tNumber of associations\n";
		foreach my $file (sort keys %$miRNAs){
			my $cont=0;
			foreach my $miRNA (sort {$miRNAs->{$file}->{$b} <=> $miRNAs->{$file}->{$a}} keys %{$miRNAs->{$file}}){
				last if($cont==5);
				$cont++;
				print SUMM "$file\t$miRNA\t" . $miRNAs->{$file}->{$miRNA} ."\n";
			}
			print SUMM "\n";
		}
	
		print SUMM "\nGenes more regulated\n";
		print SUMM "File\tGeneName\tNumber of associations\n";
		my $cont=0;
	
		foreach my $file (sort keys %$genes){
			############# Most Targeted genes ############
			my $cont=0;
			foreach my $GeneName (sort {$genes->{$file}->{$b} <=> $genes->{$file}->{$a}} keys %{$genes->{$file}}){
				last if($cont==5);
				$cont++;
				my @miRs;
				print SUMM "$file\t$GeneName\t" . $genes->{$file}->{$GeneName} ."\n";
			}
			print SUMM "\n";
		}
		close SUMM;
	}
	closedir(DIR);
	return();
}
sub get_RPKM{
	my %args=@_;
	my $file=$args{"file"};
	my $resultado;
	my $min;
	open(FILE,$file) || warn "Can't read $file\n";
	my $line=0;
	while(<FILE>){
		chomp;
		if($line !=0){
			$_=~ s/\"//g;
			my(@data)=split(/\t/);
			my $gene=shift @data;
			$resultado->{$gene}=join(",",@data);
			foreach my $rpkm(@data){
				$min->{$rpkm}++;
			}
		}
		$line++;
	}
	close FILE;
	my @mins = sort { $a <=> $b } keys %$min;
	$resultado->{MIN}=$mins[1];
	return($resultado);
}
sub correlation{

	my %args=@_;
	my $pares=$args{"pares"};
	my $all_pares=$args{"all_pares"};
	my $output_file=$args{"output"};
	my $gene_RPKM=$args{"UTRs_rpkm"};
	my $miRNA_RPKM=$args{"miRNAs_rpkm"};
	my $miARmaPath=$args{"miARmaPath"};
	my $method=$args{"method"}|| "Pearson";
	my $DE_type=$args{"DE_type"};
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	
	$output_file=~s/(.+)_\w+_\w+_results.*/$1/g;
	$output_file=$output_file."_Integrative_DifferentiallyExpressed_miRNA_mRNA_$DE_type\_$method\_correlation.xls";
	my ($tmp_file_tmp,$path)=fileparse($output_file);
	my $tmp_file="$path/.$tmp_file_tmp";

	print date()." Including an Integrative analysis based on ". ucfirst($method)." correlation ($DE_type)\n";
	
	open(TMP,">$tmp_file") || warn "Can't create integrative file $output_file\n";
	
	foreach my $miRNA (sort keys %{$pares}){
		foreach my $gene (keys %{$pares->{$miRNA}}){
			print TMP "$miRNA\t".$miRNA_RPKM->{$miRNA} ."\t$gene\t". $gene_RPKM->{$gene} ."\n";
		}
	}
	close TMP;
	
	my $R_cmd="R --no-save --silent --args $tmp_file ". lc($method) ." $output_file ". $miRNA_RPKM->{MIN} ." " . $gene_RPKM->{MIN} ." < $miARmaPath/lib/CbBio/RNASeq/R-Scripts/Cor.R >>$logfile";
	#Executing the command or if system can't be executed die showing the error.
	system($R_cmd) == 0
	    or die "TargetPrediction :: Can't create correlation file: $? ($R_cmd)";
	
		
	my $output_file2=$output_file;
	$output_file2=~s/_DifferentiallyExpressed//g;
	
	open(TMP,">$tmp_file") || warn "Can't create integrative file $output_file\n";
	foreach my $miRNA (sort keys %$all_pares){
		foreach my $gene (sort keys %{$all_pares->{$miRNA}}){
			print TMP "$miRNA\t".$miRNA_RPKM->{$miRNA} ."\t$gene\t". $gene_RPKM->{$gene} ."\n";
		}
	}
	
	close TMP;
	
	my $R_cmd="R --no-save --silent --args $tmp_file ". lc($method) ." $output_file2 ". $miRNA_RPKM->{MIN} ." " . $gene_RPKM->{MIN} ." < $miARmaPath/lib/CbBio/RNASeq/R-Scripts/Cor.R >>$logfile";
	#Executing the command or if system can't be executed die showing the error.
	system($R_cmd) == 0
	    or die "TargetPrediction :: Can't create correlation file: $? ($R_cmd)";
	
	open(COR,$output_file2) || die "Can't read $output_file2\n";
	my $line=1;
	my $int_pairs_corr;
	while(<COR>){
		chomp;
		if($line !=1){
			my($miR,$gn,$R,$pv)=split(/\t/);
			$miR=~s/\"//g;
			$gn=~s/\"//g;
			
			if($pv<=0.05 and $R<0){
				$int_pairs_corr->{$miR}->{$gn}="$R\t$pv";
			}
		}
		$line++;
	}
	return($int_pairs_corr);
	
}
sub Integration{
	my %args=@_;
	my $miRGate_pairs=$args{"pares"};
	my $corr_pairs=$args{"corr_pares"};
	my $output_file=$args{"output"};
	
	$output_file=~s/(.+)_\w+_\w+_results.*/$1/g;
	$output_file=$output_file."_Integrative_DE_miRNA_mRNA_pairs_statistical_correlation.xls";
	print STDERR "Generating $output_file\n";
	open(OUT,">$output_file") || die "Can't create $output_file\n";
	print OUT "GeneName\tmiRNA\tGene logFC\tmiRNA logFC\tPrediction Method\tTarget Site\tR coefficient\tCorrelation P-val\n";
	foreach my $miR (sort keys %{$corr_pairs}){
		if(exists $miRGate_pairs->{$miR}){
			foreach my $gn (sort {$a cmp $b} keys %{$corr_pairs->{$miR}}){
				if(exists $miRGate_pairs->{$miR}->{$gn}){
					foreach my $info ( keys %{$miRGate_pairs->{$miR}->{$gn}}){
						my $target_site;
						my @pred_meths;
						foreach my $method (keys %{$miRGate_pairs->{$miR}->{$gn}->{$info}}){
							$target_site=$miRGate_pairs->{$miR}->{$gn}->{$info}->{$method};
							push(@pred_meths,$method);
						}
						print OUT "$gn\t$miR\t$info\t".join(",",@pred_meths)."\t$target_site\t".$corr_pairs->{$miR}->{$gn} ."\n";
					}
				}
			}
		}
	}
	close OUT;
	return();
}
sub date{
	use Time::localtime;
	my $now = ctime();
	return("[$now]");
}
1;