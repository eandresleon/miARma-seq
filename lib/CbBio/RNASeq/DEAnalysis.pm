#########################################################################	
#	Differential expression analysis package							#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2014 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################

package CbBio::RNASeq::DEAnalysis;
#Export package system
require Exporter;
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(DE_noiseq DE_EdgeR DE_Analysis DE_AnalysisSummary QC_EdgeR);

use strict;
use DateTime;
use Statistics::R;
use Cwd;
use Cwd 'abs_path';

=head1 NAME

 DEAnalysis 

=head1 SYNOPSIS

DEAnalysis package is composed by three subroutines: DE_noiseq, DE_EdgeR and DE_Analysis. The aim of this package is 
perform the differential expression (DE) analysis from the tabulated files with the count of the reads 
(from htseq-count software) with Noiseq software (DE_Noiseq) or EdgeR software (DE_EdgeR). DE_Analysis is a common function
to execute one of these functions or both. 


=head1 Methods

=head2 DE_noiseq

  Example    : 
  DE_Analysis(
			projectdir=>".",
			dir=>"../4.Read_count/HtseqFormat_results/",
			file=>"../4.Read_count/HtseqFormat_results/cut_bw1-htseqresults.tab",
			targetfile=>"./targets.txt",
			label=>"cut_bw1",
			filter=>"yes",
			edger_contrastfile=>"./edger_contrast.txt",
			noiseq_contrastfile=>"./contrast.txt",
			DEsoft=>"EdgeR-Noiseq",
			filtermethod=>1,
			logfile=>"./run.log",
			verbose=>"verbose",
			cpmvalue=>2,
			repthreshold=>2,
			edger_normethod=>"rpkm",
			noiseq_normethod=>"upperquartile",
			replicates=>"yes"
		);
  
  Description: DE_Analysis performs a Differential Expression (DE) Analysis from the tabulated files with the count of the reads 
  (from htseq-count software) with Noiseq and/or EdgeR software. These analysis generate tabulated files with the results of each condition evaluated
  as well as pdf files with descriptive plots.   
  Input parameters: 
	
	Mandatory parameters:
	[projectdir] Path of the directory where the output files will be saved.
	[file] Name of the tab file which contains the number of reads from the htseq analysis
	[targetfile] Path of the target file. This is a tabulated file which contains the experimental condiction of each sample. 
	The names of the samples defined in this file will be used to the plots.
	[label] Name of the experiment to appear in the title of the plots and in the name of the results files
	[filter] This value refers to filter processing in the reads (Should be "yes" or "no")
	[logfile] Path of run.log file where execution data will be printed
	[DEsoft] Specific software to perform the Differential Expression Analysis (Allowed values: edger, noiseq or edger-noiseq)
  	 
  	 Optional parameters:
  	 [edger_contrastfile] Path of the contrast file o perform the DE analysis with EdgeR.This file has one column with the contrasts user want to evaluate. The syntax of 
  	 the contrast should be: name_of_contrast=contrast to evaluate.Any type of contrast can be done but condition name must be one (Mandatory to perform the analysis with edgeR)
  	 of the conditions present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed)
  	 [edger_normethod] Normalization method to perform the DE analysis with EdgeR. It can be one of "TMM" (default), "RLE", "upperquartile" or "none" (no normalization).
  	 [repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter with EdgeR software(2 replicates by default)
  	 [replicates] Value to indicate if replicates samples are present in the analysis to perform the DE analysis with EdgeR. It can be "yes" (by default) or "no".
  	 [bcvvalue] Value for the common BCV (square- root-dispersion) in experiments without replicates to perform the DE analysis with EdgeR. Standard values from well-controlled 
  	 experiments are 0.4 for human data (by default), 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.
  	 [cpmvalue] Cutoff for the counts per million value to be used in filter processing with methods 1 and 3 with Noiseq software and in filter processing with EdgeR (1 cpm by default).
  	 [noiseq_contrastfile] Path of the contrast file to perform the DE analysis with Noiseq.This file has one column with the contrasts user want to evaluate. 
  	 The syntax of the contrast should be: name_of_contrast=condition1-condition2. Any type of contrast can be done but condition name must be one of the conditions 
  	 present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed) (Mandatory to perform the analysis with Noiseq)
  	 [filtermethod] Method that will be used to filter proccess with Noiseq software. Method must be one of 1,2 or 3 for CPM method, Wilcoxon method and Proportion test, respectively (by default CPM method)
  	 [cutoffvalue] Cutoff for the coefficient of variation per condition to be used in filter processing with method 1 in Noiseq analysis(in percentage, 100 by default).
  	 [Rdir] Path where R software is installed
  	 [noiseq_normethod] Normalization method to perform the DE analysis with Noiseq. It can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
  	 [replicatevalue] Type of replicates to be used to perform the DE analysis with Noiseq. It can be Technical, biological or none. By default, technical replicates option is chosen.
  	 [kvalue] Counts equal to 0 are replaced by k to perform the DE analysis with Noiseq. By default, k = 0.5.
  	 [lcvalue] Length correction is done by dividing expression by length^lc to perform the DE analysis with Noiseq. By default, lc = 0.
  	 [pnrvalue] Percentage of the total reads used to simulated each sample when no replicates are available to perform the DE analysis with Noiseq. By default, pnr = 0.2.
  	 [nssvalue] Number of samples to simulate for each condition (nss>= 2) to perform the DE analysis with Noiseq. By default, nss = 5. 
  	 [vvalue] Variability in the simulated sample total reads to perform the DE analysis with Noiseq. By default, v = 0.02. Sample total reads is computed as a random value from a uniform distribution in the interval [(pnr-v)*sum(counts), (pnr+v)*sum(counts)]
  	 [qvalue] Probability of differential expression to perform the DE analysis with Noiseq. By default=0.8
  	 [verbose] Optional argument to show the execution data on screen
  
  Requeriments: DE_Analysis function requires for a correct analysis:
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

sub DE_Analysis{
	
	#Arguments provided by user are collected by %args. 
	my %args=@_;
	my $file=$args{"file"}; #Name of the tab file which contains the number of reads from the htseq analysis
	my $projectdir=$args{"projectdir"}; #Path of the directory where will be saved the results.
	my $targetfile=$args{"targetfile"}; #Path of the target file with the experimental conditions of each sample
	my $label=$args{"label"}; #Character string that will appear in the name the results file
	my $filter=$args{"filter"}; #This value refers to filter processing in the reads (Should be "yes" or "no").
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be saved
	my $DEsoft=$args{"DEsoft"}; #Specific software to perform the Differential Expression Analysis
	my $edger_normethod=$args{"edger_normethod"}; #Specific software to perform the Differential Expression Analysis
	my $rpkm=$args{"rpkm"}; #create a file with RPKM values
	my $cpm=$args{"cpm"}; #create a file with normalized CPM values
	my $miARmaPath=$args{"miARmaPath"};

  	my (undef,$dir) = fileparse($file);
	$dir=abs_path($dir);
	
	#Checking mandatory parameters
	if($DEsoft and $file and $projectdir and $targetfile and $label and $filter and $logfile){
		use File::Basename;
		#Obtaining the absolute path of the directory and the file to R execution
	  	my $projectdir = abs_path($projectdir);
	  	$file = abs_path($file);
	  	my $targetfile = abs_path($targetfile);
		
		my $addtofile=fileparse($file);
		$addtofile=~s/(.*)-ReadCount\.tab/$1/g;
		
		#Checking hidden files
		if($file =~ /^\./){
		}
		#Checking the file extension
		elsif($file =~ /.*\.tab$/){
			#Checking the software provided by the user to perform the DE analysis
			my $soft_error=0;
			my $edger_exec=0;
			my $noiseq_exec=0;
			if(lc($DEsoft) =~ /edger/){
				$edger_exec=1;
				$soft_error=1;
			}
			if(lc($DEsoft) =~ /noiseq/){
				$noiseq_exec=1;
				$soft_error=1;
			}
			if($soft_error == 0){
				warn("DE_ANALYSIS ERROR :: ".date()." Invalid value for DEsoft ($DEsoft). Allowed values are: EdgeR, Noiseq or EdgeR-Noiseq");
				help_AdapterRemoval();
			}

			#Obtaining shared variables
			my $verbose=undef;
			my $cpmvalue=undef;


			if(defined $args{"verbose"}){
				$verbose=$args{"verbose"}; #Optional argument to show the execution data on screen
			}
			if(defined $args{"cpmvalue"}){
				$cpmvalue=$args{"cpmvalue"}; #Cutoff for the counts per million value to be used in filter processing.
			}
			
			#Before DE analysis, we will perform some basic studies to gather quality of samples (in tense of data)
			
			my $repthreshold=undef;
			my $replicates=undef;
			
			if(defined $args{"repthreshold"}){
				$repthreshold=$args{"repthreshold"}; #Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)
			}

			QC_EdgeR(
				projectdir=>$projectdir,
				dir=>$dir,
				file=>fileparse($file),
				targetfile=>$targetfile,
				label=>$label ."_" . $addtofile,
				filter=>$filter,
				logfile=>$logfile,
				verbose=>$verbose,
				cpmvalue=>$cpmvalue,
				repthreshold=>$repthreshold,
				normethod=>$edger_normethod,
				miARmaPath=>$miARmaPath
			);
			
			#EdgeR analysis will be performed when user provides DEsoft with EdgeR value
			if($edger_exec == 1){
				my $edger_contrastfile=$args{"edger_contrastfile"}; #Path of the contrast file.
				$edger_contrastfile = abs_path($edger_contrastfile);
				#Optional parameters will be defined as undef variables
				my $repthreshold=undef;
				my $edger_normethod=undef;
				my $replicates=undef;
				my $bcvvalue=undef;
				
				#Obtaining optional parameters
				if(defined $args{"repthreshold"}){
					$repthreshold=$args{"repthreshold"}; #Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)
				}
				if(defined $args{"edger_normethod"}){
					$edger_normethod=$args{"edger_normethod"}; #Normalization method used in the analysis.
				}
				if(defined $args{"replicates"}){
					$replicates=$args{"replicates"}; #Value to indicate if replicates samples are present in the analysis. It can be "yes" (by default) or "no".
				}
				if(defined $args{"bcvvalue"}){
					$bcvvalue=$args{"bcvvalue"}; #Value for the common BCV (square- root-dispersion) in experiments without replicates.(0.4 by default)
				}
				my $size_file=$file;
				$size_file=~s/ReadCount/Size/g;
				if(!-e $size_file or $size_file =~ "circRNAs.tab"){
					#if size file is not found. Dont di the RPKM calculation
					$rpkm="no";
				}
				DE_EdgeR(
					projectdir=>$projectdir,
					dir=>$dir,
					size_file=>$size_file,
					file=>$file,
					targetfile=>$targetfile,
					label=>$label ."_". $addtofile,
					filter=>$filter,
					contrastfile=>$edger_contrastfile,
					logfile=>$logfile,
					verbose=>$verbose,
					cpmvalue=>$cpmvalue,
					repthreshold=>$repthreshold,
					normethod=>$edger_normethod,
					replicates=>$replicates,
					bcvvalue=>$bcvvalue,
					rpkm=>$rpkm,
					cpm=>$cpm,
					miARmaPath=>$miARmaPath
					
				);
			}

			#Noiseq analysis will be performed when user provides DEsoft with Noiseq value
			if($noiseq_exec == 1){
				my $noiseq_contrastfile=$args{"noiseq_contrastfile"}; #Path of the contrast file.
				$noiseq_contrastfile = abs_path($noiseq_contrastfile);
				#Optional parameters will be defined as undef variables
				my $lengthfile=undef;
				my $gcfile=undef;
				my $biotypefile=undef;
				my $chromsfile=undef;
				my $filtermethod=undef;
				my $cutoffvalue=undef;
				my $noiseq_normethod=undef;
				my $replicatevalue=undef;
				my $kvalue=undef;
				my $lcvalue=undef;
				my $pnrvalue=undef;
				my $nssvalue=undef;
				my $vvalue=undef;
				my $qvalue=undef;

				#Obtaining optional parameters
				if(defined $args{"lengthfile"}){
					$lengthfile=$args{"lengthfile"}; #Path of the file with the information about the length of each element.
				}
				if(defined $args{"gcfile"}){
					$gcfile=$args{"gcfile"}; #Path of the file with the information about the gc content of each element.
				}
				if(defined $args{"biotypefile"}){
					$biotypefile=$args{"biotypefile"}; #Path of the file with the information about the biotype of each element.
				}
				if(defined $args{"chromsfile"}){
					$chromsfile=$args{"chromsfile"};  #Path of the file with the genomic location information of each element.
				}
				if(defined $args{"filtermethod"}){
					$filtermethod=$args{"filtermethod"}; #Method that will be used to filter proccess.
				}
				if(defined $args{"cutoffvalue"}){
					$cutoffvalue=$args{"cutoffvalue"}; #Cutoff for the coefficient of variation per condition to be used in filter processing with method 1 
				}
				if(defined $args{"noiseq_normethod"}){
					$noiseq_normethod=$args{"noiseq_normethod"}; #Normalization method used in the analysis.
				}
				if(defined $args{"replicatevalue"}){
					$replicatevalue=$args{"replicatevalue"}; #Type of replicates to be used.
				}
				if(defined $args{"kvalue"}){
					$kvalue=$args{"kvalue"}; #Counts equal to 0 are replaced by this value.
				}
				if(defined $args{"lcvalue"}){
					$lcvalue=$args{"lcvalue"}; #Length correction is done by dividing expression by length^lc.
				}
				if(defined $args{"pnrvalue"}){
					$pnrvalue=$args{"pnrvalue"}; #Percentage of the total reads used to simulated each sample when no replicates are available.
				}
				if(defined $args{"nssvalue"}){
					$nssvalue=$args{"nssvalue"}; #Number of samples to simulate for each condition.
				}
				if(defined $args{"vvalue"}){
					$vvalue=$args{"vvalue"}; #Variability in the simulated sample total reads.
				}
				if(defined $args{"qvalue"}){
					$qvalue=$args{"qvalue"}; #Probability of differential expression.
				}

				my $size_file=$file;
				$size_file=~s/ReadCount/Size/g;
				if(!-e $size_file or $size_file =~ "circRNAs.tab"){
					#if size file is not found. Dont do the RPKM calculation
					$rpkm="no";
				}
				
				DE_noiseq(
					projectdir=>$projectdir,
					dir=>$dir,
					file=>$file,
					size_file=>$size_file,
					targetfile=>$targetfile,
					label=>$label ."_". $addtofile,
					filter=>$filter,
					contrastfile=>$noiseq_contrastfile,
					logfile=>$logfile,
					verbose=>$verbose,
					lengthfile=>$lengthfile,
					gcfile=>$gcfile,
					biotypefile=>$biotypefile,
					chromsfile=>$chromsfile,
					filtermethod=>$filtermethod,
					cpmvalue=>$cpmvalue,
					cutoffvalue=>$cutoffvalue,
					normethod=>$noiseq_normethod,
					replicatevalue=>$replicatevalue,
					kvalue=>$kvalue,
					lcvalue=>$lcvalue,
					pnrvalue=>$pnrvalue,
					nssvalue=>$nssvalue,
					vvalue=>$vvalue,
					qvalue=>$qvalue,
					rpkm=>$rpkm,
					cpm=>$cpm,
					miARmaPath=>$miARmaPath
					
				);
			}

		}else{
	  		#die ("DE_ANALYSIS ERROR :: ".date()."File($file) has an invalid format. DE_Analysis only accepts .tab files ");
			next;
	  	}
	}
	else{
		#Registering error
	   	open(LOG,">> ".$logfile) || die "DE_ANALYSIS ERROR :: ".date()."Can't open '$logfile': $!";
	    print LOG "DE_ANALYSIS ERROR :: ".date()." Directory($dir), projectdir($projectdir), file($file), DEsoft($DEsoft), targetfile($targetfile), label($label), filter($filter) and/or logfile($logfile) have not been provided";
	    close(LOG);
		#If mandatory parameters have not been provided program will die and show error message
		warn("DE_ANALYSIS ERROR :: ".date()." Directory($dir), projectdir($projectdir), file($file), DEsoft($DEsoft), targetfile($targetfile), label($label), filter($filter) and/or logfile($logfile) have not been provided");
		help_DE_Analysis();
	}
	sub help_DE_Analysis{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
			[projectdir] Path of the directory where the output files will be saved.
			[dir] Path of the directory which contains the file with the reads. This directory will be configured as working directory for R
			[file] Name of the tab file which contains the number of reads from the htseq analysis
			[targetfile] Path of the target file. This is a tabulated file which contains the experimental condiction of each sample. 
			The names of the samples defined in this file will be used to the plots.
			[label] Name of the experiment to appear in the title of the plots and in the name of the results files
			[filter] This value refers to filter processing in the reads (Should be "yes" or "no")
			[logfile] Path of run.log file where execution data will be printed
			[DEsoft] Specific software to perform the Differential Expression Analysis (Allowed values: edger, noiseq or edger-noiseq)
  	 
  	 		Optional parameters:
  	 		[edger_contrastfile] Path of the contrast file o perform the DE analysis with EdgeR.This file has one column with the contrasts user want to evaluate. The syntax of 
  	 		the contrast should be: name_of_contrast=contrast to evaluate.Any type of contrast can be done but condition name must be one (Mandatory to perform the analysis with edgeR)
  			 of the conditions present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed)
  	 		[edger_normethod] Normalization method to perform the DE analysis with EdgeR. It can be one of "TMM" (default), "RLE", "upperquartile" or "none" (no normalization).
  	 		[repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter with EdgeR software(2 replicates by default)
  	 		[replicates] Value to indicate if replicates samples are present in the analysis to perform the DE analysis with EdgeR. It can be "yes" (by default) or "no".
  	 		[bcvvalue] Value for the common BCV (square- root-dispersion) in experiments without replicates to perform the DE analysis with EdgeR. Standard values from well-controlled 
  	 		experiments are 0.4 for human data (by default), 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.
  	 		[cpmvalue] Cutoff for the counts per million value to be used in filter processing with methods 1 and 3 with Noiseq software and in filter processing with EdgeR (1 cpm by default).
  	 		[noiseq_contrastfile] Path of the contrast file to perform the DE analysis with Noiseq.This file has one column with the contrasts user want to evaluate. 
  	 		The syntax of the contrast should be: name_of_contrast=condition1-condition2. Any type of contrast can be done but condition name must be one of the conditions 
  	 		present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed) (Mandatory to perform the analysis with Noiseq)
  	 		[filtermethod] Method that will be used to filter proccess with Noiseq software. Method must be one of 1,2 or 3 for CPM method, Wilcoxon method and Proportion test, respectively (by default CPM method)
  	 		[cutoffvalue] Cutoff for the coefficient of variation per condition to be used in filter processing with method 1 in Noiseq analysis(in percentage, 100 by default).
  			[Rdir] Path where R software is installed
  			[noiseq_normethod] Normalization method to perform the DE analysis with Noiseq. It can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
  	 		[replicatevalue] Type of replicates to be used to perform the DE analysis with Noiseq. It can be Technical, biological or none. By default, technical replicates option is chosen.
  	 		[kvalue] Counts equal to 0 are replaced by k to perform the DE analysis with Noiseq. By default, k = 0.5.
  	 		[lcvalue] Length correction is done by dividing expression by length^lc to perform the DE analysis with Noiseq. By default, lc = 0.
  	 		[pnrvalue] Percentage of the total reads used to simulated each sample when no replicates are available to perform the DE analysis with Noiseq. By default, pnr = 0.2.
  	 		[nssvalue] Number of samples to simulate for each condition (nss>= 2) to perform the DE analysis with Noiseq. By default, nss = 5. 
  	 		[vvalue] Variability in the simulated sample total reads to perform the DE analysis with Noiseq. By default, v = 0.02. Sample total reads is computed as a random value from a uniform distribution in the interval [(pnr-v)*sum(counts), (pnr+v)*sum(counts)]
  	 		[qvalue] Probability of differential expression to perform the DE analysis with Noiseq. By default=0.8
  			[lengthfile] Path of the file with the information about the length of each element to perform the DE analysis with Noiseq. This file must contain the name of the element and the legnth of this element under the name.
  	 		[gcfile] Path of the file with the information about the gc content of each element to perform the DE analysis with Noiseq. This file must contain the name of the element and the gc content (%) of this element under the name.
  	 		[biotypefile] Path of the file with the information about the biotype of each element to perform the DE analysis with Noiseq. This file must contain the name of the element and the biotype of this element under the name.
  	 		[chromsfile] Path of the file with the genomic location information of each element to perform the DE analysis with Noiseq. This file must contain 4 columns corresponding to: the name of the element, the number of chromosome, and the coordinates of the start and end of the element.
  	 		[verbose] Optional argument to show the execution data on screen
  	 		               
			Examples:
			DE_Analysis(projectdir=>".", dir=>"../4.Read_count/HtseqFormat_results/", file=>"../4.Read_count/HtseqFormat_results/cut_bw1-htseqresults.tab", targetfile=>"./targets.txt", label=>"cut_bw1", filter=>"yes", 
				edger_contrastfile=>"./edger_contrast.txt", noiseq_contrastfile=>"./contrast.txt", DEsoft=>"EdgeR-Noiseq", filtermethod=>1, logfile=>"./run.log", verbose=>"verbose", cpmvalue=>2, repthreshold=>2, edger_normethod=>"rpkm",
				noiseq_normethod=>"upperquartile", replicates=>"yes");
		};

	print STDERR $usage;
	exit(); 
	}  	

}

=head1 Methods

=head2 QC_EdgeR

  Example    : 
  QC_EdgeR(
			projectdir=>".",
			dir=>"../4.Read_count/HtseqFormat_results/",
			file=>"../4.Read_count/HtseqFormat_results/cut_bw1-htseqresults.tab",
			targetfile=>"./targets.txt",
			label=>"cut_bw1",
			contrastfile=>"./contrast.txt",
			filter=>"yes",
			cpmvalue=>1,
			repthreshold=>2,
			replicates=>"yes",
		);
  
  Description: QC_EdgeR is a R function which takes a tab file with the number of the reads from htseq-count analysis and a target file 
with the experimental conditions of the samples to analyze the quality of the samples. For this purpose QC_EdgeR genrates a pdf 
report with differents plots: boxplots and density plots with the distribution of the reads in each sample with 
and without normalization, MDS and PCA plots and heatmaps of the samples with the more expressed elements(genes, miRNAs...)
and between the samples.  
  Input parameters: 
	Mandatory arguments:
  	  [projectdir] Path of the directory where will be saved the QC report.
	  [dir] the path of the directory which contains the files. This directory will be configured as working directory for R
	  [file] name of the tab file which contains the number of reads from the htseq analysis
	  [targetfile] Path of the tabulated file which contains the experimental condiction of each sample. 
	  First column must contain the names to be used to the plots, and the next columns the condition of each factor. 

  QC_EdgeR works with 1 or 2 factors.

  	[label] Name of the experiment to appear in the title of the plots and in the name of the pdf file results
  	[filter] This value refers to filter processing in the reads (Should be \"yes\" or \"no\").

Optional arguments:
  [cpmvalue] Cutoff for the counts per million value to be used in filter processing (1 cpm by default).
  [repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)", stderr())
  
  Requeriments: QC_EdgeR function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- R v3.1.0 or higher software correctly installed
  	- EdgeR package correctly installed
  	- Output file from htseq-count analysis with the number of the reads of each sample.
  	
  Returntype : pdf file with descriptive plots of the analysis and a tab file (xls extension) for condition evaluated in EdgeR_results directory.  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub QC_EdgeR{
	#First, check that R is in path:
	my @R_bin=`which R`;
	#Executing the command
	if(scalar(@R_bin)<1){
		die "QC_EdgeR ERROR :: system args failed: $? : Is R installed and exported to \$PATH ?";
	}

	#Arguments provided by user are collected by %args. 
	my %args=@_;
	my $file=$args{"file"}; #Name of the tab file which contains the number of reads from the htseq analysis
	my $dir=$args{"dir"}; #Path of the directory which contains the files. This directory will be configured as working directory for R
	my $projectdir=$args{"projectdir"}; #Path of the directory where will be saved the QC report.
	my $targetfile=$args{"targetfile"}; #Path of the target file with the experimental conditions of each sample
	my $label=$args{"label"}; #Character string that will appear in the name the results file
	my $filter=$args{"filter"}; #This value refers to filter processing in the reads (Should be "yes" or "no").
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be saved
	my $miARmaPath=$args{"miARmaPath"};
	
	
	if($file and $dir and $projectdir and $targetfile and $label and $filter and $logfile){

		#Optional parameters
		my $verbose=$args{"verbose"}; #Optional argument to show the execution data on screen

		#Defining results directory 
		my $output_dir=$projectdir;

		#Defining R command
		my $Rcommand="projectdir=\"".$output_dir."\",dir=\"".$dir."\", file=\"".$file."\", targetfile=\"".$targetfile."\", label=\"".$label."\", filter=\"".$filter."\"";

		#If user has provided any optional parameter will be added to Rcommand
		if(defined $args{"cpmvalue"}){
			my $cpmvalue=$args{"cpmvalue"}; #Cutoff for the counts per million value to be used in filter processing.
			$Rcommand.=", cpmvalue=".$cpmvalue;
		}
		if(defined $args{"repthreshold"}){
			my $repthreshold=$args{"repthreshold"}; #Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)
			$Rcommand.=", repthreshold=".$repthreshold;
		}
		if(defined $args{"normethod"}){
			my $normethod=$args{"normethod"}; #Normalization method used in the analysis.
			$Rcommand.=", normethod=\"".$normethod."\"";
		}
		if(defined $args{"replicates"}){
			my $replicates=$args{"replicates"}; #Value to indicate if replicates samples are present in the analysis. It can be "yes" (by default) or "no".
			$Rcommand.=", replicates=\"".$replicates."\"";
		}
		if(defined $args{"bcvvalue"}){
			my $bcvvalue=$args{"bcvvalue"}; #Value for the common BCV (square- root-dispersion) in experiments without replicates.(0.4 by default)
			$Rcommand.=", bcvvalue=".$bcvvalue;
		}


		# Printing the date and command execution on screen
		#print STDERR "QC_EdgeR :: ".date()." Starting Differential Expression Analysis of $file with EdgeR\n"; 

		#Creating results directory
		my $command="mkdir -p ".$output_dir;
		system($command) == 0
	   	    or die "QC_EdgeR ERROR :: ".date()." System args failed: $? ($command)";

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
		#Declaring R instructions for the differential expression analysis. DE_noiseq R function is needed 
		my $cmds = <<EOF;
		source("$miARmaPath/lib/CbBio/RNASeq/R-Scripts/QC_EdgeR.R")
		setwd("$dir")
		resultsfiles<-QC_EdgeR($Rcommand)
EOF
		
		#Printing the execution data on log file and on the screen if verbose parameter is defined 
		open (LOG,">> ".$logfile) || die "QC_EdgeR ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "QC_EdgeR :: ".date()." Executing $cmds\n";

		if($verbose){
				print STDOUT "QC_EdgeR :: ".date()." Executing $cmds\n";
		}

		#R commands execution
		my $out2 = $R->run($cmds);
		#Obtaining the files generated in the analysis
		my @media = $R->get('resultsfiles');
		foreach my $a (@media){
			print LOG "QC_EdgeR :: ".date()." The file ".$a." has been generated \n";
		}
		close LOG;
	}
	else
	{
		#Registering the error
		open(LOG,">> ".$logfile) || die "QC_EdgeR ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "QC_EdgeR ERROR :: ".date()." File($file), dir($dir), projectdir($projectdir), targetfile($targetfile), label($label), filter($filter) and/or logfile($logfile) have not been provided";
		close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn("QC_EdgeR ERROR :: ".date()." File($file), dir($dir), projectdir($projectdir), targetfile($targetfile), label($label), filter($filter) and/or logfile($logfile) have not been provided");
		help_QC_EdgeR();
	}
	sub help_QC_EdgeR{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
			[projectdir] Path of the directory where the output files will be saved.
			[dir] Path of the directory which contains the file with the reads. This directory will be configured as working directory for R
			[file] Name of the tab file which contains the number of reads from the htseq analysis
			[targetfile] Path of the target file. This is a tabulated file which contains the experimental condiction of each sample. 
			The names of the samples defined in this file will be used to the plots.
			[label] Name of the experiment to appear in the title of the plots and in the name of the results files
			[filter] This value refers to filter processing in the reads (Should be "yes" or "no")
			[contrastfile] Path of the contrast file.This file has one column with the contrasts user want to evaluate. The syntax of 
			the contrast will be: name_of_contrast=contrast to evaluate.Any type of contrast can be done but condition name must be one 
			of the conditions present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed)
			[logfile] Path of run.log file where execution data will be printed
  	 
  	 		Optional parameters:
  	 		[repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)
  			[replicates] Value to indicate if replicates samples are present in the analysis. It can be "yes" (by default) or "no".
  	 		experiments are 0.4 for human data (by default), 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.
  			[verbose] Optional argument to show the execution data on screen
  	 		               
			Examples:
			QC_EdgeR(projectdir=>".", dir=>"../4.Read_count/HtseqFormat_results/", file=>"../4.Read_count/HtseqFormat_results/cut_bw1-htseqresults.tab",targetfile=>"./targets.txt", label=>"cut_bw1", contrastfile=>"./contrast.txt", 
				filter=>"yes", cpmvalue=>1, repthreshold=>2,replicates=>"yes");
		};

	print STDERR $usage;
	exit(); 
	}  	
}

=head2 DE_noiseq

  Example    : 
  DE_noiseq(
			projectdir=>".",
			dir=>"../4.Read_count/HtseqFormat_results/",
			file=>"../4.Read_count/HtseqFormat_results/cut_bw1-htseqresults.tab",
			targetfile=>"./targets.txt",
			label=>"cut_bw1",
			filter=>"yes",
			contrastfile=>"./contrast.txt",
			lengthfile=>"./length.txt",
			gcfile=>"./gc.txt",
			biotypefile=>"./biotype.txt",
			chromsfile=>"./chroms.txt",
			filtermethod=>1,
			cpmvalue=>1,
			cutoffvalue=>100,
			normethod=>"rpkm",
			replicatevalue=>"technical"
			logfile=>"run.log",
			Rdir=>"../Applications/R".
			verbose=>"verbose"
		);
  
  Description: DE_noiseq function calls a R function called DE_noiseq which takes a tab file with the 
  number of the reads from htseq-count analysis, a target file with the experimental conditions 
  of the samples and the contrast file with the contrast to analyze the differential expression between the defined samples. 
  This function outputs a pdf file with the descriptive plots of the analysis and two tab files (xls extension to easy exploration) 
  for condition evaluated. The Noiseq_results tab file contains the whole results data and the Noiseq_DE tab file contains only the data of the 
  differentially expressed elements.   
  Input parameters: 
	Mandatory parameters:
	[projectdir] Path of the directory where the output files will be saved.
	[dir] Path of the directory which contains the file with the reads. This directory will be configured as working directory for R
	[file] Path of the tab file which contains the number of reads from the htseq analysis
	[targetfile] Path of the target file. This is a tabulated file which contains the experimental condiction of each sample. 
	The names of the samples defined in this file will be used to the plots.
	[label] Name of the experiment to appear in the title of the plots and in the name of the results files
	[filter] This value refers to filter processing in the reads (Should be "yes" or "no")
	[contrastfile] Path of the contrast file.This file has one column with the contrasts user want to evaluate. The syntax of the contrast 
	will be: name_of_contrast=condition1-condition2.Any type of contrast can be done but condition name must be one of the conditions 
	present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed)
	[logfile] Path of run.log file where execution data will be printed
  	 Optional parameters:
  	 [filtermethod] Method that will be used to filter proccess. Method must be one of 1,2 or 3 for CPM method, Wilcoxon method and Proportion test, respectively (by default CPM method)
  	 [cpmvalue] Cutoff for the counts per million value to be used in filter processing with methods 1 and 3 (1 cpm by default).
  	 [cutoffvalue] Cutoff for the coefficient of variation per condition to be used in filter processing with method 1 (in percentage, 100 by default).
  	 [Rdir] Path where R software is installed
  	 [normethod] Normalization method. It can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
  	 [replicatevalue] Type of replicates to be used. It can be Technical, biological or none. By default, technical replicates option is chosen.
  	 [kvalue] Counts equal to 0 are replaced by k. By default, k = 0.5.
  	 [lcvalue] Length correction is done by dividing expression by length^lc. By default, lc = 0.
  	 [pnrvalue] Percentage of the total reads used to simulated each sample when no replicates are available. By default, pnr = 0.2.
  	 [nssvalue] Number of samples to simulate for each condition (nss>= 2). By default, nss = 5. 
  	 [vvalue] Variability in the simulated sample total reads. By default, v = 0.02. Sample total reads is computed as a random value from a uniform distribution in the interval [(pnr-v)*sum(counts), (pnr+v)*sum(counts)]
  	 [qvalue] Probability of differential expression. By default=0.8
   	 [verbose] Optional argument to show the execution data on screen
  Requeriments: DE_noiseq function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- R v3.1.0 or higher software correctly installed
  	- Noiseq package correctly installed
  	- Output file from htseq-count analysis with the number of the reads of each sample
  	
  Returntype : pdf file with descriptive plots of the analysis and two tab files (xls extension) for condition evaluated in Noiseq_results directory.  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut



sub DE_noiseq{
	#First, check that R is in path:
	my @R_bin=`which R`;
	#Executing the command
	if(scalar(@R_bin)<1){
		die "DE_NOISEQ ERROR :: system args failed: $? : Is R installed and exported to \$PATH ?";
	}

	#Arguments provided by user are collected by %args. 
	my %args=@_;
	my $file=$args{"file"}; #Name of the tab file which contains the number of reads from the htseq analysis
	my $dir=$args{"dir"}; #Path of the directory which contains the files. This directory will be configured as working directory for R
	my $projectdir=$args{"projectdir"}; #Path of the directory where will be saved the QC report.
	my $targetfile=$args{"targetfile"}; #Path of the target file with the experimental conditions of each sample
	my $label=$args{"label"}; #Character string that will appear in the name the results file
	my $filter=$args{"filter"}; #This value refers to filter processing in the reads (Should be "yes" or "no").
	my $contrastfile=$args{"contrastfile"}; #Path of the contrast file.
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be saved
	my $rpkm=$args{"rpkm"}; #create a file with RPKM values
	my $cpm=$args{"cpm"}; #create a file with normalized CPM values
	my $file_size=$args{"size_file"}|| undef; #file with gene/tr/miRNA lengths
	my $miARmaPath=$args{"miARmaPath"};
	
	if($file and $dir and $projectdir and $targetfile and $label and $filter and $contrastfile and $logfile){
		#Optional parameters
		my $verbose=$args{"verbose"}; #Optional argument to show the execution data on screen

		#Defining results directory 
		my $output_dir=$projectdir."/Noiseq_results";

		#Defining R command
		my $Rcommand="projectdir=\"".$output_dir."\",dir=\"".$dir."\", file=\"".$file."\", targetsfile=\"".$targetfile."\", label=\"".$label."\", filter=\"".$filter."\", contrastfile=\"".$contrastfile."\"";

		#If user has provided any optional parameter will be added to Rcommand
		if(defined $args{"filtermethod"}){
			my $filtermethod=$args{"filtermethod"}; #Method that will be used to filter proccess.
			$Rcommand.=", filtermethod=\"".$filtermethod."\"";
		}
		if(defined $args{"cpmvalue"}){
			my $cpmvalue=$args{"cpmvalue"}; #Cutoff for the counts per million value to be used in filter processing.
			$Rcommand.=", cpmvalue=".$cpmvalue;
		}
		if(defined $args{"cutoffvalue"}){
			my $cutoffvalue=$args{"cutoffvalue"}; #Cutoff for the coefficient of variation per condition to be used in filter processing with method 1 
			$Rcommand.=", cutoffvalue=".$cutoffvalue;
		}
		if(defined $args{"normethod"}){
			my $normethod=$args{"normethod"}; #Normalization method used in the analysis.
			$Rcommand.=", normethod=\"".$normethod."\"";
		}
		if(defined $args{"replicatevalue"}){
			my $replicatevalue=$args{"replicatevalue"}; #Type of replicates to be used.
			$Rcommand.=", replicatevalue=\"".$replicatevalue."\"";
		}
		if(defined $args{"kvalue"}){
			my $kvalue=$args{"kvalue"}; #Counts equal to 0 are replaced by this value.
			$Rcommand.=", kvalue=".$kvalue;
		}
		if(defined $args{"lcvalue"}){
			my $lcvalue=$args{"lcvalue"}; #Length correction is done by dividing expression by length^lc.
			$Rcommand.=", lcvalue=".$lcvalue;
		}
		if(defined $args{"pnrvalue"}){
			my $pnrvalue=$args{"pnrvalue"}; #Percentage of the total reads used to simulated each sample when no replicates are available.
			$Rcommand.=", pnrvalue=".$pnrvalue;
		}
		if(defined $args{"nssvalue"}){
			my $nssvalue=$args{"nssvalue"}; #Number of samples to simulate for each condition.
			$Rcommand.=", nssvalue=".$nssvalue;
		}
		if(defined $args{"vvalue"}){
			my $vvalue=$args{"vvalue"}; #Variability in the simulated sample total reads.
			$Rcommand.=", vvalue=".$vvalue;
		}
		if(defined $args{"qvalue"}){
			my $qvalue=$args{"qvalue"}; #Probability of differential expression.
			$Rcommand.=", qvalue=".$qvalue;
		}
		if(defined $args{"rpkm"}){
			my $rpkm=$args{"rpkm"}; #RPKM
			$Rcommand.=", rpkm=\"".$rpkm."\"";
			$Rcommand.=", file_size=\"".$file_size."\"";
		}
		if(defined $args{"cpm"}){
			my $cpm=$args{"cpm"}; #normalized CPM
			$Rcommand.=", cpm=\"".$cpm."\"";
		}
		
		# Printing the date and command execution on screen
		print STDOUT "DE_NOISEQ :: ".date()." Starting Differential Expression Analysis of $file with NOISeq\n" if($verbose);

		#Creating results directory
		my $command="mkdir -p ".$output_dir;
		system($command) == 0
	   	    or die "DE_NOISEQ ERROR :: system args failed: $? ($command)";

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
		#Declaring R instructions for the differential expression analysis. DE_noiseq R function is needed 
		my $cmds = <<EOF;
		source("$miARmaPath/lib/CbBio/RNASeq/R-Scripts/DE_noiseq.R")
		setwd("$dir")
		resultsfiles<-DE_noiseq($Rcommand)
EOF
		#For testing : 	source("/Users/eandres/Proyectos/miARma/lib/CbBio/RNASeq/R-Scripts/DE_noiseq.R")
		#For testing: source("http://valkyrie.us.es/CbBio/RNASeq/R-Scripts/DE_noiseq.R")
		#Printing the execution data on log file and on the screen if verbose parameter is defined 
		open (LOG,">> ".$logfile) || die $!;
		print LOG "DE_NOISEQ :: ".date()." Executing $cmds\n";

		if($verbose){
				print STDOUT "DE_NOISEQ :: ".date()." Executing $cmds\n";
		}

		#R commands execution
		my $out2 = $R->run($cmds);
		#Obtaining the files generated in the analysis
		my @media = $R->get('resultsfiles');
		foreach my $a (@media){
			foreach my $b (@{$a}){
				print LOG "DE_NOISEQ :: ".date()." The file ".$b." has been generated \n";
			}
		}
		close LOG;
		
	}
	else 
	{
		#Registering the error
		open(LOG,">> ".$logfile) || die "DE_NOISEQ ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "DE_NOISEQ ERROR :: ".date()." File($file), dir($dir), projectdir($projectdir), targetfile($targetfile), label($label) filter($filter), contrastfile($contrastfile) and/or logfile($logfile) have not been provided";
		close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn("DE_NOISEQ ERROR :: ".date()." File($file), dir($dir), projectdir($projectdir), targetfile($targetfile), label($label) filter($filter), contrastfile($contrastfile) and/or logfile($logfile) have not been provided");
		help_DE_noiseq();
	}
	sub help_DE_noiseq{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
			[projectdir] Path of the directory where the output files will be saved.
			[dir] Path of the directory which contains the file with the reads. This directory will be configured as working directory for R
			[file] Name of the tab file which contains the number of reads from the htseq analysis
			[targetfile] Path of the target file. This is a tabulated file which contains the experimental condiction of each sample. 
			The names of the samples defined in this file will be used to the plots.
			[label] Name of the experiment to appear in the title of the plots and in the name of the results files
			[filter] This value refers to filter processing in the reads (Should be "yes" or "no")
			[contrastfile] Path of the contrast file.This file has one column with the contrasts user want to evaluate. The syntax of the contrast 
			will be: name_of_contrast=condition1-condition2.Any type of contrast can be done but condition name must be one of the conditions 
			present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed)
			[logfile] Path of run.log file where execution data will be printed
  	 
  	 		Optional parameters:
  	 		[filtermethod] Method that will be used to filter proccess. Method must be one of 1,2 or 3 for CPM method, Wilcoxon method and Proportion test, respectively (by default CPM method)
  	 		[cpmvalue] Cutoff for the counts per million value to be used in filter processing with methods 1 and 3 (1 cpm by default).
  	 		[cutoffvalue] Cutoff for the coefficient of variation per condition to be used in filter processing with method 1 (in percentage, 100 by default).
  	 		[Rdir] Path where R software is installed
  	 		[normethod] Normalization method. It can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
  	 		[replicatevalue] Type of replicates to be used. It can be Technical, biological or none. By default, technical replicates option is chosen.
  	 		[kvalue] Counts equal to 0 are replaced by k. By default, k = 0.5.
  	 		[lcvalue] Length correction is done by dividing expression by length^lc. By default, lc = 0.
  	 		[pnrvalue] Percentage of the total reads used to simulated each sample when no replicates are available. By default, pnr = 0.2.
  	 		[nssvalue] Number of samples to simulate for each condition (nss>= 2). By default, nss = 5. 
  	 		[vvalue] Variability in the simulated sample total reads. By default, v = 0.02. Sample total reads is computed as a random value from a uniform distribution in the interval [(pnr-v)*sum(counts), (pnr+v)*sum(counts)]
  	 		[qvalue] Probability of differential expression. By default=0.8
  	 		[verbose] Optional argument to show the execution data on screen
  	 		               
			Examples:
			DE_noiseq(projectdir=>".", dir=>"../4.Read_count/HtseqFormat_results/", file=>"../4.Read_count/HtseqFormat_results/cut_bw1-htseqresults.tab", targetfile=>"./targets.txt", label=>"cut_bw1", filter=>"yes", 
				contrastfile=>"./contrast.txt", lengthfile=>"./length.txt", gcfile=>"./gc.txt", biotypefile=>"./biotype.txt", chromsfile=>"./chroms.txt", filtermethod=>1, cpmvalue=>1, cutoffvalue=>100, normethod=>"rpkm", 
				replicatevalue=>"technical", logfile=>"run.log", Rdir=>"../Applications/R", verbose=>"verbose")
		};

	print STDERR $usage;
	exit(); 
	}  	
}

=head1 Methods

=head2 DE_EdgeR

  Example    : 
  DE_EdgeR(
			projectdir=>".",
			dir=>"../4.Read_count/HtseqFormat_results/",
			file=>"../4.Read_count/HtseqFormat_results/cut_bw1-htseqresults.tab",
			targetfile=>"./targets.txt",
			label=>"cut_bw1",
			contrastfile=>"./contrast.txt",
			filter=>"yes",
			cpmvalue=>1,
			repthreshold=>2,
			normethod=>"TMM",
			replicates=>"yes",
			bcvvalue=>0.4,
			logfile=>"run.log",
			Rdir=>"../Applications/R",
			verbose=>"verbose"
		);
  
  Description: DE_EdgeR function calls a R function called DE_EdgeR which takes a tab file with the 
  number of the reads from htseq-count analysis, a target file with the experimental conditions 
  of the samples and the contrast file with the contrast to analyze the differential expression between the defined samples. 
  This function outputs a pdf file with the descriptive plots of the analysis and a tab file (xls extension to easy exploration) 
  for condition evaluated. This function accepts one or two factors experiments.   
  Input parameters: 
	Mandatory parameters:
	[projectdir] Path of the directory where the output files will be saved.
	[dir] Path of the directory which contains the file with the reads. This directory will be configured as working directory for R
	[file] Name of the tab file which contains the number of reads from the htseq analysis
	[targetfile] Path of the target file. This is a tabulated file which contains the experimental condiction of each sample. 
	The names of the samples defined in this file will be used to the plots.
	[label] Name of the experiment to appear in the title of the plots and in the name of the results files
	[filter] This value refers to filter processing in the reads (Should be "yes" or "no")
	[contrastfile] Path of the contrast file.This file has one column with the contrasts user want to evaluate. The syntax of 
	the contrast will be: name_of_contrast=contrast to evaluate.Any type of contrast can be done but condition name must be one 
	of the conditions present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed)
	[logfile] Path of run.log file where execution data will be printed
  	 Optional parameters:
  	 [repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)
  	 [normethod] Normalization method. It can be one of "TMM" (default), "RLE", "upperquartile" or "none" (no normalization).
  	 [replicates] Value to indicate if replicates samples are present in the analysis. It can be "yes" (by default) or "no".
  	 [bcvvalue] Value for the common BCV (square- root-dispersion) in experiments without replicates. Standard values from well-controlled 
  	 experiments are 0.4 for human data (by default), 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.
  	 [verbose] Optional argument to show the execution data on screen
  Requeriments: DE_EdgeR function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- R v3.1.0 or higher software correctly installed
  	- EdgeR package correctly installed
  	- Output file from htseq-count analysis with the number of the reads of each sample.
  	
  Returntype : pdf file with descriptive plots of the analysis and a tab file (xls extension) for condition evaluated in EdgeR_results directory.  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub DE_EdgeR{
	#First, check that R is in path:
	my @R_bin=`which R`;
	#Executing the command
	if(scalar(@R_bin)<1){
		die "DE_EDGER ERROR :: system args failed: $? : Is R installed and exported to \$PATH ?";
	}

	#Arguments provided by user are collected by %args. 
	my %args=@_;
	my $file=$args{"file"}; #Name of the tab file which contains the number of reads from the htseq analysis
	my $dir=$args{"dir"}; #Path of the directory which contains the files. This directory will be configured as working directory for R
	my $projectdir=$args{"projectdir"}; #Path of the directory where will be saved the QC report.
	my $targetfile=$args{"targetfile"}; #Path of the target file with the experimental conditions of each sample
	my $label=$args{"label"}; #Character string that will appear in the name the results file
	my $filter=$args{"filter"}; #This value refers to filter processing in the reads (Should be "yes" or "no").
	my $contrastfile=$args{"contrastfile"}; #Path of the contrast file.
	my $logfile=$args{"logfile"}; #Path of run.log file where execution data will be saved
	my $file_size=$args{"size_file"}; #file with gene/tr/miRNA lengths
	my $miARmaPath=$args{"miARmaPath"};
	
	if($file and $dir and $projectdir and $targetfile and $label and $filter and $contrastfile and $logfile){

		#Optional parameters
		my $verbose=$args{"verbose"}; #Optional argument to show the execution data on screen

		#Defining results directory 
		my $output_dir=$projectdir."/EdgeR_results";

		#Defining R command
		my $Rcommand="projectdir=\"".$output_dir."\",dir=\"".$dir."\", file=\"".$file."\", targetfile=\"".$targetfile."\", label=\"".$label."\", contrastfile=\"".$contrastfile."\", filter=\"".$filter."\"";

		#If user has provided any optional parameter will be added to Rcommand
		if(defined $args{"cpmvalue"}){
			my $cpmvalue=$args{"cpmvalue"}; #Cutoff for the counts per million value to be used in filter processing.
			$Rcommand.=", cpmvalue=".$cpmvalue;
		}
		if(defined $args{"repthreshold"}){
			my $repthreshold=$args{"repthreshold"}; #Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)
			$Rcommand.=", repthreshold=".$repthreshold;
		}
		if(defined $args{"normethod"}){
			my $normethod=$args{"normethod"}; #Normalization method used in the analysis.
			$Rcommand.=", normethod=\"".$normethod."\"";
		}
		if(defined $args{"replicates"}){
			my $replicates=$args{"replicates"}; #Value to indicate if replicates samples are present in the analysis. It can be "yes" (by default) or "no".
			$Rcommand.=", replicates=\"".$replicates."\"";
		}
		if(defined $args{"bcvvalue"}){
			my $bcvvalue=$args{"bcvvalue"}; #Value for the common BCV (square- root-dispersion) in experiments without replicates.(0.4 by default)
			$Rcommand.=", bcvvalue=".$bcvvalue;
		}
		if(defined $args{"rpkm"}){
			my $rpkm=$args{"rpkm"}; #RPKM
			$Rcommand.=", rpkm=\"".$rpkm."\"";
			$Rcommand.=", file_size=\"".$file_size."\"";
		}
		if(defined $args{"cpm"}){
			my $cpm=$args{"cpm"}; #normalized CPM
			$Rcommand.=", cpm=\"".$cpm."\"";
		}
		# Printing the date and command execution on screen
		print STDOUT "DE_EDGER :: ".date()." Starting Differential Expression Analysis of $file with EdgeR\n" if($verbose);

		#Creating results directory
		my $command="mkdir -p ".$output_dir;
		system($command) == 0
	   	    or die "DE_EDGER ERROR :: ".date()." System args failed: $? ($command)";

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
		#Declaring R instructions for the differential expression analysis. DE_noiseq R function is needed 
		my $cmds = <<EOF;
		source("$miARmaPath/lib/CbBio/RNASeq/R-Scripts/DE_EdgeR.R")
		setwd("$dir")
		resultsfiles<-DE_EdgeR($Rcommand)
EOF
		#For testing : source("http://valkyrie.us.es/CbBio/RNASeq/R-Scripts/DE_EdgeR.R")
		#Printing the execution data on log file and on the screen if verbose parameter is defined 
		open (LOG,">> ".$logfile) || die "DE_EDGER ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "DE_EDGER :: ".date()." Executing $cmds\n";

		if($verbose){
				print STDOUT "DE_EDGER :: ".date()." Executing $cmds\n";
		}

		#R commands execution
		my $out2 = $R->run($cmds);
		#Obtaining the files generated in the analysis
		my @media = $R->get('resultsfiles');
		foreach my $a (@media){
			foreach my $b (@{$a}){
				print LOG "DE_EDGER :: ".date()." The file ".$b." has been generated \n";
			}
		}
		close LOG;
		
	}
	else
	{
		#Registering the error
		open(LOG,">> ".$logfile) || die "DE_EDGER ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "DE_EDGER ERROR :: ".date()." File($file), dir($dir), projectdir($projectdir), targetfile($targetfile), label($label), filter($filter), contrastfile($contrastfile) and/or logfile($logfile) have not been provided";
		close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn("DE_EDGER ERROR :: ".date()." File($file), dir($dir), projectdir($projectdir), targetfile($targetfile), label($label), filter($filter), contrastfile($contrastfile) and/or logfile($logfile) have not been provided");
		help_DE_EdgeR();
	}
	sub help_DE_EdgeR{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
			[projectdir] Path of the directory where the output files will be saved.
			[dir] Path of the directory which contains the file with the reads. This directory will be configured as working directory for R
			[file] Name of the tab file which contains the number of reads from the htseq analysis
			[targetfile] Path of the target file. This is a tabulated file which contains the experimental condiction of each sample. 
			The names of the samples defined in this file will be used to the plots.
			[label] Name of the experiment to appear in the title of the plots and in the name of the results files
			[filter] This value refers to filter processing in the reads (Should be "yes" or "no")
			[contrastfile] Path of the contrast file.This file has one column with the contrasts user want to evaluate. The syntax of 
			the contrast will be: name_of_contrast=contrast to evaluate.Any type of contrast can be done but condition name must be one 
			of the conditions present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed)
			[logfile] Path of run.log file where execution data will be printed
  	 
  	 		Optional parameters:
  	 		[repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)
  	 		[normethod] Normalization method. It can be one of "TMM" (default), "RLE", "upperquartile" or "none" (no normalization).
  			[replicates] Value to indicate if replicates samples are present in the analysis. It can be "yes" (by default) or "no".
  	 		[bcvvalue] Value for the common BCV (square- root-dispersion) in experiments without replicates. Standard values from well-controlled 
  	 		experiments are 0.4 for human data (by default), 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.
  			[verbose] Optional argument to show the execution data on screen
  	 		               
			Examples:
			DE_EdgeR(projectdir=>".", dir=>"../4.Read_count/HtseqFormat_results/", file=>"../4.Read_count/HtseqFormat_results/cut_bw1-htseqresults.tab",targetfile=>"./targets.txt", label=>"cut_bw1", contrastfile=>"./contrast.txt", 
				filter=>"yes", cpmvalue=>1, repthreshold=>2, normethod=>"TMM",replicates=>"yes",bcvvalue=>0.4,logfile=>"run.log",Rdir=>"../Applications/R", verbose=>"verbose");
		};

	print STDERR $usage;
	exit(); 
	}  	
}

sub DE_AnalysisSummary{
	use File::Basename;
	#Arguments provided by user are collected by %args. Dir, file, aligner, statsfile, projectdir 
	#and logfile are mandatory arguments while verbose and threads are optional.
	my %args=@_;
	my $summary_file=$args{"summary"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}; #Optional arguments to show the execution data on screen
	my $DEsoft=$args{"DEsoft"}; #Input directory where results directory will be created
	my $contrastfile=$args{"contrastfile"};
	my $comparisons;
	open(TARGET,$contrastfile) || warn date() . " WARN :: contrastfile ($contrastfile) is missing\n";
	while(<TARGET>){
		chomp;
		$_=~s/\"//g;
		if($_ !~ /Name/){
			my($comp,$form)=split(/=/);
			$comparisons->{$comp}=$form;
		}
	}
	close TARGET;
	my $summary_edger;
	my $summary_noi;
	
	if(lc($DEsoft) eq "edger"){
		opendir(EGDER, $projectdir ."/EdgeR_results/") || warn "DEAnalysis:: Folder \"$projectdir/EdgeR_results/\" is not found\n"; 
		my @edgeR_files= readdir(EGDER);
		foreach my $edge_files (sort @edgeR_files){
			if($edge_files =~ /\.xls$/){				
				foreach my $comp (sort keys %{$comparisons}){
					if($edge_files =~ /$comp/){
						open(EDGEFILE,$projectdir ."/EdgeR_results/$edge_files") || die "$! $projectdir/EdgeR_results/$edge_files";
						my $cont_pval=0;
						my $cont_fdr=0;
						while(<EDGEFILE>){
							chomp;
							$_=~s/\"//g;
							if($_ !~ /logFC/){
								my @data=split(/\t/);
								if($data[($#data-1)]<=0.05){
									$cont_pval++;
								}
								if($data[($#data)]<=0.05){
									$cont_fdr++;
								}
							}
						}
						close EDGEFILE;
						$summary_edger->{$comp}->{$edge_files}->{$cont_pval}=$cont_fdr;
					}
				}
			}
		}
	}
	if(lc($DEsoft) eq "noiseq"){
		my $summary;
		opendir(NOISEQ, $projectdir ."/Noiseq_results/") || warn "DEAnalysis:: Folder \"$projectdir/Noiseq_results/\" is not found\n"; 
		my @Noiseq_files= readdir(NOISEQ);
		foreach my $Noiseq_files (sort @Noiseq_files){
			if($Noiseq_files =~ /\.xls$/){
				foreach my $comp (sort keys %{$comparisons}){
					if($Noiseq_files =~ /$comp/){
						open(NOIFILE,$projectdir ."/Noiseq_results/$Noiseq_files");
						my $cont_pval=0;
						while(<NOIFILE>){
							chomp;
							$_=~s/\"//g;
							if($_ !~ /prob/){
								my @data=split(/\t/);
								if($data[5]>0.8){
									$cont_pval++;
								}
							}
						}
						close NOIFILE;
						$summary_noi->{$comp}->{$Noiseq_files}=$cont_pval;
					}
				}
			}
		}
	}
	if(lc($DEsoft) eq "noiseq-edger" or lc($DEsoft) eq "edger-noiseq" ){
		DE_AnalysisSummary(
		  	projectdir=>$projectdir,
			DEsoft=>"edger",
			summary=>$summary_file,
			contrastfile=>$contrastfile,
		);
		DE_AnalysisSummary(
		  	projectdir=>$projectdir,
			DEsoft=>"noiseq",
			summary=>$summary_file,
			contrastfile=>$contrastfile,
		);
	}

	if(scalar(keys %{$summary_edger})>0){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
		print SUMM "\nDifferential Expression Analysis by edgeR [".$projectdir ."/EdgeR_results/]\n";
		print SUMM "Comparison\tFile\tNumber of DE elements (Pval <=0.05)\tNumber of DE elements (FDR <=0.05)\n";
		foreach my $comp (sort keys %{$summary_edger}){
			foreach my $file (sort keys %{$summary_edger->{$comp}}){
				foreach my $pval (sort keys %{$summary_edger->{$comp}->{$file}}){
					print SUMM $comparisons->{$comp} ."\t$file\t$pval\t". $summary_edger->{$comp}->{$file}->{$pval}."\n";
				}
			}
		}
		close SUMM;
	}
	if(scalar(keys %$summary_noi)>0){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
		print SUMM "\nDifferential Expression Analysis by NoiSeq [".$projectdir ."/Noiseq_results/]\n";
		print SUMM "Comparison\tFile\tNumber of DE elements (Prob >=0.8)\n";
		foreach my $comp (sort keys %{$summary_noi}){
			foreach my $file (sort keys %{$summary_noi->{$comp}}){
				print SUMM $comparisons->{$comp}  ."\t$file\t". $summary_noi->{$comp}->{$file}."\n";
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
