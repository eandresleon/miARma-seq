#########################################################################	
#	Readcount execution package	 										#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2015 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################

package CbBio::RNASeq::Readcount;
require Exporter;
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(featureCount featureFormat CIRICount1 CIRICount2 CIRIFormat miRDeepCount miRDeepFormat featureSummary miRDeepSummary);

use strict;
use File::Basename;
use DateTime;

=head1 NAME

 Readcount

=head1 SYNOPSIS

Readcount is composed of 2 subroutines: featureCount and featureFormat. featureCount takes a bam/sam 
file and counts the number of reads of each type saving this information in a new file at 
htseq_results directory. featureFormat takes these files and creates a combined file
with the number of reads of each type in one column per file. 

=head1 Methods

=head2 htseq-count

  Example    : 
  featureCount(
  	file=>"file.bam",
  	database=>"db.gtf",
  	seqid=>"transcript_id",
  	parameters=>" -a 5", 
  	strand=>"no", 
  	featuretype=>"miRNA",
  	logfile=>"run.log",
  	inputformat=>"bam",
  	verbose=>"verbose", 
  	projectdir=>"."
  );
  Description: featureCount takes the provided file (bam/sam format) usually from bowtie analysis and 
  counts the reads of each type comparing with a provided database. The output file with 
  the number of the reads is saved in a htseq_results directory on the project directory. 
  Execution data will be saved on directory run.log file and will show on screen if verbose 
  option is selected.
  Input parameters: 
	Mandatory parameters:
  	 [file] Path of the file which is going to be processed (bam format)
  	 [logfile] Path of run.log file where execution data will be saved
  	 [database] GFF file used to calculate the number of reads in htseq-count analysis
  	 [projectdir] Directory where htseq_results directory will be created
  	 Optional parameters:
  	 [seqid] GFF attribute to be used as feature ID (default: gene_id) for htseq-count analysis
  	 [parameters] Other htseq-count parameters to perform the analysis using the htseq-count recommended syntaxis		  
  	 [strand] Whether the data is from a strand-specific assay (yes, no or reverse, yes by default) for htseq-count analysis
  	 [featuretype] Feature type (3rd column in GFF file) to be used, all features of other type are ignored (default:exon) for htseq-count analysis
  	 [inputformat] Format of the input file : sam or bam (default:sam)
  	 [verbose] Option to show the execution data on the screen   
  Returntype : File at directory htseq_results. Also returns the path of the output file.
  Requeriments: featureCount function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- htseq-count 0.5.4. or higher software correctly installed
  	- Input files on bam format 
  Caller     : web drawing code
  Status     : Stable

=cut

sub featureCount{

	#Arguments provided by user are collected by %args. File, dir, adapter, projectdir and logfile
	# are mandatory arguments while verbose is optional. 
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $arch=`uname`;
	chomp($arch);
	
	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/".$arch."/featurecounts/:". "$miARmaPath/bin/".$arch."/samtools/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/".$arch."/featurecounts/". ":$miARmaPath/bin/".$arch."/samtools/";
	}
	
	#First, check that htseq-count is in path:
	my @htseq_bin=`which featureCounts`;
	#Executing the command
	if(scalar(@htseq_bin)<1){
		die "SEQCOUNT ERROR :: system args failed: $? : Is featureCounts installed and exported to \$PATH ?";
	}
	#First, check that samtools is in path:
	my @sam_bin=`which samtools`;
	#Executing the command
	if(scalar(@sam_bin)<1){
		die "SEQCOUNT ERROR :: system args failed: $? : Is samtools installed and exported to \$PATH ?";
	}
	
	my $file=$args{"file"}; #Path of the BAM file which is going to be processed
	my $database=$args{"database"}; #gtf file used as reference to count the reads
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	
	#Variable declaration
	my $command;
	my $commanddef;
	
	#As default some htseq-count parameters has been defined, 
	my $htseqpardef= "";
	#Strand-specific analysis can be provided by the user
	if(defined($args{"strand"})){
		my $strand=$args{"strand"};
		my $feature_strand=1;
		if(lc($strand) eq "yes"){
			$feature_strand=1;
		}
		if(lc($strand) eq "no"){
			$feature_strand=0;
		}
		if(lc($strand) eq "reverse"){
			$feature_strand=2;
		}
		$htseqpardef.=" -s $feature_strand";
	}

	#feature type (3rd column in GFF file) to be used (exon by default)
	if(defined($args{"featuretype"})){
		my $featuretype=$args{"featuretype"};
		$htseqpardef.=" -t $featuretype";
	}
	#GFF attribute to be used as feature ID (gene_id as default)
	if(defined($args{"seqid"})){
		my $seqid=$args{"seqid"};
		$htseqpardef.=" -g $seqid";
	}
	#Any other htseq parameter can be provided by the user using the correct sintaxis
	if(defined($args{"parameters"})){
		my $parameters=$args{"parameters"};
		$htseqpardef.=" $parameters";
	}
	if(defined($args{"threads"})){
		my $threads=$args{"threads"};
		$htseqpardef.=" -T $threads";
	}
	if(defined($args{"quality"})){
		my $quality=$args{"quality"};
		$htseqpardef.=" -Q $quality";
	}
	if(defined($args{"Seqtype"})){
		my $Seqtype=$args{"Seqtype"};
		if(lc($Seqtype) =~ /paired/){
			$htseqpardef.=" -p ";
		}
	}	
	#Checking the mandatory parameters
	if ($file and $projectdir and $database and $logfile){
		if($file !~ /no_aligned/ and $file !~ /unmapped/){
		if($file =~ /.*\.bam$/ or $file =~ /.*\.sam$/){
			print STDOUT "SEQCOUNT :: " . date() . " Reading counts from $file\n" if($verbose);
			my $name;
			if ($file =~ /.*\.sam$/){
				#As input file contains the whole path fileparse is used to obtain the name of the file
				$name= fileparse($file, qr{\.sam.*});
			}
			elsif ($file =~ /.*\.bam$/){
				#Altough htseq can count redas in a BAM file, it needs several software, so we are going to create a sam file from a bam file
				$name= fileparse($file, qr{\.bam.*});
			}
			# else{
			# 	print STDERR "SEQCOUNT ERROR :: $inputformat has an invalid value. Allowed values are: sam or bam\n";
			# 	exit;
			# }
			#Checking the processation of the sample to classify it correctly
			$name =~ /.*_([a-z]{2,3}_bw[1-2]|his|str)/;
			my $prefix=$1;
			my $output_dir="/".$prefix."_readcount_results/";
			#htseq-count execution command
			$command="featureCounts ".$htseqpardef." -a ".$database." -o ".$projectdir.$output_dir.$name.".tab " . $file ;
			#commandef is the command will be executed by system composed of the results directory creation 
			#and the htseq_count execution. The error data will be printed on the run.log file
			$commanddef="mkdir -p ".$projectdir.$output_dir." ;".$command." 2>> ".$logfile;;
	
			#Opening the run.log and printing the execution data
			open (LOG,">> ".$logfile) || die $!;
			print LOG "SEQCOUNT :: " . date() . " Reading counts from $file\n";
			print LOG "SEQCOUNT :: ".date()." Executing $commanddef\n";
			close LOG;
			#If verbose option has been provided program will print execution data on the screen too.
			if($verbose){
				print STDOUT "SEQCOUNT :: ".date()." Executing $commanddef\n";
			}
			#Executing the command or if system can't be executed die showing the error.
			system($commanddef) == 0
			or die "SEQCOUNT ERROR :: system args failed: $? ($commanddef)";
		
			#After counting redas, we are going to conver the sam file into a sorted BAM file (a smaller file)
			#if(lc($inputformat) eq "sam"){
			#	FromSam2Bam($file);
			#} 
			#The path of the output file is returned to the main program
			if($file){
				return($projectdir.$output_dir.$name.".tab");
			}
			else{
				next;
			}
			
		}}
	}
	else{
		#Registering the error
		open(LOG,">> ".$logfile) || die "SEQCOUNT ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "SEQCOUNT ERROR :: ".date()." Projectdir ($projectdir), file ($file), logfile ($logfile) and/or database($database) have not been provided";
		close LOG;
			#If mandatory parameters have not been provided program dies and shows error message
		warn("SEQCOUNT ERROR :: Projectdir ($projectdir), file ($file), logfile ($logfile) and/or database($database) have not been provided");
		help_featureCount();
	}
	sub help_featureCount{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
  	 		[file] Path of the file which is going to be processed (bam format)
  	 		[logfile] Path of run.log file where execution data will be saved
  	 		[database] GFF file used to calculate the number of reads in htseq-count analysis
  	 		[projectdir] Directory where htseq_results directory will be created
  	 
  	 		Optional parameters:
  	 		[seqid] GFF attribute to be used as feature ID (default: gene_id) for htseq-count analysis
  	 		[parameters] Other htseq-count parameters to perform the analysis using the htseq-count recommended syntaxis		  
  	 		[strand] Whether the data is from a strand-specific assay (yes, no or reverse, yes by default) for htseq-count analysis
  	 		[featuretype] Feature type (3rd column in GFF file) to be used, all features of other type are ignored (default:exon) for htseq-count analysis
  	 		[inputformat] Format of the input file : sam or bam (default:sam)
  	 		[verbose] Option to show the execution data on the screen   
  	 		               
			Examples:
			featureCount(file=>"file.bam", database=>"db.gtf", seqid=>"transcript_id", parameters=>" -a 5", strand=>"no", featuretype=>"miRNA", logfile=>"run.log", inputformat=>"bam", verbose=>"verbose", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}  	
}

=head2 featureFormat

  Example    : 
  	featureFormat( 
  		input=>"\@results", 
  		projectdir=>".",
  		logfile=>"run.log"
  	);
  Description: featureFormat collects the paths of the files processed by featureCount by a 
  reference array and create a combined tabulated file named htseqresults.tab in the 
  htseq_results directory with the number of reads of each file in each column. The names
  of the files will be the name of the original files without extension.
  Input parameters: 
	Mandatory parameters:
  	 [projectdir] Input directory where htseqresults.tab will be saved
  	 [input] Reference array with the paths of the files which have to be joined
   	 [logfile] Path of run.log file where execution data will be saved
  Returntype : Tab file at directory htseq_results and the path of the combined file
  Requeriments: featureCount function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Output files from htseq-count anlysis on the htseq_results directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub featureFormat{
	
	#Arguments provided by user are collected by %args. Dir and input are mandatory arguments.
	my %args=@_;
	my $input=$args{"input"}; #Input is an referenced array with the paths of the files which have to be joined
	my $projectdir=$args{"projectdir"}; #Directory to create the output file
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"};
	
	#Variable declaration
	my $hashref;
	my $data;
	my $files;
	my $process;
	my @results;	
	my $size;

	if($input and $projectdir and $logfile){
		#Declaration of the output path 
		my $output_dir="/Readcount_results/";
		#Output directory creation
		my $command= "mkdir -p ".$projectdir.$output_dir;
		system($command) == 0
		or die "SEQCOUNTFORMAT ERROR :: system args failed: $? ($command)";

		#Reading each path of the input referenced array
		foreach my $file(sort {$a cmp $b} @$input){
			if($file){
				#Fileparse is used to obtain the name of the file
				my $name=fileparse($file, qr{\.tab$});
				#Obtaining the process info from the file name and the original name of the file
				$name =~ /(.*)_([a-z]{2,3}_bw[1-2]|his|str)/;
				my $originalname=$1;
				my $suffix=$2;
				#Saving the Suffix with process information in a hash
				$process->{$suffix}++;
				$files->{$originalname}++;
				#Opening each file of the array
				open (FILE, $file) or die "SEQCOUNTFORMAT ERROR :: ".date()."Can't open '$file': $!";
				#Reading the file. Each line contains the number of the miRNA, a tab and the number of reads
				while(<FILE>){
					#Deleting the last character of the line
					chomp;
					if($_ !~ /^#/ and $_ !~ /^Geneid/){
						#Using the split with the tab we can keep the name of the miRNA in $miRNA variable
						#and the number of reads in $reads variable. 
						my ($miRNAs,$chr,$start,$end,$strand,$length,$reads)=split(/\t/);
						#Creating a hash with the names of miRNAs
						$hashref->{$miRNAs}++;
						$size->{$suffix}->{$miRNAs}=$length;
						#Creating a hash with all the data, the first dimension contains the miRNA name,
						#the second the name of the file and the third is the number of reads 
						$data->{$suffix}->{$miRNAs}->{$originalname}=$reads;	
					}	
				}
				#Closing the file
				close FILE;
			}
		}
		foreach my $suffix(sort keys %$data){
			#Printing the results of the process with each combination of software
			my $fileresults= $projectdir.$output_dir.$suffix."-ReadCount.tab"; 
			#Opening the output file to write the data
			open(RESULTS,"> ".$fileresults) || die "SEQCOUNTFORMAT ERROR :: ".date()."Can't open '$fileresults': $!";
			#Printing on the first row of the results file the name of the input files sorted and 
			#separated by a tab 
			print RESULTS join("\t",(sort keys %$files))."\n";
			#The next lines will be printed going over the hash. Reading the first dimension of the hash
			open(LOG,">> ".$logfile) || die "SEQCOUNTFORMAT ERROR :: ".date()."Can't open '$logfile': $!";
			foreach my $miRNAs(sort keys %{$data->{$suffix}}){
				#At the end of the input files appeared some useless information to our fileresults  
				if (($miRNAs =~ "ambiguous") or ($miRNAs =~ "no_feature") or ($miRNAs =~ "alignment_not_unique")  or ($miRNAs =~ "not_aligned") or ($miRNAs =~ "too_low_aQual")){
					print LOG "SEQCOUNTFORMAT SUMMARY :: $miRNAs\n";
					foreach my $name(sort keys %{$data->{$suffix}->{$miRNAs}}){
			 			print LOG "SEQCOUNTFORMAT SUMMARY :: ". $data->{$suffix}->{$miRNAs}->{$name} ."\n";
					}
				}
				else{
					#Printing the name of each miRNA
					print RESULTS $miRNAs."\t";
					my @results;
					#An array is used to save the number of reads of each file for a defined miRNA.  
					foreach my $name(sort keys %{$data->{$suffix}->{$miRNAs}}){
			 			push(@results,$data->{$suffix}->{$miRNAs}->{$name});
					}
					#Printing the number of the reads of each miRNA separated by a tab 
					print RESULTS join("\t",@results)."\n";
				}
			}
			#Closing the files
			close RESULTS;
			close LOG;
			push(@results, $fileresults) if($fileresults);
		}
		
		#Saving length file for RPKM calculation
		

		foreach my $suffix (keys %{$size}){
			my $fileSize= $projectdir.$output_dir.$suffix."-Size.tab"; 
			open(SIZE,">$fileSize") || warn "Cant create $fileSize for RPKM calculation\n";
			print SIZE "Gene\tLength\n";

			foreach my $element (keys %{$size->{$suffix}}){
				print SIZE $element."\t". $size->{$suffix}->{$element} ."\n";
			}
			close SIZE;
		}
		
		#Returning the path of the results files
		open(LOG,">> ".$logfile) || die "SEQCOUNTFORMAT ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "SEQCOUNT :: ".date()." Please check the folder: $projectdir"."$output_dir.\n";
		close LOG;
		
		print STDOUT "SEQCOUNT :: ".date()." Please check the folder: $projectdir"."$output_dir.\n" if($verbose);
		
		return(@results);
	}else{
		#Registering the error
		open(LOG,">> ".$logfile) || die "SEQCOUNTFORMAT ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "SEQCOUNTFORMAT ERROR :: ".date()." Projectdir ($projectdir), input paths($input) and/or logfile($logfile) have not been provided";
		close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn("SEQCOUNTFORMAT ERROR :: ".date()." Projectdir ($projectdir), input paths($input) and/or logfile($logfile) have not been provided");
		help_featureFormat();
	}
	sub help_featureFormat{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
  	 		[projectdir] Input directory where htseqresults.tab will be saved
  	 		[input] Reference array with the paths of the files which have to be joined
   	 		[logfile] Path of run.log file where execution data will be saved
  	 
			Examples:
			featureFormat(input=>"\@results", projectdir=>".", logfile=>"run.log");
	};

	print STDERR $usage;
	exit(); 
	}  	
}

=head2 FromSam2Bam

  Example    : 
  	FromSam2Bam("File.sam");
  Description: FromSam2Bam is a function to convert sam files to bam format. This function collects the path 
  of the sam file convert it to bam saving it with the same name in the same directory, deleting the sam file 
  and returning the path of the new bam file.  
  Input parameters: 
	Mandatory parameters:
  	 [file] Path of the sam file to convert it to bam format
  Returntype : Bam file in the same directory and the path of the new bam file
  Requeriments: FromSam2Bam function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Samtools 0.1.19 or higher correctly installed
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub FromSam2Bam{
	#Obtaining the sam file
	my $sam_file=shift;
	#Declaring the new bam file with the same path and the same name but with different extension
	my $bam_file=$sam_file;
	$bam_file=~ s/.sam$//;

	#Execution command to convert sam file to bam format
	my $commanddef="samtools view -bS $sam_file | samtools sort - $bam_file";
	#Executing the command
	system($commanddef) == 0
		or die "SEQCOUNT ERROR :: system args failed: $? ($commanddef)";

	#Deleting the sam file
	unlink($sam_file);
	#Returning the path of the new bam file
	return($bam_file);
}

=head2 FromBam2Sam

  Example    : 
  	FromBam2Sam("File.sam");
  Description: FromBam2Sam is a function to convert bam files to sam format. This function collects the path 
  of the bam file, converts it to sam saving it with the same name in the same directory and returns the path 
  of the new sam file.  
  Input parameters: 
	Mandatory parameters:
  	 [file] Path of the bam file to convert it to sam format
  Returntype : Sam file in the same directory and the path of the new sam file
  Requeriments: FromBam2Sam function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Samtools 0.1.19 or higher correctly installed
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub FromBam2Sam{
	#Obtaining the bam file
	my $bam_file=shift;
	#Declaring the new sam file with the same path and the same name but with different extension
	my $sam_file=$bam_file;
	$sam_file=~ s/.bam$/.sam/;

	#Execution command to convert bam file to sam format
	my $commanddef="samtools view -h $bam_file > $sam_file";
	#Executing the command
	system($commanddef) == 0
		or die "SEQCOUNT ERROR :: system args failed: $? ($commanddef)";

	#Returning the path of the new sam file
	return($sam_file);
}


=head2 CIRICount

  Example    : 
  	featureFormat( 
  		input=>"\@results", 
  		projectdir=>".",
  		logfile=>"run.log"
  	);
  Description: featureFormat collects the paths of the files processed by featureCount by a 
  reference array and create a combined tabulated file named htseqresults.tab in the 
  htseq_results directory with the number of reads of each file in each column. The names
  of the files will be the name of the original files without extension.
  Input parameters: 
	Mandatory parameters:
  	 [projectdir] Input directory where htseqresults.tab will be saved
  	 [input] Reference array with the paths of the files which have to be joined
   	 [logfile] Path of run.log file where execution data will be saved
  Returntype : Tab file at directory htseq_results and the path of the combined file
  Requeriments: featureCount function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Output files from htseq-count anlysis on the htseq_results directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub CIRICount1{
    #Arguments provided by user are collected by %args. File, dir, adapter, projectdir and logfile
    # are mandatory arguments while verbose is optional. 
    my %args=@_;
    my $miARmaPath=$args{"miARmaPath"};
    my $arch=`uname`;
    chomp($arch);

    if ($ENV{PATH}) {
            $ENV{PATH} .= ":$miARmaPath/bin/common/CIRI/";
    }
    else {
            $ENV{PATH} = "$miARmaPath/bin/common/CIRI/";
    }

    #First, check that CIRI is in path:
    my @ciri_bin=`which CIRI_v1.2.pl`;
    #Executing the command
    if(scalar(@ciri_bin)<1){
            die "CIRICount ERROR :: system args failed: $? : Is CIRI installed and exported to \$PATH ?";
    }

    my $file=$args{"file"}; #Path of the BAM file which is going to be processed
    my $database=$args{"database"}; #gtf file used as reference to count the reads
    my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
    my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
    my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
    my $Seqtype=$args{"Seqtype"}; #Input directory where results directory will be created
    my $fasta=$args{"fasta"}; #Input directory where results directory will be created

    #Variable declaration
    my $command;
    my $commanddef;

    #Checking the mandatory parameters
    if ($file and $projectdir and $database and $logfile and $fasta and $Seqtype){
            if($file =~ /.*\.sam$/){
                    my $output_dir="/"."circRNAs_results/";
                    my($filename) = fileparse($file);
                    $filename=~s/\.sam//g;
                    if(lc($Seqtype) eq "pairedend" or lc($Seqtype) eq "paired" or lc($Seqtype) eq "paired-end"){
                            print STDOUT "CIRI :: ".date()." Checking $file for a circRNA analysis (Paired End)\n" if($verbose);
                            #CIRI execution command
                            #perl CIRI.pl -P -I test.sam -O outfile -F chr1.fa -A chr1.gtf
                            $command="CIRI_v1.2.pl -P -I $file -O ".$projectdir.$output_dir.$filename.".ciri -A " . $database ." -F " . $fasta . " -G " . $logfile;
                            #commandef is the command will be executed by system composed of the results directory creation 
                            #and the htseq_count execution. The error data will be printed on the run.log file
                            $commanddef="mkdir -p ".$projectdir.$output_dir." ;".$command." > ".$logfile ." 2>&1";
                    }
                    else{
                            #CIRI execution command
                            print STDOUT "CIRI :: ".date()." Checking $file for a circRNA analysis (Single End)\n" if($verbose);
                            $command="CIRI_v1.2.pl -S -I $file -O ".$projectdir.$output_dir.$filename.".ciri -A " . $database ." -F " . $fasta  . " -G " . $logfile;
                            #commandef is the command will be executed by system composed of the results directory creation 
                            #and the htseq_count execution. The error data will be printed on the run.log file
                            $commanddef="mkdir -p ".$projectdir.$output_dir." ;".$command." >".$logfile ." 2>&1";
                    }
                    #Opening the run.log and printing the execution data
                    open (LOG,">> ".$logfile) || die $!;
                    print LOG "CIRICount :: ".date()." Executing $commanddef\n";
                    close LOG;
                    #If verbose option has been provided program will print execution data on the screen too.
                    if($verbose){
                            print STDOUT "CIRICount :: ".date()." Executing $commanddef\n";
                    }
                    #Executing the command or if system can't be executed die showing the error.
                    system($commanddef) == 0
                    or die "CIRICount ERROR :: system args failed: $? ($commanddef)";

                    #After counting redas, we are going to conver the sam file into a sorted BAM file (a smaller file)
                    #The path of the output file is returned to the main program
                    if($file){
                            return($projectdir.$output_dir.$filename.".ciri");
                    }
                    else{
                            next;
                    }
            }
    }
	else{
		#Registering the error
		open(LOG,">> ".$logfile) || die "CIRICount ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "CIRICount ERROR :: ".date()." Projectdir ($projectdir), file ($file), logfile ($logfile) and/or database($database) have not been provided";
		close LOG;
			#If mandatory parameters have not been provided program dies and shows error message
		warn("CIRICount ERROR :: Projectdir ($projectdir), file ($file), logfile ($logfile), bwa index (fasta), sequencing type ($Seqtype) and/or database($database) have not been provided");
		help_CIRICount();
	}
	sub help_CIRICount{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
	 		[file] Path of the file which is going to be processed (bam format)
	 		[logfile] Path of run.log file where execution data will be saved
	 		[fasta] fasta file used to calculate the number of reads in CIRI analysis
	 		[projectdir] Directory where htseq_results directory will be created

	 		Optional parameters:
	 		[seqid] GFF attribute to be used as feature ID (default: gene_id) for CIRI analysis
	 		[parameters] Other CIRI parameters to perform the analysis using the CIRI recommended syntaxis		  
	 		[strand] Whether the data is from a strand-specific assay (yes, no or reverse, yes by default) for CIRI analysis
	 		[featuretype] Feature type (3rd column in GFF file) to be used, all features of other type are ignored (default:exon) for CIRI analysis
	 		[inputformat] Format of the input file : sam or bam (default:sam)
	 		[verbose] Option to show the execution data on the screen   

			Examples:
			CIRICount(file=>"file.bam", database=>"db.gtf", seqid=>"transcript_id", parameters=>" -a 5", strand=>"no", featuretype=>"miRNA", logfile=>"run.log", inputformat=>"bam", verbose=>"verbose", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}
}

sub CIRICount2{
	#Arguments provided by user are collected by %args. File, dir, adapter, projectdir and logfile
	# are mandatory arguments while verbose is optional. 
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	
	my $arch=`uname`;
	chomp($arch);

	if ($ENV{PATH}) {
		$ENV{PATH} .= ":$miARmaPath/bin/common/CIRI/";
	}
	else {
		$ENV{PATH} = "$miARmaPath/bin/common/CIRI/";
	}
	#First, check that CIRI is in path:
	my @ciri_bin=`which CIRI_v2.0.1.pl`;
	#Executing the command
	if(scalar(@ciri_bin)<1){
		die "CIRICount ERROR :: system args failed: $? : Is CIRI installed and exported to \$PATH ?";
	}
	
	my $file=$args{"file"}; #Path of the BAM file which is going to be processed
	my $database=$args{"database"}; #gtf file used as reference to count the reads
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created
	my $Seqtype=$args{"Seqtype"}; #Input directory where results directory will be created
	my $fasta=$args{"fasta"}; #Input directory where results directory will be created
	my $htseqpardef;
	#Variable declaration
	my $command;
	my $commanddef;
	
	if(defined($args{"threads"})){
		my $threads=$args{"threads"};
		$htseqpardef.=" -0 -T $threads";
	}
	
	#Checking the mandatory parameters
	if ($file and $projectdir and $database and $logfile and $fasta and $Seqtype){
		if($file =~ /.*\.sam$/){
			my $output_dir="/"."circRNAs_results/";
			my($filename) = fileparse($file);
			$filename=~s/\.sam//g;
			if(lc($Seqtype) eq "pairedend" or lc($Seqtype) eq "paired" or lc($Seqtype) eq "paired-end"){
				print STDOUT "CIRI :: ".date()." Checking $file for a circRNA analysis (Paired End)\n" if($verbose);
				#CIRI execution command
				#perl CIRI.pl -P -I test.sam -O outfile -F chr1.fa -A chr1.gtf
				$command="CIRI_v2.0.1.pl $htseqpardef -I $file -O ".$projectdir.$output_dir.$filename.".ciri -A " . $database ." -F " . $fasta . " -G " . $logfile;
				#commandef is the command will be executed by system composed of the results directory creation 
				#and the htseq_count execution. The error data will be printed on the run.log file
				$commanddef="mkdir -p ".$projectdir.$output_dir." ;".$command." > ".$logfile ." 2>&1";
			}
			else{
				#CIRI execution command
				print STDOUT "CIRI :: ".date()." Checking $file for a circRNA analysis (Single End)\n" if($verbose);
				$command="CIRI_v2.0.1.pl $htseqpardef -I $file -O ".$projectdir.$output_dir.$filename.".ciri -A " . $database ." -F " . $fasta  . " -G " . $logfile;
				#commandef is the command will be executed by system composed of the results directory creation 
				#and the htseq_count execution. The error data will be printed on the run.log file
				$commanddef="mkdir -p ".$projectdir.$output_dir." ;".$command." >".$logfile ." 2>&1";
			}
			#Opening the run.log and printing the execution data
			open (LOG,">> ".$logfile) || die $!;
			print LOG "CIRICount :: ".date()." Executing $commanddef\n";
			close LOG;
			#If verbose option has been provided program will print execution data on the screen too.
			if($verbose){
				print STDOUT "CIRICount :: ".date()." Executing $commanddef\n";
			}
			#Executing the command or if system can't be executed die showing the error.
			system($commanddef) == 0
			or die "CIRICount ERROR :: system args failed: $? ($commanddef)";
			
			#After counting redas, we are going to conver the sam file into a sorted BAM file (a smaller file)
			#The path of the output file is returned to the main program
			if($file){
				return($projectdir.$output_dir.$filename.".ciri");
			}
			else{
				next;
			}
		}
	}
	else{
		#Registering the error
		open(LOG,">> ".$logfile) || die "CIRICount ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "CIRICount ERROR :: ".date()." Projectdir ($projectdir), file ($file), logfile ($logfile) and/or database($database) have not been provided";
		close LOG;
			#If mandatory parameters have not been provided program dies and shows error message
		warn("CIRICount ERROR :: Projectdir ($projectdir), file ($file), logfile ($logfile), bwa index (fasta), sequencing type ($Seqtype) and/or database($database) have not been provided");
		help_CIRICount();
	}
	sub help_CIRICount{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
	 		[file] Path of the file which is going to be processed (bam format)
	 		[logfile] Path of run.log file where execution data will be saved
	 		[fasta] fasta file used to calculate the number of reads in CIRI analysis
	 		[projectdir] Directory where htseq_results directory will be created

	 		Optional parameters:
	 		[seqid] GFF attribute to be used as feature ID (default: gene_id) for CIRI analysis
	 		[parameters] Other CIRI parameters to perform the analysis using the CIRI recommended syntaxis		  
	 		[strand] Whether the data is from a strand-specific assay (yes, no or reverse, yes by default) for CIRI analysis
	 		[featuretype] Feature type (3rd column in GFF file) to be used, all features of other type are ignored (default:exon) for CIRI analysis
	 		[inputformat] Format of the input file : sam or bam (default:sam)
	 		[verbose] Option to show the execution data on the screen   

			Examples:
			CIRICount(file=>"file.bam", database=>"db.gtf", seqid=>"transcript_id", parameters=>" -a 5", strand=>"no", featuretype=>"miRNA", logfile=>"run.log", inputformat=>"bam", verbose=>"verbose", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}
}

=head2 CIRIFormat

  Example    : 
  	CIRIFormat( 
  		input=>"\@results", 
  		projectdir=>".",
  		logfile=>"run.log"
  	);
  Description: featureFormat collects the paths of the files processed by CIRI by a 
  reference array and create a combined tabulated file named htseqresults.tab in the 
  htseq_results directory with the number of reads of each file in each column. The names
  of the files will be the name of the original files without extension.
  Input parameters: 
	Mandatory parameters:
  	 [projectdir] Input directory where htseqresults.tab will be saved
  	 [input] Reference array with the paths of the files which have to be joined
   	 [logfile] Path of run.log file where execution data will be saved
  Returntype : Tab file at directory htseq_results and the path of the combined file
  Requeriments: featureCount function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Output files from htseq-count anlysis on the htseq_results directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub CIRIFormat{
	
	#Arguments provided by user are collected by %args. Dir and input are mandatory arguments.
	my %args=@_;
	my $input=$args{"input"}; #Input is an referenced array with the paths of the files which have to be joined
	my $projectdir=$args{"projectdir"}; #Directory to create the output file
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"};
	my $summary_file=$args{"summary"};
	
	#Variable declaration
	my $hashref;
	my $data;
	my $files;
	my @results;	

	if($input and $projectdir and $logfile){
		open(LOG,">> ".$logfile) || die "CIRI ERROR :: ".date()."Can't open '$logfile': $!";
		
		#Declaration of the output path 
		my $output_dir="/circRNAs_results/";
		#Reading each path of the input referenced array
		my $results;
		foreach my $file(sort {$a cmp $b} @$input){
			if($file){
				print LOG "LOG :: Reading $file\n";
				if($file =~ /\.ciri$/){
					#Fileparse is used to obtain the name of the file
					my $name=fileparse($file, qr{\.ciri$});
					#Obtaining the process info from the file name and the original name of the file
					$name =~ /(.*)_bwa/;
					my $originalname=$1;
					my $suffix=$2;
					if($name !~/^\./){
						print STDOUT "\tCIRI :: ". date() . " Reading $file\n" if($verbose);
						#Saving the Suffix with process information in a hash
						$files->{$originalname}++;
			
						#Opening each file of the array
						open (FILE, $file) or die "CIRI ERROR :: ".date()."Can't open '$file': $!";
						#Reading the file. Each line contains the number of the miRNA, a tab and the number of reads
						while(<FILE>){
							#Deleting the last character of the line
							chomp;
							if($_ !~ /^circRNA_ID/){
								#Using the split with the tab we can keep the name of the miRNA in $miRNA variable
								#and the number of reads in $reads variable. 
								my ($circRNAs,$chr,$start,$end,$reads)=split(/\t/);
								#Creating a hash with the names of miRNAs
								$hashref->{$circRNAs}->{$originalname}=$reads;
								#Creating a hash with all the data, the first dimension contains the miRNA name,
								#the second the name of the file and the third is the number of reads 
								$data->{$originalname}++;
								$results->{$originalname}->{$circRNAs}++;	
							}	
						}
						#Closing the file
						close FILE;
					}
				}
				else{
					print LOG "WARM :: the file [$file] is not in the correct format\n";
				}
			}
		}
		#Printing the results of the process with each combination of software
		my $fileresults= $projectdir.$output_dir."circRNAs.tab"; 
		#Opening the output file to write the data
		open(RESULTS,"> ".$fileresults) || die "SEQCOUNTFORMAT ERROR :: ".date()."Can't open '$fileresults': $!";
		#The header of the file
		print RESULTS join("\t",(sort keys %$files))."\n";
		
		#gathering the number od reads for each circRNA
		foreach my $circRNAs(sort keys %$hashref){
			#PRint the circRNAs
			print RESULTS $circRNAs ."\t";
			my @reads;
			foreach my $sample (sort keys %$files){
				#if the circRNAs has reads, print the all
				if(exists $hashref->{$circRNAs}->{$sample}){
					push (@reads,$hashref->{$circRNAs}->{$sample});
				}
				else{
					#or print 0
					push (@reads,0);
				}
			}
			print RESULTS join("\t",@reads). "\n";
		}
		close RESULTS;
		
		my $summary_path=$projectdir ."/circRNAs_results/";
		#Returning the path of the results file
		if(scalar(keys %$results)>0){
			open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
			print SUMM "\nReadCount [".$summary_path."]\n";
			print SUMM "Filename\tNumber of identified circRNAs\n";
			foreach my $processed_file (sort keys %$results){
				print SUMM $processed_file ."\t". scalar( keys %{$results->{$processed_file}})."\n";
			}
			close SUMM;
		}
		
		
		#Registering the error
		open(LOG,">> ".$logfile) || die "CIRI ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "CIRICOUNT :: ".date()." Please check the folder: $projectdir/$output_dir.\n";
		close LOG;
		
		print STDOUT "CIRICOUNT :: ".date()." Please check the folder: $projectdir/$output_dir.\n" if($verbose);
		
		return($fileresults);
		close LOG;
	}else{
		#Registering the error
		open(LOG,">> ".$logfile) || die "CIRI ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "CIRI ERROR :: ".date()." Projectdir ($projectdir), input paths($input) and/or logfile($logfile) have not been provided";
		close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn("CIRI ERROR :: ".date()." Projectdir ($projectdir), input paths($input) and/or logfile($logfile) have not been provided");
		help_CIRIFormat();
	}
	sub help_CIRIFormat{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
  	 		[projectdir] Input directory where htseqresults.tab will be saved
  	 		[input] Reference array with the paths of the files which have to be joined
   	 		[logfile] Path of run.log file where execution data will be saved
  	 
			Examples:
			CIRI(input=>"\@results", projectdir=>".", logfile=>"run.log");
	};

	print STDERR $usage;
	exit(); 
	}  	
}


=head2 miRDeepCount

  Example    : 
  	miRDeepCount( 
  		input=>"\@results", 
  		projectdir=>".",
  		logfile=>"run.log"
  	);
  Description: featureFormat collects the paths of the files processed by featureCount by a 
  reference array and create a combined tabulated file named htseqresults.tab in the 
  htseq_results directory with the number of reads of each file in each column. The names
  of the files will be the name of the original files without extension.
  Input parameters: 
	Mandatory parameters:
  	 [projectdir] Input directory where htseqresults.tab will be saved
  	 [input] Reference array with the paths of the files which have to be joined
   	 [logfile] Path of run.log file where execution data will be saved
  Returntype : Tab file at directory htseq_results and the path of the combined file
  Requeriments: featureCount function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Output files from htseq-count anlysis on the htseq_results directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub miRDeepCount{
	#Arguments provided by user are collected by %args. File, dir, adapter, projectdir and logfile
	# are mandatory arguments while verbose is optional. 
	my %args=@_;
	my $miARmaPath=$args{"miARmaPath"};
	my $file=$args{"file"}; #Path of the BAM file which is going to be processed
	my $verbose=$args{"verbose"}; #Optional arguments to show the execution data on screen
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $projectdir=$args{"projectdir"}; #Input directory where results directory will be created


	#Variable declaration
	my $command;
	my $commanddef;
	my $tab_file=undef;
	
	my $novel_miRNAs;
	my $known_miRNAs;
	
	#Checking the mandatory parameters
	if ($file and $projectdir and $logfile){
		if($file =~ /.*\.xls$/){
			my $output_dir="/"."miRDeep_results/";
			my($filename) = fileparse($file);
			$tab_file=$file;
			$tab_file=~s/\.xls//g;
			$filename=~s/\.xls//g;
			
			open(FILE,$file) || die "$! ($file)\n";
			my $starter=0;
			my $miRNAsCounts;
			my $total_miRNAs;
			
			while(<FILE>){
				chomp;
				if($_ =~/^novel/){
					$starter=1;
				}
				if($_ =~/^mature/){
					$starter=2;
				}
				if($_ =~/^\s+/){
					next;
				}
				if($_ =~/\#miRBase/){
					$starter=3;
				}
				if($starter==1){
					#Novel miRNAs
					if($_ !~/^provisional/ and $_ !~ /novel/){
						my @data=split(/\t/);
						my $read_count=$data[5];
						my $miR=$data[16];
						$read_count=~s/\.\d+//g;
						$miRNAsCounts->{$miR}+=$read_count if($miR);
						$total_miRNAs->{$miR}++;
						$novel_miRNAs->{$filename}->{$miR}++;
					}
				}
				elsif($starter==2){
					if($_ !~/^tag/){
						my @data=split(/\t/);
						my $miR=$data[9];
						my $read_count=$data[5];
						$read_count=~s/\.\d+//g;						
						$miRNAsCounts->{$miR}+=$read_count if($miR);
						$total_miRNAs->{$miR}++;
					}
				}
				elsif($starter==3){
					if($_ !~/^miRBase/){
						my @data=split(/\t/);
						my $miR=$data[0];
						my $read_count=$data[1];
						$read_count=~s/\.\d+//g;						
						#$miRNAsCounts->{$miR}+=$read_count if($miR);
						#$total_miRNAs->{$miR}++;
					}
				}
				else{
					next;
				}
			}
			close FILE;
			
			open(TAB,">$tab_file.tab") || die "$! ($tab_file.tab)\n";
			
			#miRDeep finds repeated miRNAs with diffrent number of reads, as they can come from a different precursos. So I get de averegare number:
			foreach my $miR (sort keys %$miRNAsCounts){
				#print STDERR "$miR\n";
				if($total_miRNAs->{$miR}>1){
					print TAB $miR . "\t" . sprintf("%.0f",($miRNAsCounts->{$miR}/$total_miRNAs->{$miR})) ."\n";
				}
				else{
					print TAB $miR ."\t" . $miRNAsCounts->{$miR} ."\n";
				}
			}
			close TAB;
			
			open(LOG,">> ".$logfile) || die "miRDeepCount ERROR :: ".date()."Can't open '$logfile': $!";
			print LOG "miRDeep :: ".date()." Known and novel miRNA count-file created ($tab_file.tab)\n";
			close LOG;
			
			print STDOUT "miRDeep :: ".date()." Known and novel miRNA count-file created ($tab_file.tab)\n" if($verbose);
			
			#After counting redas, we are going to conver the sam file into a sorted BAM file (a smaller file)
			#The path of the output file is returned to the main program
			if($tab_file){
				return($tab_file.".tab");
			}
			else{
				next;
			}
		}
	}
	else{
		#Registering the error
		open(LOG,">> ".$logfile) || die "miRDeepCount ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "miRDeepCount ERROR :: ".date()." Projectdir ($projectdir), file ($file), logfile ($logfile) have not been provided";
		close LOG;
			#If mandatory parameters have not been provided program dies and shows error message
		warn("miRDeepCount ERROR :: Projectdir ($projectdir), file ($file), logfile ($logfile) have not been provided");
		help_CIRICount();
	}
	sub help_miRDeepCount{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
	 		[file] Path of the file which is going to be processed (bam format)
	 		[logfile] Path of run.log file where execution data will be saved
	 		[projectdir] Directory where htseq_results directory will be created

	 		Optional parameters:
	 		[verbose] Option to show the execution data on the screen   

			Examples:
			miRDeepCount(file=>"file.bam", logfile=>"run.log", verbose=>"verbose", projectdir=>".");
	};

	print STDERR $usage;
	exit(); 
	}
}


=head2 miRDeepFormat

  Example    : 
  	miRDeepFormat( 
  		input=>"\@results", 
  		projectdir=>".",
  		logfile=>"run.log"
  	);
  Description: featureFormat collects the paths of the files processed by CIRI by a 
  reference array and create a combined tabulated file named htseqresults.tab in the 
  htseq_results directory with the number of reads of each file in each column. The names
  of the files will be the name of the original files without extension.
  Input parameters: 
	Mandatory parameters:
  	 [projectdir] Input directory where htseqresults.tab will be saved
  	 [input] Reference array with the paths of the files which have to be joined
   	 [logfile] Path of run.log file where execution data will be saved
  Returntype : Tab file at directory htseq_results and the path of the combined file
  Requeriments: featureCount function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Output files from htseq-count anlysis on the htseq_results directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub miRDeepFormat{
	
	#Arguments provided by user are collected by %args. Dir and input are mandatory arguments.
	my %args=@_;
	my $input=$args{"input"}; #Input is an referenced array with the paths of the files which have to be joined
	my $projectdir=$args{"projectdir"}; #Directory to create the output file
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $verbose=$args{"verbose"}; #write the execution data on screen
	
	
	#Variable declaration
	my $hashref;
	my $data;
	my $files;
	my @results;	

	if($input and $projectdir and $logfile){
		open(LOG,">> ".$logfile) || die "miRDeep ERROR :: ".date()."Can't open '$logfile': $!";
		
		#Declaration of the output path 
		my $output_dir="/DeNovo_ReadCount/";
		#Reading each path of the input referenced array
		foreach my $file(sort {$a cmp $b} @$input){
			print LOG "LOG :: Reading $file\n";
			if($file =~ /\.tab$/){
				#Fileparse is used to obtain the name of the file
				my $name=fileparse($file, qr{\.tab$});
				#Obtaining the process info from the file name and the original name of the file
				my $originalname=$name;
				my $suffix=$2;
				if($name !~/^\./){
					print "\tmiRDeep :: ". date() . " Reading $file\n" if($verbose);
					#Saving the Suffix with process information in a hash
					$files->{$originalname}++;
					#Opening each file of the array
					open (FILE, $file) or die "miRDeep ERROR :: ".date()."Can't open '$file': $!";
					#Reading the file. Each line contains the number of the miRNA, a tab and the number of reads
					while(<FILE>){
						#Deleting the last character of the line
						chomp;
						my ($miRNAs,$counts)=split(/\t/);
						$hashref->{$miRNAs}->{$originalname}=$counts;
						$data->{$originalname}++;	
					}
					#Closing the file
					close FILE;
				}
			}
			else{
				print LOG "WARM :: the file [$file] is not in the correct format\n";
			}

		}
		system("mkdir -p $projectdir/$output_dir");
		#Printing the results of the process with each combination of software
		my $fileresults= $projectdir.$output_dir."miRNAs.tab"; 
		#Opening the output file to write the data
		open(RESULTS,"> ".$fileresults) || die "miRDeep ERROR :: ".date()."Can't open '$fileresults': $!";
		#The header of the file
		print RESULTS join("\t",(sort keys %$data))."\n";
		#gathering the number od reads for each circRNA
		foreach my $circRNAs(sort keys %$hashref){
			#PRint the circRNAs
			print RESULTS $circRNAs ."\t";
			my @reads;
			foreach my $sample (sort keys %$data){
				#if the circRNAs has reads, print the all
				if(exists $hashref->{$circRNAs}->{$sample}){
					push (@reads,$hashref->{$circRNAs}->{$sample});
				}
				else{
					#or print 0
					push (@reads,0);
				}
			}
			print RESULTS join("\t",@reads). "\n";
		}
		close RESULTS;
		
		#Returning the path of the results file
		print STDOUT "miRDeep :: ".date()." Please check the folder: $projectdir/$output_dir.\n" if($verbose);
		
		open(LOG,">> ".$logfile) || die "miRDeep ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "miRDeep :: ".date()." Please check the folder: $projectdir/$output_dir.\n";
		close LOG;
				
		return($fileresults);
		close LOG;
	}else{
		#Registering the error
		open(LOG,">> ".$logfile) || die "miRDeep ERROR :: ".date()."Can't open '$logfile': $!";
		print LOG "miRDeep ERROR :: ".date()." Projectdir ($projectdir), input paths($input) and/or logfile($logfile) have not been provided";
		close LOG;

		#If mandatory parameters have not been provided program dies and shows error message
		warn("miRDeep ERROR :: ".date()." Projectdir ($projectdir), input paths($input) and/or logfile($logfile) have not been provided");
		help_miRDeepFormat();
	}
	sub help_miRDeepFormat{
	    my $usage = qq{
		  	$0 

			Mandatory parameters:
  	 		[projectdir] Input directory where htseqresults.tab will be saved
  	 		[input] Reference array with the paths of the files which have to be joined
   	 		[logfile] Path of run.log file where execution data will be saved
  	 
			Examples:
			miRDeepFormat(input=>"\@results", projectdir=>".", logfile=>"run.log");
	};

	print STDERR $usage;
	exit(); 
	}  	
}
sub miRDeepSummary{
	use File::Basename;
	#Arguments provided by user are collected by %args. Dir and input are mandatory arguments.
	my %args=@_;
	my $input=$args{"input"}; #Input is an referenced array with the paths of the files which have to be joined
	my $projectdir=$args{"projectdir"}; #Directory to create the output file
	my $summary_file=$args{"summary"}; #write the execution data on screen
	my $summary_path=$projectdir ."/miRDeep_results/";
	
	my $result;
	if($input and $projectdir and $summary_file){
		#Reading each path of the input referenced array
		foreach my $file(sort {$a cmp $b} @$input){
			open(TAB,$file) || warn "Can't find $file\n";
			my $filename=fileparse($file);
			$filename=~s/\.tab$//g;
			while(<TAB>){
				chomp;
				my ($miR,$read_number)=split(/\t/);
				if($miR =~ /^chr/ and $miR =~ /:/){
					#its novel
					$result->{$filename}->{"novel"}->{$miR}++;
				}
				else{
					#its Known
					$result->{$filename}->{"known"}->{$miR}++;
				}
			}
			close TAB;
		}
	}
	
	if(scalar(keys %$result)>0){
		open(SUMM,">>$summary_file") || warn "Can't create result file ($summary_file)\n";
		print SUMM "\nmiRDeep [".$summary_path."]\n";
		print SUMM "Filename\tNumber of novel miRNAs\tNumber of known miRNAs\n";
		foreach my $processed_file (sort keys %$result){
			print SUMM $processed_file ."\t". scalar( keys %{$result->{$processed_file}->{"novel"}}) ."\t". scalar( keys %{$result->{$processed_file}->{"known"}}) ."\n";
		}
		close SUMM;
	}
}

sub featureSummary{
	use File::Basename;
	
	#Arguments provided by user are collected by %args. Dir and input are mandatory arguments.
	my %args=@_;
	my $input=$args{"input"}; #Input is an referenced array with the paths of the files which have to be joined
	my $projectdir=$args{"projectdir"}; #Directory to create the output file
	my $logfile=$args{"logfile"}; #Path of the logfile to write the execution data
	my $summary_file=$args{"summary"};
	
	open(STAT,$logfile);
	my $real_file;
	my $processed;
	my $assigned;
	my $strand;
	my $summary;
	my $summary_path=$projectdir ."/Readcount_results/";
	
	while(<STAT>){
		chomp;
		if($_ =~/^SEQCOUNT :/){
			if($_ =~/Reading counts from/){
				my $file=$_;
				my @data=split(/\s+/);
				$real_file=fileparse($data[$#data]);
				$real_file=~s/\.tab$//g;
				$real_file=~s/\.bam$//g;
				$real_file=~s/\.sam$//g;
				
				$processed=0;
				$assigned=0;
				$strand="";
			}
		}
		if($_ =~ /Strand specific/){
			$strand=$_;
			$strand=~s/.*Strand specific : (.+)\s+.*/$1/g;
			$strand=~s/\s+//g;
		}
		if($_ =~ /Total reads/){
			$processed=$_;
			$processed=~s/.*Total reads : (\d+).*/$1/g;
		}
		if($_ =~ /Total fragments/){
			$processed=$_;
			$processed=~s/.*Total fragments : (\d+).*/$1/g;
		}
		if($_ =~ /Successfully assigned/){
			$assigned=$_;
			$assigned=~s/ reads //g;
			$assigned=~s/ fragments //g;
			$assigned=~s/.*Successfully assigned: (\d+) (\(\d+\.\d+%\)).*/$1 $2/g;

		}
		
		if($real_file and $assigned){
			$summary->{$real_file}="$processed\t$assigned\t$strand";
		}
	}
	close STAT;

	
	my $results;
	foreach my $file(sort {$a cmp $b} @$input){
		if($file =~ /\.tab$/){
			open(TAB,$file) || warn "Can't find $file\n";
			my $filename=fileparse($file);
			$filename=~s/\.tab$//g;
			while(<TAB>){
			 	chomp;
				if($_ !~ /^#/ and $_ !~ /^Geneid/){
					#Using the split with the tab we can keep the name of the miRNA in $miRNA variable
					#and the number of reads in $reads variable. 
					my ($miRNAs,$chr,$start,$end,$strand,$length,$reads)=split(/\t/);
					#Creating a hash with the names of miRNAs
					if($reads>0){
						$results->{$filename}->{$miRNAs}++;
					}
				}	
			}
			close TAB;
		}
	}
	if(scalar(keys %$summary)>0 and scalar(keys %$results)>0){
		open(SUMM,">>$summary_file") || warn "Can't create summary file ($summary_file)\n";
		print SUMM "\nReadCount [".$summary_path."]\n";
		print SUMM "Filename\tProcessed Reads\tAssigned reads\tStrand\tNumber of identified entities\n";
		foreach my $processed_file (sort keys %$summary){
			print SUMM $processed_file ."\t". $summary->{$processed_file} ."\t". scalar(keys %{$results->{$processed_file}}) ."\n";
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

