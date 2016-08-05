#########################################################################	
#	Benchmarking  package			 									#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2014 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################
package CbBio::RNASeq::Benchmarking;
require Exporter;
#Export package system
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(mergehtseqcount);

use strict;
use DateTime;
use File::Basename;

=head1 NAME

 Benchmarking

=head1 SYNOPSIS

Benchmarking is a module to compare the results obtained by differents RNAseq analysis methods. 
It is composed of mergehtseqcount subroutine which collects the paths of the files processed 
by htseqCount by two reference arrays and create a combined tabulated file named mergehtseqresults.tab 
in the directory with the number of reads of each file in each column.  


=head1 Methods

=head2 mergehtseqcount

  Example    : 
  	mergetseqcount( 
  		input1=>\@htseq1output,
		input2=>\@htseq2output,
		dir=>"/Users/directory/"	
  	)
  Description: mergehtseqcount collects the paths of the files processed by htseqCount by two 
  reference arrays and create a combined tabulated file named mergehtseqresults.tab in the 
  directory with the number of reads of each file in each column.  
  Input parameters: 
	Mandatory parameters:
  	 [dir] Input directory where htseqresults.tab will be saved
  	 [input1] Reference array with the paths of the files which have to be joined. These
  	 files will appear with the extension _bw1 in the results file
  	 [input2] Reference array with the paths of the files which have to be joined. These
  	 files will appear with the extension _bw2 in the results file
  Returntype : File at the provided directory and the path of the results file 
  Requeriments: htseqCount function requires for a correct analysis:
  	- Perl v5.10.0 or higher software correctly installed
  	- Input files from htseq-count anlysis  
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub mergehtseqcount{
	
	#Arguments provided by user are collected by %args. Dir, input1 and input2 are mandatory arguments.
	my %args=@_;
	my $input1=$args{"input1"}; #Input1 is an referenced array with the pathways of the files of bowtie1 which have to be joined
	my $input2=$args{"input2"}; #Input2 is an referenced array with the pathways of the files of bowtie2 which have to be joined
	my $dir=$args{"dir"}; #Directory to create the output file 
	
	#Declaration of the output path 
	my $fileresults= $dir."/mergehtseqresults.tab"; 
	#Variable declaration
	my $hashref;
	my $data;
	my $files;	

	#Reading each path of the input referenced array
	foreach my $name(@$input1){
		#Fileparse is used to obtain the name of the file
		my $file=fileparse($name);
		$files->{$file."_bw1"}++;
		#Opening each file of the array
		open (FILE, $name) or die "Can't open '$file': $!";
		#Reading the file. Each line contains the number of the miRNA, a tab and the number of reads
		while(<FILE>){
			#Deleting the last character of the line
			chomp;
			#Using the split with the tab we can keep the name of the miRNA in $miRNA variable
			#and the number of reads in $reads variable. 
			my ($miRNAs,$reads)=split(/\t/);
			#Creating a hash with the names of miRNAs
			$hashref->{$miRNAs}++;
			#Creating a hash with all the data, the first dimension contains the miRNA name,
			#the second the name of the file and the third is the number of reads 
			$data->{$miRNAs}->{$file."_bw1"}=$reads;		
		}
		#Closing the file
		close FILE;
	}
	
	#Reading each path of the input referenced array
	foreach my $name(@$input2){
		#Fileparse is used to obtain the name of the file
		my $file=fileparse($name);
		$files->{$file."_bw2"}++;
		#Opening each file of the array
		open (FILE, $name) or die "Can't open '$file': $!";
		#Reading the file. Each line contains the number of the miRNA, a tab and the number of reads
		while(<FILE>){
			#Deleting the last character of the line
			chomp;
			#Using the split with the tab we can keep the name of the miRNA in $miRNA variable
			#and the number of reads in $reads variable. 
			my ($miRNAs,$reads)=split(/\t/);
			#Creating a hash with the names of miRNAs
			$hashref->{$miRNAs}++;
			#Creating a hash with all the data, the first dimension contains the miRNA name,
			#the second the name of the file and the third is the number of reads 
			$data->{$miRNAs}->{$file."_bw2"}=$reads;		
		}
		#Closing the file
		close FILE;
	}
	
	#Opening the output file to write the data
	open(RESULTS,"> ".$fileresults) || die $!;
	#Printing on the first row of the results file the name of the input files sorted and 
	#separated by a tab 
	print RESULTS join("\t",(sort keys %$files))."\n";
	#The next lines will be printed going over the hash. Reading the first dimension of the hash
	foreach my $miRNAs(sort keys %$data){
		#At the end of the input files appeared some useless information to our fileresults  
		if (($miRNAs =~ "ambiguous") or ($miRNAs =~ "no_feature") or ($miRNAs =~ "alignment_not_unique")  or ($miRNAs =~ "not_aligned") or ($miRNAs =~ "too_low_aQual")){
		}
		else{
			#Printing the name of each miRNA
			print RESULTS $miRNAs."\t";
			my @results;
			#An array is used to save the number of reads of each file for a defined miRNA.  
			foreach my $name(sort keys %{$data->{$miRNAs}}){
		 		push(@results,$data->{$miRNAs}->{$name});
			}
			#Printing the number of the reads of each miRNA separated by a tab 
			print RESULTS join("\t",@results)."\n";
		}
	}
	#Closing the file
	close RESULTS;
	return($fileresults);
}
sub date{
	#my $dt = DateTime->now(time_zone=>'local');
	#return($dt->hms . " [" . $dt->dmy ."]");
	use Time::localtime;
	my $now = ctime();
	return("[$now]");
}	

1;
