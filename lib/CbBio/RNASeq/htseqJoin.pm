#########################################################################	
#	htseqJoin execution package 	 									#
#																		#
#	Created at Computational Biology and Bioinformatics Group (CbBio)	#
#	Institute of Biomedicine of Seville. IBIS (Spain)					#
#	Copyright (c) 2014 IBIS. All rights reserved.						#
#	mail : miARma-devel@cbbio.es 										#
#########################################################################
package CbBio::RNASeq::Join;
require Exporter;
$|=1;
@ISA=qw(Exporter);
@EXPORT=qw(htseqJoin);

use strict;
use DateTime;

=head1 NAME

 Readcount

=head1 SYNOPSIS

htseqJoin is function to crate a combined file as of the input files provided by the user. 
These inputs files have to be the output of htseq-count programme. These files contains the 
analysed element (gene_id, miRNA,...) and the number of the reads of each element. htseqJoin
read the directory provided and open each file on it (except from hidden files). Then the 
useful information is kept to write it on the results file. This file presents a first row 
with the name of the processed files and below them the number of reads of each element. 

=head1 Methods

=head2 htseqFormat

  Example    : htseqJoin(dir=>/Users/Project)
  Description: htseqJoin collects the directory where are the input files and read it 
  to open each file on it (except from hidden files).The function creates a combined 
  tabulated file named htseqcombined.tab in the provided directory with the number 
  of reads of each file in each column.  
  Returntype : File at provided directory 
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub htseqJoin{
	
	#Arguments provided by user are collected by %args. Dir and input are mandatory arguments.
	my %args=@_;
	my $dir=$args{"dir"}; #Directory to be analysed

	#Declaration of the output path 
	my $fileresults= $dir."/htseqcombined.tab"; 
	#Variable declaration
	my $hashref;
	my $data;
	my $files;	
	
	# Opening and reading the directory supplied by the user
	opendir(DIR, $dir) || die $!; 
	my @files= readdir(DIR); 

	#Reading each file of the directory
	foreach my $file( sort {$a cmp $b} @files){
		if ($file =! /^\./){
			#Saving the files name in a hash
			$files->{$file}++;
			#Opening each file of the directory
			open (FILE, $dir.$file) or die "Can't open '$file': $!";
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
				$data->{$miRNAs}->{$file}=$reads;		
			}
			#Closing the file
			close FILE;
		}
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
			foreach my $file(sort keys %{$data->{$miRNAs}}){
		 		push(@results,$data->{$miRNAs}->{$file});
			}
			#Printing the number of the reads of each miRNA separated by a tab 
			print RESULTS join("\t",@results)."\n";
		}
	}
	#Closing the file
	close RESULTS;
	
}

sub date{
	#my $dt = DateTime->now(time_zone=>'local');
	#return($dt->hms . " [" . $dt->dmy ."]");
	use Time::localtime;
	my $now = ctime();
	return("[$now]");
}

1;
