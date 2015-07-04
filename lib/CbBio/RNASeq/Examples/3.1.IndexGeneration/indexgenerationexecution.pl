#!/usr/bin/perl

use strict;
use warnings;
use CbBio::RNASeq::Aligner;

#Defining bowtie1_index variables
my $fasta="/Users/apple/Pipeline/examples/Bowtie1_index/hg19fasta/hg19.fa";
my $dir="/Users/apple/Pipeline/examples";
my $logfile="/run_".$$.".log";
my $indexname="hg19";
my $aligner="Bowtie1-Bowtie2";

# INDEXGENERATION EXECUTION
IndexGeneration(
  	aligner=>$aligner,
  	fasta=>$fasta,
  	dir=>$dir,
  	logfile=>$dir.$logfile,
  	indexname=>$indexname
  );