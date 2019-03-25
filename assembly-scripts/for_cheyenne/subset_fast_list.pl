#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# This program takes a list of desired contigs and an input fasta and subsets the fasta based on the contigs that were present in the list

my ($fasta_file, $list_file, $Out_fasta) = @ARGV;

open(my $FASTA, "<", $fasta_file) or die "Cannot open $fasta_file\n";
open(my $LIST, "<", $list_file) or die "Cannot open $list_file\n";
open(my $OUT_FASTA, ">", $Out_fasta) or die "Cannot create $Out_fasta\n";

# Store list into memory
my %list;
my $list_count=0;

while(defined(my $readline = <$LIST>)){
	chomp $readline;
	my $query = $readline;
	# Store unique
	if(!defined $list{$query}){
		$list{$query}=[];
		$list_count++;
	}
}


my $printseq=0;
my $test_num = 0;

my $num_not = 0;
while(defined(my $readline = <$FASTA>)){
	my $test = substr $readline, 0, 1;
	if($test eq ">"){
		# We're in a header
		chomp $readline;
		# Get header w/out >
		my @linearray = split(/\s/,$readline);
		my $header = substr $linearray[0], 1, length($linearray[0])-1; # Can also use s///
		if(exists $list{$header}){
			# Print this sequence to the output fasta
			#print "$header\n";
			$test_num++;
			$printseq = 1;
			print $OUT_FASTA $readline, "\n";
		}else{
			$num_not++;
			$printseq = 0;
		}
	}elsif($printseq == 1){
		print $OUT_FASTA $readline;
	}
}

print "Your list has $list_count contigs & you found $test_num in your original fasta\n";



