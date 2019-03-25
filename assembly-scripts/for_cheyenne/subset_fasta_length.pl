#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# Program to take input fasta and create new fasta with sequences above a certain length 

my ($fasta_file, $out_file, $length_threshold) = @ARGV;

my %out_fasta;
my $seqlength;
my @seqarray;
my $name;
my $filelength;
my $length = 0;
# Get file length to store last squence 
open(my $CHECK, "<", $fasta_file);
$filelength++ while <$CHECK>;
close $CHECK;

open(my $FASTA, "<", $fasta_file) or die "Cannot open $fasta_file\n";
open(my $OUT, ">", $out_file) or die "Cannot create $out_file\n";

while(defined(my $readline = <$FASTA>)){
	chomp $readline;
	if ($. == $filelength && $length >= $length_threshold){
  		# Store last sequence of file
  		$out_fasta{$name} = [@seqarray]; 	
    }
	my $test = substr $readline, 0, 1;
	if($test eq ">"){
		# If not the first sequence & above th, store the last sequence
		if(scalar(@seqarray) > 0 && $length >= $length_threshold){
			$out_fasta{$name} = [@seqarray];
			#print $name, "\n", join("", @seqarray), "\n$length\n\n\n\n";
			#print "$name\n\n$length\n\n\n\n";
		}
		# We are in a header
		$name = $readline;
		@seqarray = ();
		$length = 0;
	}else{
		# We're in sequence. Store & count
		$length = $length + length($readline); 
		push(@seqarray, $readline);
	}
}


foreach my $key (keys %out_fasta){
	print $OUT $key, "\n", join("\n",@{$out_fasta{$key}}), "\n";
}


