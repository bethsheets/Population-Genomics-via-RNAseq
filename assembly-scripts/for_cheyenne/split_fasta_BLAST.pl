#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;

# Very similar to subset_fastq_BLAST.pl - instead of subsetting by the list, it splits into two files 
# This program takes an output custom file from BLAST where the first and third column of a csv are the query and hit contig.
# It takes an input fastq and subsets the fastq based on the contigs that were ***NOT*** present in the BLAST output
# It creates a second file telling the user how many unique query & unique hits

my ($fastq_file, $BLAST_file, $Out_from_fastq, $Out_not_fastq, $Out_stats) = @ARGV;

open(my $FASTQ, "<", $fastq_file) or die "Cannot open $fastq_file\n";
open(my $BLAST, "<", $BLAST_file) or die "Cannot open $BLAST_file\n";
open(my $OUT_FROM_FASTQ, ">", $Out_from_fastq) or die "Cannot create $Out_from_fastq\n";
open(my $OUT_NOT_FASTQ, ">", $Out_not_fastq) or die "Cannot create $Out_not_fastq\n";
open(my $OUT_STATS, ">", $Out_stats) or die "Cannot create $Out_stats\n";

# Store BLAST into memory, it will always be smaller
my %queries;
my %hits;
my $hit_count=0;
my $query_count=0;

while(defined(my $readline = <$BLAST>)){
	chomp $readline;
	my @linearray = split(/,/,$readline);
	my $query = $linearray[0];
	my $hit = $linearray[2];
	# Store unique hits, unique queries
	if(!defined $queries{$query}){
		$queries{$query}=[];
		$query_count++;
	}
	if(!defined $hits{$hit}){
		$hits{$hit}=[];
		$hit_count++;
	}
}
my $printseq=0;
my $test_num = 0;

my $num_not = 0;
while(defined(my $readline = <$FASTQ>)){
	my $test = substr $readline, 0, 1;
	if($test eq ">"){
		# We're in a header
		#print $readline;
		chomp $readline;
		# Get header w/out >
		my @linearray = split(/\s/,$readline);
		my $header = substr $linearray[0], 1, length($linearray[0])-1; # Can also use s///
		if(exists $queries{$header}){
			# Print this sequence to the output fastq
			#print "$header\n";
			$test_num++;
			$printseq = 1;
			print $OUT_FROM_FASTQ $readline, "\n";
		}else{
			print $OUT_NOT_FASTQ $readline, "\n";
			$num_not++;
			$printseq = 0;
		}
	}elsif($printseq == 1){
		print $OUT_FROM_FASTQ $readline;
	}elsif($printseq == 0){
		print $OUT_NOT_FASTQ $readline;
	}
}

print "You're BLAST has $query_count queries & you found $test_num in your original fastq\n";
print $OUT_STATS "Number of unique from queries: $query_count\nNumber of unique hits: $hit_count\nNumber of queries w/out a hit: $num_not\n\nList of queries followed by 20* then list of hits\n\n\n";
print $OUT_STATS "$_\n" for keys %queries;
print $OUT_STATS "\n********************\n\n";
print $OUT_STATS "$_\n" for keys %hits;

