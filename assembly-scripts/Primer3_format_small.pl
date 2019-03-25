#!/usr/bin/perl
# Primer3_format.pl
# by Nate
# Create a primer3 input file from a fasta file of consensus seuqences...
#
# Fasta file names should be in the format "Species_XXXXXXX-##" where "XXXXXXX" identifies the locus and
# "##" identifies the sequencing target.

use strict; use warnings;

die "usage: supply a fasta file as input\n" unless @ARGV == 1;

open (FASTA, "<$ARGV[0]") or die "Error opening $ARGV[0]";

while (<FASTA>) {
	my $name_line = $_;
	my $seq_line = <FASTA>;
	chomp $name_line;
	my @nameINFO = split "-", $name_line;
	my $name = $nameINFO[0];
	$name =~ s/>//;
	my $SNP_site = $nameINFO[1];
	my $Seq_site = $SNP_site - 1;
	chomp $seq_line;

print 
"SEQUENCE_ID=$name-$SNP_site
SEQUENCE_TEMPLATE=$seq_line
SEQUENCE_TARGET=$Seq_site,1
PRIMER_TASK=pick_detection_primers
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=17
PRIMER_MAX_SIZE=25
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_PRODUCT_SIZE_RANGE=50-90
PRIMER_EXPLAIN_FLAG=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/share/PI/spalumbi/programs/primer3-2.3.7/src/primer3_config/
PRIMER_NUM_RETURN=1
PRIMER_LEFT_NUM_RETURNED=1
PRIMER_RIGHT_NUM_RETURNED=1
=\n";
}
