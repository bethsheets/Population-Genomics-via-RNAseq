#!/usr/bin/perl
# P3toIDT.pl
# by Amanda Matala
# takes native output from primer3 and formats for ordering primers from IDT
# in 96w format. Also, adds read1 seq tag to fwd sequence and read2 seq tag to rev sequence.
# Script supports no more than 384 loci in primer3 input file. Specify input file on the command line.
use strict; use warnings;
#edited by Beth Sheets Jan 2018, to call: perl P3toIDT_eas.pl Primer3out.txt Species_


my $FILENAME= $ARGV[0];
my $REPLACE= $ARGV[1];

open(FILE, "$FILENAME") or die "you forgot your file...";
my @left;
my @right;
my @locus;
my @seq;
while (<FILE>) {
  push (@locus, $&) if m/$REPLACE.*[0-9]/;
  push (@seq, $&) if m/[ACGT]+$/;
}
close FILE;

# Split the forward sequences from the reverse sequences and store in separate arrays

for (my $i = 1; $i <= @seq; $i = $i + 2) {
  push (@right, $seq[$i]);
}
for (my $i = 0; $i <= @seq; $i = $i + 2) {
  push (@left, $seq[$i]);
}
my @well = qw(A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 F1 F2 F3 F4 F5 F6 F7 F8 F9 F10 F11 F12 G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12);
my $length = @locus;
my $limit = $length - 1;

# Print the well, locus, and read1 seq. tag + the forward target sequence for the first 96 loci

print "Name\tFWD_Primer\tREV_Primer\n";
for (my $i = 0; $i <= $limit; $i = $i + 1) {
  print "$locus[$i]\tCGACAGGTTCAGAGTTCTACAGTCCGACGATC$left[$i]\tGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT$right[$i]\n";
}

