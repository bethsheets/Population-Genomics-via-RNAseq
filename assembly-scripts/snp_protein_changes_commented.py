#!/usr/bin/env python
#usage: coding.py vcf ref.fa orfpredictor_headers.txt OUT
#to get orfpredictor_headers.txt run orfpredictor on ref.fa
# then grep ">" orf.pep | sed 's/>//g' | sed 's/+//g' > orfpredictor_headers.txt


import sys
import csv
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# collect all scaffold sequences
coding=SeqIO.parse(sys.argv[2],'fasta',generic_dna)
# and place in dictionary called transcripts
transcripts=SeqIO.to_dict(coding)

OUTFILE='CHROM,POS,ID,REF,ALT,REFprot,ALTprot,EFFECT,CPOS'

frames=dict()
with open(sys.argv[3],'r') as framefile:
  reader = csv.reader(framefile,delimiter='\t')
  # place all ORF info associated with each contig/scaffold 
  # name in dictionary
  for contig,frame,start,end in reader:
    #XXX frames[contig]=(frame,start,end)
    frames[contig[1:]]=(frame,start,end) #XXX - removes '>' from scaf name
 
with open(sys.argv[1],'rU') as vcffile:
  # collect each snp (each line not headed by #) from vcf file,
  # delimit, and place in snps list
  snps = csv.reader(vcffile,delimiter='\t',quotechar='"')
  for snp in snps:
    if snp[0].startswith('#'):
      continue
    transcript=transcripts[snp[0]]
    # if the scaffold/contig name is not in the ORF dictionary,
    # specify scaffold/contig as having no ORF
    if not snp[0] in frames:
      OUTFILE=OUTFILE+'\n'+','.join(snp[0:5])+',NA,NA,noORF,NA'
      continue

    # all info that is written to out file comes from ORF file
    # the list frame holds the (frame,start,end) value for the
    # current snp in snps
    frame=frames[snp[0]]
    pos=int(snp[1])-1

    tstart=int(frame[1])-1
    tstop=int(frame[2])
    # except for alternate base, which comes from filtered snp file
    altbase=Seq(snp[4],generic_dna)
    framesign=1
    # if the frame position is a negative value, reverse compliment
    # the transcript and alternate base 
    if int(frame[0])<0:
      transcript=transcript.reverse_complement()
      pos=len(transcript)-pos-1
      altbase=altbase.reverse_complement()
      framesign=-1

    # if the position of the snp is outside the start-stop of the 
    # ORF, stop loop and label this snp as 'pUTR' (in the UTR)
    if pos < tstart or pos >= tstop:
      OUTFILE=OUTFILE+'\n'+','.join(snp[0:5])+',NA,NA,pUTR,NA'
      continue

    transcript=transcript[tstart:tstop]
    pos=pos-tstart
    cpos=pos%3
    cstart=pos-cpos
    cstop=pos-cpos+3
    codon=transcript[cstart:cstop]

    # translate ref and alt codons into proteins
    ref=codon.seq.translate()
    altcodon=codon[0:cpos]+altbase+codon[cpos+1:3]
    alt=altcodon.seq.translate()
    effect='S'
    # if the reference and alternate proteins are not the same,
    # specify as a nonsynonymous mutation (NS)
    if ref != alt:
      effect='NS'
    OUTFILE=OUTFILE+'\n'+','.join(snp[0:5])+','+str(ref)+','+str(alt)+','+effect+','+str((1+cpos)*framesign)
with open(sys.argv[4],'w') as f:
  f.write(OUTFILE)
