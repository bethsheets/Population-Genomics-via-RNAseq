#!/usr/bin/env python

''' Script to extract full contig, based on a name in header 
    usage: ./extract_contig_by_name.py seq_file.fasta contig_name 
    CYP -- 08/16/2017
'''

import sys
import re

seq_name = sys.argv[2]
found = False
with open(sys.argv[1], 'r') as f:
  for line in f: 
    if re.search(">", line):
      if found: 
        out.close()
        break
      if re.search(seq_name, line):
        found = True
        filename = seq_name+'.fasta'
        out = open(filename, 'w')
    if found: out.write(line)

