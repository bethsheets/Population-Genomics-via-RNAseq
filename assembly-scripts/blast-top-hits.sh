#!/bin/bash
#usage: bash blast-top-hits.sh <parsed_infile> outfile.txt 

awk '!x[$1]++' $1 > $2
