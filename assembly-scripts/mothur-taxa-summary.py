#!/usr/bin/env python

'''
    Uses MOTHUR output files (taxonomy and counts matrix files) to generate a
    new counts matrix organized by classifications in a specified taxonomic 
    level. For example, if Level 2 is specified, then sequences and their counts 
    will be combined in their respective Level 2 taxonomic classification.

    usage: ./mothur-taxa-summary.py taxonomy_file.csv counts_matrix.csv level_#
    ex. ./mothur-taxa-summary.py taxonomy_file.csv counts_matrix.csv 2

    output: counts-summary.csv - contains counts matrix
                                 columns: samples
                                 rows: combined level classification counts

    NOTE for running on Sherlock 2.0:
      python is not installed by default on Sherlock 2.0, so before running 
      check that python is available. If it is not loaded, run the following:
          module load python

    CYP 01/17/2017
'''

import sys
import csv

# DEBUG mode is OFF
DEBUG = False 

# global variables
seq2tax_dict = {}         # { sequence_ID : taxon_name }
counts_dict = {}          # { taxon_name : { sample-index : count } }  
tax_seq_num_dict = {}     # { taxon_name : number seqs in taxon }

'''
    Extract taxon information from taxonomy_file
    - Generates a dictionary with { sequence_ID : taxon } 
      structure, so that the taxonomic classification for
      the user specified level is accessible with seq_ID
'''
def collect_taxon_info():
  print("taxonomy file processing...")
  with open(taxonomy_infile, 'r') as tax_file:
    # create dict with taxon_name : [ sequences with that
    # taxonomic classification ]
    rdr = csv.reader(tax_file)
    counts_row_list = list(rdr)
    for row in counts_row_list:
      # sequence name is first item in row
      seq_name = row[0]
      seq_ID = seq_name.replace("-", "_")

      # the key is the taxon_name at desired level
      taxon = (row[level].split('('))[0]

      # seq2tax dictionary has format { seq_name : taxon }
      seq2tax_dict[seq_ID] = taxon
  

'''
    Collect count info from counts matrix file
    - Generates a dictionary of dictionaries with 
      { taxon : { sample-index : count } }  
      structure, so that a sample's count can be
      incremented with the counts from multiple 
      sequences, which are classified by the overall
      taxon

'''
def collect_counts_info():
  print("counts file processing...")
  with open(counts_infile, 'r') as counts_file:
    # check that the taxonomy dictionary isn't empty
    if not seq2tax_dict:
      print("ERROR: Your taxonomy dictionary is empty!")
      exit()
    rdr = csv.reader(counts_file)
    counts_row_list = list(rdr)
    sample_list = counts_row_list[0][1:]

    # using sample indices, collect counts
    for row in counts_row_list[1:]:

      # collect seq_ID from first column
      seq_name = row[0]
      seq_ID = seq_name.replace("-", "_")

      if seq_ID not in seq2tax_dict: continue
      else: tax = seq2tax_dict[seq_ID] 
       
      # keep track of number of seqs per taxon
      if tax not in tax_seq_num_dict:
        tax_seq_num_dict[tax] = 1
      else:
        tax_seq_num_dict[tax] += 1

      if DEBUG and tax == "Phycisphaerae": print(seq_ID)

      # add new taxon names to counts dictionary
      if tax not in counts_dict: counts_dict[tax] = {}

      # loop through indices in row
      for index,count in enumerate(row):
        # skip first index (sequence name)
        if index == 0: continue
        # if index not in taxon dict, insert with new count
        if index not in counts_dict[tax]:
          counts_dict[tax][index] = float(count) 
        # if index is in dict, add count to existing count value
        else:
          counts_dict[tax][index] += float(count)
        #if DEBUG and tax == "Phycisphaerae" and index == 4: print(float(count))
  return sample_list


'''
    Write all counts to outfile called 
    - counts-summary.csv
'''
def make_outfile(sample_list):
  with open("counts-summary.csv", 'w') as outfile:
    wr = csv.writer(outfile)
    # write header with sample names
    sample_list.insert(0,'Taxon')
    sample_list.append('Number_seqs_in_taxon')
    wr.writerow(sample_list)

    # add counts to new row 
    # taxon is a dictionary with format { index : count }
    for taxon in counts_dict:
      # create new row
      new_row = [taxon]
      for index in range(len(sample_list)-1):
        if index == 0: continue
        # add each count to respective index in row
        new_row.append(str(counts_dict[taxon][index]))
      # add total number of sequences count to end of row
      new_row.append(str(tax_seq_num_dict[taxon]))
        
      wr.writerow(new_row)    


### MAIN ###

'''
    Check that correct number of arguments are passed in
'''
if len(sys.argv) != 4:
  print("Please include 2 csv files and a level number.") 
  quit()
else:
  taxonomy_infile = sys.argv[1]
  counts_infile = sys.argv[2]
  level = int(sys.argv[3])
 
print("Generating a taxa count summary for level: " + str(level))
collect_taxon_info()
sample_list = collect_counts_info()
make_outfile(sample_list)
print("See you next time! :D")
