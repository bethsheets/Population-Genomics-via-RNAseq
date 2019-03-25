#!/bin/bash

# Will modify, and delete modification in, all specified files
# to bypass new Sherlock 2.0 policy of deleting files in SCRATCH
# folders if they haven't been modified in 6 months
# CYP 09/14/2017
# 
# USAGE:     ./prevent_deletion.sh 'files_to_save' [recurse]
# example 1: ./prevent_deletion.sh '*.txt'
#             will save all .txt files in the current directory
# example 2: ./prevent_deletion.sh '*'
#             will save every file in the current directory
# example 3: ./prevent_deletion.sh '*.txt' recurse
#             will save all .txt files in the current directory 
#             and in every subdirectory
# 
# NOTE: this script will not modify tar-ed or zipped files, as
#       modifying compressed files will corrupt the compression
#       To update these files, you want to copy each file into
#       new file.

### FUNCTIONS ###

# Function: main, takes $FILE_TYPE argument
# Simply saves target files in current working
# directory
function main {
  TO_SAVE=$@
  for file in $TO_SAVE; do 
    # if encounter a zipped or tarred file, leave it alone!
    if [[ $file =~ \.(gz|tar|zip)$ ]]; then
      if [[ $VERBOSE == true ]]; then echo "Did not modify $file"; fi
      continue
    elif [[ -f "$file" ]]; then
      save $file
    fi
  done
}

# Function: recurse, takes * argument
# Using specified file group, recursively saves files
# in working directory and every subdirectory within it
function recurse {
  PATTERN="*$FILE_TYPE"
  ALL=$@
  for file in $ALL; do
    # if encounter a zipped or tarred file, leave it alone!
    if [[ $file =~ \.(gz|tar|zip)$ ]]; then
      if [[ $VERBOSE == true ]]; then echo "Did not modify $file"; fi
      continue
    elif [[ -d "$file" ]]; then
      # if encounter a directory, enter the dir and call
      # this function again on the target file group
      recurse $file/*
    elif [[ $file == $PATTERN ]]; then 
      # otherwise (i.e. if it's a target file), save
      save $file
    fi
  done
}

# Function: save, takes file to save argument
# Modifies specified file, and deletes modification so that
# file appears updated without any real change to the file
function save {
  if [[ $VERBOSE == true ]]; then echo -n "Saving ${1}..."; fi
  echo "# Saving this file!" >> $1
  sed -i '$d' $file
  if [[ $VERBOSE == true ]]; then echo " Saved!"; fi
}


### BODY ###

# To display names of files being saved, set VERBOSE=true
VERBOSE=true

# check if recursive file saving is desired
if [[ $# -eq 2 ]] && [[ $2 == "recurse" ]]; then
  re=true
else 
  re=false
fi

FILE_TYPE=$1
if [[ $re == true ]]; then
  echo "Saving files in this directory and all subdirectories..."
  recurse * 
else
  echo "Saving files in this directory only..."
  main $FILE_TYPE 
fi

echo "All done! :D"
