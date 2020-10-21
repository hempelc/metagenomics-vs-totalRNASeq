#!/bin/bash

# Script to set up a SILVA SSU and LSU NR99 BLAST DB

# Usage: ./SILVA_SSU_LSU_makeblastdb_preparation.sh


# Get SILVA fasta and taxonomy files
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz
gunzip *.gz

# Get NCBI taxonomy files
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip

# Edit SILVA taxonomy file:
# SILVA taxonomy file needs to have two columns, first one with accession number
# and second one with taxonomy path
tail -n +2 taxmap_slv_ssu_ref_nr_138.1.txt > taxmap_slv_ssu_lsu_ref_nr_138.1.txt
tail -n +2 taxmap_slv_lsu_ref_nr_138.1.txt >> taxmap_slv_ssu_lsu_ref_nr_138.1.txt
sed "s/ <[a-zA-Z -,.&:'0-9]*>//g" taxmap_slv_ssu_lsu_ref_nr_138.1.txt \
| sed 's/;\t/;/g' | cut -f 1,4 \
> taxmap_slv_ssu_lsu_ref_nr_138_edited_for_NCBI_staxid_script.txt
  # 1. concatenate ssu and lsu taxmap files into one ssu_lsu file with one header
  # 2. removes <genus>, <family> etc. from taxonomic paths (otherwise the SILVA
  #    taxonomy won't match the NCBI taxonomy)
  # 3. removes the last tab of each line so that taxonomy path is in one column
  # 4. removes alignment numbers

# Edit names.dmp file into a file only containing scientific names
sed "s/ <[a-zA-Z -,.&:'0-9]*>//g" names.dmp | grep 'scientific name' \
| cut -f 1,3 | awk -F $'\t' ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' \
| grep -v 'environmental\|uncultured\|unidentified\|metagenome' \
> NCBI_staxids_scientific.txt
  # 1. Remove <genus>, <family>, and other strings in <> brackets from taxonomy
  #    (otherwise the NCBI taxonomy won't match the SILVA taxonomy)
  # 2. Extract only scientific names and staxids
  # 3. Cut out columns we need
  # 4. Invert columns for script to work
  # 5. Removes lines containing "environmental", "uncultured", "unidentified",
  #    and "metagenome"
    # Needed because
      # SILVA taxonomy can have the same taxonomic ranks (e.g., "environmental
      # sample") for different higher ranks (e.g., "nematode; environmental
      # sample" and "bacteria;environmental sample"), which would, however, be
      # assigned to the same staxid because the lower rank "environmental sample"
      # is similar
      # NCBI taxonomy can have different staxids for the same taxonomic name,
      # which will cause issues when matching
# We match SILVA taxonomy against this file (against scientific names) first

# Edit names.dmp file into a file only containing non-scientific names
sed "s/ <[a-zA-Z -,.&:'0-9]*>//g" names.dmp | grep -v 'scientific name' \
| cut -f 1,3 | awk -F $'\t' ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' \
| grep -v 'environmental\|uncultured\|unidentified\|metagenome' \
> NCBI_staxids_non_scientific.txt
	# 1. Remove <genus>, <family>, and other strings in <> brackets from taxonomy
  #    (otherwise the NCBI taxonomy won't match the SILVA taxonomy)
	# 2. Extract only non-scientific names and staxids
	# 3. Cut out columns we need
	# 4. Invert columns for script to work
	# 5. Removes lines containing "environmental", "uncultured", "unidentified",
  #    and "metagenome"
		# Needed because
			# SILVA taxonomy can have the same taxonomic ranks (e.g., "environmental
      # sample") for different higher ranks (e.g., "nematode; environmental
      # sample" and "bacteria;environmental sample"), wich would, however, be
      # assigned to the same staxid because the lower rank "environmental sample"
      # is similar
			# NCBI taxonomy can different staxids for the same taxonomic name, which
      # will cause issue when matching
# If SILVA taxonomy is not in exact scientific names, then we match against these
# to check for synonyms etc.

python - << EOF # switch to python

print("Processing files in python...")

import csv,sys

# Set variable names
NCBI_scientific_input="NCBI_staxids_scientific.txt"
NCBI_non_scientific_input="NCBI_staxids_non_scientific.txt"
SILVA_input="taxmap_slv_ssu_lsu_ref_nr_138_edited_for_NCBI_staxid_script.txt"

# Read in NCBI scientific names file as dictionary
NCBI_scientific = open(NCBI_scientific_input,'r')
# Set all file content to lowercase, so that matching of the files later is not
# depending on upper- or lowercase:
NCBI_scientific_lower = (line.lower() for line in NCBI_scientific)
reader1=csv.reader(NCBI_scientific_lower, delimiter='\t')
NCBI_scientific_dict={}
for row in reader1:
	NCBI_scientific_dict[row[0]]=row[1]

# Read in NCBI non-scientific names file as dictionary
NCBI_non_scientific = open(NCBI_non_scientific_input,'r')
# Set all file content to lowercase, so that matching of the files later is not
# depending on upper- or lowercase:
NCBI_non_scientific_lower = (line.lower() for line in NCBI_non_scientific)
reader2=csv.reader(NCBI_non_scientific_lower, delimiter='\t')
NCBI_non_scientific_dict={}
for row in reader2:
	NCBI_non_scientific_dict[row[0]]=row[1]

# Read in SILVA file as dictionary
SILVA = open(SILVA_input,'r')
# Set all file content to lowercase, so that matching of the files later is not
# depending on upper- or lowercase:
SILVA_lower = (line.lower() for line in SILVA)
reader3=csv.reader(SILVA_lower, delimiter='\t')
SILVA_dict={}
for row in reader3:
	SILVA_dict[row[0]]=row[1]

# Split up SILVA taxonomy and invert it
for key,value in SILVA_dict.items():
	SILVA_dict[key] = value.split(";")[::-1]

# Loop SILVA taxonomy dictionary over NCBI taxonomy dictionary for each line and
# match with NCBI staxid when a hit is found
output_dict={} # Make empty dictionary for matching lines
for key,value in SILVA_dict.items():
	l=0
	while l < len(value):
		if value[l] in NCBI_scientific_dict:
			output_dict[key.upper()]=NCBI_scientific_dict[value[l]]
			break
		elif value[l] in NCBI_non_scientific_dict:
			output_dict[key.upper()]=NCBI_non_scientific_dict[value[l]]
			break
		else:
			output_dict[key.upper()]='0'
			l += 1

# Saves the merged dictionary in a defined output file, delimited by tab
with open("SILVA_accession_numbers_and_NCBI_taxids.txt", 'w') as f:
    for key in output_dict.keys():
        f.write("%s\t%s\n"%(key,output_dict[key]))

EOF
# Switched back to bash


# As of 04 Sep 2020, the available SILVA LSU and SSU version contain duplicate
# sequences. I contacted the SILVA support, who said that should not have
# happened and that we can remove these duplicates. So we're going to remove
# them from the SSU file with a chunk of code that I got from the support
# (filter-fasta.awk is a separate subscript):
RE="^>(KY7649(2[1278]|58)|LNRQ01000003)\\\\."
cat SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta \
| filter-fasta.awk -v expression=$RE -v printMatch=0 \
> SILVA_138.1_SSURef_NR99_tax_silva_trunc_filtered.fasta

# Prepare SSU and LSU SILVA DB
cat SILVA_138.1_SSURef_NR99_tax_silva_trunc_filtered.fasta \
SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta \
> SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta


# Make the blast DB
#makeblastdb -dbtype 'nucl' -in SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta \
#-parse_seqids -taxid_map SILVA_accession_numbers_and_NCBI_taxids.txt

# Remove intermediate files
rm citations.dmp delnodes.dmp readme.txt SILVA_138.1_SSURef_NR99_tax_silva_trunc* \
SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta division.dmp gc.prt gencode.dmp \
SILVA_accession_numbers_and_NCBI_taxids.txt merged.dmp tax* names.dmp \
NCBI_staxids_* nodes.dmp
