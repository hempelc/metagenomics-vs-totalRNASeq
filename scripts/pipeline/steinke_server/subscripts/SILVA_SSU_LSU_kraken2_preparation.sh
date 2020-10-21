#!/bin/bash

# Script to set up SILVA SSU and LSU NR99 for kraken2

# Usage: SILVA_SSU_LSU_kraken2_preparation.sh path_to_new_DB path_and_name_for_new_generated_DB

# Set variable to path and name for new generated DB:
DB=$1

# Download files:
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz

# Unzip them:
gunzip *

# As of 04 Sep 2020, the available SILVA LSU and SSU version contain duplicate
# sequences. I contacted the SILVA support, who said that should not have
# happened and that we can remove these duplicates. So we're going to remove
# them from the SSU file with a chunk of code that I got from the support
# (filter-fasta.awk is a separate subscript):
RE="^>(KY7649(2[1278]|58)|LNRQ01000003)\\\\."
cat SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta \
| filter-fasta.awk -v expression=$RE -v printMatch=0 \
> SILVA_138.1_SSURef_NR99_tax_silva_trunc_filtered.fasta

#Concatenate the SSU und LSU fasta files:
cat SILVA_138.1_SSURef_NR99_tax_silva_trunc_filtered.fasta \
> SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta
cat SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta \
>> SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta

# As of 04 Sep 2020, the available SILVA LSU and SSU taxmap files contain
# duplicate accession IDs with different taxids in the SSU and LSU files. We're
# going to use all accession IDs from the SSU file, check which additional ones
# are in the LSU file (about 23,000 are not in the SSU file) and just take these
# extra ones from the LSU taxmap file to not overwrite SSU taxids with LSU taxids:
cat taxmap_slv_ssu_ref_nr_138.1.txt > taxmap_slv_ssu_lsu_ref_nr_138.1.txt
tail -n +2 taxmap_slv_ssu_ref_nr_138.1.txt | cut -f1 > grep_list.txt
tail -n +2 taxmap_slv_lsu_ref_nr_138.1.txt | grep -v -f grep_list.txt >> taxmap_slv_ssu_lsu_ref_nr_138.1.txt
tail -n +2 taxmap_slv_ssu_lsu_ref_nr_138.1.txt | cut -f1,6  > accids_taxids.txt

# Turn fasta into tab delimited file and edit so that accession numbers can be merged
fasta_to_tab SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta \
| sed 's/\./\t/' | sed 's/ /\t/' \
> SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_edited.tab

# Merge accession numbers of fasta file and taxid file and add kraken2-specific
# info to each sequence name (when adding fasta sequences to a kraken2 DB, they
# need to contain the taxid in the name):
mergeFilesOnColumn.pl SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_edited.tab accids_taxids.txt 1 1 \
| awk -F "\t" '{ print $1 "\t" $2 "---" $6 "\t" $3 "\t" $4}' \
| sed 's/^/>/' | sed 's/\t/\./' | sed 's/---/|kraken:taxid|/' | sed 's/\t/ /' \
| sed 's/\t/\n/' > SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_kraken2.fasta

# To make a kraken2 DB, we need files with taxonomical info in the fashion of
# NCBI (like names.dmp, nodes.dmp, and seqid2taxid from NCBI). We're going to do
# that in python:

python - << EOF # switch to python

import sys,re,pandas as pd, numpy as np, csv

print("Processing files in python...")

# Set variable names:
tax_slv_file="tax_slv_ssu_138.1.txt"


# Make nodes_SILVA_SSU_LSU.dmp:

## Read in tax file (note: we don't need the LSU tax file as well since all taxa
## from the LSU tax file are also part of the SSU tax file):
tax_slv_df1=pd.read_csv(tax_slv_file, sep='\t', header=None)

## Adding a row with root as first line, otherwise we can't make the nodes.dmp
## file as there is no parent ID for domains:
tax_slv_df1.loc[-1] = ['root;', '1', 'no rank', 'NaN', 'NaN']  # adding a row
tax_slv_df1.index = tax_slv_df1.index + 1  # shifting index
tax_slv_df1.sort_index(inplace=True)

## Make dictionary with the taxonomy paths and IDs
tax_slv_dic1={}
for index, row in tax_slv_df1.iterrows():
	tax_slv_dic1[row[0]]=row[1]

## Generate parent paths
tax_slv_df2=tax_slv_df1.copy() ### Make a copy of df1
### And remove the last rank in the paths to generate parent paths for parent
### IDs:
for index, row in tax_slv_df2.iterrows():
	row[0]=re.sub('[^;]*;$', '', row[0])

## Make list with parent paths
tax_slv_list=[]
for index, row in tax_slv_df2.iterrows():
	tax_slv_list.append(row[0])

## Make a parent ID list by matching parent paths with original IDs in dic1
parent_id_list=[]
for idx, item in enumerate(tax_slv_list):
	if item == '': # If parent path empty = no parent, insert ID 1, which stands for root
		parent_id_list.append('1')
	else:
		if tax_slv_list[idx] in tax_slv_dic1: # If parent is in the tax file, add its ID to the parent ID list
			parent_id_list.append(tax_slv_dic1[tax_slv_list[idx]])
		else: # Else: check the previous parent (had to be done because in once incident a parent was not in the tax file)
			parent_path_minus_one=re.sub('[^;]*;$', '', tax_slv_list[idx])
			parent_id_list.append(tax_slv_dic1[parent_path_minus_one])

## Columns have to be separated by "\t|\t" or "\t-\t", which doesn't work as delimiter, so
## we generate a spacer vector made of "|" and "-" and a spacer for the last
## column containing "scientific name" to fit the format:
spacer1=np.repeat("|", len(tax_slv_df1), axis=0)
spacer2=np.repeat("-", len(tax_slv_df1), axis=0)
spacer3=np.repeat("scientific name", len(tax_slv_df1), axis=0)

## Generate the final dataframe
nodes_SILVA_SSU_LSU_df=pd.DataFrame({'tax_id':tax_slv_df1[1], 'spacer1':spacer1, \
'parent_tax_id':parent_id_list, 'spacer2':spacer1, 'rank':tax_slv_df1[2], \
'spacer3':spacer1, 'spacer4':spacer2, 'spacer5':spacer1})

## Write df
nodes_SILVA_SSU_LSU_df.to_csv("nodes.dmp", sep='\t', na_rep='NA', index=False, \
header=False)


# Make names_SILVA_SSU_LSU.dmp

tax_slv_df3=tax_slv_df1.copy() # Make a copy of df1
# Remove the last ";"" and the first taxon in the taxonomy path:
tax_slv_df3[0]=tax_slv_df3[0].str.replace(r';$', '').str.replace(r'^.*;', '')

## Generate the final dataframe
names_SILVA_SSU_LSU_df=pd.DataFrame({'taxon':tax_slv_df3[1], 'spacer1':spacer1, \
'tax_id':tax_slv_df3[0], 'spacer2':spacer1, 'spacer3':spacer2, 'spacer4':spacer1, \
'spacer5':spacer3, 'spacer6':spacer1})

## Write df
names_SILVA_SSU_LSU_df.to_csv("names.dmp", sep='\t', na_rep='NA', index=False, \
header=False)

EOF
# Switched back to bash

# Make the kraken2 DB
mkdir $DB
mkdir $DB/taxonomy
mv nodes.dmp names.dmp $DB/taxonomy # Prepare the taxonomy for the kraken2 DB
kraken2-build --add-to-library SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_kraken2.fasta --db $DB
kraken2-build --build --db $DB

rm accids_taxids.txt SILVA_138.1_*SURef_NR99_tax_silva_trunc_edited.tab \
SILVA_138.1_*SURef_NR99_tax_silva_trunc.fasta taxmap_slv_*su_ref_nr_138.1.txt \
SILVA_138.1_SSURef_NR99_tax_silva_trunc_filtered.fasta tax_slv_ssu_138.1.txt \
SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_kraken2.fasta grep_list.txt
