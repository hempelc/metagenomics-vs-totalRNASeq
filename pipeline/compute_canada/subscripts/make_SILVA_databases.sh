#!/bin/bash

# Script to set up a SILVA SSU and LSU NR99 BLASTand kraken2 DB

# Usage: db)make_test.sh path_to_new_kraken2_DB path_and_name_for_new_generated_DB

# Set variable to path name for generated SILVA kraken2 DB
# Note: kraken2 won't work if the DB is moved after building
kraken2DB=$1
blastDB=$2

# Get SILVA fasta and taxonomy files
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz
gunzip *.gz

# Get NCBI taxonomy files
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
mv names.dmp names_ncbi.dmp

# Edit SILVA taxonomy file:
# SILVA taxonomy file needs to have two columns, first one with accession number
# and second one with taxonomy path
sed -e 's/\t/./' -e 's/\t/./' taxmap_slv_ssu_ref_nr_138.1.txt > taxmap_slv_ssu_lsu_ref_nr_138.1.txt
tail -n +2 taxmap_slv_ssu_lsu_ref_nr_138.1.txt | cut -f1 | sed s'/\./\\\./g' > grep_list.txt
sed -e 's/\t/./' -e 's/\t/./' taxmap_slv_lsu_ref_nr_138.1.txt | tail -n +2 \
| grep -v -f grep_list.txt >> taxmap_slv_ssu_lsu_ref_nr_138.1.txt
sed -i "s/ <[a-zA-Z -,.&:'0-9]*>//" taxmap_slv_ssu_lsu_ref_nr_138.1.txt
# 1. concatenate ssu and lsu taxmap files (with access id and start and stop merged)
#    into one ssu_lsu file with one header
#	  	 note: as of 04 Sep 2020, the available SILVA LSU and SSU taxmap files contain
#   	 duplicate accession IDs with different taxids in the SSU and LSU files. We're
#	  	 going to use all accession IDs from the SSU file, check which additional ones
# 	  are in the LSU file (about 23,000 are not in the SSU file) and just take these
# 	  extra ones from the LSU taxmap file to not overwrite SSU taxids with LSU taxids:
# 2. removes <genus>, <family> etc. from taxonomic paths (otherwise the SILVA
#    taxonomy won't match the NCBI taxonomy)


# As of 04 Sep 2020, the available SILVA LSU and SSU version contain duplicate
# sequences. I contacted the SILVA support, who said that should not have
# happened and that we can remove these duplicates. So we're going to remove
# them from the SSU file with a chunk of code that I got from the support
# (filter-fasta.awk is a separate subscript):
RE="^>(KY7649(2[1278]|58)|LNRQ01000003)\\\\."
cat SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta \
| filter-fasta.awk -v expression=$RE -v printMatch=0 \
> SILVA_138.1_SSURef_NR99_tax_silva_trunc_filtered.fasta

# Prepare the SSU and LSU SILVA DB
cat SILVA_138.1_SSURef_NR99_tax_silva_trunc_filtered.fasta \
SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta \
> SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta

# Turn fasta into tab delimited file and edit so that accession numbers can be merged later
fasta_to_tab SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta \
| sed 's/ /\t/' > SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_edited.tab

# Edit names_ncbi.dmp file into a file only containing scientific names
sed "s/ <[a-zA-Z -,.&:'0-9]*>//g" names_ncbi.dmp | grep 'scientific name' \
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

# Edit names_ncbi.dmp file into a file only containing non-scientific names
sed "s/ <[a-zA-Z -,.&:'0-9]*>//g" names_ncbi.dmp | grep -v 'scientific name' \
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


# To make a kraken2 DB, we need files with taxonomical info in the fashion of
# NCBI (like names.dmp, nodes.dmp, and seqid2taxid from NCBI). We're going to do
# that in python:

python3 - << EOF # switch to python

print("Processing files in python...")

import sys,re,pandas as pd, numpy as np, csv

# Set variable names
NCBI_scientific_input="NCBI_staxids_scientific.txt"
NCBI_non_scientific_input="NCBI_staxids_non_scientific.txt"
taxmap_input="taxmap_slv_ssu_lsu_ref_nr_138.1.txt"

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

# ======  1. makeblast

# Read in SILVA file as df
taxmap_df=pd.read_csv(taxmap_input, delimiter = "\t")
# Extract second last rank of taxonomy paths and turn them into df col "second_last_rank"
second_last_rank_lst=taxmap_df["path"].str.split(";")
taxmap_df["second_last_rank"]=second_last_rank_lst.str[-2]
# Extract organism name of taxmap df, replace unwanted characters,
# keep only first two words (idealy genus+species) and turn them into df col "species"
name=taxmap_df["organism_name"].str.replace('[\-\[\]\(\)\\/]', '').str.split()
taxmap_df['species'] = name.str[0].str.cat(name.str[1], sep=" ", na_rep='' ).str.strip()
# Replace all species names by NA if their first word (ideally genus) doesn't
# match the second last rank (which should be genus ideally)
taxmap_df['species'] = np.where(taxmap_df['second_last_rank'] == taxmap_df['species'].str.split().str[0], taxmap_df['species'], "NA")
# Merge path and species
taxmap_df["path"] = taxmap_df["path"] + taxmap_df["species"]
# Add ; at end of path to keep format consistent
taxmap_df["path"]=taxmap_df["path"] + ";"
# Keep all rows that contain species (used later)
taxmap_df_sp=taxmap_df.loc[taxmap_df['species'] != "NA"]
# Drop unwanted cols
taxmap_df=taxmap_df.drop(["organism_name", "taxid", "second_last_rank", "species"], axis=1)
# Turn all path leters into lowercase
taxmap_df["path"] = taxmap_df["path"].str.lower()
# Turn df into dict
taxmap_dict={}
for i in range(len(taxmap_df)):
	taxmap_dict[taxmap_df["primaryAccession.start.stop"][i]]=taxmap_df["path"][i]
# Split up SILVA taxonomy and invert it
for key,value in taxmap_dict.items():
	taxmap_dict[key] = value.split(";")[::-1][1:]

# Loop taxmap taxonomy dictionary over NCBI taxonomy dictionary for each line and
# match with NCBI staxid when a hit is found
output_dict={} # Make empty dictionary for matching lines
for key,value in taxmap_dict.items():
	l=0
	while l < len(value):
		if value[l] in NCBI_scientific_dict:
			output_dict[key.upper()]=NCBI_scientific_dict[value[l]]
			hit=True
			break
		elif value[l] in NCBI_non_scientific_dict:
			output_dict[key.upper()]=NCBI_non_scientific_dict[value[l]]
			hit=True
			break
		else:
			l += 1
	if hit==False:
		output_dict[key.upper()]='0'

# Saves the merged dictionary in a defined output file, delimited by tab
with open("SILVA_accession_numbers_and_NCBI_taxids.txt", 'w') as f:
    for key in output_dict.keys():
        f.write("%s\t%s\n"%(key,output_dict[key]))


# ======= 2. kraken2

# Set variable names:
tax_slv_input="tax_slv_ssu_138.1.txt"

# Make nodes.dmp:

## Read in tax file (note: we don't need the LSU tax file as well since all taxa
## from the LSU tax file are also part of the SSU tax file):
tax_slv_df_child=pd.read_csv(tax_slv_input, sep='\t', names=['path','taxID','rank','remark','release'], usecols=['path','taxID','rank'])
## Add paths from species and set their rank to species
tax_slv_df_child_sp=tax_slv_df_child.append(pd.DataFrame(taxmap_df_sp["path"]))
tax_slv_df_child_sp['rank']=tax_slv_df_child_sp['rank'].fillna("species")
## Get NCBI ID for every child and parent path, this now includes all ranks
child_ncbi_id=[]
parent_ncbi_id=[]
for child_path in tax_slv_df_child_sp["path"]:
	parent_path=re.sub('[^;]*;$', '', child_path)
	child_lst=child_path.lower().split(";")[::-1][1:]
	parent_lst=parent_path.lower().split(";")[::-1][1:]
	l=0
	hit_child=False
	while l < len(child_lst):
		if child_lst[l] in NCBI_scientific_dict:
			child_id=NCBI_scientific_dict[child_lst[l]]
			hit_child=True
			break
		elif child_lst[l] in NCBI_non_scientific_dict:
			child_id=NCBI_non_scientific_dict[child_lst[l]]
			hit_child=True
			break
		else:
			l += 1
	if hit_child==False:
		child_id='1'
	m=0
	hit_parent=False
	while m < len(parent_lst):
		if parent_lst[m] in NCBI_scientific_dict:
			if NCBI_scientific_dict[parent_lst[m]]==child_id:
				m+=1
				continue
			parent_id=NCBI_scientific_dict[parent_lst[m]]
			hit_parent=True
			break
		elif parent_lst[m] in NCBI_non_scientific_dict:
			if NCBI_non_scientific_dict[parent_lst[m]]==child_id:
				m+=1
				continue
			parent_id=NCBI_non_scientific_dict[parent_lst[m]]
			hit_parent=True
			break
		else:
			m += 1
	if hit_parent==False:
		parent_id='1'
	child_ncbi_id.append(child_id)
	parent_ncbi_id.append(parent_id)

## Add NCBI IDs to df
tax_slv_df_child_sp["child_ncbi_id"]=child_ncbi_id
tax_slv_df_child_sp["parent_ncbi_id"]=parent_ncbi_id
## Generate column indicating # of ranks
tax_slv_df_child_sp["numranks"]=tax_slv_df_child_sp["path"].str.count(";")
## Add a row with root as first line, otherwise we can't make the nodes.dmp
## file as there is no parent ID for domains:
tax_slv_df_child_sp.loc[-1] = ['root;', '1', 'no rank', '1', '1', '1']
## Only keep one path per child_ncbi_id, therefore sort df by child_ncbi_id and
## then numranks and keep only first occurence of child_ncbi_id
## (= path with the lowest # of ranks per child_ncbi_id)
tax_slv_df_child_sp=tax_slv_df_child_sp.sort_values(['child_ncbi_id', 'numranks']).drop_duplicates(subset=['child_ncbi_id']).reset_index().drop(["index"], axis=1)
## Columns have to be separated by "\t|\t" or "\t-\t", which doesn't work as delimiter, so
## we generate a spacer vector made of "|" and "-" and a spacer for the last
## column containing "scientific name" to fit the format:
spacer1=np.repeat("|", len(tax_slv_df_child_sp), axis=0)
spacer2=np.repeat("-", len(tax_slv_df_child_sp), axis=0)
spacer3=np.repeat("scientific name", len(tax_slv_df_child_sp), axis=0)

## Generate the final dataframe
nodes_SILVA_SSU_LSU_df=pd.DataFrame({'child_tax_id':tax_slv_df_child_sp["child_ncbi_id"], 'spacer1':spacer1, \
'parent_tax_id':tax_slv_df_child_sp["parent_ncbi_id"], 'spacer2':spacer1, 'rank':tax_slv_df_child_sp["rank"], \
'spacer3':spacer1, 'spacer4':spacer2, 'spacer5':spacer1})

## Write df
nodes_SILVA_SSU_LSU_df.to_csv("nodes.dmp", sep='\t', na_rep='NA', index=False, \
header=False)


# Make names.dmp
## Get last rank of every path and remove the last ";"
last_rank=tax_slv_df_child_sp["path"].str.replace(r';$', '').str.replace(r'^.*;', '')

## Generate the final dataframe
names_SILVA_SSU_LSU_df=pd.DataFrame({'tax_id':tax_slv_df_child_sp["child_ncbi_id"], 'spacer1':spacer1, \
'taxon':last_rank, 'spacer2':spacer1, 'spacer3':spacer2, 'spacer4':spacer1, \
'spacer5':spacer3, 'spacer6':spacer1})

## Write df
names_SILVA_SSU_LSU_df.to_csv("names.dmp", sep='\t', na_rep='NA', index=False, \
header=False)

print("Python done!")

EOF
# Switched back to bash


# Merge accession numbers of fasta file and taxid file and add kraken2-specific
# info to each sequence name (when adding fasta sequences to a kraken2 DB, they
# need to contain the taxid in the name):
mergeFilesOnColumn.pl SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_edited.tab \
SILVA_accession_numbers_and_NCBI_taxids.txt 1 1 \
| sed 's/^/>/g' | awk -F "\t" '{ print $1 "|kraken:taxid|" $5 " " $2 "\n" $3}' \
> SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_kraken2.fasta


# Make the kraken2 DB
mkdir $kraken2DB
mkdir $kraken2DB/taxonomy
cp nodes.dmp names.dmp $kraken2DB/taxonomy # Prepare the taxonomy for the kraken2 DB
kraken2-build --add-to-library SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_kraken2.fasta --db $kraken2DB
kraken2-build --build --db $kraken2DB

mkdir $blastDB
cd $blastDB
# Make the blast DB
makeblastdb -dbtype 'nucl' -in ../SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta \
-parse_seqids -taxid_map ../SILVA_accession_numbers_and_NCBI_taxids.txt \
-out $(basename $blastDB)
cd ..
