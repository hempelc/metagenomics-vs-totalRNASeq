#!/usr/bin/python

# Usage: ./assign_NCBI_staxids_to_CREST.py NCBI_staxids_scientific.txt NCBI_staxids_non_scientific.txt CREST_file output_file_name
# Assigns NCBI staxids to CREST taxonomy


import csv,sys,re,pandas as pd

# Set variable names
NCBI_scientific_input=sys.argv[1]
NCBI_non_scientific_input=sys.argv[2]
CREST_input=sys.argv[3]
output_name=sys.argv[4]

# Read in NCBI scientific names file as dictionary
NCBI_scientific = open(NCBI_scientific_input,'r')
NCBI_scientific_lower = (line.lower() for line in NCBI_scientific) # Sets all file content to lowercase, so that matching of the files later is not depending on upper- or lowercase)
reader1=csv.reader(NCBI_scientific_lower, delimiter='\t')
NCBI_scientific_dict={}
for row in reader1:
	NCBI_scientific_dict[row[0]]=row[1]

# Read in NCBI non-scientific names file as dictionary
NCBI_non_scientific = open(NCBI_non_scientific_input,'r')
NCBI_non_scientific_lower = (line.lower() for line in NCBI_non_scientific) # Sets all file content to lowercase, so that matching of the files later is not depending on upper- or lowercase)
reader2=csv.reader(NCBI_non_scientific_lower, delimiter='\t')
NCBI_non_scientific_dict={}
for row in reader2:
	NCBI_non_scientific_dict[row[0]]=row[1]

# Read in CREST file as dictionary
CREST_modified = open(CREST_input,'r')
reader3=csv.reader(CREST_modified, delimiter='\t')
CREST_dict={}
for row in reader3:
	CREST_dict[row[0]]=row[2]

CREST_lower_dict=dict((k, v.lower()) for k,v in CREST_dict.items())  # Converts CREST taxonomy to lowercases for matching with NCBI names
CREST_lower_no_mt_dict = {key: re.sub("\ \(mitochondrion\)", "", value) for key,value in CREST_lower_dict.items()}
CREST_lower_no_mt_chl_dict = {key: re.sub("\ \(chloroplast\)", "", value) for key,value in CREST_lower_no_mt_dict.items()}

# Split up CREST taxonomy and invert it
CREST_lower_no_mt_chl_dict_split_inverted={}
for key,value in CREST_lower_no_mt_chl_dict.items():
	CREST_lower_no_mt_chl_dict_split_inverted[key] = value.split(";")[::-1]

# Loop CREST taxonomy dictionary over NCBI taxonomy dictionary for each line and match with NCBI staxid when a hit is found
output_dict={} # Make empty dictionary for matching lines
exceptions=["environmental", "uncultured", "unidentified", "metagenome"] # when exception are found in SILVA taxonomy, the respective rank is skipped
for key,value in CREST_lower_no_mt_chl_dict_split_inverted.items():
	l=0
	while l < len(value):
                for exception in value[l].split(" "):
                        if exception in exceptions:
                                l += 1
                                continue
                if value[l] in NCBI_scientific_dict:
                        output_dict[key]=NCBI_scientific_dict[value[l]]
                        break
                elif value[l] in NCBI_non_scientific_dict:
                        output_dict[key]=NCBI_non_scientific_dict[value[l]]
                        break
                else:
                        output_dict[key]='0'
                        l += 1

# Convert dictionaries into pandas dataframes and merge them on the CREST OTU columns
CREST_df = pd.DataFrame(list(CREST_dict.items()), columns=['OTU','classification']).iloc[0:]
output_df = pd.DataFrame(list(output_dict.items()), columns=['OTU','NCBI_staxid']).iloc[0:]
output_merged = pd.merge(CREST_df, output_df, on='OTU')

# Write merged dataframe to output
output_merged.to_csv(output_name, sep='\t', na_rep='NA', index=False)
