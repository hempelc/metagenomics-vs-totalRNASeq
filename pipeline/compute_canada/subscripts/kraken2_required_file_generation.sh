#!/bin/bash

# Script to generate 3 files that are needed for the taxonomic annotation of
# kraken2 output

# Make SILVA SSU LSU taxmap:

# As of 04 Sep 2020, the available SILVA LSU and SSU taxmap files contain
# duplicate accession IDs with different taxids in the SSU and LSU files. We're
# going to use all accession IDs from the SSU file, check which additional ones
# are in the LSU file (about 23,000 are not in the SSU file) and just take these
# extra ones from the LSU taxmap file to not overwrite SSU taxids with LSU taxids:
echo -e "\nDownloading SILVA taxmap files to generate SILVA_paths_and_taxids.txt:\n"
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz
gunzip taxmap_slv_*su_ref_nr_138.1.txt.gz
cat taxmap_slv_ssu_ref_nr_138.1.txt > taxmap_slv_ssu_lsu_ref_nr_138.1.txt
tail -n +2 taxmap_slv_ssu_ref_nr_138.1.txt | cut -f1 > grep_list.txt
tail -n +2 taxmap_slv_lsu_ref_nr_138.1.txt | grep -v -f grep_list.txt \
>> taxmap_slv_ssu_lsu_ref_nr_138.1.txt
tail -n +2 taxmap_slv_ssu_lsu_ref_nr_138.1.txt | cut -f 4,6 | sort -u \
> SILVA_paths_and_taxids.txt
# Kraken2 spits out the taxid 0 when no hit is found, but 0 doesn't
# exist in the SILVA taxonomy, so manually add taxid 0 with path
# “No hits” to the SILVA path file:
echo -e "No hits;\t0" > tmp && cat SILVA_paths_and_taxids.txt >> tmp \
&& mv tmp SILVA_paths_and_taxids.txt
rm taxmap_slv_*_ref_nr_138.1.txt grep_list.txt


# Make NCBI taxonomy files:
echo -e "\nDownloading taxdmp.zip to generate NCBI_staxids_scientific.txt and NCBI_staxids_scientific.txt:\n"
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
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
rm *.dmp readme.txt taxdmp.zip gc.prt
