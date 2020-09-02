#!/bin/bash

# Script to set up SILVA SSU and LSU NR99 for kraken2

# Download files:
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz

# Unzip them:
gunzip *

# Edit both the lsu and ssu files the same way:
for i in lsu ssu; do
  # Remove header and extract accession ID and taxid:
  tail -n +2 taxmap_slv_${i}_ref_nr_138.1.txt | cut -f1,6  > ${i}_acc_id.txt
  # Turn fasta into tab delimited file and edit so that accession numbers can be merged
  fasta_to_tab SILVA_138.1_${i^^}Ref_NR99_tax_silva.fasta | sed 's/\./\t/' \
  | sed 's/ /\t/' > SILVA_138.1_${i^^}Ref_NR99_tax_silva_edited.tab
  # Merge accession numbers of fasta file and taxid file and add kraken2-specific info to each sequence name
  mergeFilesOnColumn.pl SILVA_138.1_${i^^}Ref_NR99_tax_silva_edited.tab ${i}_acc_id.txt 1 1 \
  | awk -F "\t" '{ print $1 "\t" $2 "---" $6 "\t" $3 "\t" $4}' \
  | sed 's/^/>/' | sed 's/\t/\./' | sed 's/---/|kraken:taxid|/' | sed 's/\t/ /' \
  | sed 's/\t/\n/' > SILVA_138.1_${i^^}Ref_NR99_tax_silva_kraken2.fasta
done

# Remove all intermediate files but the final kraken2 files
rm $(ls | grep -v "kraken2")

# Cat the two resulting files together and use them to build a kraken2 DB.
# Download the two taxmap files and cat them together to extract accession_number
# and ID information which is needed to translate the kraken2 results into NCBI
# taxids. A copy of that is saved on our server in the generated kraken2 DB,
# /hdd1/databases/kraken2_SILVA_SSU_LSU/taxmap_slv_lsu_ssu_ref_nr_138.1.txt
