#!/bin/bash

# Script to set up SILVA SSU and LSU NR99 for kraken2

# Usage: SILVA_SSU_LSU_kraken2_preparation.sh path_to_new_DB

DB=$1

# Download files:
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz

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

# Edit both the lsu and ssu files the same way:
for i in lsu ssu; do
  # Remove header and extract accession ID and taxid:
  tail -n +2 taxmap_slv_${i}_ref_nr_138.1.txt | cut -f1,6  > ${i}_acc_id.txt
  # Turn fasta into tab delimited file and edit so that accession numbers can be merged
  if [[ $i == 'ssu' ]] ; then
    fasta_to_tab SILVA_138.1_${i^^}Ref_NR99_tax_silva_trunc_filtered.fasta \
    | sed 's/\./\t/' | sed 's/ /\t/' \
    > SILVA_138.1_${i^^}Ref_NR99_tax_silva_trunc_edited.tab
  else
    fasta_to_tab SILVA_138.1_${i^^}Ref_NR99_tax_silva_trunc.fasta \
    | sed 's/\./\t/' | sed 's/ /\t/' \
    > SILVA_138.1_${i^^}Ref_NR99_tax_silva_trunc_edited.tab
  fi
  # Merge accession numbers of fasta file and taxid file and add kraken2-specific info to each sequence name
  mergeFilesOnColumn.pl SILVA_138.1_${i^^}Ref_NR99_tax_silva_trunc_edited.tab ${i}_acc_id.txt 1 1 \
  | awk -F "\t" '{ print $1 "\t" $2 "---" $6 "\t" $3 "\t" $4}' \
  | sed 's/^/>/' | sed 's/\t/\./' | sed 's/---/|kraken:taxid|/' | sed 's/\t/ /' \
  | sed 's/\t/\n/' > SILVA_138.1_${i^^}Ref_NR99_tax_silva_trunc_kraken2.fasta
done

cat SILVA_138.1_*SURef_NR99_tax_silva_trunc_kraken2.fasta > SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_kraken2.fasta
#cat *su_acc_id.txt > acc_ids.txt

kraken2-build --download-taxonomy --db $DB --use-ftp
kraken2-build --add-to-library SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_kraken2.fasta --db $DB
kraken2-build --build --db $DB

rm *su_acc_id.txt SILVA_138.1_*SURef_NR99_tax_silva_trunc_edited.tab \
SILVA_138.1_*SURef_NR99_tax_silva_trunc.fasta taxmap_slv_*su_ref_nr_138.1.txt \
SILVA_138.1_SSURef_NR99_tax_silva_trunc_filtered.fasta
