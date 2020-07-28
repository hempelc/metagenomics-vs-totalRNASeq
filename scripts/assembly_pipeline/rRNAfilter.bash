#!/bin/bash

wget http://hulab.ucf.edu/research/projects/rRNAFilter/software/rRNAFilter.zip
unzip rRNAFilter.zip
cd rRNAFilter/
java -jar -Xmx7g rRNAFilter_commandline.jar -i $1 -r 0
java -jar -Xmx7g rRNAFilter_commandline.jar -i $2 -r 0
cd ..
rm -r rRNAFilter
fasta_to_tab ${1}_rRNA | cut -f 1 | cut -f1 -d " " > names.txt
fasta_to_tab ${2}_rRNA | cut -f 1 | cut -f1 -d " " >> names.txt
sort -u names.txt > names_sorted.txt
seqtk subseq ${1}_rRNA names_sorted.txt > rRNAFilter_paired_R1.fa
seqtk subseq ${2}_rRNA names_sorted.txt > rRNAFilter_paired_R2.fa
