#!/bin/bash

wget http://hulab.ucf.edu/research/projects/rRNAFilter/software/rRNAFilter.zip
unzip rRNAFilter.zip
cd rRNAFilter/
fq2fa ../$1 ../${1%.*}.fa
fq2fa ../$2 ../${2%.*}.fa
java -jar -Xmx7g rRNAFilter_commandline.jar -i ../${1%.*}.fa -r 0
java -jar -Xmx7g rRNAFilter_commandline.jar -i ../${2%.*}.fa -r 0
cd ..
rm -r rRNAFilter rRNAFilter.zip
fasta_to_tab ${1%.*}.fa_rRNA | cut -f 1 | cut -f1 -d " " > names.txt
fasta_to_tab ${2%.*}.fa_rRNA | cut -f 1 | cut -f1 -d " " >> names.txt
sort -u names.txt > names_sorted.txt
seqtk subseq ${1%.*}.fa names_sorted.txt > rRNAFilter_paired_R1.fa
seqtk subseq ${2%.*}.fa names_sorted.txt > rRNAFilter_paired_R2.fa
rm names_sorted.txt names.txt
