#!/bin/bash

# Assigns taxonomy to blast searches including NCBI staxids, using the program
# ete3 and a script from David Ryder (CEFAS, Weymouth, England)
# (originally LookupTaxonoDetails.py)

# Needs ete3 being installed, folder .etetoolkit/taxa.sqlite being present in
# the home directory, and script LookupTaxonDetails3.py being in path

# If folder .etetoolkit/taxa.sqlite is not present in home directory, you need
# to run ONLY the ete3 command once first, delete the output, and then set the
# third parameter to the .etetoolit/taxa.sqlite folder that has now been
# generated in your folder

# Note that if ete3 is installed with conda, you have to activate the environment
# when running ete3


usage="$(basename "$0") -b <blast_file> -c <n> -e <path/to/.etetoolkit/taxa.sqlite>

Usage:
	-b  Blast file including column with NCBI staxids
	-c  Number of the column containing NCBI staxids
	-e  Path to ~/.etetoolkit/taxa.sqlite
	-h  Display this help and exit"


# Set specified options
while getopts ':b:c:e:h' opt; do
  case "${opt}" in
  	b) blast_file="${OPTARG}" ;;
    c) column="${OPTARG}" ;;
    e) etetoolkit="${OPTARG}" ;;
	h) echo "$usage"
	   exit ;;
	:) printf "Option -$OPTARG requires an argument."
	   echo -e "\n$usage"
	   exit ;;
   \?) printf "Invalid option: -$OPTARG"
	   echo -e "\n$usage"
	   exit
  esac
done
shift $((OPTIND - 1))


# Check if required options are set
if [[ -z "$blast_file" || -z "$column" || -z "$etetoolkit" ]]
then
   echo -e "-b, -c, and -e must be set.\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script.\n"
   exit
fi


# Making name for output variable:
blast_file_out_tmp=${blast_file%.*}_with_taxonomy.txt
blast_file_out=$(echo ${blast_file_out_tmp##*/})
# The ete3 command doesn't work on large files so we're splitting it up into
# chunks of 10000 lines:
touch matching_lineages.tsv
cut -f $column $blast_file > tmp
split -l 100000 --numeric-suffixes tmp tmp_chunk_
for id in tmp_chunk_*; do
  ete3 ncbiquery --info --search $(cat $id) >> matching_lineages.tsv
done
sed -i '1!{/^#/d;}' matching_lineages.tsv
# Running subscript:
LookupTaxonDetails3.py -b $blast_file -l matching_lineages.tsv \
-o $blast_file_out -t $column -e $etetoolkit
echo 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid lowest_rank lowest_hit superkingdom kingdom phylum subphylum class subclass order suborder infraorder family genus' \
| sed -e 's/ /\t/g' | cat - $blast_file_out > temp2 && mv temp2 $blast_file_out

rm matching_lineages.tsv tmp*
