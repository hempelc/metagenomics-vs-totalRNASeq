#!/bin/bash

# Version 0.2, made on 13 Jul 2020 by Chris Hempel (hempelc@uoguelph.ca)

# Version change: For option -t soft, if several hits have the same best bitscore
# for a sequence, they're kept and an LCA approach is applied to them

# Script to BLAST .fasta files using blastn, add taxonomy to the hits,
# and filter the hits so that each sequence gets assigned to one taxonomy

# Need to have script LookupTaxonDetails3.py in your PATH, and ete3 and blastn
# installed

# In order for this script to work - you MUST have .etetoolkit/taxa.sqlite in
# your HOME directory. If folder is not present in home directory, you need
# to run ONLY the ete3 command once first, delete the output, and then set the
# third parameter to the .etetoolit/taxa.sqlite folder that has now been
# generated in your folder

# Setting variables and importing modules
SCRIPT=$(realpath "${0}")
SCRIPTPATH=$(dirname "${SCRIPT}")
# shellcheck source=/dev/null
source "${SCRIPTPATH}"/utils.sh
cmd="${SCRIPT} $*" # Make variable containing full used command to print command in logfile

# CHRIS: Sergio wanted to remove option -e but I think it should be kept here,
# resason: this script can be run independently by itself so the option to
# change -e should be kept
usage="$(basename "$0") -i <input.fa> -f <fasta|blast> -t <soft|strict> -e <PATH/TO/.etetoolkit/taxa.sqlite> [-d <DB> -b <bitscore> -p <percentage> -c <n n n n n n> -T <threads>]

Usage:
  -i       Input file.
  -f       Format of input file:
           fasta:
              If input is a fasta file, justblast is performed first, then
              results are filtered based on option -t. Requires option -d.
           blast:
              Input is already blast output, then justblast is skipped and input
              is just filtered based on option -t. Requires blast input to be in
              the following outformat: '6 qseqid sseqid pident length mismatch
              gapopen qstart qend sstart send evalue bitscore staxids'
  -t       Type of filtering:
           soft:
              Keeps the best hit (highest bitscore) for each sequence. If multiple
              hits have the same highest bitscore, an LCA approach is applied
              (assigns the taxonomy to each sequence based on all taxonomic ranks
              that are identical in the remaining hits of each sequence)
           strict:
              Performs 3 steps:
              (1) bitscore filtering - keeps all hits with a bitscore >= (-b) and
              within (-p) % of the best bitscore per sequence.
              (2) similarity cutoff - only keeps the taxonomy of hits up to a
              certain rank, depending on the hit's blast % identity and cutoff
              values given in (-c).
              (3) LCA approach - assigns the taxonomy to each sequence based on
              all taxonomic ranks that are identical in the remaining hits of
              each sequence.
  -e       Path to .etetoolkit/taxa.sqlite (usually in home directory)
  -d       Database to use for blast.
  -b       Bitscore threshold to perform bitscore filtering on (-t) strict
           (default=155).
  -p       Percentage threshold to perform bitscore filtering on (-t) strict
           (default=0.02).
  -c       Similarity cutoff per hit based on BLAST pident values. % identify
           cutoffs have to be specified in a list divided by only spaces, for
           the ranks species, genus, family, order, class, phylum. Taxonomy is
           only kept for a rank of a hit if the hit's % identity is >= the
           respective cutoff (default=99 97 95 90 85 80).
  -T       Number of threads to use (default=16)."

# Set default options
input=''
format=''
db=''
filtering=''
bitscore='155'
percentage='0.02'
cutoff='99 97 95 90 85 80'
threads='16'

# Set specified options
while getopts ':i:f:e:d:t:b:p:c:T:h' opt; do
 	case "${opt}" in
		i) input="${OPTARG}" ;;
    f) format="${OPTARG}" ;;
    e) etetoolkit="${OPTARG}" ;;
    d) db="${OPTARG}" ;;
		t) filtering="${OPTARG}" ;;
		b) bitscore="${OPTARG}" ;;
		p) percentage="${OPTARG}" ;;
    c) cutoff="${OPTARG}" ;;
    T) threads="${OPTARG}" ;;
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
if [[ -z "$input" || -z "$format" || -z "$filtering" || -z "$etetoolkit" ]]; then
   echo -e "-i, -f, -e, and -t must be set\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script\n"
   exit
fi

if [[ $filtering != 'soft' && $filtering != 'strict' ]]; then
  echo -e "Invalid option for -t, must be set to either 'soft' or 'strict'\n"
  echo -e "$usage\n\n"
  echo -e "Exiting script\n"
  exit
fi

# CHRIS: Sergio threw this control step out but I don't know why
if [[ $format != 'fasta' && $format != 'blast' ]]; then
  echo -e "Invalid option for -f, must be set to either 'fasta' or 'blast'\n"
  echo -e "$usage\n\n"
  echo -e "Exiting script\n"
  exit
fi

if [[ $format == 'fasta' && $db == '' ]]; then
  echo -e "Option -d must be set when using -f fasta.'\n"
  echo -e "$usage\n\n"
  echo -e "Exiting script\n"
  exit
fi

# Functions
# SERGIO: I think that the assign_taxonomy_to_NCBI_staxids.sh as subscript is
# an overkill. I will put it in a function here for now
tax2ncbi() {
  blast_file="${1}"
  column="${2}"
  etetoolkit="${3}"
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

  # CHRIS: not entire content of assign_taxonomy_to_NCBI_staxids.sh script inserted here, so I'll insert the rest

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
}

##################### Write time, options etc. to output ######################

# Make open bracket to later tell script to write everything that follows into a logfile
(

  # Define starting time of script for total runtime calculation
  start=$(date +%s)
  echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
  echo -e "=================================================================\n\n"

  # Output specified options
  echo_subsection OPTIONS

  echo "Input (-i) was defined as $input"
  echo "Database (-d) was defined as $db"
  echo "Type (-t) was set to $filtering"
  echo "Bitscore threshold (-b) is $bitscore"
  echo "Percentage threshold (-p) is $percentage"
  echo "Cutoff (-c) is $cutoff"
  echo "$threads threads were used"
  echo -e "Script started with full command: $cmd\n"

  mkdir blast_filtering_results/


if [[ $format == 'fasta' ]] ; then
  echo_subsection RUNNING JUSTBLAST AGAINST DB
  blastn -query $input -db $db -out blast_filtering_results/blast_output.txt \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
  -evalue 1e-05 -num_threads $threads
  assign_taxonomy_input="blast_filtering_results/blast_output.txt"
  echo_subsection JUSTBLAST DONE
else
  assign_taxonomy_input="${input}"
fi

if [[ $filtering == 'soft' ]] ; then
  # Keeping only the best hit of each sequence:
  echo_subsection ASSIGNING TAXONOMY
  # Using a subscript:
  tax2ncbi $assign_taxonomy_input 13 $etetoolkit
  ## I think Sergio stopped here
  mv ${assign_taxonomy_input%.txt}_with_taxonomy.txt blast_filtering_results/
  sed -i '1d' blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy.txt

  echo -e "\n======== KEEPING ONLY BEST HIT PER SEQUENCE ========\n"
  # Just keep hits with the same best bitscore for each sequence, in case multiple
  # hits have the same bitscore:
  touch blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_filter.txt
  sort -u -k1,1 blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy.txt \
  | cut -f 1 | while read hit; do
    best_bitscore=$(grep "$(printf "^$hit\t")" blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy.txt \
    | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
    echo "Processing sequence $hit with bitscore $best_bitscore"
    grep $hit blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy.txt \
    | awk -v x=$best_bitscore '($12 >= x)' \
    >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_filter.txt
  done

  # Now we run an LCA approach, but it technically only runs if a sequence has
  # several hits with the same best bit score:
  touch blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_filter_and_LCA_noheader.txt
  sort -u -k1,1 blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_filter.txt \
  | cut -f 1 | while read hit ; do
    taxonomy=''
    for i in {15..27} # over all taxonomic ranks
    do
      rank_tax=$(grep $hit blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_filter.txt \
      | cut -f $i | cut -f 1-2 -d " " | uniq) # Extract taxonomy and cut down to first two words (essentially just relevant for rank "species", if more than two, e.g. subspecies info, we just want genus and species)
      if [[ $(echo "${rank_tax}" | wc -l) == 1 ]] ; then # If all taxonomy hits are the same
        if [[ ${i} == 15 ]]; then # if we look at species
          if [[ $(echo "${rank_tax}" | wc -w) == 1 ]]; then
            taxonomy=$(echo "${taxonomy}---NA")
          elif [[ "${rank_tax:0:1}" =~ [A-Z] ]]; then # if first letter is capitalized (indicates format "Genus species")
            if [[ $(echo $rank_tax | cut -f 2 -d ' ' | cut -c 1) =~ [a-z] ]]; then
              taxonomy=$(echo "${taxonomy}---${rank_tax}") # if first letter is non-capitalized (indicates format "Genus species")
            else
              taxonomy=$(echo "${taxonomy}---NA")
            fi
          else # if first letter not capitalized (indicates stuff like "uncultured bacterium" etc.)
            taxonomy=$(echo "${taxonomy}---NA")
          fi
        else # if not species rank and all taxonomy hits are the same
          taxonomy=$(echo "${taxonomy}---${rank_tax}")
        fi
      else # if taxonomy hits are not all the same
        taxonomy=$(echo "${taxonomy}---NA")
      fi
    done
    echo "${hit}${taxonomy}" | sed 's/---/\t/g' \
    >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_filter_and_LCA_noheader.txt
  done

  # Adding header and rearranging columns:
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $2}' \
  blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_filter_and_LCA_noheader.txt \
  > tmp
  echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies" \
  > blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_best_hit.txt \
  && cat tmp \
  >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_best_hit.txt
  rm tmp

  # Sort files:
  mkdir blast_filtering_results/intermediate_files/
  mv blast_filtering_results/${assign_taxonomy_input%.txt}* blast_filtering_results/intermediate_files/
  mv blast_filtering_results/intermediate_files/${assign_taxonomy_input%.txt}_with_taxonomy_and_best_hit.txt \
  blast_filtering_results/
fi

if [[ $filtering == 'strict' ]] ; then
  # Filtering reads:
  echo -e "\n======== ASSIGNING TAXONOMY ========\n"
  # Using a subscript:
  assign_taxonomy_to_NCBI_staxids.sh -b $assign_taxonomy_input -c 13 \
  -e $etetoolkit
  mv ${assign_taxonomy_input%.txt}_with_taxonomy.txt blast_filtering_results/
  sed '1d' blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy.txt \
  > blast_filtering_results/${assign_taxonomy_input%.txt}_bitscore_filtered_with_taxonomy_noheader.txt \
  && mv blast_filtering_results/${assign_taxonomy_input%.txt}_bitscore_filtered_with_taxonomy_noheader.txt \
  blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy.txt

  echo -e "\n======== PERFORMING BITSCORE FILTER ========\n"
  # Remove all hits that are below alignment length 100 (based on BASTA) and set
  # bitscore (default 155 based on CREST):
  awk -v x=$bitscore '($4 >= 100 && $12 >= x)' \
  blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy.txt \
  > blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold.txt

  # Just keep hits within first set % of best bitscore for each sequence
  touch blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt
  sort -u -k1,1 blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold.txt \
  | cut -f 1 | while read hit; do
    best_bitscore=$(grep "$(printf "^$hit\t")" blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold.txt \
    | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
    echo "Processing sequence $hit with bitscore $best_bitscore"
    bitscore_threshold=$(awk "BEGIN {print $best_bitscore - $best_bitscore * $percentage}")
    grep $hit blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold.txt \
    | awk -v x=$bitscore_threshold '($12 >= x)' \
    >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt
  done


  echo -e "\n======== PERFORMING SIMILARITY CUTOFF ========\n"
  # Setting certain columns to "NA" if they're below set threshold
  set -- $cutoff
  awk -v c1=$1 '($3 >= c1)' blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  > blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c1=$1 -v c2=$2 '(($3 < c1) && ($3 >= c2))' blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}'\
  >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c2=$2 -v c3=$3 '(($3 < c2) && ($3 >= c3))' blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}'\
  >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c3=$3 -v c4=$4 '(($3 < c3) && ($3 >= c4))' blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' \
  >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c4=$4 -v c5=$5 '(($3 < c4) && ($3 >= c5))' blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' \
  >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c5=$5 -v c6=$6 '(($3 < c5) && ($3 >= c6))' blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' \
  >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c6=$6 '($3 < c6)' blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' \
  >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt

  echo -e "\n======== PERFORMING LCA APPROACH ========\n"
  # Only keeping the taxonomy to the LCA, other ranks are set to "NA"
  touch blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt
  sort -u -k1,1 blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt \
  | cut -f 1 | while read hit ; do
    taxonomy=''
    for i in {15..27}
    do
      rank_tax=$(grep $hit blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt \
      | cut -f $i | cut -f 1-2 -d " " | uniq) # Extract taxonomy and cut down to first two words (essentially just relevant for rank "species", if more than two, e.g. subspecies info, we just want genus and species)
      if [[ $(echo "${rank_tax}" | wc -l) == 1 ]] ; then # If all taxonomy hits are the same
        if [[ ${i} == 15 ]]; then # if we look at species
          if [[ "${rank_tax:0:1}" =~ [A-Z] ]]; then # if first letter is capitalized (indicates format "Genus species")
            taxonomy=$(echo "${taxonomy}---${rank_tax}")
          else # if first letter not capitalized (indicates stuff like "uncultured bacterium" etc.)
            taxonomy=$(echo "${taxonomy}---NA")
          fi
        else # if not species rank and all taxonomy hits are the same
          taxonomy=$(echo "${taxonomy}---${rank_tax}")
        fi
      else # if taxonomy hits are not all the same
        taxonomy=$(echo "${taxonomy}---NA")
      fi
    done
    echo "${hit}${taxonomy}" | sed 's/---/\t/g' \
    >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt
  done

  # Adding header and rearranging columns:
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $2}' \
  blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt \
  > tmp
  echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies" \
  > blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt \
  && cat tmp \
  >> blast_filtering_results/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt
  rm tmp

  # Sort files:
  mkdir blast_filtering_results/intermediate_files/
  mv blast_filtering_results/${assign_taxonomy_input%.txt}* blast_filtering_results/intermediate_files/
  mv blast_filtering_results/intermediate_files/${assign_taxonomy_input%.txt}_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt \
  blast_filtering_results/
fi

# If blast was run, remove indexed files
if [[ $format == 'fasta' ]] ; then
  rm $input.fai
fi

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"

# Create log
) 2>&1 | tee blast_filtering_log.txt

mv blast_filtering_log.txt blast_filtering_results/
