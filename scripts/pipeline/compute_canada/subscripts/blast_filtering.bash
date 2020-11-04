#!/bin/bash

# Version 0.2, made on 13 Jul 2020 by Chris Hempel (hempelc@uoguelph.ca)

# Version change: For option -t soft, if several hits have the same best bitscore
# for a sequence, they're kept and an LCA approach is applied to them

# Script to BLAST .fasta files using Sergio Hleap's justblast python module,
# add taxonomy to the hits, and filter the hits so that each sequence gets
# assigned to one taxonomy

# Need to have scripts assign_taxonomy_NCBI_staxids.sh and LookupTaxonDetails3.py
# in your PATH, and ete3 and justblast installed (https://pypi.org/project/justblast/)

# In order for "assign_taxonomy_NCBI_staxids.sh" to work - you MUST have
# .etetoolkit/taxa.sqlite in your HOME directory - check the ete3 toolkit
# to see how that's set up (http://etetoolkit.org/)

cmd="$0 $@" # Make variable containing full used command to print command in logfile
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

##################### Write time, options etc. to output ######################

# Make open bracket to later tell script to write everything that follows into a logfile
(

# Define starting time of script for total runtime calculation
start=$(date +%s)
echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
echo -e "=================================================================\n\n"


# Output specified options
echo -e "======== OPTIONS ========\n"

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
  echo -e "\n======== RUNNING JUSTBLAST AGAINST DB ========\n"
  justblast $input $db --cpus $threads --evalue 1e-05 --outfmt "6 qseqid \
  sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  staxids" --out_filename blast_filtering_results/blast_output.txt
  rm -r dask-worker-space/
  assign_taxonomy_input="blast_filtering_results/blast_output.txt"
  echo -e "\n======== JUSTBLAST DONE ========\n"
else
  assign_taxonomy_input=$input
fi

if [[ $filtering == 'soft' ]] ; then
  # Keeping only the best hit of each sequence:
  echo -e "\n======== ASSIGNING TAXONOMY ========\n"
  # Using a subscript:
  assign_taxonomy_NCBI_staxids.sh -b $assign_taxonomy_input -c 13 \
  -e $etetoolkit
  mv blast_output_with_taxonomy.txt blast_filtering_results/
  sed -i '1d' blast_filtering_results/blast_output_with_taxonomy.txt

  echo -e "\n======== KEEPING ONLY BEST HIT PER SEQUENCE ========\n"
  # Just keep hits with the same bitscore for each sequence, in case multiple
  # hits have the same bitscore:
  touch blast_filtering_results/blast_output_with_taxonomy_and_bitscore_filter.txt
  sort -u -k1,1 blast_filtering_results/blast_output_with_taxonomy.txt \
  | cut -f 1 | while read hit; do
    best_bitscore=$(grep "$(printf "^$hit\t")" blast_filtering_results/blast_output_with_taxonomy.txt \
    | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
    echo "Processing sequence $hit with bitscore $best_bitscore"
    grep $hit blast_filtering_results/blast_output_with_taxonomy.txt \
    | awk -v x=$best_bitscore '($12 >= x)' \
    >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_filter.txt
  done

  # Now we run an LCA approach, but it technically only runs if a sequence has
  # several hits with the same best bit score:
  touch blast_filtering_results/blast_output_with_taxonomy_and_bitscore_filter_and_LCA_noheader.txt
  sort -u -k1,1 blast_filtering_results/blast_output_with_taxonomy_and_bitscore_filter.txt \
  | cut -f 1 | while read hit ; do
    taxonomy_hit=$(echo $hit)
    taxonomy=''
    for i in {16..27}
    do
      rank_tax=$(grep $hit blast_filtering_results/blast_output_with_taxonomy_and_bitscore_filter.txt \
      | cut -f $i | uniq)
      if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
        taxonomy=$(echo "${taxonomy}---${rank_tax}")
      else
        taxonomy=$(echo "${taxonomy}---NA")
      fi
    done
    echo "${taxonomy_hit}${taxonomy}" | sed 's/---/\t/g' \
    >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_filter_and_LCA_noheader.txt
  done

  # Adding header and rearranging columns:
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $2}' \
  blast_filtering_results/blast_output_with_taxonomy_and_bitscore_filter_and_LCA_noheader.txt \
  > tmp
  echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit" \
  > blast_filtering_results/blast_output_with_taxonomy_and_best_hit.txt \
  && cat tmp \
  >> blast_filtering_results/blast_output_with_taxonomy_and_best_hit.txt
  rm blast_filtering_results/blast_output_with_taxonomy_and_bitscore_filter_and_LCA_noheader.txt \
  tmp

  # Sort files:
  mkdir blast_filtering_results/intermediate_files/
  mv blast_filtering_results/blast_output* blast_filtering_results/intermediate_files/
  mv blast_filtering_results/intermediate_files/blast_output_with_taxonomy_and_best_hit.txt \
  blast_filtering_results/
fi

if [[ $filtering == 'strict' ]] ; then
  # Filtering reads:
  echo -e "\n======== ASSIGNING TAXONOMY ========\n"
  # Using a subscript:
  assign_taxonomy_NCBI_staxids.sh -b $assign_taxonomy_input -c 13 \
  -e $etetoolkit
  mv blast_output_with_taxonomy.txt blast_filtering_results/
  sed '1d' blast_filtering_results/blast_output_with_taxonomy.txt \
  > blast_filtering_results/blast_output_bitscore_filtered_with_taxonomy_noheader.txt \
  && mv blast_filtering_results/blast_output_bitscore_filtered_with_taxonomy_noheader.txt \
  blast_filtering_results/blast_output_with_taxonomy.txt

  echo -e "\n======== PERFORMING BITSCORE FILTER ========\n"
  # Remove all hits that are below alignment length 100 (based on BASTA) and set
  # bitscore (default 155 based on CREST):
  awk -v x=$bitscore '($4 >= 100 && $12 >= x)' \
  blast_filtering_results/blast_output_with_taxonomy.txt \
  > blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold.txt

  # Just keep hits within first set % of best bitscore for each sequence
  touch blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt
  sort -u -k1,1 blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold.txt \
  | cut -f 1 | while read hit; do
    best_bitscore=$(grep "$(printf "^$hit\t")" blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold.txt \
    | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
    echo "Processing sequence $hit with bitscore $best_bitscore"
    bitscore_threshold=$(awk "BEGIN {print $best_bitscore - $best_bitscore * $percentage}")
    grep $hit blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold.txt \
    | awk -v x=$bitscore_threshold '($12 >= x)' \
    >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt
  done


  echo -e "\n======== PERFORMING SIMILARITY CUTOFF ========\n"
  # Setting certain columns to "NA" if they're below set threshold
  set -- $cutoff
  awk -v c1=$1 '($3 >= c1)' blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  > blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c1=$1 -v c2=$2 '(($3 < c1) && ($3 >= c2))' blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}'\
  >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c2=$2 -v c3=$3 '(($3 < c2) && ($3 >= c3))' blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}'\
  >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c3=$3 -v c4=$4 '(($3 < c3) && ($3 >= c4))' blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' \
  >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c4=$4 -v c5=$5 '(($3 < c4) && ($3 >= c5))' blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' \
  >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c5=$5 -v c6=$6 '(($3 < c5) && ($3 >= c6))' blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' \
  >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c6=$6 '($3 < c6)' blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' \
  >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt

  echo -e "\n======== PERFORMING LCA APPROACH ========\n"
  # Only keeping the taxonomy to the LCA, other ranks are set to "NA"
  touch blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt
  sort -u -k1,1 blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt \
  | cut -f 1 | while read hit ; do
    taxonomy_hit=$(echo $hit)
    taxonomy=''
    for i in {16..27}
    do
      rank_tax=$(grep $hit blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt \
      | cut -f $i | uniq)
      if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
        taxonomy=$(echo "${taxonomy}---${rank_tax}")
      else
        taxonomy=$(echo "${taxonomy}---NA")
      fi
    done
    echo "${taxonomy_hit}${taxonomy}" | sed 's/---/\t/g' \
    >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt
  done

  # Adding header and rearranging columns:
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $2}' \
  blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt \
  > tmp
  echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit" \
  > blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt \
  && cat tmp \
  >> blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt
  rm blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt \
  tmp

  # Sort files:
  mkdir blast_filtering_results/intermediate_files/
  mv blast_filtering_results/blast_output* blast_filtering_results/intermediate_files/
  mv blast_filtering_results/intermediate_files/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt \
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
