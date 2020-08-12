#!/bin/bash

# Version 1.0, made on 13 Jul 2020 by Chris Hempel (hempelc@uoguelph.ca)

# Script to BLAST .fasta files using Sergio Hleap's justblast python module,
# add taxonomy to the hits, and filter the hits so that each sequence gets
# assigned to one taxonomy

# Need to have scripts assign_taxonomy_NCBI_staxids.sh and LookupTaxonDetails3.py
# in your PATH, and justblast installed (https://pypi.org/project/justblast/)
# In order for "assign_taxonomy_NCBI_staxids.sh" to work - you MUST have
# .etetoolkit/taxa.sqlite in your HOME directory - check the ete3 toolkit
# to see how that's set up (http://etetoolkit.org/)

cmd="$0 $@" # Make variable containing full used command to print command in logfile
usage="$(basename "$0") -s <sequences.fa> -d <DB> -t <soft|strict> [-b <bitscore> -p <percentage> -c <n n n n n n> -T <threads>]

Usage:
  -s       Sequences in fasta format that are to be blasted.
  -d       Database to use for blast.
  -t       Type:
           soft:
              Simply keeps the best hit (highest bitscore) for each sequence
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
sequences=''
db=''
type=''
bitscore='155'
percentage='0.02'
cutoff='99 97 95 90 85 80'
threads='16'

# Set specified options
while getopts ':s:d:t:b:p:c:T:h' opt; do
 	case "${opt}" in
		s) sequences="${OPTARG}" ;;
		d) db="${OPTARG}" ;;
		t) type="${OPTARG}" ;;
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
if [[ -z "$sequences" || -z "$db" || -z "$type" ]]
then
   echo -e "-f, -d, and -t must be set\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script\n"
   exit
fi

if [[ $type != 'soft' && $type != 'strict' ]]
then
  echo -e "Invalid option for -t, must be set to either 'soft' or 'strict'\n"
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

echo "Sequences (-s) were defined as $sequences"
echo "Database (-d) was defined as $db"
echo "Type (-t) was set to $type"
echo "Bitscore threshold (-b) is $bitscore"
echo "Percentage threshold (-p) is $percentage"
echo "Cutoff (-c) is $cutoff"
echo "$threads threads were used"
echo -e "Script started with full command: $cmd\n"


echo -e "\n======== RUNNING JUSTBLAST AGAINST DB ========\n"
mkdir blast_filtering/
justblast $(realpath $sequences) $db --cpus $threads --evalue 1e-05 --outfmt "6 qseqid \
sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
staxids" --out_filename blast_output.txt
mv blast_output.txt blast_filtering/
cd blast_filtering/
echo -e "\n======== JUSTBLAST DONE ========\n"

if [[ $type == 'soft' ]] ; then
  echo -e "\n======== ASSIGNING TAXONOMY ========\n"
  assign_taxonomy_NCBI_staxids.sh -b blast_output.txt -c 13 \
  -e ~/.etetoolkit/taxa.sqlite
  sed -i '1d' blast_output_with_taxonomy.txt

  echo -e "\n======== KEEPING ONLY BEST HIT PER SEQUENCE ========\n"
  sort -k1,1n -k12,12nr blast_output_with_taxonomy.txt \
  | sort -u -k1,1 > blast_output_with_taxonomy_sorted.txt

  # Edit file format
  cut -f1,16- blast_output_with_taxonomy_sorted.txt \
  > blast_output_with_taxonomy_sorted_cut.txt

  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $2}' \
  blast_output_with_taxonomy_sorted_cut.txt \
  > blast_output_with_taxonomy_sorted_cut_reorganized.txt


  echo -e "sequence\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies" \
  > blast_output_with_taxonomy_and_best_hit.txt \
  && cat blast_output_with_taxonomy_sorted_cut_reorganized.txt \
  >> blast_output_with_taxonomy_and_best_hit.txt


  # Sort files
  mkdir intermediate_files/
  mv blast_output.txt blast_output_with_taxonomy.txt \
  blast_output_with_taxonomy_sorted.txt blast_output_with_taxonomy_sorted_cut.txt \
  blast_output_with_taxonomy_sorted_cut_reorganized.txt intermediate_files/
fi

if [[ $type == 'strict' ]] ; then
  echo -e "\n======== ASSIGNING TAXONOMY ========\n"
  assign_taxonomy_NCBI_staxids.sh -b blast_output.txt -c 13 \
  -e ~/.etetoolkit/taxa.sqlite
  sed '1d' blast_output_with_taxonomy.txt \
  > blast_output_bitscore_filtered_with_taxonomy_noheader.txt \
  && mv blast_output_bitscore_filtered_with_taxonomy_noheader.txt \
  blast_output_with_taxonomy.txt

  echo -e "\n======== PERFORMING BITSCORE FILTER ========\n"
  # Remove all hits that are below alignment length 100 (based on BASTA) and set
  # bitscore (default 155 based on CREST)
  awk -v x=$bitscore '($4 >= 100 && $12 >= x)' \
  blast_output_with_taxonomy.txt > blast_output_with_taxonomy_and_bitscore_threshold.txt

  # Just keep hits within first set % of best bitscore for each contig
  touch blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt
  sort -u -k1,1 blast_output_with_taxonomy_and_bitscore_threshold.txt | cut -f 1 \
  | while read hit; do
    best_bitscore=$(grep ^$hit blast_output_with_taxonomy_and_bitscore_threshold.txt \
    | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
    echo "Processing sequence $hit with bitscore $best_bitscore"
    bitscore_threshold=$(awk "BEGIN {print $best_bitscore - $best_bitscore * $percentage}")
    grep $hit blast_output_with_taxonomy_and_bitscore_threshold.txt \
    | awk -v x=$bitscore_threshold '($12 >= x)' \
    >> blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt
  done


  echo -e "\n======== PERFORMING SIMILARITY CUTOFF ========\n"
  set -- $cutoff
  awk -v c1=$1 '($3 >= c1)' blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  > blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c1=$1 -v c2=$2 '(($3 < c1) && ($3 >= c2))' blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}'\
  >> blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c2=$2 -v c3=$3 '(($3 < c2) && ($3 >= c3))' blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}'\
  >> blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c3=$3 -v c4=$4 '(($3 < c3) && ($3 >= c4))' blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' \
  >> blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c4=$4 -v c5=$5 '(($3 < c4) && ($3 >= c5))' blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' \
  >> blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c5=$5 -v c6=$6 '(($3 < c5) && ($3 >= c6))' blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' \
  >> blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt
  awk -v c6=$6 '($3 < c6)' blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter.txt \
  | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' \
  >> blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt

  echo -e "\n======== PERFORMING LCA APPROACH ========\n"
  touch blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt
  sort -u -k1,1 blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt \
  | cut -f 1 | while read hit ; do
    taxonomy_hit=$(echo $hit)
    taxonomy=''
    for i in {16..27}
    do
      rank_tax=$(grep $hit blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff.txt \
      | cut -f $i | uniq)
      if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
        taxonomy=$(echo "${taxonomy}---${rank_tax}")
      else
        taxonomy=$(echo "${taxonomy}---NA")
      fi
    done
    echo "${taxonomy_hit}${taxonomy}" | sed 's/---/\t/g' \
    >> blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt
  done

  # Adding header and rearranging columns
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $2}' \
  blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt \
  > tmp
  echo -e "sequence\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit" \
  > blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt \
  && cat tmp \
  >> blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt
  rm blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA_noheader.txt \
  tmp

  # Sort files
  mkdir intermediate_files/
  mv blast_output* intermediate_files/
  mv intermediate_files/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt .
fi

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"


rm $sequences.fai

# Create log
) 2>&1 | tee blast_filtering_log.txt
