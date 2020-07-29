#!/bin/bash

# As for now, SortMeRNA is only carried out for RNA assemblers

# Need to have: assign_NCBI_staxids_to_CREST_v3.py, assign_taxonomy_NCBI_staxids.sh,
# deinterleave_fastq_reads.sh, LookupTaxonDetails3.py,
# merge_mapped_reads_and_contigs.py, mergeFilesOnColumn.pl and fasta_to_tab
# in your PATH, and justblast installed (https://pypi.org/project/justblast/)
# In order for "assign_taxonomy_NCBI_staxids.sh" to work - you MUST have
# .etetoolkit/taxa.sqlite in your HOME directory - check the ete3 toolkit
# to see how that's set up

# IMPORTANT: both kraken2 and CREST contain a script called "classify". If both
# are in your PATH, the command "classify" in the script will refer to the one
# that is LAST in your path. This needs to be CREST, because we will call
# CREST's "classify" script. So make sure CREST comes after kraken2 in your PATH.

cmd="$0 $@" # Make variable containing full used command to print command in logfile
usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq> [aDRSMmUrtTBCfsh]

Usage:
	-1 Forward reads trimmed - must state full path from root to the file
	-2 Reverse reads trimmed - must state full path from root to the file
	-a Main flag to indicate that all following flags should be used
	-D Flag to use completely pipeline only for DNA assemblers
	-R Flag to use completely pipeline only for RNA assemblers
	-S Flag to generate SPADES assembly output
	-M Flag to generate METASPADES assembly output
	-m Flag to generate MEGAHIT assembly output
	-U Flag to generate IDBA-UD assembly output
	-r Flag to generate RNASPADES assembly output
	-t Flag to generate IDBA-tran assembly output
	-T Flag to generate TRINITY assembly output
	-B Flag to use BWA and BOWTIE2 - mapping
	-C Flag to use BLAST and CREST - classification
	-f Flag to produce final output files
	-s Flag to include original assembly contig/scaffold sequences in the final output files
	-h Display this help and exit"

# Set default options
Forward_read_trimmed=''
Reverse_read_trimmed=''
SPADES='false'
METASPADES='false'
MEGAHIT='false'
IDBA_UD='false'
RNASPADES='false'
IDBA_TRAN='false'
TRINITY='false'
MAP='false'
CLASSIFICATION='false'
FINAL='false'
READS='false'

# Set specified options
while getopts ':1:2:aDRSMmUrtTBCfsh' opt; do
 	case "${opt}" in
		1) Forward_read_trimmed="${OPTARG}" ;;
		2) Reverse_read_trimmed="${OPTARG}" ;;
		a) SPADES='true'
			 METASPADES='true'
			 MEGAHIT='true'
			 IDBA_UD='true'
			 RNASPADES='true'
	 		 IDBA_TRAN='true'
	 		 TRINITY='true'
			 MAP='true'
			 CLASSIFICATION='true'
			 FINAL='true'
			 READS='true' ;;
		D) SPADES='true'
			 METASPADES='true'
			 MEGAHIT='true'
			 IDBA_UD='true'
			 MAP='true'
			 CLASSIFICATION='true'
			 FINAL='true'
			 READS='true' ;;
		R) RNASPADES='true'
		   IDBA_TRAN='true'
		   TRINITY='true'
		   MAP='true'
		   CLASSIFICATION='true'
		   FINAL='true'
		   READS='true' ;;
		S) SPADES='true' ;;
		M) METASPADES='true' ;;
		m) MEGAHIT='true' ;;
		U) IDBA_UD='true' ;;
		r) RNASPADES='true' ;;
		t) IDBA_TRAN='true' ;;
		T) TRINITY='true' ;;
		B) MAP='true' ;;
		C) CLASSIFICATION='true' ;;
		f) FINAL='true' ;;
		s) READS='true' ;;
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
if [[ -z "$Forward_read_trimmed" || -z "$Reverse_read_trimmed" ]]
then
   echo -e "-1 and -2 must be set.\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script.\n"
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

echo -e "Forward reads were defined as $Forward_read_trimmed.\n"
echo -e "Reverse reads were defined as $Reverse_read_trimmed.\n"
echo -e "Script started with full command: $cmd\n"

echo -e "======== START RUNNING SCRIPT ========\n"


######################### Beginning of DNA pipelines ########################
mkdir PIPELINE_ASSEMBLERS
cd PIPELINE_ASSEMBLERS

mkdir PIPELINE_DNA
cd PIPELINE_DNA/

echo -e "\n======== RUNNING SPADES ========\n"
if [[ $SPADES == 'true' ]] ; then
	mkdir SPADES
	spades.py -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -o SPADES
	echo -e "\n======== SPADES DONE ========\n"

	if [[ $MAP == 'true' ]] ; then

		echo -e "\n======== starting bwa index ========\n"

		bwa index -p bwa_index ./SPADES/scaffolds.fasta

		echo -e "\n======== bwa index complete. Starting bwa ========\n"

		bwa mem -t 10 bwa_index $Forward_read_trimmed $Reverse_read_trimmed > bwa_output.sam

		rm bwa_index*


		# Output file (.sam) - edit
		samtools view -F 4 bwa_output.sam > mapped_reads_bwa.sam
		samtools view -f 4 bwa_output.sam > unmapped_reads_bwa.sam
		cat mapped_reads_bwa.sam > mapped_reads_bwa.txt
		cat unmapped_reads_bwa.sam > unmapped_reads_bwa.txt
		cut -f3 mapped_reads_bwa.txt > mapped_column3_reads_bwa.txt
		cut -f3 unmapped_reads_bwa.txt > unmapped_column3_reads_bwa.txt
		sort mapped_column3_reads_bwa.txt | uniq -c > sorted_mapped_column3_reads_bwa.txt
		sort unmapped_column3_reads_bwa.txt | uniq -c > sorted_unmapped_column3_reads_bwa.txt
		column -t sorted_mapped_column3_reads_bwa.txt > aligned_mapped_bwa.txt
		column -t sorted_unmapped_column3_reads_bwa.txt > aligned_unmapped_bwa.txt
		sed 's/  */\t/g' aligned_mapped_bwa.txt > out_mapped_bwa.txt
		sed 's/  */\t/g' aligned_unmapped_bwa.txt > out_unmappped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo -e "\n======== bwa complete. Starting bowtie2 index ========\n"

		bowtie2-build -f ./SPADES/scaffolds.fasta bowtie_index

		echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"

		bowtie2 -q -x bowtie_index -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -S bowtie2_output.sam

		rm bowtie_index*


		# Output file (.sam) - edit
		samtools view -F 4 bowtie2_output.sam > mapped_reads_bowtie.sam
		samtools view -f 4 bowtie2_output.sam > unmapped_reads_bowtie.sam
		cat mapped_reads_bowtie.sam > mapped_reads_bowtie.txt
		cat unmapped_reads_bowtie.sam > unmapped_reads_bowtie.txt
		cut -f3 mapped_reads_bowtie.txt > mapped_column3_reads_bowtie.txt
		cut -f3 unmapped_reads_bowtie.txt > unmapped_column3_reads_bowtie.txt
		sort mapped_column3_reads_bowtie.txt | uniq -c > sorted_mapped_column3_reads_bowtie.txt
		sort unmapped_column3_reads_bowtie.txt | uniq -c > sorted_unmapped_column3_reads_bowtie.txt
		column -t sorted_mapped_column3_reads_bowtie.txt > aligned_mapped_bowtie.txt
		column -t sorted_unmapped_column3_reads_bowtie.txt > aligned_unmapped_bowtie.txt
		sed 's/  */\t/g' aligned_mapped_bowtie.txt > out_mapped_bowtie.txt
		sed 's/  */\t/g' aligned_unmapped_bowtie.txt > out_unmappped_bowtie.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo -e "\n======== bowtie2 complete ========\n"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ SPADES/

	else

		echo "\n======== no BWA and BOWTIE2 output generated ========\n"

	fi

	if [[ $CLASSIFICATION == 'true' ]] ; then

		# BLAST against nt

		echo -e "\n======== RUNNING JUSTBLAST AGAINST NT AND PROCESSING OUTPUT ========\n"
		justblast ./SPADES/scaffolds.fasta /hdd1/databases/nt_database_feb_2020_indexed/nt --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_nt

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_nt | sort -u -k1,1 > BLAST_output_nt_first_hit.txt

		## Blast option 2: apply bitscore threshold, best bitscore %, similarity cutoff, and an LCA approach
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_nt > BLAST_output_nt_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_nt_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_nt_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_nt_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_nt_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_nt_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST NT AND PROCESSING OUTPUT DONE ========\n"

		# BLAST against SILVA
		echo -e "\n======== RUNNING JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT ========\n"
		justblast ./SPADES/scaffolds.fasta /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_SILVA

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_SILVA | sort -u -k1,1 > BLAST_output_SILVA_first_hit.txt

		## Blast option 2: apply bitscore threshold and best bitscore %
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_SILVA > BLAST_output_SILVA_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_SILVA_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_SILVA_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_SILVA_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_SILVA_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_SILVA_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT DONE ========\n"

		mkdir SPADES/BLAST
		mv BLAST_output_nt* BLAST_output_SILVA* SPADES/BLAST/

		echo -e "\n======== RUNNING CREST ========\n"
		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./SPADES/scaffolds.fasta -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		mv BLAST_output.xml CREST
		mv CREST SPADES/

		mkdir SPADES/CLASSIFICATION
		mv SPADES/BLAST SPADES/CREST SPADES/CLASSIFICATION/

		echo -e "\n======== CREST DONE ========\n"

	else

		echo -e "\n======== NO BLAST OR CREST OUTPUT GENERATED ========\n"


	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then

			echo -e "\n======== ADD CONTIG SEQUENCES TO BOWTIE AND BWA OUTPUT AND MERGE THEM WITH BLAST+CREST ========\n"
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./SPADES/scaffolds.fasta > ./SPADES/fasta_to_tabbed.txt
			mergeFilesOnColumn.pl ./SPADES/BWA/merge_input_mapped_bwa.txt ./SPADES/fasta_to_tabbed.txt 2 1 > ./SPADES/merged_original_bwa.txt
			cut -f1,2,4 ./SPADES/merged_original_bwa.txt > ./SPADES/important_column_3_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./SPADES/important_column_3_bwa.txt > ./SPADES/final_order_bwa.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./SPADES/final_bwa_merge_ready.txt && cat ./SPADES/final_order_bwa.txt >> ./SPADES/final_bwa_merge_ready.txt

			rm ./SPADES/merged_original_bwa.txt ./SPADES/important_column_3_bwa.txt ./SPADES/final_order_bwa.txt

			mergeFilesOnColumn.pl ./SPADES/BOWTIE2/merge_input_mapped_bowtie2.txt ./SPADES/fasta_to_tabbed.txt 2 1 > ./SPADES/merged_original_bowtie2.txt
			cut -f1,2,4 ./SPADES/merged_original_bowtie2.txt > ./SPADES/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./SPADES/important_column_3_bowtie2.txt > ./SPADES/final_order_bowtie2.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./SPADES/final_bowtie2_merge_ready.txt && cat ./SPADES/final_order_bowtie2.txt >> ./SPADES/final_bowtie2_merge_ready.txt

			rm ./SPADES/fasta_to_tabbed.txt ./SPADES/merged_original_bowtie2.txt ./SPADES/important_column_3_bowtie2.txt ./SPADES/final_order_bowtie2.txt

			echo -e "\n======== contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"


		else
			echo -e "\n======== MERGE BOWIE+BWA OUTPUT WITH BLAST+CREST ========\n"
			# Prepares to merge with BLAST or CREST - does not include assembly sequence
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./SPADES/BWA/merge_input_mapped_bwa.txt > ./SPADES/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./SPADES/temp.txt > ./SPADES/temp2.txt
			sed '1d' ./SPADES/temp2.txt > ./SPADES/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./SPADES/final_bwa_merge_ready.txt && cat ./SPADES/temp3.txt >> ./SPADES/final_bwa_merge_ready.txt
 			rm ./SPADES/temp*.txt

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./SPADES/BOWTIE2/merge_input_mapped_bowtie2.txt > ./SPADES/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./SPADES/temp.txt > ./SPADES/temp2.txt
			sed '1d' ./SPADES/temp2.txt > ./SPADES/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./SPADES/final_bowtie2_merge_ready.txt && cat ./SPADES/temp3.txt >> ./SPADES/final_bowtie2_merge_ready.txt
			rm ./SPADES/temp*.txt

			echo -e "\n======== no contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"
		fi

		# Assign BLAST taxonomy - nt
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST NT OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./SPADES/CLASSIFICATION/BLAST/BLAST_output_nt_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_first_hit_with_taxonomy.txt ./SPADES/
		sed '1d' ./SPADES/BLAST_output_nt_first_hit_with_taxonomy.txt > ./SPADES/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt
		mv ./SPADES/BLAST_output_nt_first_hit_with_taxonomy.txt ./SPADES/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./SPADES/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./SPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./SPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt
		echo

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./SPADES/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_bitscore_filtered_with_taxonomy.txt ./SPADES/CLASSIFICATION/BLAST/
		sed '1d' ./SPADES/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered_with_taxonomy.txt > ./SPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./SPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt > ./SPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./SPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./SPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./SPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./SPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./SPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./SPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./SPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./SPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./SPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./SPADES/BLAST_output_nt_similarity_filtered.txt
		awk '($3 < 80)' ./SPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./SPADES/BLAST_output_nt_similarity_filtered.txt

		### Immitating BASTA
		touch ./SPADES/spades_blast_nt_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./SPADES/BLAST_output_nt_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./SPADES/BLAST_output_nt_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./SPADES/spades_blast_nt_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./SPADES/spades_blast_nt_LCA_taxonomy_merge.txt && cat ./SPADES/spades_blast_nt_LCA_taxonomy_merge_noheader.txt >> ./SPADES/spades_blast_nt_LCA_taxonomy_merge.txt


		# Assign BLAST taxonomy - SILVA
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST SILVA OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./SPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_first_hit_with_taxonomy.txt ./SPADES/
		sed '1d' ./SPADES/BLAST_output_SILVA_first_hit_with_taxonomy.txt > ./SPADES/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt
		mv ./SPADES/BLAST_output_SILVA_first_hit_with_taxonomy.txt ./SPADES/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./SPADES/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./SPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./SPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./SPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt ./SPADES/CLASSIFICATION/BLAST/
		sed '1d' ./SPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt > ./SPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./SPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt > ./SPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./SPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./SPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./SPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./SPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./SPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./SPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./SPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./SPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./SPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./SPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '($3 < 80)' ./SPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./SPADES/BLAST_output_SILVA_similarity_filtered.txt

		### Immitating BASTA
		touch ./SPADES/spades_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./SPADES/BLAST_output_SILVA_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./SPADES/BLAST_output_SILVA_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./SPADES/spades_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./SPADES/spades_blast_SILVA_LCA_taxonomy_merge.txt && cat ./SPADES/spades_blast_SILVA_LCA_taxonomy_merge_noheader.txt >> ./SPADES/spades_blast_SILVA_LCA_taxonomy_merge.txt

		echo -e "\n======== taxonomy has been assigned to BLAST files ========"
		echo -e "======== BLAST files are ready ========\n"


		# Prepare CREST file - remove 1st/3rd column, add a header
		echo -e "\n======== ADD NCBI TAXONOMY TO CREST OUTPUT ========\n"
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./SPADES/CLASSIFICATION/CREST/otus.csv ./SPADES/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./SPADES/CLASSIFICATION/CREST/CREST_output.txt > ./SPADES/CREST_seperated.txt
		cut -f2,4 ./SPADES/CREST_seperated.txt > ./SPADES/CREST_header.txt
		sed '1d' ./SPADES/CREST_header.txt > ./SPADES/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./SPADES/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./SPADES/
		sed '1d' ./SPADES/CREST_tax_ready_with_taxonomy.txt > ./SPADES/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./SPADES/CREST_merge.txt && cat ./SPADES/CREST_tax_ready_with_taxonomy_noheader.txt >> ./SPADES/CREST_merge.txt

		rm ./SPADES/CREST_seperated.txt ./SPADES/CREST_header.txt ./SPADES/CREST_tax_ready.txt ./SPADES/CREST_tax_ready_with_taxonomy.txt ./SPADES/CREST_tax_ready_with_taxonomy_noheader.txt

		echo -e "\n======== CREST file is ready ========\n"

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./SPADES/CREST_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./SPADES/CREST_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/BLAST_first_hit_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/BLAST_first_hit_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/BLAST_first_hit_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/BLAST_first_hit_nt_BWA.txt

		merge_mapped_reads_and_contigs.py ./SPADES/spades_blast_SILVA_LCA_taxonomy_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/BLAST_LCA_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./SPADES/spades_blast_nt_LCA_taxonomy_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/BLAST_LCA_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./SPADES/spades_blast_SILVA_LCA_taxonomy_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/BLAST_LCA_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./SPADES/spades_blast_nt_LCA_taxonomy_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/BLAST_LCA_nt_BWA.txt


		# Move all files - easy to find!!!
		mkdir ./SPADES/MERGE_FILES
		mv ./SPADES/*merge*.txt ./SPADES/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./SPADES/CREST_BWA.txt > ./SPADES/new.txt &&  sed 's/_/\t/3' ./SPADES/new.txt > ./SPADES/new2.txt
		sed 's/NODE_//g' ./SPADES/new2.txt > ./SPADES/new3.txt && sed 's/length_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/cov_//g' ./SPADES/new4.txt > ./SPADES/new5.txt
		sed '1d' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/CREST_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/CREST_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/CREST_BWA.txt

		sed 's/_/\t/2' ./SPADES/CREST_BOWTIE2.txt > ./SPADES/new.txt &&  sed 's/_/\t/3' ./SPADES/new.txt > ./SPADES/new2.txt
		sed 's/NODE_//g' ./SPADES/new2.txt > ./SPADES/new3.txt && sed 's/length_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/cov_//g' ./SPADES/new4.txt > ./SPADES/new5.txt
		sed '1d' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/CREST_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/CREST_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./SPADES/BLAST_first_hit_nt_BWA.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_nt_first_hit_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_nt_first_hit_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_first_hit_nt_BWA.txt

		sed '1d' ./SPADES/BLAST_first_hit_nt_BOWTIE2.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_nt_first_hit_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_nt_first_hit_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_first_hit_nt_BOWTIE2.txt

		sed '1d' ./SPADES/BLAST_LCA_nt_BWA.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_nt_LCA_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_nt_LCA_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_LCA_nt_BWA.txt

		sed '1d' ./SPADES/BLAST_LCA_nt_BOWTIE2.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_nt_LCA_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_nt_LCA_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_LCA_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./SPADES/BLAST_first_hit_SILVA_BWA.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_SILVA_first_hit_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_SILVA_first_hit_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_first_hit_SILVA_BWA.txt

		sed '1d' ./SPADES/BLAST_first_hit_SILVA_BOWTIE2.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_SILVA_first_hit_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_SILVA_first_hit_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_first_hit_SILVA_BOWTIE2.txt

		sed '1d' ./SPADES/BLAST_LCA_SILVA_BWA.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_SILVA_LCA_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_SILVA_LCA_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_LCA_SILVA_BWA.txt

		sed '1d' ./SPADES/BLAST_LCA_SILVA_BOWTIE2.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_SILVA_LCA_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_SILVA_LCA_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_LCA_SILVA_BOWTIE2.txt


		mkdir ./SPADES/FINAL_FILES_SPADES
		mv ./SPADES/*final.txt ./SPADES/FINAL_FILES_SPADES/
		mv ./SPADES/BLAST_output_nt* ./SPADES/BLAST_output_SILVA* SPADES/CLASSIFICATION/BLAST/

		echo -e "\n======== final SPADES output generated ========\n"

	else
		echo -e "\n======== no FINAL merge generated ========\n"

	fi

else
	echo -e "\n======== no SPADES output generated ========\n"

fi


echo -e "\n======== RUNNING METASPADES ========\n"
if [[ $METASPADES == 'true' ]] ; then
	mkdir METASPADES
	metaspades.py -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -o METASPADES
	echo -e "\n======== METASPADES DONE ========\n"

	if [[ $MAP == 'true' ]] ; then

		echo -e "\n======== starting bwa index ========\n"

		bwa index -p bwa_index ./METASPADES/scaffolds.fasta

		echo -e "\n======== bwa index complete. Starting bwa ========\n"

		bwa mem -t 10 bwa_index $Forward_read_trimmed $Reverse_read_trimmed > bwa_output.sam

		rm bwa_index*


		# Output file (.sam) - edit
		samtools view -F 4 bwa_output.sam > mapped_reads_bwa.sam
		samtools view -f 4 bwa_output.sam > unmapped_reads_bwa.sam
		cat mapped_reads_bwa.sam > mapped_reads_bwa.txt
		cat unmapped_reads_bwa.sam > unmapped_reads_bwa.txt
		cut -f3 mapped_reads_bwa.txt > mapped_column3_reads_bwa.txt
		cut -f3 unmapped_reads_bwa.txt > unmapped_column3_reads_bwa.txt
		sort mapped_column3_reads_bwa.txt | uniq -c > sorted_mapped_column3_reads_bwa.txt
		sort unmapped_column3_reads_bwa.txt | uniq -c > sorted_unmapped_column3_reads_bwa.txt
		column -t sorted_mapped_column3_reads_bwa.txt > aligned_mapped_bwa.txt
		column -t sorted_unmapped_column3_reads_bwa.txt > aligned_unmapped_bwa.txt
		sed 's/  */\t/g' aligned_mapped_bwa.txt > out_mapped_bwa.txt
		sed 's/  */\t/g' aligned_unmapped_bwa.txt > out_unmappped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo -e "\n======== bwa complete. Starting bowtie2 index ========\n"

		bowtie2-build -f ./METASPADES/scaffolds.fasta bowtie_index

		echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"

		bowtie2 -q -x bowtie_index -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -S bowtie2_output.sam

		rm bowtie_index*


		# Output file (.sam) - edit
		samtools view -F 4 bowtie2_output.sam > mapped_reads_bowtie.sam
		samtools view -f 4 bowtie2_output.sam > unmapped_reads_bowtie.sam
		cat mapped_reads_bowtie.sam > mapped_reads_bowtie.txt
		cat unmapped_reads_bowtie.sam > unmapped_reads_bowtie.txt
		cut -f3 mapped_reads_bowtie.txt > mapped_column3_reads_bowtie.txt
		cut -f3 unmapped_reads_bowtie.txt > unmapped_column3_reads_bowtie.txt
		sort mapped_column3_reads_bowtie.txt | uniq -c > sorted_mapped_column3_reads_bowtie.txt
		sort unmapped_column3_reads_bowtie.txt | uniq -c > sorted_unmapped_column3_reads_bowtie.txt
		column -t sorted_mapped_column3_reads_bowtie.txt > aligned_mapped_bowtie.txt
		column -t sorted_unmapped_column3_reads_bowtie.txt > aligned_unmapped_bowtie.txt
		sed 's/  */\t/g' aligned_mapped_bowtie.txt > out_mapped_bowtie.txt
		sed 's/  */\t/g' aligned_unmapped_bowtie.txt > out_unmappped_bowtie.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo -e "\n======== bowtie2 complete ========\n"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ METASPADES/

	else

		echo "\n======== no BWA and BOWTIE2 output generated ========\n"

	fi

	if [[ $CLASSIFICATION == 'true' ]] ; then

		# BLAST against nt

		echo -e "\n======== RUNNING JUSTBLAST AGAINST NT AND PROCESSING OUTPUT ========\n"
		justblast ./METASPADES/scaffolds.fasta /hdd1/databases/nt_database_feb_2020_indexed/nt --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_nt

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_nt | sort -u -k1,1 > BLAST_output_nt_first_hit.txt

		## Blast option 2: apply bitscore threshold, best bitscore %, similarity cutoff, and an LCA approach
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_nt > BLAST_output_nt_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_nt_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_nt_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_nt_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_nt_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_nt_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST NT AND PROCESSING OUTPUT DONE ========\n"

		# BLAST against SILVA
		echo -e "\n======== RUNNING JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT ========\n"
		justblast ./METASPADES/scaffolds.fasta /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_SILVA

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_SILVA | sort -u -k1,1 > BLAST_output_SILVA_first_hit.txt

		## Blast option 2: apply bitscore threshold and best bitscore %
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_SILVA > BLAST_output_SILVA_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_SILVA_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_SILVA_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_SILVA_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_SILVA_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_SILVA_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT DONE ========\n"

		mkdir METASPADES/BLAST
		mv BLAST_output_nt* BLAST_output_SILVA* METASPADES/BLAST/

		echo -e "\n======== RUNNING CREST ========\n"
		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./METASPADES/scaffolds.fasta -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		mv BLAST_output.xml CREST
		mv CREST METASPADES/

		mkdir METASPADES/CLASSIFICATION
		mv METASPADES/BLAST METASPADES/CREST METASPADES/CLASSIFICATION/

		echo -e "\n======== CREST DONE ========\n"

	else

		echo -e "\n======== NO BLAST OR CREST OUTPUT GENERATED ========\n"


	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then

			echo -e "\n======== ADD CONTIG SEQUENCES TO BOWTIE AND BWA OUTPUT AND MERGE THEM WITH BLAST+CREST ========\n"
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./METASPADES/scaffolds.fasta > ./METASPADES/fasta_to_tabbed.txt
			mergeFilesOnColumn.pl ./METASPADES/BWA/merge_input_mapped_bwa.txt ./METASPADES/fasta_to_tabbed.txt 2 1 > ./METASPADES/merged_original_bwa.txt
			cut -f1,2,4 ./METASPADES/merged_original_bwa.txt > ./METASPADES/important_column_3_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./METASPADES/important_column_3_bwa.txt > ./METASPADES/final_order_bwa.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./METASPADES/final_bwa_merge_ready.txt && cat ./METASPADES/final_order_bwa.txt >> ./METASPADES/final_bwa_merge_ready.txt

			rm ./METASPADES/merged_original_bwa.txt ./METASPADES/important_column_3_bwa.txt ./METASPADES/final_order_bwa.txt

			mergeFilesOnColumn.pl ./METASPADES/BOWTIE2/merge_input_mapped_bowtie2.txt ./METASPADES/fasta_to_tabbed.txt 2 1 > ./METASPADES/merged_original_bowtie2.txt
			cut -f1,2,4 ./METASPADES/merged_original_bowtie2.txt > ./METASPADES/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./METASPADES/important_column_3_bowtie2.txt > ./METASPADES/final_order_bowtie2.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./METASPADES/final_bowtie2_merge_ready.txt && cat ./METASPADES/final_order_bowtie2.txt >> ./METASPADES/final_bowtie2_merge_ready.txt

			rm ./METASPADES/fasta_to_tabbed.txt ./METASPADES/merged_original_bowtie2.txt ./METASPADES/important_column_3_bowtie2.txt ./METASPADES/final_order_bowtie2.txt

			echo -e "\n======== contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"


		else
			echo -e "\n======== MERGE BOWIE+BWA OUTPUT WITH BLAST+CREST ========\n"
			# Prepares to merge with BLAST or CREST - does not include assembly sequence
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./METASPADES/BWA/merge_input_mapped_bwa.txt > ./METASPADES/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./METASPADES/temp.txt > ./METASPADES/temp2.txt
			sed '1d' ./METASPADES/temp2.txt > ./METASPADES/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./METASPADES/final_bwa_merge_ready.txt && cat ./METASPADES/temp3.txt >> ./METASPADES/final_bwa_merge_ready.txt
 			rm ./METASPADES/temp*.txt

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./METASPADES/BOWTIE2/merge_input_mapped_bowtie2.txt > ./METASPADES/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./METASPADES/temp.txt > ./METASPADES/temp2.txt
			sed '1d' ./METASPADES/temp2.txt > ./METASPADES/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./METASPADES/final_bowtie2_merge_ready.txt && cat ./METASPADES/temp3.txt >> ./METASPADES/final_bowtie2_merge_ready.txt
			rm ./METASPADES/temp*.txt

			echo -e "\n======== no contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"
		fi

		# Assign BLAST taxonomy - nt
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST NT OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./METASPADES/CLASSIFICATION/BLAST/BLAST_output_nt_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_first_hit_with_taxonomy.txt ./METASPADES/
		sed '1d' ./METASPADES/BLAST_output_nt_first_hit_with_taxonomy.txt > ./METASPADES/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt
		mv ./METASPADES/BLAST_output_nt_first_hit_with_taxonomy.txt ./METASPADES/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./METASPADES/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./METASPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./METASPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt
		echo

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./METASPADES/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_bitscore_filtered_with_taxonomy.txt ./METASPADES/CLASSIFICATION/BLAST/
		sed '1d' ./METASPADES/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered_with_taxonomy.txt > ./METASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./METASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt > ./METASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./METASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./METASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./METASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./METASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./METASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./METASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./METASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./METASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./METASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./METASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '($3 < 80)' ./METASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./METASPADES/BLAST_output_nt_similarity_filtered.txt

		### Immitating BASTA
		touch ./METASPADES/METASPADES_blast_nt_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./METASPADES/BLAST_output_nt_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./METASPADES/BLAST_output_nt_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./METASPADES/METASPADES_blast_nt_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./METASPADES/METASPADES_blast_nt_LCA_taxonomy_merge.txt && cat ./METASPADES/METASPADES_blast_nt_LCA_taxonomy_merge_noheader.txt >> ./METASPADES/METASPADES_blast_nt_LCA_taxonomy_merge.txt


		# Assign BLAST taxonomy - SILVA
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST SILVA OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./METASPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_first_hit_with_taxonomy.txt ./METASPADES/
		sed '1d' ./METASPADES/BLAST_output_SILVA_first_hit_with_taxonomy.txt > ./METASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt
		mv ./METASPADES/BLAST_output_SILVA_first_hit_with_taxonomy.txt ./METASPADES/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./METASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./METASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./METASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./METASPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt ./METASPADES/CLASSIFICATION/BLAST/
		sed '1d' ./METASPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt > ./METASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./METASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt > ./METASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./METASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./METASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./METASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./METASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./METASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./METASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./METASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./METASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./METASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./METASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '($3 < 80)' ./METASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./METASPADES/BLAST_output_SILVA_similarity_filtered.txt

		### Immitating BASTA
		touch ./METASPADES/METASPADES_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./METASPADES/BLAST_output_SILVA_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./METASPADES/BLAST_output_SILVA_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./METASPADES/METASPADES_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./METASPADES/METASPADES_blast_SILVA_LCA_taxonomy_merge.txt && cat ./METASPADES/METASPADES_blast_SILVA_LCA_taxonomy_merge_noheader.txt >> ./METASPADES/METASPADES_blast_SILVA_LCA_taxonomy_merge.txt

		echo -e "\n======== taxonomy has been assigned to BLAST files ========"
		echo -e "======== BLAST files are ready ========\n"


		# Prepare CREST file - remove 1st/3rd column, add a header
		echo -e "\n======== ADD NCBI TAXONOMY TO CREST OUTPUT ========\n"
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./METASPADES/CLASSIFICATION/CREST/otus.csv ./METASPADES/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./METASPADES/CLASSIFICATION/CREST/CREST_output.txt > ./METASPADES/CREST_seperated.txt
		cut -f2,4 ./METASPADES/CREST_seperated.txt > ./METASPADES/CREST_header.txt
		sed '1d' ./METASPADES/CREST_header.txt > ./METASPADES/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./METASPADES/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./METASPADES/
		sed '1d' ./METASPADES/CREST_tax_ready_with_taxonomy.txt > ./METASPADES/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./METASPADES/CREST_merge.txt && cat ./METASPADES/CREST_tax_ready_with_taxonomy_noheader.txt >> ./METASPADES/CREST_merge.txt

		rm ./METASPADES/CREST_seperated.txt ./METASPADES/CREST_header.txt ./METASPADES/CREST_tax_ready.txt ./METASPADES/CREST_tax_ready_with_taxonomy.txt ./METASPADES/CREST_tax_ready_with_taxonomy_noheader.txt

		echo -e "\n======== CREST file is ready ========\n"

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./METASPADES/CREST_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/CREST_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/BLAST_first_hit_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/BLAST_first_hit_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/BLAST_first_hit_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/BLAST_first_hit_nt_BWA.txt

		merge_mapped_reads_and_contigs.py ./METASPADES/METASPADES_blast_SILVA_LCA_taxonomy_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/BLAST_LCA_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/METASPADES_blast_nt_LCA_taxonomy_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/BLAST_LCA_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./METASPADES/METASPADES_blast_SILVA_LCA_taxonomy_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/BLAST_LCA_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/METASPADES_blast_nt_LCA_taxonomy_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/BLAST_LCA_nt_BWA.txt


		# Move all files - easy to find!!!
		mkdir ./METASPADES/MERGE_FILES
		mv ./METASPADES/*merge*.txt ./METASPADES/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./METASPADES/CREST_BWA.txt > ./METASPADES/new.txt &&  sed 's/_/\t/3' ./METASPADES/new.txt > ./METASPADES/new2.txt
		sed 's/NODE_//g' ./METASPADES/new2.txt > ./METASPADES/new3.txt && sed 's/length_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/cov_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt
		sed '1d' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/CREST_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/CREST_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/CREST_BWA.txt

		sed 's/_/\t/2' ./METASPADES/CREST_BOWTIE2.txt > ./METASPADES/new.txt &&  sed 's/_/\t/3' ./METASPADES/new.txt > ./METASPADES/new2.txt
		sed 's/NODE_//g' ./METASPADES/new2.txt > ./METASPADES/new3.txt && sed 's/length_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/cov_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt
		sed '1d' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/CREST_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/CREST_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./METASPADES/BLAST_first_hit_nt_BWA.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_nt_first_hit_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_nt_first_hit_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_first_hit_nt_BWA.txt

		sed '1d' ./METASPADES/BLAST_first_hit_nt_BOWTIE2.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_nt_first_hit_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_nt_first_hit_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_first_hit_nt_BOWTIE2.txt

		sed '1d' ./METASPADES/BLAST_LCA_nt_BWA.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_nt_LCA_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_nt_LCA_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_LCA_nt_BWA.txt

		sed '1d' ./METASPADES/BLAST_LCA_nt_BOWTIE2.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_nt_LCA_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_nt_LCA_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_LCA_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./METASPADES/BLAST_first_hit_SILVA_BWA.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_SILVA_first_hit_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_SILVA_first_hit_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_first_hit_SILVA_BWA.txt

		sed '1d' ./METASPADES/BLAST_first_hit_SILVA_BOWTIE2.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_SILVA_first_hit_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_SILVA_first_hit_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_first_hit_SILVA_BOWTIE2.txt

		sed '1d' ./METASPADES/BLAST_LCA_SILVA_BWA.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_SILVA_LCA_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_SILVA_LCA_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_LCA_SILVA_BWA.txt

		sed '1d' ./METASPADES/BLAST_LCA_SILVA_BOWTIE2.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_SILVA_LCA_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_SILVA_LCA_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_LCA_SILVA_BOWTIE2.txt


		mkdir ./METASPADES/FINAL_FILES_METASPADES
		mv ./METASPADES/*final.txt ./METASPADES/FINAL_FILES_METASPADES/
		mv ./METASPADES/BLAST_output_nt* ./METASPADES/BLAST_output_SILVA* METASPADES/CLASSIFICATION/BLAST/

		echo -e "\n======== final METASPADES output generated ========\n"

	else
		echo -e "\n======== no FINAL merge generated ========\n"

	fi

else
	echo -e "\n======== no METASPADES output generated ========\n"

fi

echo -e "\n======== RUNNING MEGAHIT ========\n"
if [[ $MEGAHIT == 'true' ]] ; then
	mkdir MEGAHIT
	megahit --presets meta-large -t 16 -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -o ./MEGAHIT
	echo -e "\n======== MEGAHIT DONE ========\n"

	if [[ $MAP == 'true' ]] ; then

		echo -e "\n======== starting bwa index ========\n"

		bwa index -p bwa_index ./MEGAHIT/final.contigs.fa

		echo -e "\n======== bwa index complete. Starting bwa ========\n"

		bwa mem -t 10 bwa_index $Forward_read_trimmed $Reverse_read_trimmed > bwa_output.sam

		rm bwa_index*


		# Output file (.sam) - edit
		samtools view -F 4 bwa_output.sam > mapped_reads_bwa.sam
		samtools view -f 4 bwa_output.sam > unmapped_reads_bwa.sam
		cat mapped_reads_bwa.sam > mapped_reads_bwa.txt
		cat unmapped_reads_bwa.sam > unmapped_reads_bwa.txt
		cut -f3 mapped_reads_bwa.txt > mapped_column3_reads_bwa.txt
		cut -f3 unmapped_reads_bwa.txt > unmapped_column3_reads_bwa.txt
		sort mapped_column3_reads_bwa.txt | uniq -c > sorted_mapped_column3_reads_bwa.txt
		sort unmapped_column3_reads_bwa.txt | uniq -c > sorted_unmapped_column3_reads_bwa.txt
		column -t sorted_mapped_column3_reads_bwa.txt > aligned_mapped_bwa.txt
		column -t sorted_unmapped_column3_reads_bwa.txt > aligned_unmapped_bwa.txt
		sed 's/  */\t/g' aligned_mapped_bwa.txt > out_mapped_bwa.txt
		sed 's/  */\t/g' aligned_unmapped_bwa.txt > out_unmappped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo -e "\n======== bwa complete. Starting bowtie2 index ========\n"

		bowtie2-build -f ./MEGAHIT/final.contigs.fa bowtie_index

		echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"

		bowtie2 -q -x bowtie_index -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -S bowtie2_output.sam

		rm bowtie_index*


		# Output file (.sam) - edit
		samtools view -F 4 bowtie2_output.sam > mapped_reads_bowtie.sam
		samtools view -f 4 bowtie2_output.sam > unmapped_reads_bowtie.sam
		cat mapped_reads_bowtie.sam > mapped_reads_bowtie.txt
		cat unmapped_reads_bowtie.sam > unmapped_reads_bowtie.txt
		cut -f3 mapped_reads_bowtie.txt > mapped_column3_reads_bowtie.txt
		cut -f3 unmapped_reads_bowtie.txt > unmapped_column3_reads_bowtie.txt
		sort mapped_column3_reads_bowtie.txt | uniq -c > sorted_mapped_column3_reads_bowtie.txt
		sort unmapped_column3_reads_bowtie.txt | uniq -c > sorted_unmapped_column3_reads_bowtie.txt
		column -t sorted_mapped_column3_reads_bowtie.txt > aligned_mapped_bowtie.txt
		column -t sorted_unmapped_column3_reads_bowtie.txt > aligned_unmapped_bowtie.txt
		sed 's/  */\t/g' aligned_mapped_bowtie.txt > out_mapped_bowtie.txt
		sed 's/  */\t/g' aligned_unmapped_bowtie.txt > out_unmappped_bowtie.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo -e "\n======== bowtie2 complete ========\n"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ MEGAHIT/

	else

		echo "\n======== no BWA and BOWTIE2 output generated ========\n"

	fi

	if [[ $CLASSIFICATION == 'true' ]] ; then

		# BLAST against nt

		echo -e "\n======== RUNNING JUSTBLAST AGAINST NT AND PROCESSING OUTPUT ========\n"
		justblast ./MEGAHIT/final.contigs.fa/hdd1/databases/nt_database_feb_2020_indexed/nt --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_nt

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_nt | sort -u -k1,1 > BLAST_output_nt_first_hit.txt

		## Blast option 2: apply bitscore threshold, best bitscore %, similarity cutoff, and an LCA approach
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_nt > BLAST_output_nt_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_nt_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_nt_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_nt_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_nt_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_nt_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST NT AND PROCESSING OUTPUT DONE ========\n"

		# BLAST against SILVA
		echo -e "\n======== RUNNING JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT ========\n"
		justblast ./MEGAHIT/final.contigs.fa /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_SILVA

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_SILVA | sort -u -k1,1 > BLAST_output_SILVA_first_hit.txt

		## Blast option 2: apply bitscore threshold and best bitscore %
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_SILVA > BLAST_output_SILVA_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_SILVA_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_SILVA_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_SILVA_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_SILVA_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_SILVA_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT DONE ========\n"

		mkdir MEGAHIT/BLAST
		mv BLAST_output_nt* BLAST_output_SILVA* MEGAHIT/BLAST/

		echo -e "\n======== RUNNING CREST ========\n"
		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./MEGAHIT/final.contigs.fa -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		mv BLAST_output.xml CREST
		mv CREST MEGAHIT/

		mkdir MEGAHIT/CLASSIFICATION
		mv MEGAHIT/BLAST MEGAHIT/CREST MEGAHIT/CLASSIFICATION/

		echo -e "\n======== CREST DONE ========\n"

	else

		echo -e "\n======== NO BLAST OR CREST OUTPUT GENERATED ========\n"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then

			echo -e "\n======== ADD CONTIG SEQUENCES TO BOWTIE AND BWA OUTPUT AND MERGE THEM WITH BLAST+CREST ========\n"
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./MEGAHIT/final.contigs.fa > ./MEGAHIT/fasta_to_tabbed.txt
			mergeFilesOnColumn.pl ./MEGAHIT/BWA/merge_input_mapped_bwa.txt ./MEGAHIT/fasta_to_tabbed.txt 2 1 > ./MEGAHIT/merged_original_bwa.txt
			cut -f1,2,4 ./MEGAHIT/merged_original_bwa.txt > ./MEGAHIT/important_column_3_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./MEGAHIT/important_column_3_bwa.txt > ./MEGAHIT/final_order_bwa.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./MEGAHIT/final_bwa_merge_ready.txt && cat ./MEGAHIT/final_order_bwa.txt >> ./MEGAHIT/final_bwa_merge_ready.txt

			rm ./MEGAHIT/merged_original_bwa.txt ./MEGAHIT/important_column_3_bwa.txt ./MEGAHIT/final_order_bwa.txt

			mergeFilesOnColumn.pl ./MEGAHIT/BOWTIE2/merge_input_mapped_bowtie2.txt ./MEGAHIT/fasta_to_tabbed.txt 2 1 > ./MEGAHIT/merged_original_bowtie2.txt
			cut -f1,2,4 ./MEGAHIT/merged_original_bowtie2.txt > ./MEGAHIT/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./MEGAHIT/important_column_3_bowtie2.txt > ./MEGAHIT/final_order_bowtie2.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./MEGAHIT/final_bowtie2_merge_ready.txt && cat ./MEGAHIT/final_order_bowtie2.txt >> ./MEGAHIT/final_bowtie2_merge_ready.txt

			rm ./MEGAHIT/fasta_to_tabbed.txt ./MEGAHIT/merged_original_bowtie2.txt ./MEGAHIT/important_column_3_bowtie2.txt ./MEGAHIT/final_order_bowtie2.txt

			echo -e "\n======== contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"


		else
			echo -e "\n======== MERGE BOWIE+BWA OUTPUT WITH BLAST+CREST ========\n"
			# Prepares to merge with BLAST or CREST - does not include assembly sequence
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./MEGAHIT/BWA/merge_input_mapped_bwa.txt > ./MEGAHIT/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./MEGAHIT/temp.txt > ./MEGAHIT/temp2.txt
			sed '1d' ./MEGAHIT/temp2.txt > ./MEGAHIT/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./MEGAHIT/final_bwa_merge_ready.txt && cat ./MEGAHIT/temp3.txt >> ./MEGAHIT/final_bwa_merge_ready.txt
 			rm ./MEGAHIT/temp*.txt

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./MEGAHIT/BOWTIE2/merge_input_mapped_bowtie2.txt > ./MEGAHIT/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./MEGAHIT/temp.txt > ./MEGAHIT/temp2.txt
			sed '1d' ./MEGAHIT/temp2.txt > ./MEGAHIT/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./MEGAHIT/final_bowtie2_merge_ready.txt && cat ./MEGAHIT/temp3.txt >> ./MEGAHIT/final_bowtie2_merge_ready.txt
			rm ./MEGAHIT/temp*.txt

			echo -e "\n======== no contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"
		fi

		# Assign BLAST taxonomy - nt
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST NT OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_nt_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_first_hit_with_taxonomy.txt ./MEGAHIT/
		sed '1d' ./MEGAHIT/BLAST_output_nt_first_hit_with_taxonomy.txt > ./MEGAHIT/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt
		mv ./MEGAHIT/BLAST_output_nt_first_hit_with_taxonomy.txt ./MEGAHIT/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./MEGAHIT/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/BLAST_output_nt_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./MEGAHIT/BLAST_output_nt_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt
		echo

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_bitscore_filtered_with_taxonomy.txt ./MEGAHIT/CLASSIFICATION/BLAST/
		sed '1d' ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered_with_taxonomy.txt > ./MEGAHIT/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./MEGAHIT/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt > ./MEGAHIT/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./MEGAHIT/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./MEGAHIT/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./MEGAHIT/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./MEGAHIT/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./MEGAHIT/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./MEGAHIT/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./MEGAHIT/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./MEGAHIT/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./MEGAHIT/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./MEGAHIT/BLAST_output_nt_similarity_filtered.txt
		awk '($3 < 80)' ./MEGAHIT/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./MEGAHIT/BLAST_output_nt_similarity_filtered.txt

		### Immitating BASTA
		touch ./MEGAHIT/MEGAHIT_blast_nt_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./MEGAHIT/BLAST_output_nt_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./MEGAHIT/BLAST_output_nt_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./MEGAHIT/MEGAHIT_blast_nt_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/MEGAHIT_blast_nt_LCA_taxonomy_merge.txt && cat ./MEGAHIT/MEGAHIT_blast_nt_LCA_taxonomy_merge_noheader.txt >> ./MEGAHIT/MEGAHIT_blast_nt_LCA_taxonomy_merge.txt


		# Assign BLAST taxonomy - SILVA
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST SILVA OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_SILVA_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_first_hit_with_taxonomy.txt ./MEGAHIT/
		sed '1d' ./MEGAHIT/BLAST_output_SILVA_first_hit_with_taxonomy.txt > ./MEGAHIT/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt
		mv ./MEGAHIT/BLAST_output_SILVA_first_hit_with_taxonomy.txt ./MEGAHIT/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./MEGAHIT/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./MEGAHIT/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt ./MEGAHIT/CLASSIFICATION/BLAST/
		sed '1d' ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt > ./MEGAHIT/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./MEGAHIT/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt > ./MEGAHIT/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./MEGAHIT/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./MEGAHIT/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./MEGAHIT/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./MEGAHIT/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./MEGAHIT/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./MEGAHIT/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./MEGAHIT/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./MEGAHIT/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./MEGAHIT/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./MEGAHIT/BLAST_output_SILVA_similarity_filtered.txt
		awk '($3 < 80)' ./MEGAHIT/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./MEGAHIT/BLAST_output_SILVA_similarity_filtered.txt

		### Immitating BASTA
		touch ./MEGAHIT/MEGAHIT_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./MEGAHIT/BLAST_output_SILVA_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./MEGAHIT/BLAST_output_SILVA_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./MEGAHIT/MEGAHIT_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/MEGAHIT_blast_SILVA_LCA_taxonomy_merge.txt && cat ./MEGAHIT/MEGAHIT_blast_SILVA_LCA_taxonomy_merge_noheader.txt >> ./MEGAHIT/MEGAHIT_blast_SILVA_LCA_taxonomy_merge.txt

		echo -e "\n======== taxonomy has been assigned to BLAST files ========"
		echo -e "======== BLAST files are ready ========\n"


		# Prepare CREST file - remove 1st/3rd column, add a header
		echo -e "\n======== ADD NCBI TAXONOMY TO CREST OUTPUT ========\n"
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./MEGAHIT/CLASSIFICATION/CREST/otus.csv ./MEGAHIT/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./MEGAHIT/CLASSIFICATION/CREST/CREST_output.txt > ./MEGAHIT/CREST_seperated.txt
		cut -f2,4 ./MEGAHIT/CREST_seperated.txt > ./MEGAHIT/CREST_header.txt
		sed '1d' ./MEGAHIT/CREST_header.txt > ./MEGAHIT/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/
		sed '1d' ./MEGAHIT/CREST_tax_ready_with_taxonomy.txt > ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/CREST_merge.txt && cat ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt >> ./MEGAHIT/CREST_merge.txt

		rm ./MEGAHIT/CREST_seperated.txt ./MEGAHIT/CREST_header.txt ./MEGAHIT/CREST_tax_ready.txt ./MEGAHIT/CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt

		echo -e "\n======== CREST file is ready ========\n"

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./MEGAHIT/CREST_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./MEGAHIT/CREST_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_first_hit_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_first_hit_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_first_hit_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_first_hit_nt_BWA.txt

		merge_mapped_reads_and_contigs.py ./MEGAHIT/MEGAHIT_blast_SILVA_LCA_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_LCA_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./MEGAHIT/MEGAHIT_blast_nt_LCA_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_LCA_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./MEGAHIT/MEGAHIT_blast_SILVA_LCA_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_LCA_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./MEGAHIT/MEGAHIT_blast_nt_LCA_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_LCA_nt_BWA.txt


		# Move all files - easy to find!!!
		mkdir ./MEGAHIT/MERGE_FILES
		mv ./MEGAHIT/*merge*.txt ./MEGAHIT/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./MEGAHIT/CREST_BWA.txt > ./MEGAHIT/new.txt &&  sed 's/_/\t/3' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
		sed 's/NODE_//g' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt && sed 's/length_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/cov_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt
		sed '1d' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/CREST_BWA_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/CREST_BWA_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/CREST_BWA.txt

		sed 's/_/\t/2' ./MEGAHIT/CREST_BOWTIE2.txt > ./MEGAHIT/new.txt &&  sed 's/_/\t/3' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
		sed 's/NODE_//g' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt && sed 's/length_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/cov_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt
		sed '1d' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/CREST_BOWTIE2_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/CREST_BOWTIE2_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./MEGAHIT/BLAST_first_hit_nt_BWA.txt > ./MEGAHIT/new.txt
		sed 's/_/\t/2' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt &&  sed 's/_/\t/3' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
		sed 's/NODE_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/length_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt && sed 's/cov_//g' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_first_hit_BWA_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/BLAST_nt_first_hit_BWA_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_first_hit_nt_BWA.txt

		sed '1d' ./MEGAHIT/BLAST_first_hit_nt_BOWTIE2.txt > ./MEGAHIT/new.txt
		sed 's/_/\t/2' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt &&  sed 's/_/\t/3' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
		sed 's/NODE_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/length_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt && sed 's/cov_//g' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_first_hit_BOWTIE2_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/BLAST_nt_first_hit_BOWTIE2_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_first_hit_nt_BOWTIE2.txt

		sed '1d' ./MEGAHIT/BLAST_LCA_nt_BWA.txt > ./MEGAHIT/new.txt
		sed 's/_/\t/2' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt &&  sed 's/_/\t/3' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
		sed 's/NODE_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/length_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt && sed 's/cov_//g' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_LCA_BWA_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/BLAST_nt_LCA_BWA_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_LCA_nt_BWA.txt

		sed '1d' ./MEGAHIT/BLAST_LCA_nt_BOWTIE2.txt > ./MEGAHIT/new.txt
		sed 's/_/\t/2' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt &&  sed 's/_/\t/3' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
		sed 's/NODE_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/length_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt && sed 's/cov_//g' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_LCA_BOWTIE2_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/BLAST_nt_LCA_BOWTIE2_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_LCA_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./MEGAHIT/BLAST_first_hit_SILVA_BWA.txt > ./MEGAHIT/new.txt
		sed 's/_/\t/2' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt &&  sed 's/_/\t/3' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
		sed 's/NODE_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/length_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt && sed 's/cov_//g' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_first_hit_BWA_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/BLAST_SILVA_first_hit_BWA_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_first_hit_SILVA_BWA.txt

		sed '1d' ./MEGAHIT/BLAST_first_hit_SILVA_BOWTIE2.txt > ./MEGAHIT/new.txt
		sed 's/_/\t/2' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt &&  sed 's/_/\t/3' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
		sed 's/NODE_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/length_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt && sed 's/cov_//g' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_first_hit_BOWTIE2_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/BLAST_SILVA_first_hit_BOWTIE2_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_first_hit_SILVA_BOWTIE2.txt

		sed '1d' ./MEGAHIT/BLAST_LCA_SILVA_BWA.txt > ./MEGAHIT/new.txt
		sed 's/_/\t/2' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt &&  sed 's/_/\t/3' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
		sed 's/NODE_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/length_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt && sed 's/cov_//g' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_LCA_BWA_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/BLAST_SILVA_LCA_BWA_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_LCA_SILVA_BWA.txt

		sed '1d' ./MEGAHIT/BLAST_LCA_SILVA_BOWTIE2.txt > ./MEGAHIT/new.txt
		sed 's/_/\t/2' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt &&  sed 's/_/\t/3' ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
		sed 's/NODE_//g' ./MEGAHIT/new3.txt > ./MEGAHIT/new4.txt && sed 's/length_//g' ./MEGAHIT/new4.txt > ./MEGAHIT/new5.txt && sed 's/cov_//g' ./MEGAHIT/new5.txt > ./MEGAHIT/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_LCA_BOWTIE2_final.txt && cat ./MEGAHIT/new6.txt >> ./MEGAHIT/BLAST_SILVA_LCA_BOWTIE2_final.txt
		rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_LCA_SILVA_BOWTIE2.txt


		mkdir ./MEGAHIT/FINAL_FILES_MEGAHIT
		mv ./MEGAHIT/*final.txt ./MEGAHIT/FINAL_FILES_MEGAHIT/
		mv ./MEGAHIT/BLAST_output_nt* ./MEGAHIT/BLAST_output_SILVA* MEGAHIT/CLASSIFICATION/BLAST/

		echo -e "\n======== final MEGAHIT output generated ========\n"

	else
		echo -e "\n======== no FINAL merge generated ========\n"

	fi

else
	echo -e "\n======== no MEGAHIT output generated ========\n"

fi

echo -e "\n======== RUNNING IDBA_UD ========\n"
if [[ $IDBA_UD == 'true' ]] ; then
	fq2fa --merge --filter $Forward_read_trimmed $Reverse_read_trimmed new_input.fa
	idba_ud --num_threads 16 --pre_correction -r ./new_input.fa -o IDBA_UD
	echo -e "\n======== IDBA_UD DONE ========\n"

	if [[ $MAP == 'true' ]] ; then

		echo -e "\n======== starting bwa index ========\n"

		bwa index -p bwa_index ./IDBA_UD/contig.fa

		echo -e "\n======== bwa index complete. Starting bwa ========\n"

		bwa mem -t 10 bwa_index $Forward_read_trimmed $Reverse_read_trimmed > bwa_output.sam

		rm bwa_index*


		# Output file (.sam) - edit
		samtools view -F 4 bwa_output.sam > mapped_reads_bwa.sam
		samtools view -f 4 bwa_output.sam > unmapped_reads_bwa.sam
		cat mapped_reads_bwa.sam > mapped_reads_bwa.txt
		cat unmapped_reads_bwa.sam > unmapped_reads_bwa.txt
		cut -f3 mapped_reads_bwa.txt > mapped_column3_reads_bwa.txt
		cut -f3 unmapped_reads_bwa.txt > unmapped_column3_reads_bwa.txt
		sort mapped_column3_reads_bwa.txt | uniq -c > sorted_mapped_column3_reads_bwa.txt
		sort unmapped_column3_reads_bwa.txt | uniq -c > sorted_unmapped_column3_reads_bwa.txt
		column -t sorted_mapped_column3_reads_bwa.txt > aligned_mapped_bwa.txt
		column -t sorted_unmapped_column3_reads_bwa.txt > aligned_unmapped_bwa.txt
		sed 's/  */\t/g' aligned_mapped_bwa.txt > out_mapped_bwa.txt
		sed 's/  */\t/g' aligned_unmapped_bwa.txt > out_unmappped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo -e "\n======== bwa complete. Starting bowtie2 index ========\n"

		bowtie2-build -f ./IDBA_UD/contig.fa bowtie_index

		echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"

		bowtie2 -q -x bowtie_index -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -S bowtie2_output.sam

		rm bowtie_index*


		# Output file (.sam) - edit
		samtools view -F 4 bowtie2_output.sam > mapped_reads_bowtie.sam
		samtools view -f 4 bowtie2_output.sam > unmapped_reads_bowtie.sam
		cat mapped_reads_bowtie.sam > mapped_reads_bowtie.txt
		cat unmapped_reads_bowtie.sam > unmapped_reads_bowtie.txt
		cut -f3 mapped_reads_bowtie.txt > mapped_column3_reads_bowtie.txt
		cut -f3 unmapped_reads_bowtie.txt > unmapped_column3_reads_bowtie.txt
		sort mapped_column3_reads_bowtie.txt | uniq -c > sorted_mapped_column3_reads_bowtie.txt
		sort unmapped_column3_reads_bowtie.txt | uniq -c > sorted_unmapped_column3_reads_bowtie.txt
		column -t sorted_mapped_column3_reads_bowtie.txt > aligned_mapped_bowtie.txt
		column -t sorted_unmapped_column3_reads_bowtie.txt > aligned_unmapped_bowtie.txt
		sed 's/  */\t/g' aligned_mapped_bowtie.txt > out_mapped_bowtie.txt
		sed 's/  */\t/g' aligned_unmapped_bowtie.txt > out_unmappped_bowtie.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo -e "\n======== bowtie2 complete ========\n"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ IDBA_UD/

	else

		echo "\n======== no BWA and BOWTIE2 output generated ========\n"

	fi

	if [[ $CLASSIFICATION == 'true' ]] ; then

		# BLAST against nt
		echo -e "\n======== RUNNING JUSTBLAST AGAINST NT AND PROCESSING OUTPUT ========\n"
		justblast ./IDBA_UD/contig.fa/hdd1/databases/nt_database_feb_2020_indexed/nt --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_nt

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_nt | sort -u -k1,1 > BLAST_output_nt_first_hit.txt

		## Blast option 2: apply bitscore threshold, best bitscore %, similarity cutoff, and an LCA approach
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_nt > BLAST_output_nt_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_nt_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_nt_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_nt_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_nt_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_nt_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST NT AND PROCESSING OUTPUT DONE ========\n"

		# BLAST against SILVA
		echo -e "\n======== RUNNING JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT ========\n"
		justblast ./IDBA_UD/contig.fa /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_SILVA

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_SILVA | sort -u -k1,1 > BLAST_output_SILVA_first_hit.txt

		## Blast option 2: apply bitscore threshold and best bitscore %
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_SILVA > BLAST_output_SILVA_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_SILVA_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_SILVA_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_SILVA_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_SILVA_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_SILVA_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT DONE ========\n"

		mkdir IDBA_UD/BLAST
		mv BLAST_output_nt* BLAST_output_SILVA* IDBA_UD/BLAST/

		echo -e "\n======== RUNNING CREST ========\n"
		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./IDBA_UD/contig.fa -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		mv BLAST_output.xml CREST
		mv CREST IDBA_UD/

		mkdir IDBA_UD/CLASSIFICATION
		mv IDBA_UD/BLAST IDBA_UD/CREST IDBA_UD/CLASSIFICATION/

		echo -e "\n======== CREST DONE ========\n"

	else

		echo -e "\n======== NO BLAST OR CREST OUTPUT GENERATED ========\n"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then

			echo -e "\n======== ADD CONTIG SEQUENCES TO BOWTIE AND BWA OUTPUT AND MERGE THEM WITH BLAST+CREST ========\n"
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./IDBA_UD/contig.fa > ./IDBA_UD/fasta_to_tabbed.txt
			mergeFilesOnColumn.pl ./IDBA_UD/BWA/merge_input_mapped_bwa.txt ./IDBA_UD/fasta_to_tabbed.txt 2 1 > ./IDBA_UD/merged_original_bwa.txt
			cut -f1,2,4 ./IDBA_UD/merged_original_bwa.txt > ./IDBA_UD/important_column_3_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./IDBA_UD/important_column_3_bwa.txt > ./IDBA_UD/final_order_bwa.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA_UD/final_bwa_merge_ready.txt && cat ./IDBA_UD/final_order_bwa.txt >> ./IDBA_UD/final_bwa_merge_ready.txt

			rm ./IDBA_UD/merged_original_bwa.txt ./IDBA_UD/important_column_3_bwa.txt ./IDBA_UD/final_order_bwa.txt

			mergeFilesOnColumn.pl ./IDBA_UD/BOWTIE2/merge_input_mapped_bowtie2.txt ./IDBA_UD/fasta_to_tabbed.txt 2 1 > ./IDBA_UD/merged_original_bowtie2.txt
			cut -f1,2,4 ./IDBA_UD/merged_original_bowtie2.txt > ./IDBA_UD/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./IDBA_UD/important_column_3_bowtie2.txt > ./IDBA_UD/final_order_bowtie2.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA_UD/final_bowtie2_merge_ready.txt && cat ./IDBA_UD/final_order_bowtie2.txt >> ./IDBA_UD/final_bowtie2_merge_ready.txt

			rm ./IDBA_UD/fasta_to_tabbed.txt ./IDBA_UD/merged_original_bowtie2.txt ./IDBA_UD/important_column_3_bowtie2.txt ./IDBA_UD/final_order_bowtie2.txt

			echo -e "\n======== contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"


		else
			echo -e "\n======== MERGE BOWIE+BWA OUTPUT WITH BLAST+CREST ========\n"
			# Prepares to merge with BLAST or CREST - does not include assembly sequence
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./IDBA_UD/BWA/merge_input_mapped_bwa.txt > ./IDBA_UD/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./IDBA_UD/temp.txt > ./IDBA_UD/temp2.txt
			sed '1d' ./IDBA_UD/temp2.txt > ./IDBA_UD/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA_UD/final_bwa_merge_ready.txt && cat ./IDBA_UD/temp3.txt >> ./IDBA_UD/final_bwa_merge_ready.txt
 			rm ./IDBA_UD/temp*.txt

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./IDBA_UD/BOWTIE2/merge_input_mapped_bowtie2.txt > ./IDBA_UD/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./IDBA_UD/temp.txt > ./IDBA_UD/temp2.txt
			sed '1d' ./IDBA_UD/temp2.txt > ./IDBA_UD/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA_UD/final_bowtie2_merge_ready.txt && cat ./IDBA_UD/temp3.txt >> ./IDBA_UD/final_bowtie2_merge_ready.txt
			rm ./IDBA_UD/temp*.txt

			echo -e "\n======== no contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"
		fi

		# Assign BLAST taxonomy - nt
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST NT OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_nt_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_first_hit_with_taxonomy.txt ./IDBA_UD/
		sed '1d' ./IDBA_UD/BLAST_output_nt_first_hit_with_taxonomy.txt > ./IDBA_UD/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt
		mv ./IDBA_UD/BLAST_output_nt_first_hit_with_taxonomy.txt ./IDBA_UD/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./IDBA_UD/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/BLAST_output_nt_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./IDBA_UD/BLAST_output_nt_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt
		echo

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_bitscore_filtered_with_taxonomy.txt ./IDBA_UD/CLASSIFICATION/BLAST/
		sed '1d' ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered_with_taxonomy.txt > ./IDBA_UD/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./IDBA_UD/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt > ./IDBA_UD/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./IDBA_UD/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./IDBA_UD/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./IDBA_UD/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./IDBA_UD/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./IDBA_UD/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./IDBA_UD/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./IDBA_UD/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./IDBA_UD/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./IDBA_UD/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./IDBA_UD/BLAST_output_nt_similarity_filtered.txt
		awk '($3 < 80)' ./IDBA_UD/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./IDBA_UD/BLAST_output_nt_similarity_filtered.txt

		### Immitating BASTA
		touch ./IDBA_UD/IDBA_UD_blast_nt_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./IDBA_UD/BLAST_output_nt_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./IDBA_UD/BLAST_output_nt_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./IDBA_UD/IDBA_UD_blast_nt_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/IDBA_UD_blast_nt_LCA_taxonomy_merge.txt && cat ./IDBA_UD/IDBA_UD_blast_nt_LCA_taxonomy_merge_noheader.txt >> ./IDBA_UD/IDBA_UD_blast_nt_LCA_taxonomy_merge.txt


		# Assign BLAST taxonomy - SILVA
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST SILVA OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_SILVA_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_first_hit_with_taxonomy.txt ./IDBA_UD/
		sed '1d' ./IDBA_UD/BLAST_output_SILVA_first_hit_with_taxonomy.txt > ./IDBA_UD/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt
		mv ./IDBA_UD/BLAST_output_SILVA_first_hit_with_taxonomy.txt ./IDBA_UD/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./IDBA_UD/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./IDBA_UD/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt ./IDBA_UD/CLASSIFICATION/BLAST/
		sed '1d' ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt > ./IDBA_UD/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./IDBA_UD/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt > ./IDBA_UD/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./IDBA_UD/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./IDBA_UD/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./IDBA_UD/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./IDBA_UD/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./IDBA_UD/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./IDBA_UD/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./IDBA_UD/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./IDBA_UD/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./IDBA_UD/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./IDBA_UD/BLAST_output_SILVA_similarity_filtered.txt
		awk '($3 < 80)' ./IDBA_UD/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./IDBA_UD/BLAST_output_SILVA_similarity_filtered.txt

		### Immitating BASTA
		touch ./IDBA_UD/IDBA_UD_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./IDBA_UD/BLAST_output_SILVA_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./IDBA_UD/BLAST_output_SILVA_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./IDBA_UD/IDBA_UD_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/IDBA_UD_blast_SILVA_LCA_taxonomy_merge.txt && cat ./IDBA_UD/IDBA_UD_blast_SILVA_LCA_taxonomy_merge_noheader.txt >> ./IDBA_UD/IDBA_UD_blast_SILVA_LCA_taxonomy_merge.txt

		echo -e "\n======== taxonomy has been assigned to BLAST files ========"
		echo -e "======== BLAST files are ready ========\n"


		# Prepare CREST file - remove 1st/3rd column, add a header
		echo -e "\n======== ADD NCBI TAXONOMY TO CREST OUTPUT ========\n"
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA_UD/CLASSIFICATION/CREST/otus.csv ./IDBA_UD/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./IDBA_UD/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA_UD/CREST_seperated.txt
		cut -f2,4 ./IDBA_UD/CREST_seperated.txt > ./IDBA_UD/CREST_header.txt
		sed '1d' ./IDBA_UD/CREST_header.txt > ./IDBA_UD/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./IDBA_UD/
		sed '1d' ./IDBA_UD/CREST_tax_ready_with_taxonomy.txt > ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/CREST_merge.txt && cat ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA_UD/CREST_merge.txt

		rm ./IDBA_UD/CREST_seperated.txt ./IDBA_UD/CREST_header.txt ./IDBA_UD/CREST_tax_ready.txt ./IDBA_UD/CREST_tax_ready_with_taxonomy.txt ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt

		echo -e "\n======== CREST file is ready ========\n"

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./IDBA_UD/CREST_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./IDBA_UD/CREST_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_first_hit_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_first_hit_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_first_hit_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_first_hit_nt_BWA.txt

		merge_mapped_reads_and_contigs.py ./IDBA_UD/IDBA_UD_blast_SILVA_LCA_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_LCA_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./IDBA_UD/IDBA_UD_blast_nt_LCA_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_LCA_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./IDBA_UD/IDBA_UD_blast_SILVA_LCA_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_LCA_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./IDBA_UD/IDBA_UD_blast_nt_LCA_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_LCA_nt_BWA.txt


		# Move all files - easy to find!!!
		mkdir ./IDBA_UD/MERGE_FILES
		mv ./IDBA_UD/*merge*.txt ./IDBA_UD/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./IDBA_UD/CREST_BWA.txt > ./IDBA_UD/new.txt &&  sed 's/_/\t/3' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
		sed 's/NODE_//g' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt && sed 's/length_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/cov_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt
		sed '1d' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/CREST_BWA_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/CREST_BWA_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/CREST_BWA.txt

		sed 's/_/\t/2' ./IDBA_UD/CREST_BOWTIE2.txt > ./IDBA_UD/new.txt &&  sed 's/_/\t/3' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
		sed 's/NODE_//g' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt && sed 's/length_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/cov_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt
		sed '1d' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/CREST_BOWTIE2_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/CREST_BOWTIE2_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./IDBA_UD/BLAST_first_hit_nt_BWA.txt > ./IDBA_UD/new.txt
		sed 's/_/\t/2' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt &&  sed 's/_/\t/3' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
		sed 's/NODE_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/length_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt && sed 's/cov_//g' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_first_hit_BWA_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/BLAST_nt_first_hit_BWA_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_first_hit_nt_BWA.txt

		sed '1d' ./IDBA_UD/BLAST_first_hit_nt_BOWTIE2.txt > ./IDBA_UD/new.txt
		sed 's/_/\t/2' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt &&  sed 's/_/\t/3' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
		sed 's/NODE_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/length_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt && sed 's/cov_//g' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_first_hit_BOWTIE2_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/BLAST_nt_first_hit_BOWTIE2_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_first_hit_nt_BOWTIE2.txt

		sed '1d' ./IDBA_UD/BLAST_LCA_nt_BWA.txt > ./IDBA_UD/new.txt
		sed 's/_/\t/2' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt &&  sed 's/_/\t/3' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
		sed 's/NODE_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/length_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt && sed 's/cov_//g' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_LCA_BWA_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/BLAST_nt_LCA_BWA_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_LCA_nt_BWA.txt

		sed '1d' ./IDBA_UD/BLAST_LCA_nt_BOWTIE2.txt > ./IDBA_UD/new.txt
		sed 's/_/\t/2' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt &&  sed 's/_/\t/3' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
		sed 's/NODE_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/length_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt && sed 's/cov_//g' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_LCA_BOWTIE2_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/BLAST_nt_LCA_BOWTIE2_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_LCA_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./IDBA_UD/BLAST_first_hit_SILVA_BWA.txt > ./IDBA_UD/new.txt
		sed 's/_/\t/2' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt &&  sed 's/_/\t/3' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
		sed 's/NODE_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/length_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt && sed 's/cov_//g' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_first_hit_BWA_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/BLAST_SILVA_first_hit_BWA_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_first_hit_SILVA_BWA.txt

		sed '1d' ./IDBA_UD/BLAST_first_hit_SILVA_BOWTIE2.txt > ./IDBA_UD/new.txt
		sed 's/_/\t/2' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt &&  sed 's/_/\t/3' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
		sed 's/NODE_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/length_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt && sed 's/cov_//g' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_first_hit_BOWTIE2_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/BLAST_SILVA_first_hit_BOWTIE2_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_first_hit_SILVA_BOWTIE2.txt

		sed '1d' ./IDBA_UD/BLAST_LCA_SILVA_BWA.txt > ./IDBA_UD/new.txt
		sed 's/_/\t/2' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt &&  sed 's/_/\t/3' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
		sed 's/NODE_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/length_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt && sed 's/cov_//g' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_LCA_BWA_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/BLAST_SILVA_LCA_BWA_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_LCA_SILVA_BWA.txt

		sed '1d' ./IDBA_UD/BLAST_LCA_SILVA_BOWTIE2.txt > ./IDBA_UD/new.txt
		sed 's/_/\t/2' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt &&  sed 's/_/\t/3' ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
		sed 's/NODE_//g' ./IDBA_UD/new3.txt > ./IDBA_UD/new4.txt && sed 's/length_//g' ./IDBA_UD/new4.txt > ./IDBA_UD/new5.txt && sed 's/cov_//g' ./IDBA_UD/new5.txt > ./IDBA_UD/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_LCA_BOWTIE2_final.txt && cat ./IDBA_UD/new6.txt >> ./IDBA_UD/BLAST_SILVA_LCA_BOWTIE2_final.txt
		rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_LCA_SILVA_BOWTIE2.txt


		mkdir ./IDBA_UD/FINAL_FILES_IDBA_UD
		mv ./IDBA_UD/*final.txt ./IDBA_UD/FINAL_FILES_IDBA_UD/
		mv ./IDBA_UD/BLAST_output_nt* ./IDBA_UD/BLAST_output_SILVA* IDBA_UD/CLASSIFICATION/BLAST/

		echo -e "\n======== final IDBA_UD output generated ========\n"

	else
		echo -e "\n======== no FINAL merge generated ========\n"

	fi

else
	echo -e "\n======== no IDBA_UD output generated ========\n"

fi

mkdir ASSEMBLERS
mv SPADES METASPADES MEGAHIT IDBA_UD ASSEMBLERS/


######################### Beginning of RNA pipelines ########################
cd ..
mkdir PIPELINE_RNA
cd PIPELINE_RNA/

echo -e "\n======== RUNNING SORTMERNA ========\n"
mkdir SORTMERNA
sortmerna --ref /hdd1/databases/sortmerna_silva_databases/silva-bac-16s-id90.fasta --ref /hdd1/databases/sortmerna_silva_databases/silva-arc-16s-id95.fasta --ref /hdd1/databases/sortmerna_silva_databases/silva-euk-18s-id95.fasta --reads $Forward_read_trimmed --reads $Reverse_read_trimmed --paired_in -other -fastx 1 -num_alignments 1 -v -workdir SORTMERNA
echo -e "\n======== SORTMERNA DONE ========\n"

echo -e "\n======== RUNNING RNASPADES ========\n"
if [[ $RNASPADES == 'true' ]] ; then
	mkdir RNASPADES
	rnaspades.py --12 ./SORTMERNA/out/aligned.fq -o RNASPADES
	echo -e "\n======== RNASPADES DONE ========\n"

	if [[ $MAP == 'true' ]] ; then

		echo -e "\n======== starting bwa index ========\n"

		bwa index -p bwa_index ./RNASPADES/transcripts.fasta

		echo -e "\n======== bwa index complete. Starting bwa ========\n"

		bwa mem -t 10 bwa_index ./SORTMERNA/out/aligned.fq > bwa_output.sam

		rm bwa_index*


		# Output file (.sam) - edit
		samtools view -F 4 bwa_output.sam > mapped_reads_bwa.sam
		samtools view -f 4 bwa_output.sam > unmapped_reads_bwa.sam
		cat mapped_reads_bwa.sam > mapped_reads_bwa.txt
		cat unmapped_reads_bwa.sam > unmapped_reads_bwa.txt
		cut -f3 mapped_reads_bwa.txt > mapped_column3_reads_bwa.txt
		cut -f3 unmapped_reads_bwa.txt > unmapped_column3_reads_bwa.txt
		sort mapped_column3_reads_bwa.txt | uniq -c > sorted_mapped_column3_reads_bwa.txt
		sort unmapped_column3_reads_bwa.txt | uniq -c > sorted_unmapped_column3_reads_bwa.txt
		column -t sorted_mapped_column3_reads_bwa.txt > aligned_mapped_bwa.txt
		column -t sorted_unmapped_column3_reads_bwa.txt > aligned_unmapped_bwa.txt
		sed 's/  */\t/g' aligned_mapped_bwa.txt > out_mapped_bwa.txt
		sed 's/  */\t/g' aligned_unmapped_bwa.txt > out_unmappped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo -e "\n======== bwa complete. Starting bowtie2 index ========\n"

		bowtie2-build -f ./RNASPADES/transcripts.fasta bowtie_index

		echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"

		bowtie2 -q -x bowtie_index --interleaved ./SORTMERNA/out/aligned.fq -S bowtie2_output.sam

		rm bowtie_index*


		# Output file (.sam) - edit
		samtools view -F 4 bowtie2_output.sam > mapped_reads_bowtie.sam
		samtools view -f 4 bowtie2_output.sam > unmapped_reads_bowtie.sam
		cat mapped_reads_bowtie.sam > mapped_reads_bowtie.txt
		cat unmapped_reads_bowtie.sam > unmapped_reads_bowtie.txt
		cut -f3 mapped_reads_bowtie.txt > mapped_column3_reads_bowtie.txt
		cut -f3 unmapped_reads_bowtie.txt > unmapped_column3_reads_bowtie.txt
		sort mapped_column3_reads_bowtie.txt | uniq -c > sorted_mapped_column3_reads_bowtie.txt
		sort unmapped_column3_reads_bowtie.txt | uniq -c > sorted_unmapped_column3_reads_bowtie.txt
		column -t sorted_mapped_column3_reads_bowtie.txt > aligned_mapped_bowtie.txt
		column -t sorted_unmapped_column3_reads_bowtie.txt > aligned_unmapped_bowtie.txt
		sed 's/  */\t/g' aligned_mapped_bowtie.txt > out_mapped_bowtie.txt
		sed 's/  */\t/g' aligned_unmapped_bowtie.txt > out_unmappped_bowtie.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo -e "\n======== bowtie2 complete ========\n"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ RNASPADES/

	else

		echo "\n======== no BWA and BOWTIE2 output generated ========\n"

	fi

	if [[ $CLASSIFICATION == 'true' ]] ; then

		# BLAST against nt
		echo -e "\n======== RUNNING JUSTBLAST AGAINST NT AND PROCESSING OUTPUT ========\n"
		justblast ./RNASPADES/transcripts.fasta/hdd1/databases/nt_database_feb_2020_indexed/nt --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_nt

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_nt | sort -u -k1,1 > BLAST_output_nt_first_hit.txt

		## Blast option 2: apply bitscore threshold, best bitscore %, similarity cutoff, and an LCA approach
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_nt > BLAST_output_nt_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_nt_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_nt_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_nt_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_nt_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_nt_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST NT AND PROCESSING OUTPUT DONE ========\n"

		# BLAST against SILVA
		echo -e "\n======== RUNNING JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT ========\n"
		justblast ./RNASPADES/transcripts.fasta /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_SILVA

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_SILVA | sort -u -k1,1 > BLAST_output_SILVA_first_hit.txt

		## Blast option 2: apply bitscore threshold and best bitscore %
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_SILVA > BLAST_output_SILVA_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_SILVA_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_SILVA_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_SILVA_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_SILVA_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_SILVA_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT DONE ========\n"

		mkdir RNASPADES/BLAST
		mv BLAST_output_nt* BLAST_output_SILVA* RNASPADES/BLAST/

		echo -e "\n======== RUNNING CREST ========\n"
		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./RNASPADES/transcripts.fasta -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		mv BLAST_output.xml CREST
		mv CREST RNASPADES/

		mkdir RNASPADES/CLASSIFICATION
		mv RNASPADES/BLAST RNASPADES/CREST RNASPADES/CLASSIFICATION/

		echo -e "\n======== CREST DONE ========\n"

	else

		echo -e "\n======== NO BLAST OR CREST OUTPUT GENERATED ========\n"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then

			echo -e "\n======== ADD CONTIG SEQUENCES TO BOWTIE AND BWA OUTPUT AND MERGE THEM WITH BLAST+CREST ========\n"
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./RNASPADES/transcripts.fasta > ./RNASPADES/fasta_to_tabbed.txt
			mergeFilesOnColumn.pl ./RNASPADES/BWA/merge_input_mapped_bwa.txt ./RNASPADES/fasta_to_tabbed.txt 2 1 > ./RNASPADES/merged_original_bwa.txt
			cut -f1,2,4 ./RNASPADES/merged_original_bwa.txt > ./RNASPADES/important_column_3_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./RNASPADES/important_column_3_bwa.txt > ./RNASPADES/final_order_bwa.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./RNASPADES/final_bwa_merge_ready.txt && cat ./RNASPADES/final_order_bwa.txt >> ./RNASPADES/final_bwa_merge_ready.txt

			rm ./RNASPADES/merged_original_bwa.txt ./RNASPADES/important_column_3_bwa.txt ./RNASPADES/final_order_bwa.txt

			mergeFilesOnColumn.pl ./RNASPADES/BOWTIE2/merge_input_mapped_bowtie2.txt ./RNASPADES/fasta_to_tabbed.txt 2 1 > ./RNASPADES/merged_original_bowtie2.txt
			cut -f1,2,4 ./RNASPADES/merged_original_bowtie2.txt > ./RNASPADES/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./RNASPADES/important_column_3_bowtie2.txt > ./RNASPADES/final_order_bowtie2.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./RNASPADES/final_bowtie2_merge_ready.txt && cat ./RNASPADES/final_order_bowtie2.txt >> ./RNASPADES/final_bowtie2_merge_ready.txt

			rm ./RNASPADES/fasta_to_tabbed.txt ./RNASPADES/merged_original_bowtie2.txt ./RNASPADES/important_column_3_bowtie2.txt ./RNASPADES/final_order_bowtie2.txt

			echo -e "\n======== contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"


		else
			echo -e "\n======== MERGE BOWIE+BWA OUTPUT WITH BLAST+CREST ========\n"
			# Prepares to merge with BLAST or CREST - does not include assembly sequence
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./RNASPADES/BWA/merge_input_mapped_bwa.txt > ./RNASPADES/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./RNASPADES/temp.txt > ./RNASPADES/temp2.txt
			sed '1d' ./RNASPADES/temp2.txt > ./RNASPADES/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./RNASPADES/final_bwa_merge_ready.txt && cat ./RNASPADES/temp3.txt >> ./RNASPADES/final_bwa_merge_ready.txt
 			rm ./RNASPADES/temp*.txt

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./RNASPADES/BOWTIE2/merge_input_mapped_bowtie2.txt > ./RNASPADES/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./RNASPADES/temp.txt > ./RNASPADES/temp2.txt
			sed '1d' ./RNASPADES/temp2.txt > ./RNASPADES/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./RNASPADES/final_bowtie2_merge_ready.txt && cat ./RNASPADES/temp3.txt >> ./RNASPADES/final_bowtie2_merge_ready.txt
			rm ./RNASPADES/temp*.txt

			echo -e "\n======== no contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"
		fi

		# Assign BLAST taxonomy - nt
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST NT OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./RNASPADES/CLASSIFICATION/BLAST/BLAST_output_nt_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_first_hit_with_taxonomy.txt ./RNASPADES/
		sed '1d' ./RNASPADES/BLAST_output_nt_first_hit_with_taxonomy.txt > ./RNASPADES/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt
		mv ./RNASPADES/BLAST_output_nt_first_hit_with_taxonomy.txt ./RNASPADES/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./RNASPADES/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./RNASPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./RNASPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt
		echo

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./RNASPADES/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_bitscore_filtered_with_taxonomy.txt ./RNASPADES/CLASSIFICATION/BLAST/
		sed '1d' ./RNASPADES/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered_with_taxonomy.txt > ./RNASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./RNASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt > ./RNASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./RNASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./RNASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./RNASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./RNASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./RNASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./RNASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./RNASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./RNASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./RNASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./RNASPADES/BLAST_output_nt_similarity_filtered.txt
		awk '($3 < 80)' ./RNASPADES/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./RNASPADES/BLAST_output_nt_similarity_filtered.txt

		### Immitating BASTA
		touch ./RNASPADES/RNASPADES_blast_nt_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./RNASPADES/BLAST_output_nt_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./RNASPADES/BLAST_output_nt_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./RNASPADES/RNASPADES_blast_nt_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./RNASPADES/RNASPADES_blast_nt_LCA_taxonomy_merge.txt && cat ./RNASPADES/RNASPADES_blast_nt_LCA_taxonomy_merge_noheader.txt >> ./RNASPADES/RNASPADES_blast_nt_LCA_taxonomy_merge.txt


		# Assign BLAST taxonomy - SILVA
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST SILVA OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./RNASPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_first_hit_with_taxonomy.txt ./RNASPADES/
		sed '1d' ./RNASPADES/BLAST_output_SILVA_first_hit_with_taxonomy.txt > ./RNASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt
		mv ./RNASPADES/BLAST_output_SILVA_first_hit_with_taxonomy.txt ./RNASPADES/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./RNASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./RNASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./RNASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./RNASPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt ./RNASPADES/CLASSIFICATION/BLAST/
		sed '1d' ./RNASPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt > ./RNASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./RNASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt > ./RNASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./RNASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./RNASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./RNASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./RNASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./RNASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./RNASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./RNASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./RNASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./RNASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./RNASPADES/BLAST_output_SILVA_similarity_filtered.txt
		awk '($3 < 80)' ./RNASPADES/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./RNASPADES/BLAST_output_SILVA_similarity_filtered.txt

		### Immitating BASTA
		touch ./RNASPADES/RNASPADES_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./RNASPADES/BLAST_output_SILVA_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./RNASPADES/BLAST_output_SILVA_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./RNASPADES/RNASPADES_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./RNASPADES/RNASPADES_blast_SILVA_LCA_taxonomy_merge.txt && cat ./RNASPADES/RNASPADES_blast_SILVA_LCA_taxonomy_merge_noheader.txt >> ./RNASPADES/RNASPADES_blast_SILVA_LCA_taxonomy_merge.txt

		echo -e "\n======== taxonomy has been assigned to BLAST files ========"
		echo -e "======== BLAST files are ready ========\n"


		# Prepare CREST file - remove 1st/3rd column, add a header
		echo -e "\n======== ADD NCBI TAXONOMY TO CREST OUTPUT ========\n"
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./RNASPADES/CLASSIFICATION/CREST/otus.csv ./RNASPADES/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./RNASPADES/CLASSIFICATION/CREST/CREST_output.txt > ./RNASPADES/CREST_seperated.txt
		cut -f2,4 ./RNASPADES/CREST_seperated.txt > ./RNASPADES/CREST_header.txt
		sed '1d' ./RNASPADES/CREST_header.txt > ./RNASPADES/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./RNASPADES/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./RNASPADES/
		sed '1d' ./RNASPADES/CREST_tax_ready_with_taxonomy.txt > ./RNASPADES/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./RNASPADES/CREST_merge.txt && cat ./RNASPADES/CREST_tax_ready_with_taxonomy_noheader.txt >> ./RNASPADES/CREST_merge.txt

		rm ./RNASPADES/CREST_seperated.txt ./RNASPADES/CREST_header.txt ./RNASPADES/CREST_tax_ready.txt ./RNASPADES/CREST_tax_ready_with_taxonomy.txt ./RNASPADES/CREST_tax_ready_with_taxonomy_noheader.txt

		echo -e "\n======== CREST file is ready ========\n"

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./RNASPADES/CREST_merge.txt ./RNASPADES/final_bwa_merge_ready.txt ./RNASPADES/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./RNASPADES/CREST_merge.txt ./RNASPADES/final_bowtie2_merge_ready.txt ./RNASPADES/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./RNASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./RNASPADES/final_bowtie2_merge_ready.txt ./RNASPADES/BLAST_first_hit_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./RNASPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./RNASPADES/final_bowtie2_merge_ready.txt ./RNASPADES/BLAST_first_hit_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./RNASPADES/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./RNASPADES/final_bwa_merge_ready.txt ./RNASPADES/BLAST_first_hit_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./RNASPADES/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./RNASPADES/final_bwa_merge_ready.txt ./RNASPADES/BLAST_first_hit_nt_BWA.txt

		merge_mapped_reads_and_contigs.py ./RNASPADES/RNASPADES_blast_SILVA_LCA_taxonomy_merge.txt ./RNASPADES/final_bowtie2_merge_ready.txt ./RNASPADES/BLAST_LCA_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./RNASPADES/RNASPADES_blast_nt_LCA_taxonomy_merge.txt ./RNASPADES/final_bowtie2_merge_ready.txt ./RNASPADES/BLAST_LCA_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./RNASPADES/RNASPADES_blast_SILVA_LCA_taxonomy_merge.txt ./RNASPADES/final_bwa_merge_ready.txt ./RNASPADES/BLAST_LCA_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./RNASPADES/RNASPADES_blast_nt_LCA_taxonomy_merge.txt ./RNASPADES/final_bwa_merge_ready.txt ./RNASPADES/BLAST_LCA_nt_BWA.txt


		# Move all files - easy to find!!!
		mkdir ./RNASPADES/MERGE_FILES
		mv ./RNASPADES/*merge*.txt ./RNASPADES/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./RNASPADES/CREST_BWA.txt > ./RNASPADES/new.txt &&  sed 's/_/\t/3' ./RNASPADES/new.txt > ./RNASPADES/new2.txt
		sed 's/NODE_//g' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt && sed 's/length_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/cov_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt
		sed '1d' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/CREST_BWA_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/CREST_BWA_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/CREST_BWA.txt

		sed 's/_/\t/2' ./RNASPADES/CREST_BOWTIE2.txt > ./RNASPADES/new.txt &&  sed 's/_/\t/3' ./RNASPADES/new.txt > ./RNASPADES/new2.txt
		sed 's/NODE_//g' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt && sed 's/length_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/cov_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt
		sed '1d' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/CREST_BOWTIE2_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/CREST_BOWTIE2_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./RNASPADES/BLAST_first_hit_nt_BWA.txt > ./RNASPADES/new.txt
		sed 's/_/\t/2' ./RNASPADES/new.txt > ./RNASPADES/new2.txt &&  sed 's/_/\t/3' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt
		sed 's/NODE_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/length_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt && sed 's/cov_//g' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/BLAST_nt_first_hit_BWA_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/BLAST_nt_first_hit_BWA_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/BLAST_first_hit_nt_BWA.txt

		sed '1d' ./RNASPADES/BLAST_first_hit_nt_BOWTIE2.txt > ./RNASPADES/new.txt
		sed 's/_/\t/2' ./RNASPADES/new.txt > ./RNASPADES/new2.txt &&  sed 's/_/\t/3' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt
		sed 's/NODE_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/length_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt && sed 's/cov_//g' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/BLAST_nt_first_hit_BOWTIE2_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/BLAST_nt_first_hit_BOWTIE2_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/BLAST_first_hit_nt_BOWTIE2.txt

		sed '1d' ./RNASPADES/BLAST_LCA_nt_BWA.txt > ./RNASPADES/new.txt
		sed 's/_/\t/2' ./RNASPADES/new.txt > ./RNASPADES/new2.txt &&  sed 's/_/\t/3' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt
		sed 's/NODE_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/length_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt && sed 's/cov_//g' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/BLAST_nt_LCA_BWA_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/BLAST_nt_LCA_BWA_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/BLAST_LCA_nt_BWA.txt

		sed '1d' ./RNASPADES/BLAST_LCA_nt_BOWTIE2.txt > ./RNASPADES/new.txt
		sed 's/_/\t/2' ./RNASPADES/new.txt > ./RNASPADES/new2.txt &&  sed 's/_/\t/3' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt
		sed 's/NODE_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/length_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt && sed 's/cov_//g' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/BLAST_nt_LCA_BOWTIE2_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/BLAST_nt_LCA_BOWTIE2_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/BLAST_LCA_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./RNASPADES/BLAST_first_hit_SILVA_BWA.txt > ./RNASPADES/new.txt
		sed 's/_/\t/2' ./RNASPADES/new.txt > ./RNASPADES/new2.txt &&  sed 's/_/\t/3' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt
		sed 's/NODE_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/length_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt && sed 's/cov_//g' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/BLAST_SILVA_first_hit_BWA_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/BLAST_SILVA_first_hit_BWA_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/BLAST_first_hit_SILVA_BWA.txt

		sed '1d' ./RNASPADES/BLAST_first_hit_SILVA_BOWTIE2.txt > ./RNASPADES/new.txt
		sed 's/_/\t/2' ./RNASPADES/new.txt > ./RNASPADES/new2.txt &&  sed 's/_/\t/3' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt
		sed 's/NODE_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/length_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt && sed 's/cov_//g' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/BLAST_SILVA_first_hit_BOWTIE2_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/BLAST_SILVA_first_hit_BOWTIE2_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/BLAST_first_hit_SILVA_BOWTIE2.txt

		sed '1d' ./RNASPADES/BLAST_LCA_SILVA_BWA.txt > ./RNASPADES/new.txt
		sed 's/_/\t/2' ./RNASPADES/new.txt > ./RNASPADES/new2.txt &&  sed 's/_/\t/3' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt
		sed 's/NODE_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/length_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt && sed 's/cov_//g' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/BLAST_SILVA_LCA_BWA_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/BLAST_SILVA_LCA_BWA_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/BLAST_LCA_SILVA_BWA.txt

		sed '1d' ./RNASPADES/BLAST_LCA_SILVA_BOWTIE2.txt > ./RNASPADES/new.txt
		sed 's/_/\t/2' ./RNASPADES/new.txt > ./RNASPADES/new2.txt &&  sed 's/_/\t/3' ./RNASPADES/new2.txt > ./RNASPADES/new3.txt
		sed 's/NODE_//g' ./RNASPADES/new3.txt > ./RNASPADES/new4.txt && sed 's/length_//g' ./RNASPADES/new4.txt > ./RNASPADES/new5.txt && sed 's/cov_//g' ./RNASPADES/new5.txt > ./RNASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./RNASPADES/BLAST_SILVA_LCA_BOWTIE2_final.txt && cat ./RNASPADES/new6.txt >> ./RNASPADES/BLAST_SILVA_LCA_BOWTIE2_final.txt
		rm ./RNASPADES/new*.txt ./RNASPADES/BLAST_LCA_SILVA_BOWTIE2.txt


		mkdir ./RNASPADES/FINAL_FILES_RNASPADES
		mv ./RNASPADES/*final.txt ./RNASPADES/FINAL_FILES_RNASPADES/
		mv ./RNASPADES/BLAST_output_nt* ./RNASPADES/BLAST_output_SILVA* RNASPADES/CLASSIFICATION/BLAST/

		echo -e "\n======== final RNASPADES output generated ========\n"

	else
		echo -e "\n======== no FINAL merge generated ========\n"

	fi

else
	echo -e "\n======== no RNASPADES output generated ========\n"

fi


echo -e "\n======== RUNNING IDBA_TRAN ========\n"
if [[ $IDBA_TRAN == 'true' ]] ; then
	fq2fa ./SORTMERNA/out/aligned.fq ./SORTMERNA/out/aligned.fa
	idba_tran --num_threads 16 --pre_correction -r ./SORTMERNA/out/aligned.fa -o IDBA_TRAN
	echo -e "\n======== IDBA_TRAN DONE ========\n"

	if [[ $MAP == 'true' ]] ; then

		echo -e "\n======== starting bwa index ========\n"

		bwa index -p bwa_index ./IDBA_TRAN/contig.fa

		echo -e "\n======== bwa index complete. Starting bwa ========\n"

		bwa mem -t 10 bwa_index ./SORTMERNA/out/aligned.fq > bwa_output.sam

		rm bwa_index*


		# Output file (.sam) - edit
		samtools view -F 4 bwa_output.sam > mapped_reads_bwa.sam
		samtools view -f 4 bwa_output.sam > unmapped_reads_bwa.sam
		cat mapped_reads_bwa.sam > mapped_reads_bwa.txt
		cat unmapped_reads_bwa.sam > unmapped_reads_bwa.txt
		cut -f3 mapped_reads_bwa.txt > mapped_column3_reads_bwa.txt
		cut -f3 unmapped_reads_bwa.txt > unmapped_column3_reads_bwa.txt
		sort mapped_column3_reads_bwa.txt | uniq -c > sorted_mapped_column3_reads_bwa.txt
		sort unmapped_column3_reads_bwa.txt | uniq -c > sorted_unmapped_column3_reads_bwa.txt
		column -t sorted_mapped_column3_reads_bwa.txt > aligned_mapped_bwa.txt
		column -t sorted_unmapped_column3_reads_bwa.txt > aligned_unmapped_bwa.txt
		sed 's/  */\t/g' aligned_mapped_bwa.txt > out_mapped_bwa.txt
		sed 's/  */\t/g' aligned_unmapped_bwa.txt > out_unmappped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo -e "\n======== bwa complete. Starting bowtie2 index ========\n"

		bowtie2-build -f ./IDBA_TRAN/contig.fa bowtie_index

		echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"

		bowtie2 -q -x bowtie_index --interleaved ./SORTMERNA/out/aligned.fq -S bowtie2_output.sam

		rm bowtie_index*


		# Output file (.sam) - edit
		samtools view -F 4 bowtie2_output.sam > mapped_reads_bowtie.sam
		samtools view -f 4 bowtie2_output.sam > unmapped_reads_bowtie.sam
		cat mapped_reads_bowtie.sam > mapped_reads_bowtie.txt
		cat unmapped_reads_bowtie.sam > unmapped_reads_bowtie.txt
		cut -f3 mapped_reads_bowtie.txt > mapped_column3_reads_bowtie.txt
		cut -f3 unmapped_reads_bowtie.txt > unmapped_column3_reads_bowtie.txt
		sort mapped_column3_reads_bowtie.txt | uniq -c > sorted_mapped_column3_reads_bowtie.txt
		sort unmapped_column3_reads_bowtie.txt | uniq -c > sorted_unmapped_column3_reads_bowtie.txt
		column -t sorted_mapped_column3_reads_bowtie.txt > aligned_mapped_bowtie.txt
		column -t sorted_unmapped_column3_reads_bowtie.txt > aligned_unmapped_bowtie.txt
		sed 's/  */\t/g' aligned_mapped_bowtie.txt > out_mapped_bowtie.txt
		sed 's/  */\t/g' aligned_unmapped_bowtie.txt > out_unmappped_bowtie.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo -e "\n======== bowtie2 complete ========\n"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ IDBA_TRAN/

	else

		echo "\n======== no BWA and BOWTIE2 output generated ========\n"

	fi

	if [[ $CLASSIFICATION == 'true' ]] ; then

		# BLAST against nt
		echo -e "\n======== RUNNING JUSTBLAST AGAINST NT AND PROCESSING OUTPUT ========\n"
		justblast ./IDBA_TRAN/contig.fa/hdd1/databases/nt_database_feb_2020_indexed/nt --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_nt

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_nt | sort -u -k1,1 > BLAST_output_nt_first_hit.txt

		## Blast option 2: apply bitscore threshold, best bitscore %, similarity cutoff, and an LCA approach
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_nt > BLAST_output_nt_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_nt_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_nt_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_nt_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_nt_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_nt_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST NT AND PROCESSING OUTPUT DONE ========\n"

		# BLAST against SILVA
		echo -e "\n======== RUNNING JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT ========\n"
		justblast ./IDBA_TRAN/contig.fa /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_SILVA

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_SILVA | sort -u -k1,1 > BLAST_output_SILVA_first_hit.txt

		## Blast option 2: apply bitscore threshold and best bitscore %
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_SILVA > BLAST_output_SILVA_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_SILVA_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_SILVA_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_SILVA_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_SILVA_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_SILVA_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT DONE ========\n"

		mkdir IDBA_TRAN/BLAST
		mv BLAST_output_nt* BLAST_output_SILVA* IDBA_TRAN/BLAST/

		echo -e "\n======== RUNNING CREST ========\n"
		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./IDBA_TRAN/contig.fa -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		mv BLAST_output.xml CREST
		mv CREST IDBA_TRAN/

		mkdir IDBA_TRAN/CLASSIFICATION
		mv IDBA_TRAN/BLAST IDBA_TRAN/CREST IDBA_TRAN/CLASSIFICATION/

		echo -e "\n======== CREST DONE ========\n"

	else

		echo -e "\n======== NO BLAST OR CREST OUTPUT GENERATED ========\n"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then

			echo -e "\n======== ADD CONTIG SEQUENCES TO BOWTIE AND BWA OUTPUT AND MERGE THEM WITH BLAST+CREST ========\n"
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./IDBA_TRAN/contig.fa > ./IDBA_TRAN/fasta_to_tabbed.txt
			mergeFilesOnColumn.pl ./IDBA_TRAN/BWA/merge_input_mapped_bwa.txt ./IDBA_TRAN/fasta_to_tabbed.txt 2 1 > ./IDBA_TRAN/merged_original_bwa.txt
			cut -f1,2,4 ./IDBA_TRAN/merged_original_bwa.txt > ./IDBA_TRAN/important_column_3_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./IDBA_TRAN/important_column_3_bwa.txt > ./IDBA_TRAN/final_order_bwa.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA_TRAN/final_bwa_merge_ready.txt && cat ./IDBA_TRAN/final_order_bwa.txt >> ./IDBA_TRAN/final_bwa_merge_ready.txt

			rm ./IDBA_TRAN/merged_original_bwa.txt ./IDBA_TRAN/important_column_3_bwa.txt ./IDBA_TRAN/final_order_bwa.txt

			mergeFilesOnColumn.pl ./IDBA_TRAN/BOWTIE2/merge_input_mapped_bowtie2.txt ./IDBA_TRAN/fasta_to_tabbed.txt 2 1 > ./IDBA_TRAN/merged_original_bowtie2.txt
			cut -f1,2,4 ./IDBA_TRAN/merged_original_bowtie2.txt > ./IDBA_TRAN/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./IDBA_TRAN/important_column_3_bowtie2.txt > ./IDBA_TRAN/final_order_bowtie2.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA_TRAN/final_bowtie2_merge_ready.txt && cat ./IDBA_TRAN/final_order_bowtie2.txt >> ./IDBA_TRAN/final_bowtie2_merge_ready.txt

			rm ./IDBA_TRAN/fasta_to_tabbed.txt ./IDBA_TRAN/merged_original_bowtie2.txt ./IDBA_TRAN/important_column_3_bowtie2.txt ./IDBA_TRAN/final_order_bowtie2.txt

			echo -e "\n======== contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"


		else
			echo -e "\n======== MERGE BOWIE+BWA OUTPUT WITH BLAST+CREST ========\n"
			# Prepares to merge with BLAST or CREST - does not include assembly sequence
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./IDBA_TRAN/BWA/merge_input_mapped_bwa.txt > ./IDBA_TRAN/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./IDBA_TRAN/temp.txt > ./IDBA_TRAN/temp2.txt
			sed '1d' ./IDBA_TRAN/temp2.txt > ./IDBA_TRAN/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA_TRAN/final_bwa_merge_ready.txt && cat ./IDBA_TRAN/temp3.txt >> ./IDBA_TRAN/final_bwa_merge_ready.txt
 			rm ./IDBA_TRAN/temp*.txt

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./IDBA_TRAN/BOWTIE2/merge_input_mapped_bowtie2.txt > ./IDBA_TRAN/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./IDBA_TRAN/temp.txt > ./IDBA_TRAN/temp2.txt
			sed '1d' ./IDBA_TRAN/temp2.txt > ./IDBA_TRAN/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA_TRAN/final_bowtie2_merge_ready.txt && cat ./IDBA_TRAN/temp3.txt >> ./IDBA_TRAN/final_bowtie2_merge_ready.txt
			rm ./IDBA_TRAN/temp*.txt

			echo -e "\n======== no contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"
		fi

		# Assign BLAST taxonomy - nt
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST NT OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_nt_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_first_hit_with_taxonomy.txt ./IDBA_TRAN/
		sed '1d' ./IDBA_TRAN/BLAST_output_nt_first_hit_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt
		mv ./RNASPADES/BLAST_output_nt_first_hit_with_taxonomy.txt ./RNASPADES/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./IDBA_TRAN/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/BLAST_output_nt_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./IDBA_TRAN/BLAST_output_nt_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt
		echo

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_bitscore_filtered_with_taxonomy.txt ./IDBA_TRAN/CLASSIFICATION/BLAST/
		sed '1d' ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./IDBA_TRAN/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt > ./IDBA_TRAN/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./IDBA_TRAN/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./IDBA_TRAN/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./IDBA_TRAN/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./IDBA_TRAN/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./IDBA_TRAN/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./IDBA_TRAN/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./IDBA_TRAN/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./IDBA_TRAN/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./IDBA_TRAN/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./IDBA_TRAN/BLAST_output_nt_similarity_filtered.txt
		awk '($3 < 80)' ./IDBA_TRAN/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./IDBA_TRAN/BLAST_output_nt_similarity_filtered.txt

		### Immitating BASTA
		touch ./IDBA_TRAN/IDBA_TRAN_blast_nt_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./IDBA_TRAN/BLAST_output_nt_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./IDBA_TRAN/BLAST_output_nt_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./IDBA_TRAN/IDBA_TRAN_blast_nt_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/IDBA_TRAN_blast_nt_LCA_taxonomy_merge.txt && cat ./IDBA_TRAN/IDBA_TRAN_blast_nt_LCA_taxonomy_merge_noheader.txt >> ./IDBA_TRAN/IDBA_TRAN_blast_nt_LCA_taxonomy_merge.txt


		# Assign BLAST taxonomy - SILVA
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST SILVA OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_SILVA_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_first_hit_with_taxonomy.txt ./IDBA_TRAN/
		sed '1d' ./IDBA_TRAN/BLAST_output_SILVA_first_hit_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt
		mv ./IDBA_TRAN/BLAST_output_SILVA_first_hit_with_taxonomy.txt ./IDBA_TRAN/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./IDBA_TRAN/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./IDBA_TRAN/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_first_hit_with_taxonomy.txt ./IDBA_TRAN/CLASSIFICATION/BLAST/
		sed '1d' ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./IDBA_TRAN/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt > ./IDBA_TRAN/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./IDBA_TRAN/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./IDBA_TRAN/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./IDBA_TRAN/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./IDBA_TRAN/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./IDBA_TRAN/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./IDBA_TRAN/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./IDBA_TRAN/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./IDBA_TRAN/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./IDBA_TRAN/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./IDBA_TRAN/BLAST_output_SILVA_similarity_filtered.txt
		awk '($3 < 80)' ./IDBA_TRAN/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./IDBA_TRAN/BLAST_output_SILVA_similarity_filtered.txt

		### Immitating BASTA
		touch ./IDBA_TRAN/IDBA_TRAN_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./IDBA_TRAN/BLAST_output_SILVA_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./IDBA_TRAN/BLAST_output_SILVA_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./IDBA_TRAN/IDBA_TRAN_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/IDBA_TRAN_blast_SILVA_LCA_taxonomy_merge.txt && cat ./IDBA_TRAN/IDBA_TRAN_blast_SILVA_LCA_taxonomy_merge_noheader.txt >> ./IDBA_TRAN/IDBA_TRAN_blast_SILVA_LCA_taxonomy_merge.txt

		echo -e "\n======== taxonomy has been assigned to BLAST files ========"
		echo -e "======== BLAST files are ready ========\n"


		# Prepare CREST file - remove 1st/3rd column, add a header
		echo -e "\n======== ADD NCBI TAXONOMY TO CREST OUTPUT ========\n"
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA_TRAN/CLASSIFICATION/CREST/otus.csv ./IDBA_TRAN/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./IDBA_TRAN/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA_TRAN/CREST_seperated.txt
		cut -f2,4 ./IDBA_TRAN/CREST_seperated.txt > ./IDBA_TRAN/CREST_header.txt
		sed '1d' ./IDBA_TRAN/CREST_header.txt > ./IDBA_TRAN/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./IDBA_TRAN/
		sed '1d' ./IDBA_TRAN/CREST_tax_ready_with_taxonomy.txt > ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/CREST_merge.txt && cat ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA_TRAN/CREST_merge.txt

		rm ./IDBA_TRAN/CREST_seperated.txt ./IDBA_TRAN/CREST_header.txt ./IDBA_TRAN/CREST_tax_ready.txt ./IDBA_TRAN/CREST_tax_ready_with_taxonomy.txt ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt

		echo -e "\n======== CREST file is ready ========\n"

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_first_hit_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_first_hit_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_first_hit_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_first_hit_nt_BWA.txt

		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/IDBA_TRAN_blast_SILVA_LCA_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_LCA_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/IDBA_TRAN_blast_nt_LCA_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_LCA_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/IDBA_TRAN_blast_SILVA_LCA_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_LCA_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./IDBA_TRAN/IDBA_TRAN_blast_nt_LCA_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_LCA_nt_BWA.txt


		# Move all files - easy to find!!!
		mkdir ./IDBA_TRAN/MERGE_FILES
		mv ./IDBA_TRAN/*merge*.txt ./IDBA_TRAN/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./IDBA_TRAN/CREST_BWA.txt > ./IDBA_TRAN/new.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt && sed 's/length_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/cov_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt
		sed '1d' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/CREST_BWA_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/CREST_BWA_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/CREST_BWA.txt

		sed 's/_/\t/2' ./IDBA_TRAN/CREST_BOWTIE2.txt > ./IDBA_TRAN/new.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt && sed 's/length_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/cov_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt
		sed '1d' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/CREST_BOWTIE2_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/CREST_BOWTIE2_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./IDBA_TRAN/BLAST_first_hit_nt_BWA.txt > ./IDBA_TRAN/new.txt
		sed 's/_/\t/2' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/length_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt && sed 's/cov_//g' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_first_hit_BWA_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/BLAST_nt_first_hit_BWA_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_first_hit_nt_BWA.txt

		sed '1d' ./IDBA_TRAN/BLAST_first_hit_nt_BOWTIE2.txt > ./IDBA_TRAN/new.txt
		sed 's/_/\t/2' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/length_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt && sed 's/cov_//g' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_first_hit_BOWTIE2_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/BLAST_nt_first_hit_BOWTIE2_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_first_hit_nt_BOWTIE2.txt

		sed '1d' ./IDBA_TRAN/BLAST_LCA_nt_BWA.txt > ./IDBA_TRAN/new.txt
		sed 's/_/\t/2' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/length_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt && sed 's/cov_//g' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_LCA_BWA_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/BLAST_nt_LCA_BWA_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_LCA_nt_BWA.txt

		sed '1d' ./IDBA_TRAN/BLAST_LCA_nt_BOWTIE2.txt > ./IDBA_TRAN/new.txt
		sed 's/_/\t/2' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/length_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt && sed 's/cov_//g' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_LCA_BOWTIE2_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/BLAST_nt_LCA_BOWTIE2_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_LCA_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./IDBA_TRAN/BLAST_first_hit_SILVA_BWA.txt > ./IDBA_TRAN/new.txt
		sed 's/_/\t/2' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/length_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt && sed 's/cov_//g' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_first_hit_BWA_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/BLAST_SILVA_first_hit_BWA_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_first_hit_SILVA_BWA.txt

		sed '1d' ./IDBA_TRAN/BLAST_first_hit_SILVA_BOWTIE2.txt > ./IDBA_TRAN/new.txt
		sed 's/_/\t/2' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/length_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt && sed 's/cov_//g' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_first_hit_BOWTIE2_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/BLAST_SILVA_first_hit_BOWTIE2_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_first_hit_SILVA_BOWTIE2.txt

		sed '1d' ./IDBA_TRAN/BLAST_LCA_SILVA_BWA.txt > ./IDBA_TRAN/new.txt
		sed 's/_/\t/2' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/length_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt && sed 's/cov_//g' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_LCA_BWA_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/BLAST_SILVA_LCA_BWA_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_LCA_SILVA_BWA.txt

		sed '1d' ./IDBA_TRAN/BLAST_LCA_SILVA_BOWTIE2.txt > ./IDBA_TRAN/new.txt
		sed 's/_/\t/2' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt &&  sed 's/_/\t/3' ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
		sed 's/NODE_//g' ./IDBA_TRAN/new3.txt > ./IDBA_TRAN/new4.txt && sed 's/length_//g' ./IDBA_TRAN/new4.txt > ./IDBA_TRAN/new5.txt && sed 's/cov_//g' ./IDBA_TRAN/new5.txt > ./IDBA_TRAN/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_LCA_BOWTIE2_final.txt && cat ./IDBA_TRAN/new6.txt >> ./IDBA_TRAN/BLAST_SILVA_LCA_BOWTIE2_final.txt
		rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_LCA_SILVA_BOWTIE2.txt


		mkdir ./IDBA_TRAN/FINAL_FILES_IDBA_TRAN
		mv ./IDBA_TRAN/*final.txt ./IDBA_TRAN/FINAL_FILES_IDBA_TRAN/
		mv ./IDBA_TRAN/BLAST_output_nt* ./IDBA_TRAN/BLAST_output_SILVA* IDBA_TRAN/CLASSIFICATION/BLAST/

		echo -e "\n======== final IDBA_TRAN output generated ========\n"

	else
		echo -e "\n======== no FINAL merge generated ========\n"

	fi

else
	echo -e "\n======== no IDBA_TRAN output generated ========\n"

fi


echo -e "\n======== RUNNING TRINITY ========\n"
if [[ $TRINITY == 'true' ]] ; then
	deinterleave_fastq_reads.sh < ./SORTMERNA/out/aligned.fq ./SORTMERNA/out/aligned_1.fq ./SORTMERNA/out/aligned_2.fq
	Trinity --seqType fq --max_memory 64G --left ./SORTMERNA/out/aligned_1.fq --right ./SORTMERNA/out/aligned_2.fq --CPU 16 --output TRINITY
	cat ./TRINITY/Trinity.fasta | sed 's/ len/_len/g' > ./TRINITY/Trinity_with_length.fasta
	echo -e "\n======== TRINITY DONE ========\n"

	if [[ $MAP == 'true' ]] ; then

		echo -e "\n======== starting bwa index ========\n"

		bwa index -p bwa_index ./TRINITY/Trinity_with_length.fasta

		echo -e "\n======== bwa index complete. Starting bwa ========\n"

		bwa mem -t 10 bwa_index ./SORTMERNA/out/aligned.fq > bwa_output.sam

		rm bwa_index*


		# Output file (.sam) - edit
		samtools view -F 4 bwa_output.sam > mapped_reads_bwa.sam
		samtools view -f 4 bwa_output.sam > unmapped_reads_bwa.sam
		cat mapped_reads_bwa.sam > mapped_reads_bwa.txt
		cat unmapped_reads_bwa.sam > unmapped_reads_bwa.txt
		cut -f3 mapped_reads_bwa.txt > mapped_column3_reads_bwa.txt
		cut -f3 unmapped_reads_bwa.txt > unmapped_column3_reads_bwa.txt
		sort mapped_column3_reads_bwa.txt | uniq -c > sorted_mapped_column3_reads_bwa.txt
		sort unmapped_column3_reads_bwa.txt | uniq -c > sorted_unmapped_column3_reads_bwa.txt
		column -t sorted_mapped_column3_reads_bwa.txt > aligned_mapped_bwa.txt
		column -t sorted_unmapped_column3_reads_bwa.txt > aligned_unmapped_bwa.txt
		sed 's/  */\t/g' aligned_mapped_bwa.txt > out_mapped_bwa.txt
		sed 's/  */\t/g' aligned_unmapped_bwa.txt > out_unmappped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo -e "\n======== bwa complete. Starting bowtie2 index ========\n"

		bowtie2-build -f ./TRINITY/Trinity_with_length.fasta bowtie_index

		echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"

		bowtie2 -q -x bowtie_index --interleaved ./SORTMERNA/out/aligned.fq -S bowtie2_output.sam

		rm bowtie_index*


		# Output file (.sam) - edit
		samtools view -F 4 bowtie2_output.sam > mapped_reads_bowtie.sam
		samtools view -f 4 bowtie2_output.sam > unmapped_reads_bowtie.sam
		cat mapped_reads_bowtie.sam > mapped_reads_bowtie.txt
		cat unmapped_reads_bowtie.sam > unmapped_reads_bowtie.txt
		cut -f3 mapped_reads_bowtie.txt > mapped_column3_reads_bowtie.txt
		cut -f3 unmapped_reads_bowtie.txt > unmapped_column3_reads_bowtie.txt
		sort mapped_column3_reads_bowtie.txt | uniq -c > sorted_mapped_column3_reads_bowtie.txt
		sort unmapped_column3_reads_bowtie.txt | uniq -c > sorted_unmapped_column3_reads_bowtie.txt
		column -t sorted_mapped_column3_reads_bowtie.txt > aligned_mapped_bowtie.txt
		column -t sorted_unmapped_column3_reads_bowtie.txt > aligned_unmapped_bowtie.txt
		sed 's/  */\t/g' aligned_mapped_bowtie.txt > out_mapped_bowtie.txt
		sed 's/  */\t/g' aligned_unmapped_bowtie.txt > out_unmappped_bowtie.txt
		echo -e "counts\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "counts\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo -e "\n======== bowtie2 complete ========\n"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ TRINITY/

	else

		echo "\n======== no BWA and BOWTIE2 output generated ========\n"

	fi

	if [[ $CLASSIFICATION == 'true' ]] ; then

		# BLAST against nt
		echo -e "\n======== RUNNING JUSTBLAST AGAINST NT AND PROCESSING OUTPUT ========\n"
		justblast ./TRINITY/Trinity_with_length.fasta /hdd1/databases/nt_database_feb_2020_indexed/nt --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_nt

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_nt | sort -u -k1,1 > BLAST_output_nt_first_hit.txt

		## Blast option 2: apply bitscore threshold, best bitscore %, similarity cutoff, and an LCA approach
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_nt > BLAST_output_nt_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_nt_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_nt_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_nt_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_nt_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_nt_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST NT AND PROCESSING OUTPUT DONE ========\n"

		# BLAST against SILVA
		echo -e "\n======== RUNNING JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT ========\n"
		justblast ./TRINITY/Trinity_with_length.fasta /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_SILVA

		## BLAST option 1: just keep first hit per contig
		sort -k1,1n -k12,12nr BLAST_output_SILVA | sort -u -k1,1 > BLAST_output_SILVA_first_hit.txt

		## Blast option 2: apply bitscore threshold and best bitscore %
		### Remove all hits that are below alignment length 100 (based on BASTA) and bitscore 155 (based on CREST)
		awk '($4 >= 100 && $12 >= 155 )' BLAST_output_SILVA > BLAST_output_SILVA_bitscore_threshold.txt

		### Just keep hits within first 2% of best bitscore for each contig
		touch BLAST_output_SILVA_bitscore_filtered.txt
		sort -u -k1,1 BLAST_output_SILVA_bitscore_threshold.txt | cut -f 1 | while read contig ; do
		  bitscore=$(grep $contig BLAST_output_SILVA_bitscore_threshold.txt | sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f 12)
		  echo "Processing contig $contig with bitscore $bitscore"
		  bitscore_threshold=$(awk "BEGIN {print $bitscore - $bitscore * 2/100}") #change 2/100 to % you like for bitscore inclusion
		  grep $contig BLAST_output_SILVA_bitscore_threshold.txt | awk -v x=$bitscore_threshold '($12 >= x)' >> BLAST_output_SILVA_bitscore_filtered.txt
		done

		echo -e "\n======== JUSTBLAST AGAINST SILVA AND PROCESSING OUTPUT DONE ========\n"

		mkdir TRINITY/BLAST
		mv BLAST_output_nt* BLAST_output_SILVA* TRINITY/BLAST/

		echo -e "\n======== RUNNING CREST ========\n"
		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./TRINITY/Trinity_with_length.fasta -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		mv BLAST_output.xml CREST
		mv CREST TRINITY/

		mkdir TRINITY/CLASSIFICATION
		mv TRINITY/BLAST TRINITY/CREST TRINITY/CLASSIFICATION/

		echo -e "\n======== CREST DONE ========\n"

	else

		echo -e "\n======== NO BLAST OR CREST OUTPUT GENERATED ========\n"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then

			echo -e "\n======== ADD CONTIG SEQUENCES TO BOWTIE AND BWA OUTPUT AND MERGE THEM WITH BLAST+CREST ========\n"
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./TRINITY/Trinity_with_length.fasta > ./TRINITY/fasta_to_tabbed.txt
			mergeFilesOnColumn.pl ./TRINITY/BWA/merge_input_mapped_bwa.txt ./TRINITY/fasta_to_tabbed.txt 2 1 > ./TRINITY/merged_original_bwa.txt
			cut -f1,2,4 ./TRINITY/merged_original_bwa.txt > ./TRINITY/important_column_3_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./TRINITY/important_column_3_bwa.txt > ./TRINITY/final_order_bwa.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./TRINITY/final_bwa_merge_ready.txt && cat ./TRINITY/final_order_bwa.txt >> ./TRINITY/final_bwa_merge_ready.txt

			rm ./TRINITY/merged_original_bwa.txt ./TRINITY/important_column_3_bwa.txt ./TRINITY/final_order_bwa.txt

			mergeFilesOnColumn.pl ./TRINITY/BOWTIE2/merge_input_mapped_bowtie2.txt ./TRINITY/fasta_to_tabbed.txt 2 1 > ./TRINITY/merged_original_bowtie2.txt
			cut -f1,2,4 ./TRINITY/merged_original_bowtie2.txt > ./TRINITY/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./TRINITY/important_column_3_bowtie2.txt > ./TRINITY/final_order_bowtie2.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./TRINITY/final_bowtie2_merge_ready.txt && cat ./TRINITY/final_order_bowtie2.txt >> ./TRINITY/final_bowtie2_merge_ready.txt

			rm ./TRINITY/fasta_to_tabbed.txt ./TRINITY/merged_original_bowtie2.txt ./TRINITY/important_column_3_bowtie2.txt ./TRINITY/final_order_bowtie2.txt

			echo -e "\n======== contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"


		else
			echo -e "\n======== MERGE BOWIE+BWA OUTPUT WITH BLAST+CREST ========\n"
			# Prepares to merge with BLAST or CREST - does not include assembly sequence
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./TRINITY/BWA/merge_input_mapped_bwa.txt > ./TRINITY/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./TRINITY/temp.txt > ./TRINITY/temp2.txt
			sed '1d' ./TRINITY/temp2.txt > ./TRINITY/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./TRINITY/final_bwa_merge_ready.txt && cat ./TRINITY/temp3.txt >> ./TRINITY/final_bwa_merge_ready.txt
 			rm ./TRINITY/temp*.txt

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./TRINITY/BOWTIE2/merge_input_mapped_bowtie2.txt > ./TRINITY/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./TRINITY/temp.txt > ./TRINITY/temp2.txt
			sed '1d' ./TRINITY/temp2.txt > ./TRINITY/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./TRINITY/final_bowtie2_merge_ready.txt && cat ./TRINITY/temp3.txt >> ./TRINITY/final_bowtie2_merge_ready.txt
			rm ./TRINITY/temp*.txt

			echo -e "\n======== no contig sequences were included ========"
			echo -e "======== BWA and BOWTIE2 files are ready ========\n"
		fi

		# Assign BLAST taxonomy - nt
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST NT OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./TRINITY/CLASSIFICATION/BLAST/BLAST_output_nt_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_first_hit_with_taxonomy.txt ./TRINITY/
		sed '1d' ./TRINITY/BLAST_output_nt_first_hit_with_taxonomy.txt > ./TRINITY/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt
		mv ./TRINITY/BLAST_output_nt_first_hit_with_taxonomy.txt ./TRINITY/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./TRINITY/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRINITY/BLAST_output_nt_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./TRINITY/BLAST_output_nt_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt
		echo

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./TRINITY/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_nt_bitscore_filtered_with_taxonomy.txt ./TRINITY/CLASSIFICATION/BLAST/
		sed '1d' ./TRINITY/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered_with_taxonomy.txt > ./TRINITY/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./TRINITY/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt > ./TRINITY/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./TRINITY/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./TRINITY/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./TRINITY/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./TRINITY/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./TRINITY/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./TRINITY/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./TRINITY/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./TRINITY/BLAST_output_nt_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./TRINITY/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./TRINITY/BLAST_output_nt_similarity_filtered.txt
		awk '($3 < 80)' ./TRINITY/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./TRINITY/BLAST_output_nt_similarity_filtered.txt

		### Immitating BASTA
		touch ./TRINITY/TRINITY_blast_nt_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./TRINITY/BLAST_output_nt_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./TRINITY/BLAST_output_nt_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./TRINITY/TRINITY_blast_nt_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRINITY/TRINITY_blast_nt_LCA_taxonomy_merge.txt && cat ./TRINITY/TRINITY_blast_nt_LCA_taxonomy_merge_noheader.txt >> ./TRINITY/TRINITY_blast_nt_LCA_taxonomy_merge.txt


		# Assign BLAST taxonomy - SILVA
		echo -e "\n======== ADD NCBI TAXONOMY TO BLAST SILVA OUTPUTS ========\n"

		## BLAST option 1: assign taxonomy to first hit per contig
		assign_taxonomy_NCBI_staxids.sh -b ./TRINITY/CLASSIFICATION/BLAST/BLAST_output_SILVA_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_first_hit_with_taxonomy.txt ./TRINITY/
		sed '1d' ./TRINITY/BLAST_output_SILVA_first_hit_with_taxonomy.txt > ./TRINITY/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt
		mv ./TRINITY/BLAST_output_SILVA_first_hit_with_taxonomy.txt ./TRINITY/CLASSIFICATION/BLAST/
		cut -f 1-13,15-27 ./TRINITY/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRINITY/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./TRINITY/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt
		rm tmp_cut.txt

		## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
		assign_taxonomy_NCBI_staxids.sh -b ./TRINITY/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt ./TRINITY/CLASSIFICATION/BLAST/
		sed '1d' ./TRINITY/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt > ./TRINITY/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt

		### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
		awk '($3 >= 99)' ./TRINITY/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt > ./TRINITY/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 99) && ($3 >= 97))' ./TRINITY/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./TRINITY/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 97) && ($3 >= 95))' ./TRINITY/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./TRINITY/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 95) && ($3 >= 90))' ./TRINITY/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./TRINITY/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 90) && ($3 >= 85))' ./TRINITY/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./TRINITY/BLAST_output_SILVA_similarity_filtered.txt
		awk '(($3 < 85) && ($3 >= 80))' ./TRINITY/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./TRINITY/BLAST_output_SILVA_similarity_filtered.txt
		awk '($3 < 80)' ./TRINITY/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./TRINITY/BLAST_output_SILVA_similarity_filtered.txt

		### Immitating BASTA
		touch ./TRINITY/TRINITY_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		sort -u -k1,1 ./TRINITY/BLAST_output_SILVA_similarity_filtered.txt | cut -f 1 | while read contig ; do
		  taxonomy_contig=$(echo $contig)
		  taxonomy=''
		  for i in {16..27}
		  do
		    rank_tax=$(grep $contig ./TRINITY/BLAST_output_SILVA_similarity_filtered.txt | cut -f $i | uniq)
		    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
		      taxonomy=$(echo "${taxonomy}---${rank_tax}")
		    else
		      taxonomy=$(echo "${taxonomy}---NA")
		    fi
		  done
		  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./TRINITY/TRINITY_blast_SILVA_LCA_taxonomy_merge_noheader.txt
		done
		echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRINITY/TRINITY_blast_SILVA_LCA_taxonomy_merge.txt && cat ./TRINITY/TRINITY_blast_SILVA_LCA_taxonomy_merge_noheader.txt >> ./TRINITY/TRINITY_blast_SILVA_LCA_taxonomy_merge.txt

		echo -e "\n======== taxonomy has been assigned to BLAST files ========"
		echo -e "======== BLAST files are ready ========\n"


		# Prepare CREST file - remove 1st/3rd column, add a header
		echo -e "\n======== ADD NCBI TAXONOMY TO CREST OUTPUT ========\n"
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./TRINITY/CLASSIFICATION/CREST/otus.csv ./TRINITY/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./TRINITY/CLASSIFICATION/CREST/CREST_output.txt > ./TRINITY/CREST_seperated.txt
		cut -f2,4 ./TRINITY/CREST_seperated.txt > ./TRINITY/CREST_header.txt
		sed '1d' ./TRINITY/CREST_header.txt > ./TRINITY/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./TRINITY/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./TRINITY/
		sed '1d' ./TRINITY/CREST_tax_ready_with_taxonomy.txt > ./TRINITY/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRINITY/CREST_merge.txt && cat ./TRINITY/CREST_tax_ready_with_taxonomy_noheader.txt >> ./TRINITY/CREST_merge.txt

		rm ./TRINITY/CREST_seperated.txt ./TRINITY/CREST_header.txt ./TRINITY/CREST_tax_ready.txt ./TRINITY/CREST_tax_ready_with_taxonomy.txt ./TRINITY/CREST_tax_ready_with_taxonomy_noheader.txt

		echo -e "\n======== CREST file is ready ========\n"

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./TRINITY/CREST_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/CREST_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/BLAST_first_hit_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/BLAST_first_hit_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/BLAST_first_hit_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/BLAST_first_hit_nt_BWA.txt

		merge_mapped_reads_and_contigs.py ./TRINITY/TRINITY_blast_SILVA_LCA_taxonomy_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/BLAST_LCA_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/TRINITY_blast_nt_LCA_taxonomy_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/BLAST_LCA_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./TRINITY/TRINITY_blast_SILVA_LCA_taxonomy_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/BLAST_LCA_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/TRINITY_blast_nt_LCA_taxonomy_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/BLAST_LCA_nt_BWA.txt


		# Move all files - easy to find!!!
		mkdir ./TRINITY/MERGE_FILES
		mv ./TRINITY/*merge*.txt ./TRINITY/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./TRINITY/CREST_BWA.txt > ./TRINITY/new.txt &&  sed 's/_/\t/3' ./TRINITY/new.txt > ./TRINITY/new2.txt
		sed 's/NODE_//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt && sed 's/length_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/cov_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt
		sed '1d' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/CREST_BWA_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/CREST_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/CREST_BWA.txt

		sed 's/_/\t/2' ./TRINITY/CREST_BOWTIE2.txt > ./TRINITY/new.txt &&  sed 's/_/\t/3' ./TRINITY/new.txt > ./TRINITY/new2.txt
		sed 's/NODE_//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt && sed 's/length_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/cov_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt
		sed '1d' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/CREST_BOWTIE2_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/CREST_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./TRINITY/BLAST_first_hit_nt_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/2' ./TRINITY/new.txt > ./TRINITY/new2.txt &&  sed 's/_/\t/3' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/NODE_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/length_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt && sed 's/cov_//g' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_nt_first_hit_BWA_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/BLAST_nt_first_hit_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_first_hit_nt_BWA.txt

		sed '1d' ./TRINITY/BLAST_first_hit_nt_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/2' ./TRINITY/new.txt > ./TRINITY/new2.txt &&  sed 's/_/\t/3' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/NODE_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/length_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt && sed 's/cov_//g' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_nt_first_hit_BOWTIE2_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/BLAST_nt_first_hit_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_first_hit_nt_BOWTIE2.txt

		sed '1d' ./TRINITY/BLAST_LCA_nt_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/2' ./TRINITY/new.txt > ./TRINITY/new2.txt &&  sed 's/_/\t/3' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/NODE_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/length_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt && sed 's/cov_//g' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_nt_LCA_BWA_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/BLAST_nt_LCA_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_LCA_nt_BWA.txt

		sed '1d' ./TRINITY/BLAST_LCA_nt_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/2' ./TRINITY/new.txt > ./TRINITY/new2.txt &&  sed 's/_/\t/3' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/NODE_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/length_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt && sed 's/cov_//g' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_nt_LCA_BOWTIE2_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/BLAST_nt_LCA_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_LCA_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./TRINITY/BLAST_first_hit_SILVA_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/2' ./TRINITY/new.txt > ./TRINITY/new2.txt &&  sed 's/_/\t/3' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/NODE_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/length_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt && sed 's/cov_//g' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_SILVA_first_hit_BWA_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/BLAST_SILVA_first_hit_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_first_hit_SILVA_BWA.txt

		sed '1d' ./TRINITY/BLAST_first_hit_SILVA_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/2' ./TRINITY/new.txt > ./TRINITY/new2.txt &&  sed 's/_/\t/3' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/NODE_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/length_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt && sed 's/cov_//g' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_SILVA_first_hit_BOWTIE2_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/BLAST_SILVA_first_hit_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_first_hit_SILVA_BOWTIE2.txt

		sed '1d' ./TRINITY/BLAST_LCA_SILVA_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/2' ./TRINITY/new.txt > ./TRINITY/new2.txt &&  sed 's/_/\t/3' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/NODE_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/length_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt && sed 's/cov_//g' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_SILVA_LCA_BWA_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/BLAST_SILVA_LCA_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_LCA_SILVA_BWA.txt

		sed '1d' ./TRINITY/BLAST_LCA_SILVA_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/2' ./TRINITY/new.txt > ./TRINITY/new2.txt &&  sed 's/_/\t/3' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/NODE_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt && sed 's/length_//g' ./TRINITY/new4.txt > ./TRINITY/new5.txt && sed 's/cov_//g' ./TRINITY/new5.txt > ./TRINITY/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_SILVA_LCA_BOWTIE2_final.txt && cat ./TRINITY/new6.txt >> ./TRINITY/BLAST_SILVA_LCA_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_LCA_SILVA_BOWTIE2.txt


		mkdir ./TRINITY/FINAL_FILES_TRINITY
		mv ./TRINITY/*final.txt ./TRINITY/FINAL_FILES_TRINITY/
		mv ./TRINITY/BLAST_output_nt* ./TRINITY/BLAST_output_SILVA* TRINITY/CLASSIFICATION/BLAST/

		echo -e "\n======== final TRINITY output generated ========\n"

	else
		echo -e "\n======== no FINAL merge generated ========\n"

	fi

else
	echo -e "\n======== no TRINITY output generated ========\n"

fi

mkdir ASSEMBLERS
mv RNASPADES IDBA_TRAN TRINITY ASSEMBLERS/


######################### Final output directory generation ####################

# Extract all final output files out of PIPELINE_DNA/ASSEMBLERS and
# PIPELINE_RNA/ASSEMBLERS folders into one folder
cd ..

for assembler_folder in PIPELINE_DNA/ASSEMBLERS/*
do
  for file in $assembler_folder/FINAL_FILES_${assembler_folder##*/}/*
  do
    cp $file ${assembler_folder##*/}_${file##*/}
  done
done

for assembler_folder in PIPELINE_RNA/ASSEMBLERS/*
do
  for file in $assembler_folder/FINAL_FILES_${assembler_folder##*/}/*
  do
    cp $file ${assembler_folder##*/}_${file##*/}
  done
done

mkdir final_files
mv *.txt final_files


# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"


# Create log
) 2>&1 | tee PIPELINE_ASSEMBLERS_LOG.txt