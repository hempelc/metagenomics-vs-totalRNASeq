#!/bin/bash

# As for now, SortMeRNA is only carried out for RNA assemblers
# Need to have: assign_NCBI_staxids_to_CREST_v3.py, assign_taxonomy_NCBI_staxids.sh,
# deinterleave_fastq_reads.sh, LookupTaxonDetails3.py,
# merge_mapped_reads_and_contigs.py, mergeFilesOnColumn.pl and fasta_to_tab
# in your PATH, and justblast installed (https://pypi.org/project/justblast/)
# In order for "assign_taxonomy_NCBI_staxids.sh" to work - you MUST have
# .etetoolkit/taxa.sqlite in your HOME directory - check the ete3 toolkit
# to see how that's set up

cmd="$0 $@" # Make variable containing full used command to print command in logfile

usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq> [aDRSMmUrtTBCfsh]

Usage:
	-1 Forward reads trimmed - must state full path from root to the file
	-2 Reverse reads trimmed - must state full path from root to the file
	-a Main flag to indicate that all following flags should be used
	-D Flag to use completely pipeline only for DNA assemblers
	-R Flag to use completely pipeline only for RNA assemblers
	-S Flag to generate TRANSABYSS assembly output
	-M Flag to generate METATRANSABYSS assembly output
	-m Flag to generate MEGAHIT assembly output
	-U Flag to generate IDBA-UD assembly output
	-r Flag to generate RNATRANSABYSS assembly output
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
TRANSABYSS='false'
METATRANSABYSS='false'
MEGAHIT='false'
IDBA_UD='false'
RNATRANSABYSS='false'
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
		a) TRANSABYSS='true'
			 METATRANSABYSS='true'
			 MEGAHIT='true'
			 IDBA_UD='true'
			 RNATRANSABYSS='true'
	 		 IDBA_TRAN='true'
	 		 TRINITY='true'
			 MAP='true'
			 CLASSIFICATION='true'
			 FINAL='true'
			 READS='true' ;;
		D) TRANSABYSS='true'
			 METATRANSABYSS='true'
			 MEGAHIT='true'
			 IDBA_UD='true'
			 MAP='true'
			 CLASSIFICATION='true'
			 FINAL='true'
			 READS='true' ;;
		R) RNATRANSABYSS='true'
		   IDBA_TRAN='true'
		   TRINITY='true'
		   MAP='true'
		   CLASSIFICATION='true'
		   FINAL='true'
		   READS='true' ;;
		S) TRANSABYSS='true' ;;
		M) METATRANSABYSS='true' ;;
		m) MEGAHIT='true' ;;
		U) IDBA_UD='true' ;;
		r) RNATRANSABYSS='true' ;;
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

##################### Write time, options etc. to output ######################

# Make open bracket to later tell script to write everything that follows into a logfile
(

# Define starting time of script for total runtime calculation
start=$(date +%s)

echo -e "\n======== RUNNING TRANSABYSS ========\n"
if [[ $TRANSABYSS == 'true' ]] ; then
	transabyss --pe $Forward_read_trimmed $Reverse_read_trimmed --threads 16 --outdir TRANSABYSS
	sed 's/ /_/g' TRANSABYSS/transabyss-final.fa > TRANSABYSS/transabyss-final_edited.fa
	echo -e "\n======== TRANSABYSS DONE ========\n"

		if [[ $MAP == 'true' ]] ; then

			echo -e "\n======== starting bwa index ========\n"

			bwa index -p bwa_index TRANSABYSS/transabyss-final_edited.fa

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

			bowtie2-build -f TRANSABYSS/transabyss-final_edited.fa bowtie_index

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

			mv BWA/ BOWTIE2/ TRANSABYSS/

		else

			echo "\n======== no BWA and BOWTIE2 output generated ========\n"

		fi

		if [[ $CLASSIFICATION == 'true' ]] ; then

			# BLAST against nt

			echo -e "\n======== RUNNING JUSTBLAST AGAINST NT AND PROCESSING OUTPUT ========\n"
			justblast TRANSABYSS/transabyss-final_edited.fa /hdd1/databases/nt_database_feb_2020_indexed/nt --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_nt

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
			justblast TRANSABYSS/transabyss-final_edited.fa /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta --cpus 16 --evalue 1e-05 --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" --out_filename BLAST_output_SILVA

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

			mkdir TRANSABYSS/BLAST
			mv BLAST_output_nt* BLAST_output_SILVA* TRANSABYSS/BLAST/

			echo -e "\n======== RUNNING CREST ========\n"
			blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query TRANSABYSS/transabyss-final_edited.fa -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
			classify BLAST_output.xml -o CREST

			mv BLAST_output.xml CREST
			mv CREST TRANSABYSS/

			mkdir TRANSABYSS/CLASSIFICATION
			mv TRANSABYSS/BLAST TRANSABYSS/CREST TRANSABYSS/CLASSIFICATION/

			echo -e "\n======== CREST DONE ========\n"

		else

			echo -e "\n======== NO BLAST OR CREST OUTPUT GENERATED ========\n"


		fi

		if [[ $FINAL == 'true' ]]; then
			if [[ $READS == 'true' ]]; then

				echo -e "\n======== ADD CONTIG SEQUENCES TO BOWTIE AND BWA OUTPUT AND MERGE THEM WITH BLAST+CREST ========\n"
				# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
				fasta_to_tab TRANSABYSS/transabyss-final_edited.fa > ./TRANSABYSS/fasta_to_tabbed.txt
				mergeFilesOnColumn.pl ./TRANSABYSS/BWA/merge_input_mapped_bwa.txt ./TRANSABYSS/fasta_to_tabbed.txt 2 1 > ./TRANSABYSS/merged_original_bwa.txt
				cut -f1,2,4 ./TRANSABYSS/merged_original_bwa.txt > ./TRANSABYSS/important_column_3_bwa.txt
				awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./TRANSABYSS/important_column_3_bwa.txt > ./TRANSABYSS/final_order_bwa.txt
				echo -e "contig_number\tcounts\tassembly_sequence" > ./TRANSABYSS/final_bwa_merge_ready.txt && cat ./TRANSABYSS/final_order_bwa.txt >> ./TRANSABYSS/final_bwa_merge_ready.txt

				rm ./TRANSABYSS/merged_original_bwa.txt ./TRANSABYSS/important_column_3_bwa.txt ./TRANSABYSS/final_order_bwa.txt

				mergeFilesOnColumn.pl ./TRANSABYSS/BOWTIE2/merge_input_mapped_bowtie2.txt ./TRANSABYSS/fasta_to_tabbed.txt 2 1 > ./TRANSABYSS/merged_original_bowtie2.txt
				cut -f1,2,4 ./TRANSABYSS/merged_original_bowtie2.txt > ./TRANSABYSS/important_column_3_bowtie2.txt
				awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./TRANSABYSS/important_column_3_bowtie2.txt > ./TRANSABYSS/final_order_bowtie2.txt
				echo -e "contig_number\tcounts\tassembly_sequence" > ./TRANSABYSS/final_bowtie2_merge_ready.txt && cat ./TRANSABYSS/final_order_bowtie2.txt >> ./TRANSABYSS/final_bowtie2_merge_ready.txt

				rm ./TRANSABYSS/fasta_to_tabbed.txt ./TRANSABYSS/merged_original_bowtie2.txt ./TRANSABYSS/important_column_3_bowtie2.txt ./TRANSABYSS/final_order_bowtie2.txt

				echo -e "\n======== contig sequences were included ========"
				echo -e "======== BWA and BOWTIE2 files are ready ========\n"


			else
				echo -e "\n======== MERGE BOWIE+BWA OUTPUT WITH BLAST+CREST ========\n"
				# Prepares to merge with BLAST or CREST - does not include assembly sequence
				awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./TRANSABYSS/BWA/merge_input_mapped_bwa.txt > ./TRANSABYSS/temp.txt
				sed 's/$/\tno assembly sequence was included/' ./TRANSABYSS/temp.txt > ./TRANSABYSS/temp2.txt
				sed '1d' ./TRANSABYSS/temp2.txt > ./TRANSABYSS/temp3.txt
				echo -e "contig_number\tcounts\tassembly_sequence" > ./TRANSABYSS/final_bwa_merge_ready.txt && cat ./TRANSABYSS/temp3.txt >> ./TRANSABYSS/final_bwa_merge_ready.txt
	 			rm ./TRANSABYSS/temp*.txt

				awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./TRANSABYSS/BOWTIE2/merge_input_mapped_bowtie2.txt > ./TRANSABYSS/temp.txt
				sed 's/$/\tno assembly sequence was included/' ./TRANSABYSS/temp.txt > ./TRANSABYSS/temp2.txt
				sed '1d' ./TRANSABYSS/temp2.txt > ./TRANSABYSS/temp3.txt
				echo -e "contig_number\tcounts\tassembly_sequence" > ./TRANSABYSS/final_bowtie2_merge_ready.txt && cat ./TRANSABYSS/temp3.txt >> ./TRANSABYSS/final_bowtie2_merge_ready.txt
				rm ./TRANSABYSS/temp*.txt

				echo -e "\n======== no contig sequences were included ========"
				echo -e "======== BWA and BOWTIE2 files are ready ========\n"
			fi

			# Assign BLAST taxonomy - nt
			echo -e "\n======== ADD NCBI TAXONOMY TO BLAST NT OUTPUTS ========\n"

			## BLAST option 1: assign taxonomy to first hit per contig
			assign_taxonomy_NCBI_staxids.sh -b ./TRANSABYSS/CLASSIFICATION/BLAST/BLAST_output_nt_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv BLAST_output_nt_first_hit_with_taxonomy.txt ./TRANSABYSS/
			sed '1d' ./TRANSABYSS/BLAST_output_nt_first_hit_with_taxonomy.txt > ./TRANSABYSS/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt
			mv ./TRANSABYSS/BLAST_output_nt_first_hit_with_taxonomy.txt ./TRANSABYSS/CLASSIFICATION/BLAST/
			cut -f 1-13,15-27 ./TRANSABYSS/BLAST_output_nt_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRANSABYSS/BLAST_output_nt_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./TRANSABYSS/BLAST_output_nt_first_hit_with_taxonomy_merge.txt
			rm tmp_cut.txt
			echo

			## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
			assign_taxonomy_NCBI_staxids.sh -b ./TRANSABYSS/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv BLAST_output_nt_bitscore_filtered_with_taxonomy.txt ./TRANSABYSS/CLASSIFICATION/BLAST/
			sed '1d' ./TRANSABYSS/CLASSIFICATION/BLAST/BLAST_output_nt_bitscore_filtered_with_taxonomy.txt > ./TRANSABYSS/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt

			### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
			awk '($3 >= 99)' ./TRANSABYSS/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt > ./TRANSABYSS/BLAST_output_nt_similarity_filtered.txt
			awk '(($3 < 99) && ($3 >= 97))' ./TRANSABYSS/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./TRANSABYSS/BLAST_output_nt_similarity_filtered.txt
			awk '(($3 < 97) && ($3 >= 95))' ./TRANSABYSS/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./TRANSABYSS/BLAST_output_nt_similarity_filtered.txt
			awk '(($3 < 95) && ($3 >= 90))' ./TRANSABYSS/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./TRANSABYSS/BLAST_output_nt_similarity_filtered.txt
			awk '(($3 < 90) && ($3 >= 85))' ./TRANSABYSS/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./TRANSABYSS/BLAST_output_nt_similarity_filtered.txt
			awk '(($3 < 85) && ($3 >= 80))' ./TRANSABYSS/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./TRANSABYSS/BLAST_output_nt_similarity_filtered.txt
			awk '($3 < 80)' ./TRANSABYSS/BLAST_output_nt_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./TRANSABYSS/BLAST_output_nt_similarity_filtered.txt

			### Immitating BASTA
			touch ./TRANSABYSS/TRANSABYSS_blast_nt_LCA_taxonomy_merge_noheader.txt
			sort -u -k1,1 ./TRANSABYSS/BLAST_output_nt_similarity_filtered.txt | cut -f 1 | while read contig ; do
			  taxonomy_contig=$(echo $contig)
			  taxonomy=''
			  for i in {16..27}
			  do
			    rank_tax=$(grep $contig ./TRANSABYSS/BLAST_output_nt_similarity_filtered.txt | cut -f $i | uniq)
			    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
			      taxonomy=$(echo "${taxonomy}---${rank_tax}")
			    else
			      taxonomy=$(echo "${taxonomy}---NA")
			    fi
			  done
			  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./TRANSABYSS/TRANSABYSS_blast_nt_LCA_taxonomy_merge_noheader.txt
			done
			echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRANSABYSS/TRANSABYSS_blast_nt_LCA_taxonomy_merge.txt && cat ./TRANSABYSS/TRANSABYSS_blast_nt_LCA_taxonomy_merge_noheader.txt >> ./TRANSABYSS/TRANSABYSS_blast_nt_LCA_taxonomy_merge.txt


			# Assign BLAST taxonomy - SILVA
			echo -e "\n======== ADD NCBI TAXONOMY TO BLAST SILVA OUTPUTS ========\n"

			## BLAST option 1: assign taxonomy to first hit per contig
			assign_taxonomy_NCBI_staxids.sh -b ./TRANSABYSS/CLASSIFICATION/BLAST/BLAST_output_SILVA_first_hit.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv BLAST_output_SILVA_first_hit_with_taxonomy.txt ./TRANSABYSS/
			sed '1d' ./TRANSABYSS/BLAST_output_SILVA_first_hit_with_taxonomy.txt > ./TRANSABYSS/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt
			mv ./TRANSABYSS/BLAST_output_SILVA_first_hit_with_taxonomy.txt ./TRANSABYSS/CLASSIFICATION/BLAST/
			cut -f 1-13,15-27 ./TRANSABYSS/BLAST_output_SILVA_first_hit_with_taxonomy_noheader.txt > tmp_cut.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRANSABYSS/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt && cat tmp_cut.txt >> ./TRANSABYSS/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt
			rm tmp_cut.txt
			echo

			## Blast option 2: apply similarity cutoff and LCA approach immitating BASTA
			assign_taxonomy_NCBI_staxids.sh -b ./TRANSABYSS/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered.txt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt ./TRANSABYSS/CLASSIFICATION/BLAST/
			sed '1d' ./TRANSABYSS/CLASSIFICATION/BLAST/BLAST_output_SILVA_bitscore_filtered_with_taxonomy.txt > ./TRANSABYSS/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt

			### For similarity cutoff: from CREST: For the species, genus, family, order, class and phylum ranks the respective default cut-offs are 99%, 97%, 95%, 90%, 85% and 80%.
			awk '($3 >= 99)' ./TRANSABYSS/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt > ./TRANSABYSS/BLAST_output_SILVA_similarity_filtered.txt
			awk '(($3 < 99) && ($3 >= 97))' ./TRANSABYSS/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; print}' >> ./TRANSABYSS/BLAST_output_SILVA_similarity_filtered.txt
			awk '(($3 < 97) && ($3 >= 95))' ./TRANSABYSS/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; print}' >> ./TRANSABYSS/BLAST_output_SILVA_similarity_filtered.txt
			awk '(($3 < 95) && ($3 >= 90))' ./TRANSABYSS/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; print}' >> ./TRANSABYSS/BLAST_output_SILVA_similarity_filtered.txt
			awk '(($3 < 90) && ($3 >= 85))' ./TRANSABYSS/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; print}' >> ./TRANSABYSS/BLAST_output_SILVA_similarity_filtered.txt
			awk '(($3 < 85) && ($3 >= 80))' ./TRANSABYSS/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; print}' >> ./TRANSABYSS/BLAST_output_SILVA_similarity_filtered.txt
			awk '($3 < 80)' ./TRANSABYSS/BLAST_output_SILVA_bitscore_filtered_with_taxonomy_noheader.txt | awk '{ FS="\t"; $0=$0; OFS="\t"; $16="NA"; $27="NA"; $26="NA"; $25="NA"; $24="NA"; $23="NA"; $22="NA"; $21="NA"; $20="NA"; $19="NA"; $18="NA"; print}' >> ./TRANSABYSS/BLAST_output_SILVA_similarity_filtered.txt

			### Immitating BASTA
			touch ./TRANSABYSS/TRANSABYSS_blast_SILVA_LCA_taxonomy_merge_noheader.txt
			sort -u -k1,1 ./TRANSABYSS/BLAST_output_SILVA_similarity_filtered.txt | cut -f 1 | while read contig ; do
			  taxonomy_contig=$(echo $contig)
			  taxonomy=''
			  for i in {16..27}
			  do
			    rank_tax=$(grep $contig ./TRANSABYSS/BLAST_output_SILVA_similarity_filtered.txt | cut -f $i | uniq)
			    if [[ $(echo "{$rank_tax}" | wc -l) == 1 ]] ; then
			      taxonomy=$(echo "${taxonomy}---${rank_tax}")
			    else
			      taxonomy=$(echo "${taxonomy}---NA")
			    fi
			  done
			  echo "${taxonomy_contig}${taxonomy}" | sed 's/---/\t/g' >> ./TRANSABYSS/TRANSABYSS_blast_SILVA_LCA_taxonomy_merge_noheader.txt
			done
			echo -e "contig_number\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRANSABYSS/TRANSABYSS_blast_SILVA_LCA_taxonomy_merge.txt && cat ./TRANSABYSS/TRANSABYSS_blast_SILVA_LCA_taxonomy_merge_noheader.txt >> ./TRANSABYSS/TRANSABYSS_blast_SILVA_LCA_taxonomy_merge.txt

			echo -e "\n======== taxonomy has been assigned to BLAST files ========"
			echo -e "======== BLAST files are ready ========\n"


			# Prepare CREST file - remove 1st/3rd column, add a header
			echo -e "\n======== ADD NCBI TAXONOMY TO CREST OUTPUT ========\n"
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./TRANSABYSS/CLASSIFICATION/CREST/otus.csv ./TRANSABYSS/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./TRANSABYSS/CLASSIFICATION/CREST/CREST_output.txt > ./TRANSABYSS/CREST_seperated.txt
			cut -f2,4 ./TRANSABYSS/CREST_seperated.txt > ./TRANSABYSS/CREST_header.txt
			sed '1d' ./TRANSABYSS/CREST_header.txt > ./TRANSABYSS/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./TRANSABYSS/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./TRANSABYSS/
			sed '1d' ./TRANSABYSS/CREST_tax_ready_with_taxonomy.txt > ./TRANSABYSS/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRANSABYSS/CREST_merge.txt && cat ./TRANSABYSS/CREST_tax_ready_with_taxonomy_noheader.txt >> ./TRANSABYSS/CREST_merge.txt

			rm ./TRANSABYSS/CREST_seperated.txt ./TRANSABYSS/CREST_header.txt ./TRANSABYSS/CREST_tax_ready.txt ./TRANSABYSS/CREST_tax_ready_with_taxonomy.txt ./TRANSABYSS/CREST_tax_ready_with_taxonomy_noheader.txt

			echo -e "\n======== CREST file is ready ========\n"

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./TRANSABYSS/CREST_merge.txt ./TRANSABYSS/final_bwa_merge_ready.txt ./TRANSABYSS/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./TRANSABYSS/CREST_merge.txt ./TRANSABYSS/final_bowtie2_merge_ready.txt ./TRANSABYSS/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./TRANSABYSS/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./TRANSABYSS/final_bowtie2_merge_ready.txt ./TRANSABYSS/BLAST_first_hit_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./TRANSABYSS/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./TRANSABYSS/final_bowtie2_merge_ready.txt ./TRANSABYSS/BLAST_first_hit_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./TRANSABYSS/BLAST_output_SILVA_first_hit_with_taxonomy_merge.txt ./TRANSABYSS/final_bwa_merge_ready.txt ./TRANSABYSS/BLAST_first_hit_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./TRANSABYSS/BLAST_output_nt_first_hit_with_taxonomy_merge.txt ./TRANSABYSS/final_bwa_merge_ready.txt ./TRANSABYSS/BLAST_first_hit_nt_BWA.txt

			merge_mapped_reads_and_contigs.py ./TRANSABYSS/TRANSABYSS_blast_SILVA_LCA_taxonomy_merge.txt ./TRANSABYSS/final_bowtie2_merge_ready.txt ./TRANSABYSS/BLAST_LCA_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./TRANSABYSS/TRANSABYSS_blast_nt_LCA_taxonomy_merge.txt ./TRANSABYSS/final_bowtie2_merge_ready.txt ./TRANSABYSS/BLAST_LCA_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./TRANSABYSS/TRANSABYSS_blast_SILVA_LCA_taxonomy_merge.txt ./TRANSABYSS/final_bwa_merge_ready.txt ./TRANSABYSS/BLAST_LCA_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./TRANSABYSS/TRANSABYSS_blast_nt_LCA_taxonomy_merge.txt ./TRANSABYSS/final_bwa_merge_ready.txt ./TRANSABYSS/BLAST_LCA_nt_BWA.txt


			# Move all files - easy to find!!!
			mkdir ./TRANSABYSS/MERGE_FILES
			mv ./TRANSABYSS/*merge*.txt ./TRANSABYSS/MERGE_FILES/

			# Edit NODE_#_Length_#_Coverage_# in all the final files
			# CREST files
			sed 's/_/\t/2' ./TRANSABYSS/CREST_BWA.txt > ./TRANSABYSS/new.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt
			sed 's/NODE_//g' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt && sed 's/length_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/cov_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt
			sed '1d' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/CREST_BWA_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/CREST_BWA_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/CREST_BWA.txt

			sed 's/_/\t/2' ./TRANSABYSS/CREST_BOWTIE2.txt > ./TRANSABYSS/new.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt
			sed 's/NODE_//g' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt && sed 's/length_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/cov_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt
			sed '1d' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/CREST_BOWTIE2_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/CREST_BOWTIE2_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./TRANSABYSS/BLAST_first_hit_nt_BWA.txt > ./TRANSABYSS/new.txt
			sed 's/_/\t/2' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt
			sed 's/NODE_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/length_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt && sed 's/cov_//g' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/BLAST_nt_first_hit_BWA_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/BLAST_nt_first_hit_BWA_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/BLAST_first_hit_nt_BWA.txt

			sed '1d' ./TRANSABYSS/BLAST_first_hit_nt_BOWTIE2.txt > ./TRANSABYSS/new.txt
			sed 's/_/\t/2' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt
			sed 's/NODE_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/length_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt && sed 's/cov_//g' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/BLAST_nt_first_hit_BOWTIE2_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/BLAST_nt_first_hit_BOWTIE2_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/BLAST_first_hit_nt_BOWTIE2.txt

			sed '1d' ./TRANSABYSS/BLAST_LCA_nt_BWA.txt > ./TRANSABYSS/new.txt
			sed 's/_/\t/2' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt
			sed 's/NODE_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/length_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt && sed 's/cov_//g' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/BLAST_nt_LCA_BWA_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/BLAST_nt_LCA_BWA_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/BLAST_LCA_nt_BWA.txt

			sed '1d' ./TRANSABYSS/BLAST_LCA_nt_BOWTIE2.txt > ./TRANSABYSS/new.txt
			sed 's/_/\t/2' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt
			sed 's/NODE_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/length_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt && sed 's/cov_//g' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/BLAST_nt_LCA_BOWTIE2_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/BLAST_nt_LCA_BOWTIE2_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/BLAST_LCA_nt_BOWTIE2.txt

			# BLAST SILVA files
			sed '1d' ./TRANSABYSS/BLAST_first_hit_SILVA_BWA.txt > ./TRANSABYSS/new.txt
			sed 's/_/\t/2' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt
			sed 's/NODE_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/length_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt && sed 's/cov_//g' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/BLAST_SILVA_first_hit_BWA_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/BLAST_SILVA_first_hit_BWA_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/BLAST_first_hit_SILVA_BWA.txt

			sed '1d' ./TRANSABYSS/BLAST_first_hit_SILVA_BOWTIE2.txt > ./TRANSABYSS/new.txt
			sed 's/_/\t/2' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt
			sed 's/NODE_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/length_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt && sed 's/cov_//g' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/BLAST_SILVA_first_hit_BOWTIE2_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/BLAST_SILVA_first_hit_BOWTIE2_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/BLAST_first_hit_SILVA_BOWTIE2.txt

			sed '1d' ./TRANSABYSS/BLAST_LCA_SILVA_BWA.txt > ./TRANSABYSS/new.txt
			sed 's/_/\t/2' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt
			sed 's/NODE_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/length_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt && sed 's/cov_//g' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/BLAST_SILVA_LCA_BWA_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/BLAST_SILVA_LCA_BWA_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/BLAST_LCA_SILVA_BWA.txt

			sed '1d' ./TRANSABYSS/BLAST_LCA_SILVA_BOWTIE2.txt > ./TRANSABYSS/new.txt
			sed 's/_/\t/2' ./TRANSABYSS/new.txt > ./TRANSABYSS/new2.txt &&  sed 's/_/\t/3' ./TRANSABYSS/new2.txt > ./TRANSABYSS/new3.txt
			sed 's/NODE_//g' ./TRANSABYSS/new3.txt > ./TRANSABYSS/new4.txt && sed 's/length_//g' ./TRANSABYSS/new4.txt > ./TRANSABYSS/new5.txt && sed 's/cov_//g' ./TRANSABYSS/new5.txt > ./TRANSABYSS/new6.txt
			echo -e "contig_number\tcontig_length\tcontig_coverage\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRANSABYSS/BLAST_SILVA_LCA_BOWTIE2_final.txt && cat ./TRANSABYSS/new6.txt >> ./TRANSABYSS/BLAST_SILVA_LCA_BOWTIE2_final.txt
			rm ./TRANSABYSS/new*.txt ./TRANSABYSS/BLAST_LCA_SILVA_BOWTIE2.txt

			mkdir ./TRANSABYSS/FINAL_FILES_TRANSABYSS
			mv ./TRANSABYSS/*final.txt ./TRANSABYSS/FINAL_FILES_TRANSABYSS/
			mv ./TRANSABYSS/BLAST_output_nt* ./TRANSABYSS/BLAST_output_SILVA* TRANSABYSS/CLASSIFICATION/BLAST/

		echo -e "\n======== final TRANSABYSS output generated ========\n"

	else
		echo -e "\n======== no FINAL merge generated ========\n"

	fi

else
	echo -e "\n======== no TRANSABYSS output generated ========\n"

fi

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"


# Create log
) 2>&1 | tee PIPELINE_ASSEMBLERS_LOG.txt
