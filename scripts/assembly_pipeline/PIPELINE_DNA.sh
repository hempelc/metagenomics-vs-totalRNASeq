#!/bin/bash

# Need to have: assign_NCBI_staxids_to_CREST_v3.py, assign_taxonomy_NCBI_staxids.sh, deinterleave_fastq_reads.sh, LookupTaxonDetails3.py, merge_mapped_reads_and_contigs.py, mergeFilesOnColumn.pl and fasta_to_tab in your PATH
# In order for "assign_taxonomy_NCBI_staxids.sh" to work - MUST have .etetoolkit/taxa.sqlite in your HOME directory

usage="$(basename "$0") -1 Forward_read_trimmed -2 Reverse_read_trimmed [-S] [-M] [-m] [-I] [-B] [-C] [-f] [-r] [-h]

Usage:
	-1 Forward_read_trimmed - must state PATH to the file
	-2 Reverse_read_trimmed - must state PATH to the file
	-S Flag to generate SPADES assembly output 
	-M Flag to generate METASPADES assembly output
	-m Flag to generate MEGAHIT assembly output
	-I Flag to generate IDBA assembly output
	-B Flag to use BWA and BOWTIE2 - mapping
	-C Flag to use BLAST and CREST - classification
	-f Flag to produce final output files
	-r Flag to include original assembly contig/scaffold sequences in the final output files
	-h Display this help and exit"

# Set default options
Forward_read_trimmed=''
Reverse_read_trimmed=''
SPADES='false'
METASPADES='false'
MEGAHIT='false'
IDBA='false'
MAP='false'
CLASSIFICATION='false'
FINAL='false'
READS='false'

# Set specified options
while getopts ':1:2:SMmIBCfrh' opt; do
 	case "${opt}" in
		1) Forward_read_trimmed="${OPTARG}" ;;
		2) Reverse_read_trimmed="${OPTARG}" ;;
		S) SPADES='true' ;;
		M) METASPADES='true' ;;
		m) MEGAHIT='true' ;;
		I) IDBA='true' ;;
		B) MAP='true' ;;
		C) CLASSIFICATION='true' ;;
		f) FINAL='true' ;;
		r) READS='true' ;;
		h) echo "$usage"
			exit ;;
		:) printf "Option -$OPTARG requires an argument." >&2 ;;
		\?) printf "Invalid option: -$OPTARG" >&2 ;;
	esac
done
shift $((OPTIND - 1))


# Define starting time of script for total runtime calculation
start=$(date +%s)
echo -e "\nSTART RUNNING SCRIPT AT $(date)\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"

######################### Beginning of actual pipeline ######################## 

(

mkdir PIPELINE_DNA
cd PIPELINE_DNA/
	
if [[ $SPADES == 'true' ]] ; then
	mkdir SPADES 
	spades.py -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -o SPADES

	if [[ $MAP == 'true' ]] ; then

		echo "starting bwa index"

		bwa index -p bwa_index ./SPADES/scaffolds.fasta

		echo "bwa index complete."
		echo "starting bwa."

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

		echo "bwa complete."

		echo "starting bowtie2 index"

		bowtie2-build -f ./SPADES/scaffolds.fasta bowtie_index

		echo "bowtie2 index complete."
		echo "starting bowtie2."

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

		echo "bowtie2 complete"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ SPADES/

	else 

		echo "no BWA and BOWTIE2 output generated"

	fi

	if [[ $CLASSIFICATION == 'true' ]] ; then

		blastn -db /hdd1/databases/nt_database_feb_2020_indexed/nt -query ./SPADES/scaffolds.fasta -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_nt

		echo "BLAST nt output complete."

		blastn -db /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta -query ./SPADES/scaffolds.fasta -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_SILVA

		echo "BLAST SILVA output complete."

		mkdir SPADES/BLAST
		mv BLAST_output_nt BLAST_output_SILVA SPADES/BLAST/

		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./SPADES/scaffolds.fasta -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		rm BLAST_output.xml
		mv CREST SPADES/

		mkdir SPADES/CLASSIFICATION
		mv SPADES/BLAST SPADES/CREST SPADES/CLASSIFICATION/

		echo "CREST output complete."

	else 

		echo "no BLAST or CREST output generated"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then
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

			echo "assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."

		else 
			# Prepares to merge with BLAST of CREST - does not include assembly sequence
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

			echo "no assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."
		fi

		# Assign BLAST taxonomy - nt
		assign_taxonomy_NCBI_staxids.sh -b ./SPADES/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv _with_taxonomy.txt ./SPADES/BLAST_output_nt_with_taxonomy.txt
		sed '1d' ./SPADES/BLAST_output_nt_with_taxonomy.txt > ./SPADES/BLAST_output_nt_with_taxonomy_noheader.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./SPADES/BLAST_output_nt_with_taxonomy_merge.txt && cat ./SPADES/BLAST_output_nt_with_taxonomy_noheader.txt >> ./SPADES/BLAST_output_nt_with_taxonomy_merge.txt

		# Assign BLAST taxonomy - SILVA 
		assign_taxonomy_NCBI_staxids.sh -b ./SPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv _with_taxonomy.txt ./SPADES/BLAST_output_SILVA_with_taxonomy.txt
		sed '1d' ./SPADES/BLAST_output_SILVA_with_taxonomy.txt > ./SPADES/BLAST_output_SILVA_with_taxonomy_noheader.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./SPADES/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./SPADES/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./SPADES/BLAST_output_SILVA_with_taxonomy_merge.txt

		rm ./SPADES/BLAST_output_nt_with_taxonomy.txt ./SPADES/BLAST_output_nt_with_taxonomy_noheader.txt ./SPADES/BLAST_output_SILVA_with_taxonomy.txt ./SPADES/BLAST_output_SILVA_with_taxonomy_noheader.txt

		echo "taxonomy has been assigned to BLAST files."
		echo "BLAST files are ready."

		# Prepare CREST file - remove 1st/3rd column, add a header
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./SPADES/CLASSIFICATION/CREST/otus.csv ./SPADES/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./SPADES/CLASSIFICATION/CREST/CREST_output.txt > ./SPADES/CREST_seperated.txt
		cut -f2,4 ./SPADES/CREST_seperated.txt > ./SPADES/CREST_header.txt
		sed '1d' ./SPADES/CREST_header.txt > ./SPADES/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./SPADES/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./SPADES/
		sed '1d' ./SPADES/CREST_tax_ready_with_taxonomy.txt > ./SPADES/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./SPADES/CREST_merge.txt && cat ./SPADES/CREST_tax_ready_with_taxonomy_noheader.txt >> ./SPADES/CREST_merge.txt

		rm ./SPADES/CREST_seperated.txt ./SPADES/CREST_header.txt ./SPADES/CREST_tax_ready.txt ./SPADES/CREST_tax_ready_with_taxonomy.txt ./SPADES/CREST_tax_ready_with_taxonomy_noheader.txt 

		echo "CREST file is ready."

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./SPADES/CREST_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./SPADES/CREST_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/BLAST_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_nt_with_taxonomy_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/BLAST_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/BLAST_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_nt_with_taxonomy_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/BLAST_nt_BWA.txt

		# Move all files - easy to find!!!
		mkdir ./SPADES/MERGE_FILES
		mv ./SPADES/final_bwa_merge_ready.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/CREST_merge.txt ./SPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./SPADES/BLAST_output_nt_with_taxonomy_merge.txt ./SPADES/MERGE_FILES/

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
		sed '1d' ./SPADES/BLAST_nt_BWA.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_nt_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_nt_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_nt_BWA.txt

		sed '1d' ./SPADES/BLAST_nt_BOWTIE2.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_nt_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_nt_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./SPADES/BLAST_SILVA_BWA.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_SILVA_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_SILVA_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_SILVA_BWA.txt

		sed '1d' ./SPADES/BLAST_SILVA_BOWTIE2.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_SILVA_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_SILVA_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_SILVA_BOWTIE2.txt

		mkdir ./SPADES/FINAL_FILES_SPADES
		mv ./SPADES/CREST_BWA_final.txt ./SPADES/CREST_BOWTIE2_final.txt ./SPADES/BLAST_SILVA_BOWTIE2_final.txt ./SPADES/BLAST_nt_BOWTIE2_final.txt ./SPADES/BLAST_SILVA_BWA_final.txt ./SPADES/BLAST_nt_BWA_final.txt ./SPADES/FINAL_FILES_SPADES/

		echo "final SPADES output generated."

	else 

		echo "no FINAL merge generated."

	fi

else
	echo "no SPADES output generated."
fi 


if [[ $METASPADES == 'true' ]] ; then
	mkdir METASPADES 
	metaspades.py -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -o METASPADES

	if [[ $MAP == 'true' ]] ; then

		echo "starting bwa index"

		bwa index -p bwa_index ./METASPADES/scaffolds.fasta

		echo "bwa index complete."
		echo "starting bwa."

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

		echo "bwa complete."

		echo "starting bowtie2 index"

		bowtie2-build -f ./METASPADES/scaffolds.fasta bowtie_index

		echo "bowtie2 index complete."
		echo "starting bowtie2."

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

		echo "bowtie2 complete"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ METASPADES/

	else 

		echo "no BWA and BOWTIE2 output generated"

	fi


	if [[ $CLASSIFICATION == 'true' ]] ; then

		blastn -db /hdd1/databases/nt_database_feb_2020_indexed/nt -query ./METASPADES/scaffolds.fasta -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_nt

		echo "BLAST nt output complete."

		blastn -db /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta -query ./METASPADES/scaffolds.fasta -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_SILVA

		echo "BLAST SILVA output complete."

		mkdir METASPADES/BLAST
		mv BLAST_output_nt BLAST_output_SILVA METASPADES/BLAST/

		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./METASPADES/scaffolds.fasta -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		rm BLAST_output.xml
		mv CREST METASPADES/

		mkdir METASPADES/CLASSIFICATION
		mv METASPADES/BLAST METASPADES/CREST METASPADES/CLASSIFICATION/

		echo "CREST output complete."

	else 

		echo "no BLAST or CREST output generated"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then
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

			echo "assembly sequence was added."
			echo "BWA and BOWTIE2 files are ready."

		else 
			# Prepares to merge with BLAST of CREST - does not include assembly sequence 
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

			echo "no assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."
		fi

		# Assign BLAST taxonomy - nt
		assign_taxonomy_NCBI_staxids.sh -b ./METASPADES/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv _with_taxonomy.txt ./METASPADES/BLAST_output_nt_with_taxonomy.txt
		sed '1d' ./METASPADES/BLAST_output_nt_with_taxonomy.txt > ./METASPADES/BLAST_output_nt_with_taxonomy_noheader.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./METASPADES/BLAST_output_nt_with_taxonomy_merge.txt && cat ./METASPADES/BLAST_output_nt_with_taxonomy_noheader.txt >> ./METASPADES/BLAST_output_nt_with_taxonomy_merge.txt

		# Assign BLAST taxonomy - SILVA 
		assign_taxonomy_NCBI_staxids.sh -b ./METASPADES/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv _with_taxonomy.txt ./METASPADES/BLAST_output_SILVA_with_taxonomy.txt
		sed '1d' ./METASPADES/BLAST_output_SILVA_with_taxonomy.txt > ./METASPADES/BLAST_output_SILVA_with_taxonomy_noheader.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./METASPADES/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./METASPADES/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./METASPADES/BLAST_output_SILVA_with_taxonomy_merge.txt

		rm ./METASPADES/BLAST_output_nt_with_taxonomy.txt ./METASPADES/BLAST_output_nt_with_taxonomy_noheader.txt ./METASPADES/BLAST_output_SILVA_with_taxonomy.txt ./METASPADES/BLAST_output_SILVA_with_taxonomy_noheader.txt

		echo "taxonomy has been assigned to BLAST files."
		echo "BLAST files are ready."

		# Prepare CREST file - remove 1st/3rd column, add a header
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./METASPADES/CLASSIFICATION/CREST/otus.csv ./METASPADES/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./METASPADES/CLASSIFICATION/CREST/CREST_output.txt > ./METASPADES/CREST_seperated.txt
		cut -f2,4 ./METASPADES/CREST_seperated.txt > ./METASPADES/CREST_header.txt
		sed '1d' ./METASPADES/CREST_header.txt > ./METASPADES/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./METASPADES/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./METASPADES/
		sed '1d' ./METASPADES/CREST_tax_ready_with_taxonomy.txt > ./METASPADES/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./METASPADES/CREST_merge.txt && cat ./METASPADES/CREST_tax_ready_with_taxonomy_noheader.txt >> ./METASPADES/CREST_merge.txt

		rm ./METASPADES/CREST_seperated.txt ./METASPADES/CREST_header.txt ./METASPADES/CREST_tax_ready.txt ./METASPADES/CREST_tax_ready_with_taxonomy.txt ./METASPADES/CREST_tax_ready_with_taxonomy_noheader.txt 

		echo "CREST file is ready."

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./METASPADES/CREST_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/CREST_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/BLAST_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_nt_with_taxonomy_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/BLAST_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/BLAST_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_nt_with_taxonomy_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/BLAST_nt_BWA.txt

		# Move all files - easy to find!!!
		mkdir ./METASPADES/MERGE_FILES
		mv ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/CREST_merge.txt ./METASPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./METASPADES/BLAST_output_nt_with_taxonomy_merge.txt ./METASPADES/MERGE_FILES/

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
		sed '1d' ./METASPADES/BLAST_nt_BWA.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_nt_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_nt_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_nt_BWA.txt

		sed '1d' ./METASPADES/BLAST_nt_BOWTIE2.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_nt_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_nt_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./METASPADES/BLAST_SILVA_BWA.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_SILVA_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_SILVA_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_SILVA_BWA.txt

		sed '1d' ./METASPADES/BLAST_SILVA_BOWTIE2.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_SILVA_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_SILVA_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_SILVA_BOWTIE2.txt

		mkdir ./METASPADES/FINAL_FILES_METASPADES
		mv ./METASPADES/CREST_BWA_final.txt ./METASPADES/CREST_BOWTIE2_final.txt ./METASPADES/BLAST_SILVA_BOWTIE2_final.txt ./METASPADES/BLAST_nt_BOWTIE2_final.txt ./METASPADES/BLAST_SILVA_BWA_final.txt ./METASPADES/BLAST_nt_BWA_final.txt ./METASPADES/FINAL_FILES_METASPADES/

		echo "final METASPADES output generated."

	else 

		echo "no FINAL merge generated."

	fi

else
	echo "no METASPADES output generated."
fi 


if [[ $MEGAHIT == 'true' ]] ; then
	mkdir MEGAHIT 
	megahit --presets meta-large -t 16 -1 $Forward_read_trimmed -2 $Reverse_read_trimmed
	mv megahit_out MEGAHIT/

	if [[ $MAP == 'true' ]] ; then

		echo "starting bwa index"

		bwa index -p bwa_index ./MEGAHIT/megahit_out/final.contigs.fa

		echo "bwa index complete."
		echo "starting bwa."

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

		echo "bwa complete."

		echo "starting bowtie2 index"

		bowtie2-build -f ./MEGAHIT/megahit_out/final.contigs.fa bowtie_index

		echo "bowtie2 index complete."
		echo "starting bowtie2."

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

		echo "bowtie2 complete"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ MEGAHIT/megahit_out/

	else 

		echo "no BWA and BOWTIE2 output generated"

	fi


	if [[ $CLASSIFICATION == 'true' ]] ; then

		blastn -db /hdd1/databases/nt_database_feb_2020_indexed/nt -query ./MEGAHIT/megahit_out/final.contigs.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_nt

		echo "BLAST nt output complete."

		blastn -db /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta -query ./MEGAHIT/megahit_out/final.contigs.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_SILVA

		echo "BLAST SILVA output complete."

		mkdir MEGAHIT/megahit_out/BLAST
		mv BLAST_output_nt BLAST_output_SILVA MEGAHIT/megahit_out/BLAST/

		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./MEGAHIT/megahit_out/final.contigs.fa -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		rm BLAST_output.xml
		mv CREST MEGAHIT/megahit_out/

		mkdir MEGAHIT/megahit_out/CLASSIFICATION
		mv MEGAHIT/megahit_out/BLAST MEGAHIT/megahit_out/CREST MEGAHIT/megahit_out/CLASSIFICATION/

		echo "CREST output complete."

	else 

		echo "no BLAST or CREST output generated"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			# Prepare assembly output - fasta to tabular
			fasta_to_tab ./MEGAHIT/megahit_out/final.contigs.fa > ./MEGAHIT/megahit_out/fasta_to_tabbed.txt
			sed 's/ /\t/g' ./MEGAHIT/megahit_out/fasta_to_tabbed.txt > ./MEGAHIT/megahit_out/fasta_to_tabbed_tab.txt
			cut -f1,4,5 ./MEGAHIT/megahit_out/fasta_to_tabbed_tab.txt > ./MEGAHIT/megahit_out/important_columns.txt
			sed 's/len=//g' ./MEGAHIT/megahit_out/important_columns.txt > ./MEGAHIT/megahit_out/tab_to_merge.txt

			# Merge BWA and assembly sequences 
			mergeFilesOnColumn.pl ./MEGAHIT/megahit_out/BWA/merge_input_mapped_bwa.txt ./MEGAHIT/megahit_out/tab_to_merge.txt 2 1 > ./MEGAHIT/megahit_out/merged_original_bwa.txt
			cut -f1,2,4,5 ./MEGAHIT/megahit_out/merged_original_bwa.txt > ./MEGAHIT/megahit_out/important_column_4_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $1, $4}' ./MEGAHIT/megahit_out/important_column_4_bwa.txt > ./MEGAHIT/megahit_out/final_order_bwa.txt
			echo -e "contig_number\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt && cat ./MEGAHIT/megahit_out/final_order_bwa.txt >> ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt

			rm ./MEGAHIT/megahit_out/fasta_to_tabbed.txt ./MEGAHIT/megahit_out/fasta_to_tabbed_tab.txt ./MEGAHIT/megahit_out/important_columns.txt ./MEGAHIT/megahit_out/merged_original_bwa.txt ./MEGAHIT/megahit_out/important_column_4_bwa.txt ./MEGAHIT/megahit_out/final_order_bwa.txt

			# Merge BOWTIE2 and assembly sequences
			mergeFilesOnColumn.pl ./MEGAHIT/megahit_out/BOWTIE2/merge_input_mapped_bowtie2.txt ./MEGAHIT/megahit_out/tab_to_merge.txt 2 1 > ./MEGAHIT/megahit_out/merged_original_bowtie2.txt
			cut -f1,2,4,5 ./MEGAHIT/megahit_out/merged_original_bowtie2.txt > ./MEGAHIT/megahit_out/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $1, $4}' ./MEGAHIT/megahit_out/important_column_3_bowtie2.txt > ./MEGAHIT/megahit_out/final_order_bowtie2.txt
			echo -e "contig_number\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt && cat ./MEGAHIT/megahit_out/final_order_bowtie2.txt >> ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt

			rm ./MEGAHIT/megahit_out/tab_to_merge.txt ./MEGAHIT/megahit_out/merged_original_bowtie2.txt ./MEGAHIT/megahit_out/important_column_3_bowtie2.txt ./MEGAHIT/megahit_out/final_order_bowtie2.txt

			echo "assembly sequence was added."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/megahit_out/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy.txt > ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt && cat ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_noheader.txt >> ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA 
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/megahit_out/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy.txt > ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy.txt ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_noheader.txt ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy.txt ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."


			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./MEGAHIT/megahit_out/CLASSIFICATION/CREST/otus.csv ./MEGAHIT/megahit_out/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./MEGAHIT/megahit_out/CLASSIFICATION/CREST/CREST_output.txt > ./MEGAHIT/megahit_out/CREST_seperated.txt
			cut -f2,4 ./MEGAHIT/megahit_out/CREST_seperated.txt > ./MEGAHIT/megahit_out/CREST_header.txt
			sed '1d' ./MEGAHIT/megahit_out/CREST_header.txt > ./MEGAHIT/megahit_out/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/megahit_out/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/megahit_out/
			sed '1d' ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy.txt > ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/megahit_out/CREST_merge.txt && cat ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy_noheader.txt >> ./MEGAHIT/megahit_out/CREST_merge.txt

			rm ./MEGAHIT/megahit_out/CREST_seperated.txt ./MEGAHIT/megahit_out/CREST_header.txt ./MEGAHIT/megahit_out/CREST_tax_ready.txt ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/CREST_merge.txt ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt ./MEGAHIT/megahit_out/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/CREST_merge.txt ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt ./MEGAHIT/megahit_out/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt ./MEGAHIT/megahit_out/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./MEGAHIT/megahit_out/MERGE_FILES
			mv ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt ./MEGAHIT/megahit_out/CREST_merge.txt ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/MERGE_FILES/

			# Edit k*_ - in all files 
			# CREST files
			sed '1d' ./MEGAHIT/megahit_out/CREST_BWA.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt
			cut -f2-19 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/CREST_BWA_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/CREST_BWA_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/CREST_BWA.txt

			sed '1d' ./MEGAHIT/megahit_out/CREST_BOWTIE2.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt
			cut -f2-19 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/CREST_BOWTIE2_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/CREST_BOWTIE2_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./MEGAHIT/megahit_out/BLAST_nt_BWA.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt
			cut -f2-30 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/BLAST_nt_BWA_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/BLAST_nt_BWA_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/BLAST_nt_BWA.txt

			sed '1d' ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt 
			cut -f2-30 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2.txt


			# BLAST SILVA files
			sed '1d' ./MEGAHIT/megahit_out/BLAST_SILVA_BWA.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt 
			cut -f2-30 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/BLAST_SILVA_BWA_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/BLAST_SILVA_BWA_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BWA.txt

			sed '1d' ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt 
			cut -f2-30 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2.txt

			mkdir ./MEGAHIT/megahit_out/FINAL_FILES_MEGAHIT
			mv ./MEGAHIT/megahit_out/CREST_BWA_final.txt ./MEGAHIT/megahit_out/CREST_BOWTIE2_final.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2_final.txt ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2_final.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BWA_final.txt ./MEGAHIT/megahit_out/BLAST_nt_BWA_final.txt ./MEGAHIT/megahit_out/FINAL_FILES_MEGAHIT/

		else 

			# Prepares to merge with BLAST of CREST - does not include assembly sequence 
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./MEGAHIT/megahit_out/BWA/merge_input_mapped_bwa.txt > ./MEGAHIT/megahit_out/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./MEGAHIT/megahit_out/temp.txt > ./MEGAHIT/megahit_out/temp2.txt
			sed '1d' ./MEGAHIT/megahit_out/temp2.txt > ./MEGAHIT/megahit_out/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt && cat ./MEGAHIT/megahit_out/temp3.txt >> ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt
 			rm ./MEGAHIT/megahit_out/temp*.txt 

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./MEGAHIT/megahit_out/BOWTIE2/merge_input_mapped_bowtie2.txt > ./MEGAHIT/megahit_out/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./MEGAHIT/megahit_out/temp.txt > ./MEGAHIT/megahit_out/temp2.txt
			sed '1d' ./MEGAHIT/megahit_out/temp2.txt > ./MEGAHIT/megahit_out/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt && cat ./MEGAHIT/megahit_out/temp3.txt >> ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt
			rm ./MEGAHIT/megahit_out/temp*.txt

			echo "no assembly sequence added."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/megahit_out/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy.txt > ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt && cat ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_noheader.txt >> ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA 
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/megahit_out/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy.txt > ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy.txt ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_noheader.txt ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy.txt ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."


			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./MEGAHIT/megahit_out/CLASSIFICATION/CREST/otus.csv ./MEGAHIT/megahit_out/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./MEGAHIT/megahit_out/CLASSIFICATION/CREST/CREST_output.txt > ./MEGAHIT/megahit_out/CREST_seperated.txt
			cut -f2,4 ./MEGAHIT/megahit_out/CREST_seperated.txt > ./MEGAHIT/megahit_out/CREST_header.txt
			sed '1d' ./MEGAHIT/megahit_out/CREST_header.txt > ./MEGAHIT/megahit_out/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/megahit_out/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/megahit_out/
			sed '1d' ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy.txt > ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/megahit_out/CREST_merge.txt && cat ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy_noheader.txt >> ./MEGAHIT/megahit_out/CREST_merge.txt

			rm ./MEGAHIT/megahit_out/CREST_seperated.txt ./MEGAHIT/megahit_out/CREST_header.txt ./MEGAHIT/megahit_out/CREST_tax_ready.txt ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/megahit_out/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."


			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/CREST_merge.txt ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt ./MEGAHIT/megahit_out/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/CREST_merge.txt ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt ./MEGAHIT/megahit_out/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt ./MEGAHIT/megahit_out/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./MEGAHIT/megahit_out/MERGE_FILES
			mv ./MEGAHIT/megahit_out/final_bwa_merge_ready.txt ./MEGAHIT/megahit_out/final_bowtie2_merge_ready.txt ./MEGAHIT/megahit_out/CREST_merge.txt ./MEGAHIT/megahit_out/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/megahit_out/MERGE_FILES/
		

			# CREST files
			sed '1d' ./MEGAHIT/megahit_out/CREST_BWA.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt
			cut -f2-18 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/CREST_BWA_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/CREST_BWA_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/CREST_BWA.txt

			sed '1d' ./MEGAHIT/megahit_out/CREST_BOWTIE2.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt
			cut -f2-18 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/CREST_BOWTIE2_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/CREST_BOWTIE2_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/CREST_BOWTIE2.txt


			# BLAST nt files
			sed '1d' ./MEGAHIT/megahit_out/BLAST_nt_BWA.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt
			cut -f2-29 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/BLAST_nt_BWA_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/BLAST_nt_BWA_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/BLAST_nt_BWA.txt

			sed '1d' ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt 
			cut -f2-29 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2.txt


			# BLAST SILVA files
			sed '1d' ./MEGAHIT/megahit_out/BLAST_SILVA_BWA.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt 
			cut -f2-29 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/BLAST_SILVA_BWA_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/BLAST_SILVA_BWA_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BWA.txt

			sed '1d' ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2.txt > ./MEGAHIT/megahit_out/new.txt
			sed 's/_/\t/1' ./MEGAHIT/megahit_out/new.txt > ./MEGAHIT/megahit_out/new2.txt 
			cut -f2-29 ./MEGAHIT/megahit_out/new2.txt > ./MEGAHIT/megahit_out/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2_final.txt && cat ./MEGAHIT/megahit_out/new3.txt >> ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2_final.txt
			rm ./MEGAHIT/megahit_out/new*.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2.txt

			mkdir ./MEGAHIT/megahit_out/FINAL_FILES_MEGAHIT
			mv ./MEGAHIT/megahit_out/CREST_BWA_final.txt ./MEGAHIT/megahit_out/CREST_BOWTIE2_final.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BOWTIE2_final.txt ./MEGAHIT/megahit_out/BLAST_nt_BOWTIE2_final.txt ./MEGAHIT/megahit_out/BLAST_SILVA_BWA_final.txt ./MEGAHIT/megahit_out/BLAST_nt_BWA_final.txt ./MEGAHIT/megahit_out/FINAL_FILES_MEGAHIT/

			echo "final MEGAHIT output generated."

		fi

	else 

		echo "no FINAL merge generated."

	fi	

else
	echo "no MEGAHIT output generated."
fi


if [[ $IDBA == 'true' ]] ; then
	fq2fa --merge --filter $Forward_read_trimmed $Reverse_read_trimmed new_input.fa
	idba_ud --num_threads 16 --pre_correction -r ./new_input.fa -o IDBA

		if [[ $MAP == 'true' ]] ; then

		echo "starting bwa index"

		bwa index -p bwa_index ./IDBA/contig.fa

		echo "bwa index complete."
		echo "starting bwa."

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

		echo "bwa complete."

		echo "starting bowtie2 index"

		bowtie2-build -f ./IDBA/contig.fa bowtie_index

		echo "bowtie2 index complete."
		echo "starting bowtie2."

		bowtie2 -q -x bowtie_index -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -S bowtie2_output.sam

		rm bowtie_index* new_input.fa


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

		echo "bowtie2 complete"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ IDBA/

	else 

		echo "no BWA and BOWTIE2 output generated"

	fi


	if [[ $CLASSIFICATION == 'true' ]] ; then

		blastn -db /hdd1/databases/nt_database_feb_2020_indexed/nt -query ./IDBA/contig.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_nt

		echo "BLAST nt output complete."

		blastn -db /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta -query ./IDBA/contig.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_SILVA

		echo "BLAST SILVA output complete."

		mkdir IDBA/BLAST
		mv BLAST_output_nt BLAST_output_SILVA IDBA/BLAST/

		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./IDBA/contig.fa -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		rm BLAST_output.xml
		mv CREST IDBA/

		mkdir IDBA/CLASSIFICATION
		mv IDBA/BLAST IDBA/CREST IDBA/CLASSIFICATION/

		echo "CREST output complete."

	else 

		echo "no BLAST or CREST output generated"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./IDBA/contig.fa > ./IDBA/fasta_to_tabbed.txt
			sed 's/_/\t/2' ./IDBA/fasta_to_tabbed.txt > ./IDBA/fasta_to_tabbed_tab.txt
			sed 's/_/\t/3' ./IDBA/fasta_to_tabbed_tab.txt > ./IDBA/fasta_to_tabbed_tab_2.txt
			sed 's/ /\t/g' ./IDBA/fasta_to_tabbed_tab_2.txt > ./IDBA/cut_tab_ready.txt
			cut -f1,3,5,6 ./IDBA/cut_tab_ready.txt > ./IDBA/important_columns.txt
			rm ./IDBA/fasta_to_tabbed.txt ./IDBA/fasta_to_tabbed_tab.txt ./IDBA/fasta_to_tabbed_tab_2.txt ./IDBA/cut_tab_ready.txt

			mergeFilesOnColumn.pl ./IDBA/BWA/merge_input_mapped_bwa.txt ./IDBA/important_columns.txt 2 1 > ./IDBA/merged_original_bwa.txt
			cut -f1,2,4,5,6 ./IDBA/merged_original_bwa.txt > ./IDBA/important_columns_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $4, $1, $5}' ./IDBA/important_columns_bwa.txt > ./IDBA/final_order_bwa.txt
			echo -e "contig_number\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA/final_bwa_merge_ready.txt && cat ./IDBA/final_order_bwa.txt >> ./IDBA/final_bwa_merge_ready.txt
			rm ./IDBA/merged_original_bwa.txt ./IDBA/important_columns_bwa.txt ./IDBA/final_order_bwa.txt


			mergeFilesOnColumn.pl ./IDBA/BOWTIE2/merge_input_mapped_bowtie2.txt ./IDBA/important_columns.txt 2 1 > ./IDBA/merged_original_bowtie2.txt
			cut -f1,2,4,5,6 ./IDBA/merged_original_bowtie2.txt > ./IDBA/important_columns_bowtie.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $4, $1, $5}' ./IDBA/important_columns_bowtie.txt > ./IDBA/final_order_bowtie2.txt
			echo -e "contig_number\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA/final_bowtie2_merge_ready.txt && cat ./IDBA/final_order_bowtie2.txt >> ./IDBA/final_bowtie2_merge_ready.txt
			rm ./IDBA/merged_original_bowtie2.txt ./IDBA/important_columns_bowtie.txt ./IDBA/final_order_bowtie2.txt ./IDBA/important_columns.txt

			echo "assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./IDBA/BLAST_output_nt_with_taxonomy.txt > ./IDBA/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt && cat ./IDBA/BLAST_output_nt_with_taxonomy_noheader.txt >> ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA 
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./IDBA/BLAST_output_SILVA_with_taxonomy.txt > ./IDBA/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./IDBA/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./IDBA/BLAST_output_nt_with_taxonomy.txt ./IDBA/BLAST_output_nt_with_taxonomy_noheader.txt ./IDBA/BLAST_output_SILVA_with_taxonomy.txt ./IDBA/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."

			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA/CLASSIFICATION/CREST/otus.csv ./IDBA/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./IDBA/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA/CREST_seperated.txt
			cut -f2,4 ./IDBA/CREST_seperated.txt > ./IDBA/CREST_header.txt
			sed '1d' ./IDBA/CREST_header.txt > ./IDBA/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./IDBA/
			sed '1d' ./IDBA/CREST_tax_ready_with_taxonomy.txt > ./IDBA/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA/CREST_merge.txt && cat ./IDBA/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA/CREST_merge.txt

			rm ./IDBA/CREST_seperated.txt ./IDBA/CREST_header.txt ./IDBA/CREST_tax_ready.txt ./IDBA/CREST_tax_ready_with_taxonomy.txt ./IDBA/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./IDBA/CREST_merge.txt ./IDBA/final_bwa_merge_ready.txt ./IDBA/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA/CREST_merge.txt ./IDBA/final_bowtie2_merge_ready.txt ./IDBA/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA/final_bowtie2_merge_ready.txt ./IDBA/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA/final_bowtie2_merge_ready.txt ./IDBA/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA/final_bwa_merge_ready.txt ./IDBA/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA/final_bwa_merge_ready.txt ./IDBA/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./IDBA/MERGE_FILES
			mv ./IDBA/final_bwa_merge_ready.txt ./IDBA/final_bowtie2_merge_ready.txt ./IDBA/CREST_merge.txt ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA/MERGE_FILES/

			# Edit contig-* in all the final files
			# CREST files
			sed '1d' ./IDBA/CREST_BWA.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-20 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA/CREST_BWA_final.txt && cat ./IDBA/new3.txt >> ./IDBA/CREST_BWA_final.txt
			rm ./IDBA/new*.txt ./IDBA/CREST_BWA.txt

			sed '1d' ./IDBA/CREST_BOWTIE2.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-20 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA/CREST_BOWTIE2_final.txt && cat ./IDBA/new3.txt >> ./IDBA/CREST_BOWTIE2_final.txt
			rm ./IDBA/new*.txt ./IDBA/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./IDBA/BLAST_nt_BWA.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-31 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA/BLAST_nt_BWA_final.txt && cat ./IDBA/new3.txt >> ./IDBA/BLAST_nt_BWA_final.txt
			rm ./IDBA/new*.txt ./IDBA/BLAST_nt_BWA.txt

			sed '1d' ./IDBA/BLAST_nt_BOWTIE2.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-31 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA/BLAST_nt_BOWTIE2_final.txt && cat ./IDBA/new3.txt >> ./IDBA/BLAST_nt_BOWTIE2_final.txt
			rm ./IDBA/new*.txt ./IDBA/BLAST_nt_BOWTIE2.txt

			# BLAST SILVA files
			sed '1d' ./IDBA/BLAST_SILVA_BWA.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-31 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA/BLAST_SILVA_BWA_final.txt && cat ./IDBA/new3.txt >> ./IDBA/BLAST_SILVA_BWA_final.txt
			rm ./IDBA/new*.txt ./IDBA/BLAST_SILVA_BWA.txt

			sed '1d' ./IDBA/BLAST_SILVA_BOWTIE2.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-31 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA/BLAST_SILVA_BOWTIE2_final.txt && cat ./IDBA/new3.txt >> ./IDBA/BLAST_SILVA_BOWTIE2_final.txt
			rm ./IDBA/new*.txt ./IDBA/BLAST_SILVA_BOWTIE2.txt

			mkdir ./IDBA/FINAL_FILES_IDBA
			mv ./IDBA/CREST_BWA_final.txt ./IDBA/CREST_BOWTIE2_final.txt ./IDBA/BLAST_SILVA_BOWTIE2_final.txt ./IDBA/BLAST_nt_BOWTIE2_final.txt ./IDBA/BLAST_SILVA_BWA_final.txt ./IDBA/BLAST_nt_BWA_final.txt ./IDBA/FINAL_FILES_IDBA/

			echo "final IDBA output generated."

		else 
			# Prepares to merge with BLAST of CREST - does not include assembly sequence
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./IDBA/BWA/merge_input_mapped_bwa.txt > ./IDBA/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./IDBA/temp.txt > ./IDBA/temp2.txt
			sed '1d' ./IDBA/temp2.txt > ./IDBA/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA/final_bwa_merge_ready.txt && cat ./IDBA/temp3.txt >> ./IDBA/final_bwa_merge_ready.txt
 			rm ./IDBA/temp*.txt 

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./IDBA/BOWTIE2/merge_input_mapped_bowtie2.txt > ./IDBA/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./IDBA/temp.txt > ./IDBA/temp2.txt
			sed '1d' ./IDBA/temp2.txt > ./IDBA/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./IDBA/final_bowtie2_merge_ready.txt && cat ./IDBA/temp3.txt >> ./IDBA/final_bowtie2_merge_ready.txt
			rm ./IDBA/temp*.txt 

			echo "no assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./IDBA/BLAST_output_nt_with_taxonomy.txt > ./IDBA/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt && cat ./IDBA/BLAST_output_nt_with_taxonomy_noheader.txt >> ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA 
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./IDBA/BLAST_output_SILVA_with_taxonomy.txt > ./IDBA/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./IDBA/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./IDBA/BLAST_output_nt_with_taxonomy.txt ./IDBA/BLAST_output_nt_with_taxonomy_noheader.txt ./IDBA/BLAST_output_SILVA_with_taxonomy.txt ./IDBA/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."


			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA/CLASSIFICATION/CREST/otus.csv ./IDBA/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./IDBA/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA/CREST_seperated.txt
			cut -f2,4 ./IDBA/CREST_seperated.txt > ./IDBA/CREST_header.txt
			sed '1d' ./IDBA/CREST_header.txt > ./IDBA/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./IDBA/
			sed '1d' ./IDBA/CREST_tax_ready_with_taxonomy.txt > ./IDBA/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA/CREST_merge.txt && cat ./IDBA/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA/CREST_merge.txt

			rm ./IDBA/CREST_seperated.txt ./IDBA/CREST_header.txt ./IDBA/CREST_tax_ready.txt ./IDBA/CREST_tax_ready_with_taxonomy.txt ./IDBA/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./IDBA/CREST_merge.txt ./IDBA/final_bwa_merge_ready.txt ./IDBA/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA/CREST_merge.txt ./IDBA/final_bowtie2_merge_ready.txt ./IDBA/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA/final_bowtie2_merge_ready.txt ./IDBA/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA/final_bowtie2_merge_ready.txt ./IDBA/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA/final_bwa_merge_ready.txt ./IDBA/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA/final_bwa_merge_ready.txt ./IDBA/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./IDBA/MERGE_FILES
			mv ./IDBA/final_bwa_merge_ready.txt ./IDBA/final_bowtie2_merge_ready.txt ./IDBA/CREST_merge.txt ./IDBA/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA/MERGE_FILES/


			# Edit contig-* in all the final files
			# CREST files
			sed '1d' ./IDBA/CREST_BWA.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-18 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA/CREST_BWA_final.txt && cat ./IDBA/new3.txt >> ./IDBA/CREST_BWA_final.txt
			rm ./IDBA/new*.txt ./IDBA/CREST_BWA.txt

			sed '1d' ./IDBA/CREST_BOWTIE2.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-18 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA/CREST_BOWTIE2_final.txt && cat ./IDBA/new3.txt >> ./IDBA/CREST_BOWTIE2_final.txt
			rm ./IDBA/new*.txt ./IDBA/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./IDBA/BLAST_nt_BWA.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-29 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA/BLAST_nt_BWA_final.txt && cat ./IDBA/new3.txt >> ./IDBA/BLAST_nt_BWA_final.txt
			rm ./IDBA/new*.txt ./IDBA/BLAST_nt_BWA.txt

			sed '1d' ./IDBA/BLAST_nt_BOWTIE2.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-29 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA/BLAST_nt_BOWTIE2_final.txt && cat ./IDBA/new3.txt >> ./IDBA/BLAST_nt_BOWTIE2_final.txt
			rm ./IDBA/new*.txt ./IDBA/BLAST_nt_BOWTIE2.txt

			# BLAST SILVA files
			sed '1d' ./IDBA/BLAST_SILVA_BWA.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-29 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA/BLAST_SILVA_BWA_final.txt && cat ./IDBA/new3.txt >> ./IDBA/BLAST_SILVA_BWA_final.txt
			rm ./IDBA/new*.txt ./IDBA/BLAST_SILVA_BWA.txt

			sed '1d' ./IDBA/BLAST_SILVA_BOWTIE2.txt > ./IDBA/new.txt
			sed 's/_/\t/1' ./IDBA/new.txt > ./IDBA/new2.txt
			cut -f2-29 ./IDBA/new2.txt > ./IDBA/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA/BLAST_SILVA_BOWTIE2_final.txt && cat ./IDBA/new3.txt >> ./IDBA/BLAST_SILVA_BOWTIE2_final.txt
			rm ./IDBA/new*.txt ./IDBA/BLAST_SILVA_BOWTIE2.txt

			mkdir ./IDBA/FINAL_FILES_IDBA
			mv ./IDBA/CREST_BWA_final.txt ./IDBA/CREST_BOWTIE2_final.txt ./IDBA/BLAST_SILVA_BOWTIE2_final.txt ./IDBA/BLAST_nt_BOWTIE2_final.txt ./IDBA/BLAST_SILVA_BWA_final.txt ./IDBA/BLAST_nt_BWA_final.txt ./IDBA/FINAL_FILES_IDBA/

			echo "final IDBA output generated."

		fi

	else 

		echo "no FINAL merge generated."

	fi

else

	echo "no IDBA output generated."

fi 

mkdir ASSEMBLERS
mv SPADES METASPADES MEGAHIT IDBA ASSEMBLERS/


# Display runtime
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nSCRIPT RUNTIME: $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"

# Create log
) 2>&1 | tee PIPELINE_DNA_LOG.txt
mv PIPELINE_DNA_LOG.txt PIPELINE_DNA/

echo "all done."