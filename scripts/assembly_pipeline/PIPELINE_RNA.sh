#!/bin/bash

# Need to have: assign_NCBI_staxids_to_CREST_v3.py, assign_taxonomy_NCBI_staxids.sh, deinterleave_fastq_reads.sh, LookupTaxonDetails3.py, merge_mapped_reads_and_contigs.py, mergeFilesOnColumn.pl and fasta_to_tab in your PATH
# In order for "assign_taxonomy_NCBI_staxids.sh" to work - MUST have .etetoolkit/taxa.sqlite in your HOME directory

usage="$(basename "$0") -1 Forward_read_trimmed -2 Reverse_read_trimmed [-S] [-I] [-T] [-B] [-C] [-f] [-r] [-h]

Usage:
	-1 Forward_read_trimmed - must state PATH to the file
	-2 Reverse_read_trimmed - must state PATH to the file
	-r Flag to generate rnaSPAdes assembly output
	-t Flag to generate IDBA-tran assembly output
	-T Flag to generate TRINITY assembly output
	-B Flag to use BWA and BOWTIE2 - mapping
	-C Flag to use BLAST and CREST - classification
	-f Flag to produce final output files
	-s Flag to include assembly contig/scaffold sequences in the final output files
	-h Display this help and exit"

# Set default options
Forward_read_trimmed=''
Reverse_read_trimmed=''
rnaSPAdes='false'
IDBA_TRAN='false'
TRINITY='false'
MAP='false'
CLASSIFICATION='false'
FINAL='false'
READS='false'

# Set specified options
while getopts ':1:2:SITBCfrh' opt; do
 	case "${opt}" in
		1) Forward_read_trimmed="${OPTARG}" ;;
		2) Reverse_read_trimmed="${OPTARG}" ;;
		r) rnaSPAdes='true' ;;
		t) IDBA_TRAN='true' ;;
		T) TRINITY='true' ;;
		B) MAP='true' ;;
		C) CLASSIFICATION='true' ;;
		f) FINAL='true' ;;
		s) READS='true' ;;
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
mkdir PIPELINE_RNA
cd PIPELINE_RNA/

mkdir SORTMERNA
sortmerna --ref /hdd1/databases/sortmerna_silva_databases/silva-bac-16s-id90.fasta --ref /hdd1/databases/sortmerna_silva_databases/silva-arc-16s-id95.fasta --ref /hdd1/databases/sortmerna_silva_databases/silva-euk-18s-id95.fasta --reads $Forward_read_trimmed --reads $Reverse_read_trimmed --paired_in -other -fastx 1 -num_alignments 1 -v -workdir SORTMERNA

if [[ $rnaSPAdes == 'true' ]] ; then
	mkdir rnaSPAdes
	rnaspades.py --12 ./SORTMERNA/out/aligned.fq -o rnaSPAdes

	if [[ $MAP == 'true' ]] ; then

		echo "starting bwa index"

		bwa index -p bwa_index ./rnaSPAdes/transcripts.fasta

		echo "bwa index complete."
		echo "starting bwa."

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
		echo -e "mapped\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "mapped\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo "bwa complete."

		echo "starting bowtie2 index"

		bowtie2-build -f ./rnaSPAdes/transcripts.fasta bowtie_index

		echo "bowtie2 index complete."
		echo "starting bowtie2."

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
		echo -e "mapped\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "mapped\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo "bowtie2 complete"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ rnaSPAdes/

	else

		echo "no BWA and BOWTIE2 output generated"

	fi

	if [[ $CLASSIFICATION == 'true' ]] ; then

		blastn -db /hdd1/databases/nt_database_feb_2020_indexed/nt -query ./rnaSPAdes/transcripts.fasta -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_nt

		echo "BLAST nt output complete."

		blastn -db /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta -query ./rnaSPAdes/transcripts.fasta -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_SILVA

		echo "BLAST SILVA output complete."

		mkdir rnaSPAdes/BLAST
		mv BLAST_output_nt BLAST_output_SILVA rnaSPAdes/BLAST/

		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./rnaSPAdes/transcripts.fasta -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		rm BLAST_output.xml
		mv CREST rnaSPAdes/

		mkdir rnaSPAdes/CLASSIFICATION
		mv rnaSPAdes/BLAST rnaSPAdes/CREST rnaSPAdes/CLASSIFICATION/

		echo "CREST output complete."

	else

		echo "no BLAST or CREST output generated"

	fi


	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./rnaSPAdes/transcripts.fasta > ./rnaSPAdes/fasta_to_tabbed.txt
			mergeFilesOnColumn.pl ./rnaSPAdes/BWA/merge_input_mapped_bwa.txt ./rnaSPAdes/fasta_to_tabbed.txt 2 1 > ./rnaSPAdes/merged_original_bwa.txt
			cut -f1,2,4 ./rnaSPAdes/merged_original_bwa.txt > ./rnaSPAdes/important_column_3_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./rnaSPAdes/important_column_3_bwa.txt > ./rnaSPAdes/final_order_bwa.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./rnaSPAdes/final_bwa_merge_ready.txt && cat ./rnaSPAdes/final_order_bwa.txt >> ./rnaSPAdes/final_bwa_merge_ready.txt

			rm ./rnaSPAdes/merged_original_bwa.txt ./rnaSPAdes/important_column_3_bwa.txt ./rnaSPAdes/final_order_bwa.txt

			mergeFilesOnColumn.pl ./rnaSPAdes/BOWTIE2/merge_input_mapped_bowtie2.txt ./rnaSPAdes/fasta_to_tabbed.txt 2 1 > ./rnaSPAdes/merged_original_bowtie2.txt
			cut -f1,2,4 ./rnaSPAdes/merged_original_bowtie2.txt > ./rnaSPAdes/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./rnaSPAdes/important_column_3_bowtie2.txt > ./rnaSPAdes/final_order_bowtie2.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./rnaSPAdes/final_bowtie2_merge_ready.txt && cat ./rnaSPAdes/final_order_bowtie2.txt >> ./rnaSPAdes/final_bowtie2_merge_ready.txt

			rm ./rnaSPAdes/fasta_to_tabbed.txt ./rnaSPAdes/merged_original_bowtie2.txt ./rnaSPAdes/important_column_3_bowtie2.txt ./rnaSPAdes/final_order_bowtie2.txt

			echo "assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."

		else
			# Prepares to merge with BLAST of CREST - does not include assembly sequence
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./rnaSPAdes/BWA/merge_input_mapped_bwa.txt > ./rnaSPAdes/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./rnaSPAdes/temp.txt > ./rnaSPAdes/temp2.txt
			sed '1d' ./rnaSPAdes/temp2.txt > ./rnaSPAdes/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./rnaSPAdes/final_bwa_merge_ready.txt && cat ./rnaSPAdes/temp3.txt >> ./rnaSPAdes/final_bwa_merge_ready.txt
 			rm ./rnaSPAdes/temp*.txt

			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1}' ./rnaSPAdes/BOWTIE2/merge_input_mapped_bowtie2.txt > ./rnaSPAdes/temp.txt
			sed 's/$/\tno assembly sequence was included/' ./rnaSPAdes/temp.txt > ./rnaSPAdes/temp2.txt
			sed '1d' ./rnaSPAdes/temp2.txt > ./rnaSPAdes/temp3.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./rnaSPAdes/final_bowtie2_merge_ready.txt && cat ./rnaSPAdes/temp3.txt >> ./rnaSPAdes/final_bowtie2_merge_ready.txt
			rm ./rnaSPAdes/temp*.txt

			echo "no assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."
		fi


		# Assign BLAST taxonomy - nt
		assign_taxonomy_NCBI_staxids.sh -b ./rnaSPAdes/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv _with_taxonomy.txt ./rnaSPAdes/BLAST_output_nt_with_taxonomy.txt
		sed '1d' ./rnaSPAdes/BLAST_output_nt_with_taxonomy.txt > ./rnaSPAdes/BLAST_output_nt_with_taxonomy_noheader.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./rnaSPAdes/BLAST_output_nt_with_taxonomy_merge.txt && cat ./rnaSPAdes/BLAST_output_nt_with_taxonomy_noheader.txt >> ./rnaSPAdes/BLAST_output_nt_with_taxonomy_merge.txt

		# Assign BLAST taxonomy - SILVA
		assign_taxonomy_NCBI_staxids.sh -b ./rnaSPAdes/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv _with_taxonomy.txt ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy.txt
		sed '1d' ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy.txt > ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_noheader.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_merge.txt

		rm ./rnaSPAdes/BLAST_output_nt_with_taxonomy.txt ./rnaSPAdes/BLAST_output_nt_with_taxonomy_noheader.txt ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy.txt ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_noheader.txt

		echo "taxonomy has been assigned to BLAST files."
		echo "BLAST files are ready."


		# Prepare CREST file - remove 1st/3rd column, add a header
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./rnaSPAdes/CLASSIFICATION/CREST/otus.csv ./rnaSPAdes/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./rnaSPAdes/CLASSIFICATION/CREST/CREST_output.txt > ./rnaSPAdes/CREST_seperated.txt
		cut -f2,4 ./rnaSPAdes/CREST_seperated.txt > ./rnaSPAdes/CREST_header.txt
		sed '1d' ./rnaSPAdes/CREST_header.txt > ./rnaSPAdes/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./rnaSPAdes/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./rnaSPAdes/
		sed '1d' ./rnaSPAdes/CREST_tax_ready_with_taxonomy.txt > ./rnaSPAdes/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./rnaSPAdes/CREST_merge.txt && cat ./rnaSPAdes/CREST_tax_ready_with_taxonomy_noheader.txt >> ./rnaSPAdes/CREST_merge.txt

		rm ./rnaSPAdes/CREST_seperated.txt ./rnaSPAdes/CREST_header.txt ./rnaSPAdes/CREST_tax_ready.txt ./rnaSPAdes/CREST_tax_ready_with_taxonomy.txt ./rnaSPAdes/CREST_tax_ready_with_taxonomy_noheader.txt

		echo "CREST file is ready."


		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./rnaSPAdes/CREST_merge.txt ./rnaSPAdes/final_bwa_merge_ready.txt ./rnaSPAdes/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./rnaSPAdes/CREST_merge.txt ./rnaSPAdes/final_bowtie2_merge_ready.txt ./rnaSPAdes/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_merge.txt ./rnaSPAdes/final_bowtie2_merge_ready.txt ./rnaSPAdes/BLAST_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./rnaSPAdes/BLAST_output_nt_with_taxonomy_merge.txt ./rnaSPAdes/final_bowtie2_merge_ready.txt ./rnaSPAdes/BLAST_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_merge.txt ./rnaSPAdes/final_bwa_merge_ready.txt ./rnaSPAdes/BLAST_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./rnaSPAdes/BLAST_output_nt_with_taxonomy_merge.txt ./rnaSPAdes/final_bwa_merge_ready.txt ./rnaSPAdes/BLAST_nt_BWA.txt


		# Move all files - easy to find!!!
		mkdir ./rnaSPAdes/MERGE_FILES
		mv ./rnaSPAdes/final_bwa_merge_ready.txt ./rnaSPAdes/final_bowtie2_merge_ready.txt ./rnaSPAdes/CREST_merge.txt ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_merge.txt ./rnaSPAdes/BLAST_output_nt_with_taxonomy_merge.txt ./rnaSPAdes/MERGE_FILES/


		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./rnaSPAdes/CREST_BWA.txt > ./rnaSPAdes/new.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt
		sed 's/NODE_//g' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt && sed 's/length_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/cov_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt
		sed '1d' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-20 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/CREST_BWA_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/CREST_BWA_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/CREST_BWA.txt

		sed 's/_/\t/2' ./rnaSPAdes/CREST_BOWTIE2.txt > ./rnaSPAdes/new.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt
		sed 's/NODE_//g' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt && sed 's/length_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/cov_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt
		sed '1d' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-20 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/CREST_BOWTIE2_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/CREST_BOWTIE2_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/CREST_BOWTIE2.txt


		# BLAST nt files
		sed '1d' ./rnaSPAdes/BLAST_nt_BWA.txt > ./rnaSPAdes/new.txt
		sed 's/_/\t/2' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt
		sed 's/NODE_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/length_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt && sed 's/cov_//g' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-31 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/BLAST_nt_BWA_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/BLAST_nt_BWA_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/BLAST_nt_BWA.txt

		sed '1d' ./rnaSPAdes/BLAST_nt_BOWTIE2.txt > ./rnaSPAdes/new.txt
		sed 's/_/\t/2' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt
		sed 's/NODE_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/length_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt && sed 's/cov_//g' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-31 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/BLAST_nt_BOWTIE2_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/BLAST_nt_BOWTIE2_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/BLAST_nt_BOWTIE2.txt


		# BLAST SILVA files
		sed '1d' ./rnaSPAdes/BLAST_SILVA_BWA.txt > ./rnaSPAdes/new.txt
		sed 's/_/\t/2' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt
		sed 's/NODE_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/length_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt && sed 's/cov_//g' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-31 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/BLAST_SILVA_BWA_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/BLAST_SILVA_BWA_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/BLAST_SILVA_BWA.txt

		sed '1d' ./rnaSPAdes/BLAST_SILVA_BOWTIE2.txt > ./rnaSPAdes/new.txt
		sed 's/_/\t/2' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt
		sed 's/NODE_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/length_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt && sed 's/cov_//g' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-31 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/BLAST_SILVA_BOWTIE2_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/BLAST_SILVA_BOWTIE2_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/BLAST_SILVA_BOWTIE2.txt

		mkdir ./rnaSPAdes/FINAL_FILES_rnaSPAdes
		mv ./rnaSPAdes/CREST_BWA_final.txt ./rnaSPAdes/CREST_BOWTIE2_final.txt ./rnaSPAdes/BLAST_SILVA_BOWTIE2_final.txt ./rnaSPAdes/BLAST_nt_BOWTIE2_final.txt ./rnaSPAdes/BLAST_SILVA_BWA_final.txt ./rnaSPAdes/BLAST_nt_BWA_final.txt ./rnaSPAdes/FINAL_FILES_rnaSPAdes/

		echo "final rnaSPAdes output generated."

	else

		echo "no FINAL merge generated."

	fi


else
	echo "no rnaSPAdes output generated."
fi



if [[ $IDBA_TRAN == 'true' ]] ; then
	fq2fa ./SORTMERNA/out/aligned.fq ./SORTMERNA/out/aligned.fa
	idba_tran --num_threads 16 --pre_correction -r ./SORTMERNA/out/aligned.fa -o IDBA_TRAN

		if [[ $MAP == 'true' ]] ; then

		echo "starting bwa index"

		bwa index -p bwa_index ./IDBA_TRAN/contig.fa

		echo "bwa index complete."
		echo "starting bwa."

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
		echo -e "mapped\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "mapped\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo "bwa complete."

		echo "starting bowtie2 index"

		bowtie2-build -f ./IDBA_TRAN/contig.fa bowtie_index

		echo "bowtie2 index complete."
		echo "starting bowtie2."

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
		echo -e "mapped\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "mapped\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo "bowtie2 complete"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ IDBA_TRAN/

	else

		echo "no BWA and BOWTIE2 output generated"

	fi


	if [[ $CLASSIFICATION == 'true' ]] ; then

		blastn -db /hdd1/databases/nt_database_feb_2020_indexed/nt -query ./IDBA_TRAN/contig.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_nt

		echo "BLAST nt output complete."

		blastn -db /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta -query ./IDBA_TRAN/contig.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_SILVA

		echo "BLAST SILVA output complete."

		mkdir IDBA_TRAN/BLAST
		mv BLAST_output_nt BLAST_output_SILVA IDBA_TRAN/BLAST/

		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./IDBA_TRAN/contig.fa -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		rm BLAST_output.xml
		mv CREST IDBA_TRAN/

		mkdir IDBA_TRAN/CLASSIFICATION
		mv IDBA_TRAN/BLAST IDBA_TRAN/CREST IDBA_TRAN/CLASSIFICATION/

	else

		echo "no BLAST or CREST output generated"

	fi


	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./IDBA_TRAN/contig.fa > ./IDBA_TRAN/fasta_to_tabbed.txt
			sed 's/_/\t/2' ./IDBA_TRAN/fasta_to_tabbed.txt > ./IDBA_TRAN/fasta_to_tabbed_tab.txt
			sed 's/_/\t/3' ./IDBA_TRAN/fasta_to_tabbed_tab.txt > ./IDBA_TRAN/fasta_to_tabbed_tab_2.txt
			sed 's/ /\t/g' ./IDBA_TRAN/fasta_to_tabbed_tab_2.txt > ./IDBA_TRAN/cut_tab_ready.txt
			cut -f1,3,5,6 ./IDBA_TRAN/cut_tab_ready.txt > ./IDBA_TRAN/important_columns.txt
			rm ./IDBA_TRAN/fasta_to_tabbed.txt ./IDBA_TRAN/fasta_to_tabbed_tab.txt ./IDBA_TRAN/fasta_to_tabbed_tab_2.txt ./IDBA_TRAN/cut_tab_ready.txt

			mergeFilesOnColumn.pl ./IDBA_TRAN/BWA/merge_input_mapped_bwa.txt ./IDBA_TRAN/important_columns.txt 2 1 > ./IDBA_TRAN/merged_original_bwa.txt
			cut -f1,2,4,5,6 ./IDBA_TRAN/merged_original_bwa.txt > ./IDBA_TRAN/important_columns_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $4, $1, $5}' ./IDBA_TRAN/important_columns_bwa.txt > ./IDBA_TRAN/final_order_bwa.txt
			echo -e "contig_number\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/final_bwa_merge_ready.txt && cat ./IDBA_TRAN/final_order_bwa.txt >> ./IDBA_TRAN/final_bwa_merge_ready.txt
			rm ./IDBA_TRAN/merged_original_bwa.txt ./IDBA_TRAN/important_columns_bwa.txt ./IDBA_TRAN/final_order_bwa.txt


			mergeFilesOnColumn.pl ./IDBA_TRAN/BOWTIE2/merge_input_mapped_bowtie2.txt ./IDBA_TRAN/important_columns.txt 2 1 > ./IDBA_TRAN/merged_original_bowtie2.txt
			cut -f1,2,4,5,6 ./IDBA_TRAN/merged_original_bowtie2.txt > ./IDBA_TRAN/important_columns_bowtie.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $4, $1, $5}' ./IDBA_TRAN/important_columns_bowtie.txt > ./IDBA_TRAN/final_order_bowtie2.txt
			echo -e "contig_number\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/final_bowtie2_merge_ready.txt && cat ./IDBA_TRAN/final_order_bowtie2.txt >> ./IDBA_TRAN/final_bowtie2_merge_ready.txt
			rm ./IDBA_TRAN/merged_original_bowtie2.txt ./IDBA_TRAN/important_columns_bowtie.txt ./IDBA_TRAN/final_order_bowtie2.txt ./IDBA_TRAN/important_columns.txt

			echo "assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_TRAN/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./IDBA_TRAN/BLAST_output_nt_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt && cat ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_noheader.txt >> ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./IDBA_TRAN/BLAST_output_nt_with_taxonomy.txt ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_noheader.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."

			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA_TRAN/CLASSIFICATION/CREST/otus.csv ./IDBA_TRAN/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./IDBA_TRAN/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA_TRAN/CREST_seperated.txt
			cut -f2,4 ./IDBA_TRAN/CREST_seperated.txt > ./IDBA_TRAN/CREST_header.txt
			sed '1d' ./IDBA_TRAN/CREST_header.txt > ./IDBA_TRAN/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./IDBA_TRAN/
			sed '1d' ./IDBA_TRAN/CREST_tax_ready_with_taxonomy.txt > ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt

			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/CREST_merge.txt && cat ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA_TRAN/CREST_merge.txt

			rm ./IDBA_TRAN/CREST_seperated.txt ./IDBA_TRAN/CREST_header.txt ./IDBA_TRAN/CREST_tax_ready.txt ./IDBA_TRAN/CREST_tax_ready_with_taxonomy.txt ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./IDBA_TRAN/MERGE_FILES
			mv ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_TRAN/MERGE_FILES/

			# Edit contig-* in all the final files
			# CREST files
			sed '1d' ./IDBA_TRAN/CREST_BWA.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-20 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/CREST_BWA_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/CREST_BWA_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/CREST_BWA.txt

			sed '1d' ./IDBA_TRAN/CREST_BOWTIE2.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-20 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence"  > ./IDBA_TRAN/CREST_BOWTIE2_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/CREST_BOWTIE2_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./IDBA_TRAN/BLAST_nt_BWA.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-31 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_BWA_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_nt_BWA_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_nt_BWA.txt

			sed '1d' ./IDBA_TRAN/BLAST_nt_BOWTIE2.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-31 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_BOWTIE2_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_nt_BOWTIE2_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_nt_BOWTIE2.txt

			# BLAST SILVA files
			sed '1d' ./IDBA_TRAN/BLAST_SILVA_BWA.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-31 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_BWA_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_SILVA_BWA_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_SILVA_BWA.txt

			sed '1d' ./IDBA_TRAN/BLAST_SILVA_BOWTIE2.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-31 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_BOWTIE2_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_SILVA_BOWTIE2_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_SILVA_BOWTIE2.txt

			mkdir ./IDBA_TRAN/FINAL_FILES_IDBA_TRAN
			mv ./IDBA_TRAN/CREST_BWA_final.txt ./IDBA_TRAN/CREST_BOWTIE2_final.txt ./IDBA_TRAN/BLAST_SILVA_BOWTIE2_final.txt ./IDBA_TRAN/BLAST_nt_BOWTIE2_final.txt ./IDBA_TRAN/BLAST_SILVA_BWA_final.txt ./IDBA_TRAN/BLAST_nt_BWA_final.txt ./IDBA_TRAN/FINAL_FILES_IDBA_TRAN/

			echo "final IDBA_TRAN output generated."

		else
			# Prepares to merge with BLAST of CREST - does not include assembly sequence
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

			echo "no assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_TRAN/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./IDBA_TRAN/BLAST_output_nt_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt && cat ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_noheader.txt >> ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./IDBA_TRAN/BLAST_output_nt_with_taxonomy.txt ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_noheader.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."


			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA_TRAN/CLASSIFICATION/CREST/otus.csv ./IDBA_TRAN/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./IDBA_TRAN/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA_TRAN/CREST_seperated.txt
			cut -f2,4 ./IDBA_TRAN/CREST_seperated.txt > ./IDBA_TRAN/CREST_header.txt
			sed '1d' ./IDBA_TRAN/CREST_header.txt > ./IDBA_TRAN/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./IDBA_TRAN/
			sed '1d' ./IDBA_TRAN/CREST_tax_ready_with_taxonomy.txt > ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt

			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/CREST_merge.txt && cat ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA_TRAN/CREST_merge.txt

			rm ./IDBA_TRAN/CREST_seperated.txt ./IDBA_TRAN/CREST_header.txt ./IDBA_TRAN/CREST_tax_ready.txt ./IDBA_TRAN/CREST_tax_ready_with_taxonomy.txt ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./IDBA_TRAN/MERGE_FILES
			mv ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_TRAN/MERGE_FILES/


			# Edit contig-* in all the final files
			# CREST files
			sed '1d' ./IDBA_TRAN/CREST_BWA.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-18 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/CREST_BWA_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/CREST_BWA_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/CREST_BWA.txt

			sed '1d' ./IDBA_TRAN/CREST_BOWTIE2.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-18 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/CREST_BOWTIE2_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/CREST_BOWTIE2_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./IDBA_TRAN/BLAST_nt_BWA.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-29 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_BWA_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_nt_BWA_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_nt_BWA.txt

			sed '1d' ./IDBA_TRAN/BLAST_nt_BOWTIE2.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-29 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_BOWTIE2_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_nt_BOWTIE2_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_nt_BOWTIE2.txt

			# BLAST SILVA files
			sed '1d' ./IDBA_TRAN/BLAST_SILVA_BWA.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-29 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_BWA_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_SILVA_BWA_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_SILVA_BWA.txt

			sed '1d' ./IDBA_TRAN/BLAST_SILVA_BOWTIE2.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-29 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_BOWTIE2_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_SILVA_BOWTIE2_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_SILVA_BOWTIE2.txt

			mkdir ./IDBA_TRAN/FINAL_FILES_IDBA_TRAN
			mv ./IDBA_TRAN/CREST_BWA_final.txt ./IDBA_TRAN/CREST_BOWTIE2_final.txt ./IDBA_TRAN/BLAST_SILVA_BOWTIE2_final.txt ./IDBA_TRAN/BLAST_nt_BOWTIE2_final.txt ./IDBA_TRAN/BLAST_SILVA_BWA_final.txt ./IDBA_TRAN/BLAST_nt_BWA_final.txt ./IDBA_TRAN/FINAL_FILES_IDBA_TRAN/

			echo "final IDBA_TRAN output generated."

		fi

	else

		echo "no FINAL merge generated."

	fi

else
	echo "no IDBA_TRAN output generated."
fi



if [[ $TRINITY == 'true' ]] ; then
	deinterleave_fastq_reads.sh < ./SORTMERNA/out/aligned.fq ./SORTMERNA/out/aligned_1.fq ./SORTMERNA/out/aligned_2.fq
	Trinity --seqType fq --max_memory 64G --left ./SORTMERNA/out/aligned_1.fq --right ./SORTMERNA/out/aligned_2.fq --CPU 16 --output TRINITY
	cat ./TRINITY/Trinity.fasta | sed 's/ len/_len/g' > ./TRINITY/Trinity_with_length.fasta

	if [[ $MAP == 'true' ]] ; then

		echo "starting bwa index"

		bwa index -p bwa_index ./TRINITY/Trinity_with_length.fasta

		echo "bwa index complete."
		echo "starting bwa."

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
		echo -e "mapped\tcontig_number" > merge_input_mapped_bwa.txt && cat out_mapped_bwa.txt >> merge_input_mapped_bwa.txt
		echo -e "mapped\tcontig_number" > merge_input_unmapped_bwa.txt && cat out_unmappped_bwa.txt >> merge_input_unmapped_bwa.txt

		rm mapped_reads_bwa.sam unmapped_reads_bwa.sam mapped_reads_bwa.txt unmapped_reads_bwa.txt mapped_column3_reads_bwa.txt unmapped_column3_reads_bwa.txt sorted_mapped_column3_reads_bwa.txt sorted_unmapped_column3_reads_bwa.txt aligned_mapped_bwa.txt aligned_unmapped_bwa.txt out_mapped_bwa.txt out_unmappped_bwa.txt

		echo "bwa complete."

		echo "starting bowtie2 index"

		bowtie2-build -f ./TRINITY/Trinity_with_length.fasta bowtie_index

		echo "bowtie2 index complete."
		echo "starting bowtie2."

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
		echo -e "mapped\tcontig_number" > merge_input_mapped_bowtie2.txt && cat out_mapped_bowtie.txt >> merge_input_mapped_bowtie2.txt
		echo -e "mapped\tcontig_number" > merge_input_unmapped_bowtie2.txt && cat out_unmappped_bowtie.txt >> merge_input_unmapped_bowtie2.txt

		rm mapped_reads_bowtie.sam unmapped_reads_bowtie.sam mapped_reads_bowtie.txt unmapped_reads_bowtie.txt mapped_column3_reads_bowtie.txt unmapped_column3_reads_bowtie.txt sorted_mapped_column3_reads_bowtie.txt sorted_unmapped_column3_reads_bowtie.txt aligned_mapped_bowtie.txt aligned_unmapped_bowtie.txt out_mapped_bowtie.txt out_unmappped_bowtie.txt

		echo "bowtie2 complete"

		mkdir BWA/
		mkdir BOWTIE2/

		mv bwa_output.sam merge_input_mapped_bwa.txt merge_input_unmapped_bwa.txt BWA/
		mv bowtie2_output.sam merge_input_mapped_bowtie2.txt merge_input_unmapped_bowtie2.txt BOWTIE2/

		mv BWA/ BOWTIE2/ TRINITY/

	else

		echo "no BWA and BOWTIE2 output generated"

	fi


	if [[ $CLASSIFICATION == 'true' ]] ; then

		blastn -db /hdd1/databases/nt_database_feb_2020_indexed/nt -query ./TRINITY/Trinity_with_length.fasta -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_nt

		echo "BLAST nt output complete."

		blastn -db /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta -query ./TRINITY/Trinity_with_length.fasta -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_SILVA

		echo "BLAST SILVA output complete."

		mkdir TRINITY/BLAST
		mv BLAST_output_nt BLAST_output_SILVA TRINITY/BLAST/

		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./TRINITY/Trinity_with_length.fasta -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		rm BLAST_output.xml
		mv CREST TRINITY/

		mkdir TRINITY/CLASSIFICATION
		mv TRINITY/BLAST TRINITY/CREST TRINITY/CLASSIFICATION/

	else

		echo "no BLAST or CREST output generated"

	fi


if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./TRINITY/Trinity_with_length.fasta > ./TRINITY/fasta_to_tabbed.txt
			sed 's/ /\t/1' ./TRINITY/fasta_to_tabbed.txt > ./TRINITY/tabbed.txt
			cut -f1,3 ./TRINITY/tabbed.txt > ./TRINITY/ready.txt


			mergeFilesOnColumn.pl ./TRINITY/BWA/merge_input_mapped_bwa.txt ./TRINITY/ready.txt 2 1 > ./TRINITY/merged_original_bwa.txt
			cut -f1,2,4 ./TRINITY/merged_original_bwa.txt > ./TRINITY/important_column_3_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./TRINITY/important_column_3_bwa.txt > ./TRINITY/final_order_bwa.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./TRINITY/final_bwa_merge_ready.txt && cat ./TRINITY/final_order_bwa.txt >> ./TRINITY/final_bwa_merge_ready.txt

			rm ./TRINITY/merged_original_bwa.txt ./TRINITY/important_column_3_bwa.txt ./TRINITY/final_order_bwa.txt

			mergeFilesOnColumn.pl ./TRINITY/BOWTIE2/merge_input_mapped_bowtie2.txt ./TRINITY/ready.txt 2 1 > ./TRINITY/merged_original_bowtie2.txt
			cut -f1,2,4 ./TRINITY/merged_original_bowtie2.txt > ./TRINITY/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./TRINITY/important_column_3_bowtie2.txt > ./TRINITY/final_order_bowtie2.txt
			echo -e "contig_number\tcounts\tassembly_sequence" > ./TRINITY/final_bowtie2_merge_ready.txt && cat ./TRINITY/final_order_bowtie2.txt >> ./TRINITY/final_bowtie2_merge_ready.txt

			rm ./TRINITY/fasta_to_tabbed.txt ./TRINITY/tabbed.txt ./TRINITY/ready.txt ./TRINITY/merged_original_bowtie2.txt ./TRINITY/important_column_3_bowtie2.txt ./TRINITY/final_order_bowtie2.txt

			echo "assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."

		else
			# Prepares to merge with BLAST of CREST - does not include assembly sequence
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

			echo "no assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."
		fi

		# Assign BLAST taxonomy - nt
		assign_taxonomy_NCBI_staxids.sh -b ./TRINITY/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv _with_taxonomy.txt ./TRINITY/BLAST_output_nt_with_taxonomy.txt
		sed '1d' ./TRINITY/BLAST_output_nt_with_taxonomy.txt > ./TRINITY/BLAST_output_nt_with_taxonomy_noheader.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRINITY/BLAST_output_nt_with_taxonomy_merge.txt && cat ./TRINITY/BLAST_output_nt_with_taxonomy_noheader.txt >> ./TRINITY/BLAST_output_nt_with_taxonomy_merge.txt

		# Assign BLAST taxonomy - SILVA
		assign_taxonomy_NCBI_staxids.sh -b ./TRINITY/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
		mv _with_taxonomy.txt ./TRINITY/BLAST_output_SILVA_with_taxonomy.txt
		sed '1d' ./TRINITY/BLAST_output_SILVA_with_taxonomy.txt > ./TRINITY/BLAST_output_SILVA_with_taxonomy_noheader.txt
		echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRINITY/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./TRINITY/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./TRINITY/BLAST_output_SILVA_with_taxonomy_merge.txt

		rm ./TRINITY/BLAST_output_nt_with_taxonomy.txt ./TRINITY/BLAST_output_nt_with_taxonomy_noheader.txt ./TRINITY/BLAST_output_SILVA_with_taxonomy.txt ./TRINITY/BLAST_output_SILVA_with_taxonomy_noheader.txt

		echo "taxonomy has been assigned to BLAST files."
		echo "BLAST files are ready."

		# Prepare CREST file - remove 1st/3rd column, add a header
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./TRINITY/CLASSIFICATION/CREST/otus.csv ./TRINITY/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./TRINITY/CLASSIFICATION/CREST/CREST_output.txt > ./TRINITY/CREST_seperated.txt
		cut -f2,4 ./TRINITY/CREST_seperated.txt > ./TRINITY/CREST_header.txt
		sed '1d' ./TRINITY/CREST_header.txt > ./TRINITY/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./TRINITY/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./TRINITY/
		sed '1d' ./TRINITY/CREST_tax_ready_with_taxonomy.txt > ./TRINITY/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRINITY/CREST_merge.txt && cat ./TRINITY/CREST_tax_ready_with_taxonomy_noheader.txt >> ./TRINITY/CREST_merge.txt

		rm ./TRINITY/CREST_seperated.txt ./TRINITY/CREST_header.txt ./TRINITY/CREST_tax_ready.txt ./TRINITY/CREST_tax_ready_with_taxonomy.txt ./TRINITY/CREST_tax_ready_with_taxonomy_noheader.txt

		echo "CREST file is ready."

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./TRINITY/CREST_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/CREST_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_SILVA_with_taxonomy_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/BLAST_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_nt_with_taxonomy_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/BLAST_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_SILVA_with_taxonomy_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/BLAST_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_nt_with_taxonomy_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/BLAST_nt_BWA.txt

		# Move all files - easy to find!!!
		mkdir ./TRINITY/MERGE_FILES
		mv ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/CREST_merge.txt ./TRINITY/BLAST_output_SILVA_with_taxonomy_merge.txt ./TRINITY/BLAST_output_nt_with_taxonomy_merge.txt ./TRINITY/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed '1d' ./TRINITY/CREST_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/CREST_BWA_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/CREST_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/CREST_BWA.txt

		sed '1d' ./TRINITY/CREST_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/CREST_BOWTIE2_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/CREST_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./TRINITY/BLAST_nt_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_nt_BWA_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/BLAST_nt_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_nt_BWA.txt

		sed '1d' ./TRINITY/BLAST_nt_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_nt_BOWTIE2_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/BLAST_nt_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./TRINITY/BLAST_SILVA_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_SILVA_BWA_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/BLAST_SILVA_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_SILVA_BWA.txt

		sed '1d' ./TRINITY/BLAST_SILVA_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_SILVA_BOWTIE2_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/BLAST_SILVA_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_SILVA_BOWTIE2.txt

		mkdir ./TRINITY/FINAL_FILES_TRINITY
		mv ./TRINITY/CREST_BWA_final.txt ./TRINITY/CREST_BOWTIE2_final.txt ./TRINITY/BLAST_SILVA_BOWTIE2_final.txt ./TRINITY/BLAST_nt_BOWTIE2_final.txt ./TRINITY/BLAST_SILVA_BWA_final.txt ./TRINITY/BLAST_nt_BWA_final.txt ./TRINITY/FINAL_FILES_TRINITY/

		echo "final TRINITY output generated."

	else

		echo "no FINAL merge generated."

	fi

else
	echo "no TRINITY output generated."
fi

mkdir ASSEMBLERS
mv rnaSPAdes IDBA_TRAN TRINITY ASSEMBLERS/

# Display runtime
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nSCRIPT RUNTIME: $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"

# Create log
) 2>&1 | tee PIPELINE_RNA_LOG.txt
mv PIPELINE_RNA_LOG.txt PIPELINE_RNA/

echo "all done."
