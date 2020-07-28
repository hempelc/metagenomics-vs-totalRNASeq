#!/bin/bash

# Need to have: assign_NCBI_staxids_to_CREST_v3.py, assign_taxonomy_NCBI_staxids.sh, deinterleave_fastq_reads.sh, LookupTaxonDetails3.py, merge_mapped_reads_and_contigs.py, mergeFilesOnColumn.pl and fasta_to_tab in your PATH
# In order for "assign_taxonomy_NCBI_staxids.sh" to work - MUST have .etetoolkit/taxa.sqlite in your HOME directory
cmd="$0 $@" # Make variable containing the used command to print command later
usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq> [aSMmIBCfrh]

Usage:
	-1 Forward reads trimmed - must state whole path to the file
	-2 Reverse reads trimmed - must state whole path to the file
	-a Main flag to indicate that all following flags should be used
	-S Flag to generate SPADES assembly output
	-M Flag to generate METASPADES assembly output
	-m Flag to generate MEGAHIT assembly output
	-I Flag to generate IDBA-UD assembly output
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
IDBA_UD='false'
MAP='false'
CLASSIFICATION='false'
FINAL='false'
READS='false'

# Set specified options
while getopts ':1:2:aSMmIBCfrh' opt; do
 	case "${opt}" in
		1) Forward_read_trimmed="${OPTARG}" ;;
		2) Reverse_read_trimmed="${OPTARG}" ;;
		a) SPADES='true'
		   METASPADES='true'
		   MEGAHIT='true'
		   IDBA_UD='true'
		   MAP='true'
		   CLASSIFICATION='true'
		   FINAL='true'
		   READS='true' ;;
		S) SPADES='true' ;;
		M) METASPADES='true' ;;
		m) MEGAHIT='true' ;;
		I) IDBA_UD='true' ;;
		B) MAP='true' ;;
		C) CLASSIFICATION='true' ;;
		f) FINAL='true' ;;
		r) READS='true' ;;
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

# Making output directory for this script and open bracket to later tell script
# to write everything that follows into a logfile
mkdir PIPELINE_DNA
cd PIPELINE_DNA/
(

# Define starting time of script for total runtime calculation
start=$(date +%s)
echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
echo -e "=================================================================\n\n"


# Output specified options
echo -e "======== OPTIONS ========\n"

echo -e "Forward reads were defined as $Forward_read_trimmed.\n"
echo -e "Reverse reads were defined as $Reverse_read_trimmed.\n"
echo -e "Script started as $cmd\n"

echo -e "======== START RUNNING SCRIPT ========\n"


######################### Beginning of actual pipeline ########################

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
	megahit --presets meta-large -t 16 -1 $Forward_read_trimmed -2 $Reverse_read_trimmed -o ./MEGAHIT

	if [[ $MAP == 'true' ]] ; then

		echo "starting bwa index"

		bwa index -p bwa_index ./MEGAHIT/final.contigs.fa

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

		bowtie2-build -f ./MEGAHIT/final.contigs.fa bowtie_index

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

		mv BWA/ BOWTIE2/ MEGAHIT/

	else

		echo "no BWA and BOWTIE2 output generated"

	fi


	if [[ $CLASSIFICATION == 'true' ]] ; then

		blastn -db /hdd1/databases/nt_database_feb_2020_indexed/nt -query ./MEGAHIT/final.contigs.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_nt

		echo "BLAST nt output complete."

		blastn -db /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta -query ./MEGAHIT/final.contigs.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_SILVA

		echo "BLAST SILVA output complete."

		mkdir MEGAHIT/BLAST
		mv BLAST_output_nt BLAST_output_SILVA MEGAHIT/BLAST/

		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./MEGAHIT/final.contigs.fa -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		rm BLAST_output.xml
		mv CREST MEGAHIT/

		mkdir MEGAHIT/CLASSIFICATION
		mv MEGAHIT/BLAST MEGAHIT/CREST MEGAHIT/CLASSIFICATION/

		echo "CREST output complete."

	else

		echo "no BLAST or CREST output generated"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			# Prepare assembly output - fasta to tabular
			fasta_to_tab ./MEGAHIT/final.contigs.fa > ./MEGAHIT/fasta_to_tabbed.txt
			sed 's/ /\t/g' ./MEGAHIT/fasta_to_tabbed.txt > ./MEGAHIT/fasta_to_tabbed_tab.txt
			cut -f1,4,5 ./MEGAHIT/fasta_to_tabbed_tab.txt > ./MEGAHIT/important_columns.txt
			sed 's/len=//g' ./MEGAHIT/important_columns.txt > ./MEGAHIT/tab_to_merge.txt

			# Merge BWA and assembly sequences
			mergeFilesOnColumn.pl ./MEGAHIT/BWA/merge_input_mapped_bwa.txt ./MEGAHIT/tab_to_merge.txt 2 1 > ./MEGAHIT/merged_original_bwa.txt
			cut -f1,2,4,5 ./MEGAHIT/merged_original_bwa.txt > ./MEGAHIT/important_column_4_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $1, $4}' ./MEGAHIT/important_column_4_bwa.txt > ./MEGAHIT/final_order_bwa.txt
			echo -e "contig_number\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/final_bwa_merge_ready.txt && cat ./MEGAHIT/final_order_bwa.txt >> ./MEGAHIT/final_bwa_merge_ready.txt

			rm ./MEGAHIT/fasta_to_tabbed.txt ./MEGAHIT/fasta_to_tabbed_tab.txt ./MEGAHIT/important_columns.txt ./MEGAHIT/merged_original_bwa.txt ./MEGAHIT/important_column_4_bwa.txt ./MEGAHIT/final_order_bwa.txt

			# Merge BOWTIE2 and assembly sequences
			mergeFilesOnColumn.pl ./MEGAHIT/BOWTIE2/merge_input_mapped_bowtie2.txt ./MEGAHIT/tab_to_merge.txt 2 1 > ./MEGAHIT/merged_original_bowtie2.txt
			cut -f1,2,4,5 ./MEGAHIT/merged_original_bowtie2.txt > ./MEGAHIT/important_column_3_bowtie2.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $1, $4}' ./MEGAHIT/important_column_3_bowtie2.txt > ./MEGAHIT/final_order_bowtie2.txt
			echo -e "contig_number\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/final_bowtie2_merge_ready.txt && cat ./MEGAHIT/final_order_bowtie2.txt >> ./MEGAHIT/final_bowtie2_merge_ready.txt

			rm ./MEGAHIT/tab_to_merge.txt ./MEGAHIT/merged_original_bowtie2.txt ./MEGAHIT/important_column_3_bowtie2.txt ./MEGAHIT/final_order_bowtie2.txt

			echo "assembly sequence was added."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./MEGAHIT/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./MEGAHIT/BLAST_output_nt_with_taxonomy.txt > ./MEGAHIT/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt && cat ./MEGAHIT/BLAST_output_nt_with_taxonomy_noheader.txt >> ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./MEGAHIT/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./MEGAHIT/BLAST_output_SILVA_with_taxonomy.txt > ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./MEGAHIT/BLAST_output_nt_with_taxonomy.txt ./MEGAHIT/BLAST_output_nt_with_taxonomy_noheader.txt ./MEGAHIT/BLAST_output_SILVA_with_taxonomy.txt ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."


			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./MEGAHIT/CLASSIFICATION/CREST/otus.csv ./MEGAHIT/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./MEGAHIT/CLASSIFICATION/CREST/CREST_output.txt > ./MEGAHIT/CREST_seperated.txt
			cut -f2,4 ./MEGAHIT/CREST_seperated.txt > ./MEGAHIT/CREST_header.txt
			sed '1d' ./MEGAHIT/CREST_header.txt > ./MEGAHIT/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/
			sed '1d' ./MEGAHIT/CREST_tax_ready_with_taxonomy.txt > ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/CREST_merge.txt && cat ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt >> ./MEGAHIT/CREST_merge.txt

			rm ./MEGAHIT/CREST_seperated.txt ./MEGAHIT/CREST_header.txt ./MEGAHIT/CREST_tax_ready.txt ./MEGAHIT/CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./MEGAHIT/CREST_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/CREST_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./MEGAHIT/MERGE_FILES
			mv ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/CREST_merge.txt ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/MERGE_FILES/

			# Edit k*_ - in all files
			# CREST files
			sed '1d' ./MEGAHIT/CREST_BWA.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-19 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/CREST_BWA_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/CREST_BWA_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/CREST_BWA.txt

			sed '1d' ./MEGAHIT/CREST_BOWTIE2.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-19 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/CREST_BOWTIE2_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/CREST_BOWTIE2_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./MEGAHIT/BLAST_nt_BWA.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-30 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_BWA_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_nt_BWA_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_nt_BWA.txt

			sed '1d' ./MEGAHIT/BLAST_nt_BOWTIE2.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-30 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_BOWTIE2_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_nt_BOWTIE2_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_nt_BOWTIE2.txt


			# BLAST SILVA files
			sed '1d' ./MEGAHIT/BLAST_SILVA_BWA.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-30 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_BWA_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_SILVA_BWA_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_SILVA_BWA.txt

			sed '1d' ./MEGAHIT/BLAST_SILVA_BOWTIE2.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-30 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_BOWTIE2_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_SILVA_BOWTIE2_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_SILVA_BOWTIE2.txt

			mkdir ./MEGAHIT/FINAL_FILES_MEGAHIT
			mv ./MEGAHIT/CREST_BWA_final.txt ./MEGAHIT/CREST_BOWTIE2_final.txt ./MEGAHIT/BLAST_SILVA_BOWTIE2_final.txt ./MEGAHIT/BLAST_nt_BOWTIE2_final.txt ./MEGAHIT/BLAST_SILVA_BWA_final.txt ./MEGAHIT/BLAST_nt_BWA_final.txt ./MEGAHIT/FINAL_FILES_MEGAHIT/

		else

			# Prepares to merge with BLAST of CREST - does not include assembly sequence
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

			echo "no assembly sequence added."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./MEGAHIT/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./MEGAHIT/BLAST_output_nt_with_taxonomy.txt > ./MEGAHIT/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt && cat ./MEGAHIT/BLAST_output_nt_with_taxonomy_noheader.txt >> ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./MEGAHIT/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./MEGAHIT/BLAST_output_SILVA_with_taxonomy.txt > ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./MEGAHIT/BLAST_output_nt_with_taxonomy.txt ./MEGAHIT/BLAST_output_nt_with_taxonomy_noheader.txt ./MEGAHIT/BLAST_output_SILVA_with_taxonomy.txt ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."


			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./MEGAHIT/CLASSIFICATION/CREST/otus.csv ./MEGAHIT/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./MEGAHIT/CLASSIFICATION/CREST/CREST_output.txt > ./MEGAHIT/CREST_seperated.txt
			cut -f2,4 ./MEGAHIT/CREST_seperated.txt > ./MEGAHIT/CREST_header.txt
			sed '1d' ./MEGAHIT/CREST_header.txt > ./MEGAHIT/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/
			sed '1d' ./MEGAHIT/CREST_tax_ready_with_taxonomy.txt > ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/CREST_merge.txt && cat ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt >> ./MEGAHIT/CREST_merge.txt

			rm ./MEGAHIT/CREST_seperated.txt ./MEGAHIT/CREST_header.txt ./MEGAHIT/CREST_tax_ready.txt ./MEGAHIT/CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."


			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./MEGAHIT/CREST_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/CREST_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./MEGAHIT/MERGE_FILES
			mv ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/CREST_merge.txt ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/MERGE_FILES/


			# CREST files
			sed '1d' ./MEGAHIT/CREST_BWA.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-18 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/CREST_BWA_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/CREST_BWA_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/CREST_BWA.txt

			sed '1d' ./MEGAHIT/CREST_BOWTIE2.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-18 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/CREST_BOWTIE2_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/CREST_BOWTIE2_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/CREST_BOWTIE2.txt


			# BLAST nt files
			sed '1d' ./MEGAHIT/BLAST_nt_BWA.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-29 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_BWA_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_nt_BWA_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_nt_BWA.txt

			sed '1d' ./MEGAHIT/BLAST_nt_BOWTIE2.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-29 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_BOWTIE2_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_nt_BOWTIE2_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_nt_BOWTIE2.txt


			# BLAST SILVA files
			sed '1d' ./MEGAHIT/BLAST_SILVA_BWA.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-29 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_BWA_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_SILVA_BWA_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_SILVA_BWA.txt

			sed '1d' ./MEGAHIT/BLAST_SILVA_BOWTIE2.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-29 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_BOWTIE2_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_SILVA_BOWTIE2_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_SILVA_BOWTIE2.txt

			mkdir ./MEGAHIT/FINAL_FILES_MEGAHIT
			mv ./MEGAHIT/CREST_BWA_final.txt ./MEGAHIT/CREST_BOWTIE2_final.txt ./MEGAHIT/BLAST_SILVA_BOWTIE2_final.txt ./MEGAHIT/BLAST_nt_BOWTIE2_final.txt ./MEGAHIT/BLAST_SILVA_BWA_final.txt ./MEGAHIT/BLAST_nt_BWA_final.txt ./MEGAHIT/megahit_out/FINAL_FILES_MEGAHIT/

			echo "final MEGAHIT output generated."

		fi

	else

		echo "no FINAL merge generated."

	fi

else
	echo "no MEGAHIT output generated."
fi


if [[ $IDBA_UD == 'true' ]] ; then
	fq2fa --merge --filter $Forward_read_trimmed $Reverse_read_trimmed new_input.fa
	idba_ud --num_threads 16 --pre_correction -r ./new_input.fa -o IDBA_UD

		if [[ $MAP == 'true' ]] ; then

		echo "starting bwa index"

		bwa index -p bwa_index ./IDBA_UD/contig.fa

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

		bowtie2-build -f ./IDBA_UD/contig.fa bowtie_index

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

		mv BWA/ BOWTIE2/ IDBA_UD/

	else

		echo "no BWA and BOWTIE2 output generated"

	fi


	if [[ $CLASSIFICATION == 'true' ]] ; then

		blastn -db /hdd1/databases/nt_database_feb_2020_indexed/nt -query ./IDBA_UD/contig.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_nt

		echo "BLAST nt output complete."

		blastn -db /hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta -query ./IDBA_UD/contig.fa -outfmt "6 std staxids" -max_target_seqs 5 -num_threads 16 -out BLAST_output_SILVA

		echo "BLAST SILVA output complete."

		mkdir IDBA_UD/BLAST
		mv BLAST_output_nt BLAST_output_SILVA IDBA_UD/BLAST/

		blastn -db ~/programs/CREST/LCAClassifier/parts/flatdb/silvamod/silvamod128.fasta -query ./IDBA_UD/contig.fa -outfmt "5" -max_target_seqs 5 -num_threads 16 -out BLAST_output.xml
		classify BLAST_output.xml -o CREST

		rm BLAST_output.xml
		mv CREST IDBA_UD/

		mkdir IDBA_UD/CLASSIFICATION
		mv IDBA_UD/BLAST IDBA_UD/CREST IDBA_UD/CLASSIFICATION/

		echo "CREST output complete."

	else

		echo "no BLAST or CREST output generated"

	fi

	if [[ $FINAL == 'true' ]]; then
		if [[ $READS == 'true' ]]; then
			# Adding full assembly sequence to BWA or BOWTIE2 and prepares to merge with BLAST or CREST
			fasta_to_tab ./IDBA_UD/contig.fa > ./IDBA_UD/fasta_to_tabbed.txt
			sed 's/_/\t/2' ./IDBA_UD/fasta_to_tabbed.txt > ./IDBA_UD/fasta_to_tabbed_tab.txt
			sed 's/_/\t/3' ./IDBA_UD/fasta_to_tabbed_tab.txt > ./IDBA_UD/fasta_to_tabbed_tab_2.txt
			sed 's/ /\t/g' ./IDBA_UD/fasta_to_tabbed_tab_2.txt > ./IDBA_UD/cut_tab_ready.txt
			cut -f1,3,5,6 ./IDBA_UD/cut_tab_ready.txt > ./IDBA_UD/important_columns.txt
			rm ./IDBA_UD/fasta_to_tabbed.txt ./IDBA_UD/fasta_to_tabbed_tab.txt ./IDBA_UD/fasta_to_tabbed_tab_2.txt ./IDBA_UD/cut_tab_ready.txt

			mergeFilesOnColumn.pl ./IDBA_UD/BWA/merge_input_mapped_bwa.txt ./IDBA_UD/important_columns.txt 2 1 > ./IDBA_UD/merged_original_bwa.txt
			cut -f1,2,4,5,6 ./IDBA_UD/merged_original_bwa.txt > ./IDBA_UD/important_columns_bwa.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $4, $1, $5}' ./IDBA_UD/important_columns_bwa.txt > ./IDBA_UD/final_order_bwa.txt
			echo -e "contig_number\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/final_bwa_merge_ready.txt && cat ./IDBA_UD/final_order_bwa.txt >> ./IDBA_UD/final_bwa_merge_ready.txt
			rm ./IDBA_UD/merged_original_bwa.txt ./IDBA_UD/important_columns_bwa.txt ./IDBA_UD/final_order_bwa.txt


			mergeFilesOnColumn.pl ./IDBA_UD/BOWTIE2/merge_input_mapped_bowtie2.txt ./IDBA_UD/important_columns.txt 2 1 > ./IDBA_UD/merged_original_bowtie2.txt
			cut -f1,2,4,5,6 ./IDBA_UD/merged_original_bowtie2.txt > ./IDBA_UD/important_columns_bowtie.txt
			awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $4, $1, $5}' ./IDBA_UD/important_columns_bowtie.txt > ./IDBA_UD/final_order_bowtie2.txt
			echo -e "contig_number\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/final_bowtie2_merge_ready.txt && cat ./IDBA_UD/final_order_bowtie2.txt >> ./IDBA_UD/final_bowtie2_merge_ready.txt
			rm ./IDBA_UD/merged_original_bowtie2.txt ./IDBA_UD/important_columns_bowtie.txt ./IDBA_UD/final_order_bowtie2.txt ./IDBA_UD/important_columns.txt

			echo "assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_UD/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./IDBA_UD/BLAST_output_nt_with_taxonomy.txt > ./IDBA_UD/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt && cat ./IDBA_UD/BLAST_output_nt_with_taxonomy_noheader.txt >> ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./IDBA_UD/BLAST_output_SILVA_with_taxonomy.txt > ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./IDBA_UD/BLAST_output_nt_with_taxonomy.txt ./IDBA_UD/BLAST_output_nt_with_taxonomy_noheader.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."

			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA_UD/CLASSIFICATION/CREST/otus.csv ./IDBA_UD/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./IDBA_UD/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA_UD/CREST_seperated.txt
			cut -f2,4 ./IDBA_UD/CREST_seperated.txt > ./IDBA_UD/CREST_header.txt
			sed '1d' ./IDBA_UD/CREST_header.txt > ./IDBA_UD/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./IDBA_UD/
			sed '1d' ./IDBA_UD/CREST_tax_ready_with_taxonomy.txt > ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/CREST_merge.txt && cat ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA_UD/CREST_merge.txt

			rm ./IDBA_UD/CREST_seperated.txt ./IDBA_UD/CREST_header.txt ./IDBA_UD/CREST_tax_ready.txt ./IDBA_UD/CREST_tax_ready_with_taxonomy.txt ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./IDBA_UD/CREST_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_UD/CREST_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./IDBA_UD/MERGE_FILES
			mv ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/CREST_merge.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_UD/MERGE_FILES/

			# Edit contig-* in all the final files
			# CREST files
			sed '1d' ./IDBA_UD/CREST_BWA.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-20 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/CREST_BWA_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/CREST_BWA_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/CREST_BWA.txt

			sed '1d' ./IDBA_UD/CREST_BOWTIE2.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-20 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/CREST_BOWTIE2_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/CREST_BOWTIE2_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./IDBA_UD/BLAST_nt_BWA.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-31 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_BWA_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_nt_BWA_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_nt_BWA.txt

			sed '1d' ./IDBA_UD/BLAST_nt_BOWTIE2.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-31 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_BOWTIE2_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_nt_BOWTIE2_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_nt_BOWTIE2.txt

			# BLAST SILVA files
			sed '1d' ./IDBA_UD/BLAST_SILVA_BWA.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-31 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_BWA_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_SILVA_BWA_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_SILVA_BWA.txt

			sed '1d' ./IDBA_UD/BLAST_SILVA_BOWTIE2.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-31 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_BOWTIE2_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_SILVA_BOWTIE2_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_SILVA_BOWTIE2.txt

			mkdir ./IDBA_UD/FINAL_FILES_IDBA_UD
			mv ./IDBA_UD/CREST_BWA_final.txt ./IDBA_UD/CREST_BOWTIE2_final.txt ./IDBA_UD/BLAST_SILVA_BOWTIE2_final.txt ./IDBA_UD/BLAST_nt_BOWTIE2_final.txt ./IDBA_UD/BLAST_SILVA_BWA_final.txt ./IDBA_UD/BLAST_nt_BWA_final.txt ./IDBA_UD/FINAL_FILES_IDBA_UD/

			echo "final IDBA_UD output generated."

		else
			# Prepares to merge with BLAST of CREST - does not include assembly sequence
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

			echo "no assembly sequence was included."
			echo "BWA and BOWTIE2 files are ready."

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_UD/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./IDBA_UD/BLAST_output_nt_with_taxonomy.txt > ./IDBA_UD/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt && cat ./IDBA_UD/BLAST_output_nt_with_taxonomy_noheader.txt >> ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./IDBA_UD/BLAST_output_SILVA_with_taxonomy.txt > ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./IDBA_UD/BLAST_output_nt_with_taxonomy.txt ./IDBA_UD/BLAST_output_nt_with_taxonomy_noheader.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."


			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA_UD/CLASSIFICATION/CREST/otus.csv ./IDBA_UD/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./IDBA_UD/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA_UD/CREST_seperated.txt
			cut -f2,4 ./IDBA_UD/CREST_seperated.txt > ./IDBA_UD/CREST_header.txt
			sed '1d' ./IDBA_UD/CREST_header.txt > ./IDBA_UD/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./IDBA_UD/
			sed '1d' ./IDBA_UD/CREST_tax_ready_with_taxonomy.txt > ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/CREST_merge.txt && cat ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA_UD/CREST_merge.txt

			rm ./IDBA_UD/CREST_seperated.txt ./IDBA_UD/CREST_header.txt ./IDBA_UD/CREST_tax_ready.txt ./IDBA_UD/CREST_tax_ready_with_taxonomy.txt ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./IDBA_UD/CREST_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_UD/CREST_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./IDBA_UD/MERGE_FILES
			mv ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/CREST_merge.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_UD/MERGE_FILES/


			# Edit contig-* in all the final files
			# CREST files
			sed '1d' ./IDBA_UD/CREST_BWA.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-18 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/CREST_BWA_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/CREST_BWA_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/CREST_BWA.txt

			sed '1d' ./IDBA_UD/CREST_BOWTIE2.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-18 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/CREST_BOWTIE2_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/CREST_BOWTIE2_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./IDBA_UD/BLAST_nt_BWA.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-29 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_BWA_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_nt_BWA_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_nt_BWA.txt

			sed '1d' ./IDBA_UD/BLAST_nt_BOWTIE2.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-29 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_BOWTIE2_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_nt_BOWTIE2_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_nt_BOWTIE2.txt

			# BLAST SILVA files
			sed '1d' ./IDBA_UD/BLAST_SILVA_BWA.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-29 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_BWA_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_SILVA_BWA_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_SILVA_BWA.txt

			sed '1d' ./IDBA_UD/BLAST_SILVA_BOWTIE2.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-29 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_BOWTIE2_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_SILVA_BOWTIE2_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_SILVA_BOWTIE2.txt

			mkdir ./IDBA_UD/FINAL_FILES_IDBA_UD
			mv ./IDBA_UD/CREST_BWA_final.txt ./IDBA_UD/CREST_BOWTIE2_final.txt ./IDBA_UD/BLAST_SILVA_BOWTIE2_final.txt ./IDBA_UD/BLAST_nt_BOWTIE2_final.txt ./IDBA_UD/BLAST_SILVA_BWA_final.txt ./IDBA_UD/BLAST_nt_BWA_final.txt ./IDBA_UD/FINAL_FILES_IDBA_UD/

			echo "final IDBA_UD output generated."

		fi

	else

		echo "no FINAL merge generated."

	fi

else

	echo "no IDBA_UD output generated."

fi

mkdir ASSEMBLERS
mv SPADES METASPADES MEGAHIT IDBA_UD ASSEMBLERS/

# Extract all final output files out of ASSEMBLER follders into one folder
cd ASSEMBLERS
for assembler_folder in *
do
  for file in $assembler_folder/FINAL_FILES_$assembler_folder/*
  do
    cp $file ${assembler_folder}_${file##*/}
  done
done
cd ..
mkdir final_files
mv ASSEMBLERS/*.txt final_files


# Display runtime
echo -e "=================================================================\n\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"


# Create log
) 2>&1 | tee PIPELINE_DNA_LOG.txt
mv PIPELINE_DNA_LOG.txt PIPELINE_DNA/
