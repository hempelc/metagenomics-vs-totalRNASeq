#!/bin/bash

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
forward_reads=''
reverse_reads=''
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
		1) forward_reads="${OPTARG}" ;;
		2) reverse_reads="${OPTARG}" ;;
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
 if [[ -z "$forward_reads" || -z "$reverse_reads" ]]
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

echo -e "Forward reads were defined as $forward_reads.\n"
echo -e "Reverse reads were defined as $reverse_reads.\n"
echo -e "Script started with full command: $cmd\n"

echo -e "======== START RUNNING SCRIPT ========\n"


######################### Step 1: trimming ################################
echo -e "======== START STEP 1: TRIMMING ========\n"

# Trimming is done with separate script:
trimming_with_phred_scores_and_fastqc_report.sh -T /hdd1/programs_for_pilot/Trimmomatic-0.39/trimmomatic-0.39.jar -1 $forward_reads -2 $reverse_reads
mv trimming_with_phred_scores_and_fastqc_report_output/ step_1_trimming/

# trimmed reads are then in ./step_1_trimming/trimmomatic/*_trimmed_at_phred_*/*P.fastq <-- *1P.fastq for forward, *2P.fastq for reverse

echo -e "======== FINISHED STEP 1: TRIMMING ========\n"

######################### Step 2: rRNA sorting ################################

echo -e "======== START STEP 2: rRNA SORTING ========\n"
for trimming_results in step_1_trimming/trimmomatic/*/; do
	mkdir ${trimming_results}step_2_rrna_sorting/
	cd ${trimming_results}step_2_rrna_sorting/

	echo -e "\n======== SORTING TRIMMED READS IN FOLDER $trimming_results ========\n"

	echo -e "\n======== RUNNING SORTMERNA ========\n"
	mkdir SORTMERNA/
	sortmerna --ref /hdd1/databases/sortmerna_silva_databases/silva-bac-16s-id90.fasta \
	--ref /hdd1/databases/sortmerna_silva_databases/silva-arc-16s-id95.fasta \
	--ref /hdd1/databases/sortmerna_silva_databases/silva-euk-18s-id95.fasta \
	--reads ../*1P.fastq --reads ../*2P.fastq --paired_in -other -fastx 1 \
	-num_alignments 1 -v \
	-workdir SORTMERNA/
	deinterleave_fastq_reads.sh < SORTMERNA/out/aligned.fq SORTMERNA/out/aligned_R1.fq SORTMERNA/out/aligned_R2.fq

	echo -e "\n======== SORTMERNA DONE ========\n"

	echo -e "\n======== RUNNING rRNAFILTER ========\n"
	mkdir rRNAFILTER/
	cd rRNAFILTER/
	fq2fa ../../*1P.fastq R1.fa
	fq2fa ../../*2P.fastq R2.fa
	wget http://hulab.ucf.edu/research/projects/rRNAFilter/software/rRNAFilter.zip
	unzip rRNAFilter.zip
	cd rRNAFilter/
	java -jar -Xmx7g rRNAFilter_commandline.jar -i ../R1.fa -r 0
	java -jar -Xmx7g rRNAFilter_commandline.jar -i ../R2.fa -r 0
	cd ..
	rm -r rRNAFilter rRNAFilter.zip
	fasta_to_tab R1.fa_rRNA | cut -f 1 | cut -f1 -d " " > names.txt
	fasta_to_tab R2.fa_rRNA | cut -f 1 | cut -f1 -d " " >> names.txt
	sort -u names.txt > names_sorted.txt
	seqtk subseq R1.fa names_sorted.txt > rRNAFilter_paired_R1.fa
	seqtk subseq R2.fa names_sorted.txt > rRNAFilter_paired_R2.fa
	rm names_sorted.txt names.txt
	cd ..
	echo -e "\n======== rRNAFILTER DONE ========\n"

	echo -e "\n======== RUNNING BARRNAP ========\n"
# code for barrnap
	echo -e "\n======== BARRNAP DONE ========\n"

	echo -e "\n======== MAKING FOLDER UNSORTED/ AND COPYING UNSORTED READS IN THERE TO KEEP THE FOLDER STRUCTURE CONSTANT ========\n"
	mkdir UNSORTED/
	cp ../*1P.fastq ../*2P.fastq UNSORTED/

	echo -e "\n======== FINISHED SORTING TRIMMED READS IN FOLDER $trimming_results ========\n"
#	cd ../../../../

	echo -e "======== FINISHED STEP 2: rRNA SORTING ========\n"

	######################### Step 3: Assembly ################################

	echo -e "======== START STEP 3: ASSEMBLY ========\n"

  rrna_filter_results_list=(rRNAFILTER SORTMERNA UNSORTED) # <-- MISSES BARRNAP YET, NEEDS TO BE ADDED
	for rrna_filter_results in $rrna_filter_results_list; do
		if [[ $rrna_filter_results == 'rRNAFILTER' ]] ; then
			R1_sorted='rRNAFILTER/rRNAFilter_paired_R1.fa'
			R2_sorted='rRNAFILTER/rRNAFilter_paired_R2.fa'
  	elif [[ $rrna_filter_results == 'SORTMERNA' ]] ; then
			R1_sorted='SORTMERNA/out/aligned_R1.fq'
			R2_sorted='SORTMERNA/out/aligned_R2.fq'
# 	elif [[ $rrna_filter_results == 'BARRNAP' ]] ; then
#			R1_sorted=
#			R2_sorted=
	  else
			R1_sorted='UNSORTED/*1P.fastq'
			R2_sorted='UNSORTED/*2P.fastq'
		fi

		mkdir $rrna_filter_results/step_3_assembly/
    cd $rrna_filter_results/step_3_assembly/

		echo -e "\n======== RUNNING SPADES ========\n"
		mkdir SPADES/
		spades.py -1 $R1_sorted -2 $R2_sorted -o SPADES/
		echo -e "\n======== SPADES DONE ========\n"

		echo -e "\n======== RUNNING METASPADES ========\n"
		mkdir METASPADES/
		metaspades.py -1 $R1_sorted -2 $R2_sorted -o METASPADES/
		echo -e "\n======== METASPADES DONE ========\n"

		echo -e "\n======== RUNNING MEGAHIT ========\n"
		mkdir MEGAHIT/
		megahit --presets meta-large -t 16 -1 $R1_sorted -2 $R2_sorted -o MEGAHIT/
		echo -e "\n======== MEGAHIT DONE ========\n"

		echo -e "\n======== RUNNING IDBA_UD ========\n"
		fq2fa --merge --filter $R1_sorted $R2_sorted idba_ud_input.fa
		idba_ud --num_threads 16 --pre_correction -r idba_ud_input.fa -o IDBA_UD/
		mv idba_ud_input.fa IDBA_UD/
		echo -e "\n======== IDBA_UD DONE ========\n"

		echo -e "\n======== RUNNING RNASPADES ========\n"
		mkdir RNASPADES/
		rnaspades.py -1 $R1_sorted -2 $R2_sorted -o RNASPADES/
		echo -e "\n======== RNASPADES DONE ========\n"

		echo -e "\n======== RUNNING IDBA_TRAN ========\n"
		fq2fa --merge --filter $R1_sorted $R2_sorted idba_tran_input.fa
		idba_tran --num_threads 16 --pre_correction -r idba_tran_input.fa -o IDBA_TRAN/
		mv idba_ud_input.fa IDBA_TRAN/
		echo -e "\n======== IDBA_TRAN DONE ========\n"

		echo -e "\n======== RUNNING TRINITY ========\n"
		Trinity --seqType fq --max_memory 64G --left $R1_sorted --right $R2_sorted --CPU 16 --output TRINITY/
		cat TRINITY/Trinity.fasta | sed 's/ len/_len/g' > TRINITY/Trinity_with_length.fasta
		echo -e "\n======== TRINITY DONE ========\n"

		echo -e "\n======== RUNNING TRANSABYSS ========\n"
		# TRANSABYSS CODE
		echo -e "\n======== TRANSABYSS DONE ========\n"

  done

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"


# Create log
) 2>&1 | tee PIPELINE_ASSEMBLERS_LOG.txt
