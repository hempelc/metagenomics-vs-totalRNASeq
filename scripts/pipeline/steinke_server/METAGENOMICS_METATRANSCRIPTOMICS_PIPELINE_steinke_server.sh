#!/bin/bash

# Version 0.1
# Written by Natalie Wright (nwrigh06@uoguelph.ca) and Chris Hempel (hempelc@uoguelph.ca)

# This is a pipeline for Chris Hempel's first PhD chapter

# It trims raw paired-end input sequences at 4 PHRED scores, filters rRNA
# with 4 approaches, uses 8 assemblers, maps trimmed reads back to scaffolds
# using 2 mappers, and assigns taxonomy to scaffolds using 2 databases with
# 3 classification approaches.

# The output is a folder called METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/
# that contains tab-separated, taxonomically annotated scaffolds and read counts
# for all combinations of the above described steps.

# The pipeline requires the following subscripts, which are all located in the
# subscripts/ directory:
	# assign_NCBI_staxids_to_CREST_v4.py, fasta_to_tab, mergeFilesOnColumn.pl,
	# assign_taxonomy_to_NCBI_staxids.sh  fastqc_on_R1_R2_and_optional_trimming.sh,
	# merge_on_outer.py, blast_filtering.bash, filter-fasta.awk,
	# SILVA_SSU_LSU_kraken2_preparation.sh, deinterleave_fastq_reads.sh,
	# LookupTaxonDetails3.py, SILVA_SSU_LSU_makeblastdb_preparation.sh

# The pipeline requires the following programs/python packages (versions we used
# when writing this script are indicated in brackets):
	# FastQC (0.11.5), Trimmomatic (0.33), sortmeRNA (4.0.0), barrnap (0.9),
	# rRNAFILTER (1.1), SPADES (3.14.0), METASPADES (3.14.0), RNASPADES (3.14.0),
	# MEGAHIT (1.2.9), IDBA-UD (1.1.1), IDBA-TRAN (1.1.1), Trinity (2.10.0),
	# bowtie2 (2.3.3.1), bwa (0.7.17), blast (2.10.0+), seqtk (1.2-r94),
	# samtools (1.10), python module justblast (2020.0.3), python module ete3 (3.1.2)
	
	# Note 1: we had to edit IDBA prior to compiling it because it didn't work
	# using long reads and the -l option. This seems to be a common problem and
	# can be circumvented following for example the instructions in
	# http://seqanswers.com/forums/showthread.php?t=29109, and see also
	# https://github.com/loneknightpy/idba/issues/26

	# Note 2: we had to install ete3 via conda in a conda environment, so we have
	# to activate that environment when running this script

cmd="$0 $@" # Make variable containing the entire entered command to print command to logfile
usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq> [-t <n>]
Usage:
	-1 Full path to forward reads in .fastq/.fq format
	-2 Full path to reverse reads in .fastq/.fq format
	-t Number of threads (default:16)
	-h Display this help and exit"

# Set default options:
threads='16'

# Set specified options:
while getopts ':1:2:t:h' opt; do
 	case "${opt}" in
		1) forward_reads="${OPTARG}" ;;
		2) reverse_reads="${OPTARG}" ;;
		t) threads="${OPTARG}" ;;
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

# Check if required options are set:
if [[ -z "$forward_reads" || -z "$reverse_reads" ]]
then
   echo -e "-1 and -2 must be set.\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script.\n"
   exit
fi


##################### Write start time and options to output ######################

# Make open bracket to later tell script to write everything that follows into a logfile:
(

# Define starting time of script for total runtime calculation:
start=$(date +%s)
echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
echo -e "=================================================================\n\n"


# Output specified options:
echo -e "======== OPTIONS ========\n"

echo -e "Forward reads were defined as $forward_reads.\n"
echo -e "Reverse reads were defined as $reverse_reads.\n"
echo -e "Number of threads was set to $threads.\n"
echo -e "Script started with full command: $cmd\n"



######################### Start of the actual script ################################
echo -e "++++++++ START RUNNING SCRIPT ++++++++\n"

# Activate the conda ete3 environment within this script to be able to run ete3.
# I found this solution # to activate conda environments in scripts here:
# https://github.com/conda/conda/issues/7980.
eval "$(conda shell.bash hook)" # Without this, the conda environment cannot be
# activated within the script
conda activate ete3 # ete3 is our conda environemnt in which we installed ete3

# Make output directory and directory for final files:
mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/
mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/
cd METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/

# Save full current path in variable to make navigation between directories easier:
base_directory=$(pwd)

######################### Step 1: trimming ################################
echo -e "++++++++ START STEP 1: TRIMMING AND ERROR CORRECTION ++++++++\n"

# Trimming is done with a separate subscript:
fastqc_on_R1_R2_and_optional_trimming.sh \
-T /hdd1/programs_for_pilot/Trimmomatic-0.39/trimmomatic-0.39.jar \
-1 $forward_reads -2 $reverse_reads -t yes -p $threads
mv trimming_with_phred_scores_and_fastqc_report_output/ step_1_trimming/

# Running error correction module of SPAdes on all trimmed reads
for trimming_results in step_1_trimming/trimmomatic/*; do
	echo -e "\n======== ERROR-CORRECTING READS IN FOLDER $trimming_results ========\n"
	spades.py -1 $trimming_results/*1P.fastq -2 $trimming_results/*2P.fastq \
	--only-error-correction --disable-gzip-output -o $trimming_results/error_correction \
	-t $threads
	mv $trimming_results/error_correction/corrected/*1P*.fastq \
	$trimming_results/error_correction/corrected/*2P*.fastq $trimming_results
	# Rename weird name of error-corrected reads:
	R1=$(echo $trimming_results/*1P.00.0_0.cor.fastq) \
  && 	mv $trimming_results/*1P.00.0_0.cor.fastq ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
  sed -r -i 's/ BH:.{2,6}//g' ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
	R2=$(echo $trimming_results/*2P.00.0_0.cor.fastq) \
  && mv $trimming_results/*2P.00.0_0.cor.fastq ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
  sed -r -i 's/ BH:.{2,6}//g' ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
	rm -r $trimming_results/error_correction/
	echo -e "\n======== FINISHED ERROR-CORRECTING READS IN FOLDER $trimming_results ========\n"
done

echo -e "++++++++ FINISHED STEP 1: TRIMMING AND ERROR CORRECTION ++++++++\n"


######################### Step 2: rRNA sorting ################################

# NOTE: To enable each combination of programs in each step, we run nested for
# loops and close them all at the very end of the script.

# For loop 1: loop over the folders for trimmed data at different PHRED scores
# for rRNA filtration:
for trimming_results in step_1_trimming/trimmomatic/*; do
	mkdir $trimming_results/step_2_rrna_sorting/
	cd $trimming_results/step_2_rrna_sorting/

	echo -e "++++++++ START STEP 2: rRNA SORTING OF TRIMMED READS IN FOLDER $trimming_results ++++++++\n"

	echo -e "\n======== CONVERT READS IN FASTA FORMAT FOR rRNAFILTER AND BARRNAP ========\n"
	mkdir reads_in_fasta_format/
	fq2fa ../*1P_error_corrected.fastq reads_in_fasta_format/R1.fa
	fq2fa ../*2P_error_corrected.fastq reads_in_fasta_format/R2.fa
	echo -e "\n======== READS TO FASTA CONVERSION DONE ========\n"

	echo -e "\n======== RUNNING SORTMERNA ========\n"
	mkdir SORTMERNA/
	sortmerna --ref /hdd1/databases/sortmerna_silva_databases/silva-bac-16s-id90.fasta \
	--ref /hdd1/databases/sortmerna_silva_databases/silva-arc-16s-id95.fasta \
	--ref /hdd1/databases/sortmerna_silva_databases/silva-euk-18s-id95.fasta \
	--reads ../*1P_error_corrected.fastq --reads ../*2P_error_corrected.fastq \
	--paired_in -other -fastx 1 -num_alignments 1 -v -workdir SORTMERNA/ \
	--threads 1:1:$threads
	# SortMeRNA interleaves reads, which we don't want, so we deinterleave them:
	deinterleave_fastq_reads.sh < SORTMERNA/out/aligned.fastq \
	SORTMERNA/out/aligned_R1.fq SORTMERNA/out/aligned_R2.fq
	echo -e "\n======== SORTMERNA DONE ========\n"

	echo -e "\n======== RUNNING rRNAFILTER ========\n"
	mkdir rRNAFILTER/
	cd rRNAFILTER/
	# rRNAFilter only worked for us when we started it within the directory
	# containing the .jar file. To simplify switching to that directory, we simply
	# download the small program within the script and delete it after usage:
	wget http://hulab.ucf.edu/research/projects/rRNAFilter/software/rRNAFilter.zip
	unzip rRNAFilter.zip
	cd rRNAFilter/
	# We use 7GB for the rRNAFilter .jar, as shown in the rRNAFilter manual:
	java -jar -Xmx7g rRNAFilter_commandline.jar \
	-i ../../reads_in_fasta_format/R1.fa -r 0
	java -jar -Xmx7g rRNAFilter_commandline.jar \
	-i ../../reads_in_fasta_format/R2.fa -r 0
	mv ../../reads_in_fasta_format/R*.fa_rRNA ..
	cd ..
	rm -r rRNAFilter rRNAFilter.zip
	# We want to keep paired reads, so we extract all rRNA read names that were
	# found in R1 and R2, save them as one list, and extract all reads from both
	# R1 and R2 reads. That way, even if only one read from a pair was identified
	# as rRNA, we keep the pair of reads:
	fasta_to_tab R1.fa_rRNA | cut -f 1 | cut -f1 -d " " > names.txt
	fasta_to_tab R2.fa_rRNA | cut -f 1 | cut -f1 -d " " >> names.txt
	sort -u names.txt > names_sorted.txt
	seqtk subseq ../reads_in_fasta_format/R1.fa names_sorted.txt \
	> rRNAFilter_paired_R1.fa
	seqtk subseq ../reads_in_fasta_format/R2.fa names_sorted.txt \
	> rRNAFilter_paired_R2.fa
	rm names_sorted.txt names.txt
	cd ..
	echo -e "\n======== rRNAFILTER DONE ========\n"

	echo -e "\n======== RUNNING BARRNAP ========\n"
	mkdir BARRNAP/
	for kingdom in euk bac arc; do # barrnap needs to be run on kingdoms separately
		echo -e "\n======== RUNNING BARRNAP ON KINGDOM $kingdom AND R1 READS ========\n"
		barrnap --quiet --lencutoff 0.000001 --reject 0.000001 --kingdom $kingdom \
		--threads $threads --outseq BARRNAP/${kingdom}_reads1.fa \
		reads_in_fasta_format/R1.fa
		echo -e "\n======== RUNNING BARRNAP ON KINGDOM $kingdom AND R2 READS ========\n"
		barrnap --quiet --lencutoff 0.000001 --reject 0.000001 --kingdom $kingdom \
		--threads $threads --outseq BARRNAP/${kingdom}_reads2.fa \
		reads_in_fasta_format/R2.fa
		rm reads_in_fasta_format/*.fai
		sed 's/.*::/>/g' BARRNAP/${kingdom}_reads1.fa | sed 's/:[^:]*$//g' \
		> BARRNAP/${kingdom}_reads1_edited.fa
		sed 's/.*::/>/g' BARRNAP/${kingdom}_reads2.fa | sed 's/:[^:]*$//g' \
		> BARRNAP/${kingdom}_reads2_edited.fa
	done
	# Concatenating results from the three kingdoms and R1 and R2 files
	cat BARRNAP/*edited.fa > BARRNAP/all_reads.fa
	# We want to keep paired reads, so we extract all rRNA read names that were
	# found in R1 and R2 for the three kingdoms (in all_reads.fa), save them as
	# one list, and extract all reads from both R1 and R2 reads. That way, even if
	# only one read from a pair was identified as rRNA, we keep the pair of reads:
	fasta_to_tab BARRNAP/all_reads.fa | cut -f 1 | cut -f1 -d " " | sort -u \
	> BARRNAP/names_sorted.txt
	seqtk subseq reads_in_fasta_format/R1.fa BARRNAP/names_sorted.txt \
	> BARRNAP/barrnap_paired_R1.fa
	seqtk subseq reads_in_fasta_format/R2.fa BARRNAP/names_sorted.txt \
	> BARRNAP/barrnap_paired_R2.fa
	rm BARRNAP/names_sorted.txt
	echo -e "\n======== BARRNAP DONE ========\n"

	echo -e "\n======== MAKING FOLDER UNSORTED/ AND COPYING UNSORTED READS IN THERE TO KEEP THE FOLDER STRUCTURE CONSTANT ========\n"
	mkdir UNSORTED/
	cp ../*1P_error_corrected.fastq ../*2P_error_corrected.fastq UNSORTED/

	echo -e "\n++++++++ FINISHED STEP 2: rRNA SORTING OF TRIMMED READS IN FOLDER $trimming_results ++++++++\n"


	######################### Step 3: Assembly ################################

	# For loop 2: loop over the folders for rRNA filtered data for assembly:
  for rrna_filter_results in rRNAFILTER SORTMERNA BARRNAP UNSORTED; do
		if [[ $rrna_filter_results == 'rRNAFILTER' ]]; then
			R1_sorted='rRNAFILTER/rRNAFilter_paired_R1.fa'
			R2_sorted='rRNAFILTER/rRNAFilter_paired_R2.fa'
  	elif [[ $rrna_filter_results == 'SORTMERNA' ]]; then
			R1_sorted='SORTMERNA/out/aligned_R1.fq'
			R2_sorted='SORTMERNA/out/aligned_R2.fq'
		elif [[ $rrna_filter_results == 'BARRNAP' ]]; then
			R1_sorted='BARRNAP/barrnap_paired_R1.fa'
			R2_sorted='BARRNAP/barrnap_paired_R2.fa'
	  else # UNSORTED
			R1_sorted='UNSORTED/*1P_error_corrected.fastq'
			R2_sorted='UNSORTED/*2P_error_corrected.fastq'
		fi

		echo -e "++++++++ START STEP 3: ASSEMBLY OF TRIMMED READS IN FOLDER $trimming_results/ AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ ++++++++\n"
		mkdir $rrna_filter_results/step_3_assembly/
		cd $rrna_filter_results/step_3_assembly/

		echo -e "\n======== RUNNING SPADES ========\n"
		mkdir SPADES/
		spades.py -1 ../../$R1_sorted -2 ../../$R2_sorted --only-assembler \
		-o SPADES/ -t $threads
		echo -e "\n======== SPADES DONE ========\n"

		echo -e "\n======== RUNNING METASPADES ========\n"
		mkdir METASPADES/
		metaspades.py -1 ../../$R1_sorted -2 ../../$R2_sorted --only-assembler \
		-o METASPADES/ -t $threads
		echo -e "\n======== METASPADES DONE ========\n"

		echo -e "\n======== RUNNING MEGAHIT ========\n"
		megahit --presets meta-large -t $threads -1 ../../$R1_sorted \
		-2 ../../$R2_sorted -o MEGAHIT/
		echo -e "\n======== MEGAHIT DONE ========\n"

		echo -e "\n======== RUNNING IDBA_UD ========\n"
		# Note: we had to edit IDBA prior to compiling it because it didn't work
		# using long reads and the -l option. This seems to be a common problem and
		# can be circumvented following for example the instructions in
		# http://seqanswers.com/forums/showthread.php?t=29109, and see also
		# https://github.com/loneknightpy/idba/issues/26
		# IDBA_UD only takes interleaved fasta files
		fq2fa --merge --filter ../../$R1_sorted ../../$R2_sorted idba_ud_input.fa
		idba_ud --num_threads $threads --pre_correction -r idba_ud_input.fa \
    -o IDBA_UD/
		mv idba_ud_input.fa IDBA_UD/
		echo -e "\n======== IDBA_UD DONE ========\n"

		echo -e "\n======== RUNNING RNASPADES ========\n"
		mkdir RNASPADES/
		rnaspades.py -1 ../../$R1_sorted -2 ../../$R2_sorted --only-assembler \
		-o RNASPADES/ -t $threads
		echo -e "\n======== RNASPADES DONE ========\n"

		echo -e "\n======== RUNNING IDBA_TRAN ========\n"
		# IDBA_TRAN only takes interleaved fasta files
		fq2fa --merge ../../$R1_sorted ../../$R2_sorted idba_tran_input.fa
		idba_tran --num_threads $threads --pre_correction -l idba_tran_input.fa \
    -o IDBA_TRAN/
		mv idba_tran_input.fa IDBA_TRAN/
		echo -e "\n======== IDBA_TRAN DONE ========\n"

		echo -e "\n======== RUNNING TRINITY ========\n"
		# Barrnap and rRNAFilter output fasta files which has to be indicated to Trinity:
    if [[ $rrna_filter_results == "rRNAFILTER" \
		|| $rrna_filter_results == "BARRNAP" ]]; then
  		Trinity --seqType fa --max_memory $(echo $(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024)-5)))G --left ../../$R1_sorted --right \
      ../../$R2_sorted --CPU $threads --output TRINITY/ # The max_memory command simply takes the maximum RAM size in GB and subtracts 5GB
    else
      Trinity --seqType fq --max_memory $(echo $(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024)-5)))G --left ../../$R1_sorted --right \
      ../../$R2_sorted --CPU $threads --output TRINITY/ # The max_memory command simply takes the maximum RAM size in GB and subtracts 5GB
    fi
    cat TRINITY/Trinity.fasta | sed 's/ len/_len/g' \
		> TRINITY/Trinity_with_length.fasta  # Edit for universal format
		echo -e "\n======== TRINITY DONE ========\n"

		echo -e "\n======== RUNNING TRANSABYSS ========\n"
    transabyss --pe ../../$R1_sorted ../../$R2_sorted --threads $threads \
    --outdir TRANSABYSS/
    sed 's/ /_/g' TRANSABYSS/transabyss-final.fa \
		> TRANSABYSS/transabyss-final_edited.fa # Edit for universal format
		echo -e "\n======== TRANSABYSS DONE ========\n"

		echo -e "\n++++++++ FINISHED STEP 3: ASSEMBLY OF TRIMMED READS IN FOLDER $trimming_results/ AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ ++++++++\n"


		######################### Step 4: Mapping ################################

		# For loop 3: loop over the folders for assembled data for mapping:
		for assembly_results in SPADES METASPADES MEGAHIT IDBA_UD RNASPADES IDBA_TRAN TRINITY TRANSABYSS; do
			if [[ $assembly_results == 'SPADES' ]]; then
				scaffolds='SPADES/scaffolds.fasta'
			elif [[ $assembly_results == 'METASPADES' ]]; then
				scaffolds='METASPADES/scaffolds.fasta'
			elif [[ $assembly_results == 'MEGAHIT' ]]; then
				scaffolds='MEGAHIT/final.contigs.fa'
			elif [[ $assembly_results == 'IDBA_UD' ]]; then
				scaffolds='IDBA_UD/scaffold.fa'
			elif [[ $assembly_results == 'RNASPADES' ]]; then
				scaffolds='RNASPADES/transcripts.fasta'
			elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
				scaffolds='IDBA_TRAN/contig.fa'
			elif [[ $assembly_results == 'TRINITY' ]]; then
				scaffolds='TRINITY/Trinity_with_length.fasta'
			else # TRANSABYSS
				scaffolds='TRANSABYSS/transabyss-final_edited.fa'
			fi

			echo -e "++++++++ START STEP 4: MAPPING OF TRIMMED READS IN FOLDER $trimming_results/ AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ AND ASSEMBLY IN FOLDER $assembly_results/ ++++++++\n"
			mkdir $assembly_results/step_4_mapping/
			cd $assembly_results/step_4_mapping/

			# We insert a loop here for the mappers because their output is edited
			# with identical commands:
      for mapper in BWA BOWTIE2; do

        mkdir $mapper
        cd $mapper

        if [[ $mapper == 'BWA' ]]; then
          echo -e "\n======== Starting bwa index ========\n"
          bwa index -p bwa_index ../../../$scaffolds
          echo -e "\n======== bwa index complete. Starting bwa mem ========\n"
          bwa mem -t $threads bwa_index ../../../../../../*1P_error_corrected.fastq \
					../../../../../../*2P_error_corrected.fastq > ${mapper}_output.sam
          rm bwa_index*
					echo -e "\n======== bwa mem complete ========\n"
  			else
          echo -e "\n======== Starting bowtie2 index ========\n"
    			bowtie2-build -f ../../../$scaffolds bowtie_index
    			echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"
    			bowtie2 -q -x bowtie_index -1 ../../../../../../*1P_error_corrected.fastq \
					-2 ../../../../../../*2P_error_corrected.fastq -S ${mapper}_output.sam \
					-p $threads
    			rm bowtie_index*
          echo -e "\n======== bowtie2 complete ========\n"
        fi

  			# Editing the mapper outputs:
  			samtools view -F 4 ${mapper}_output.sam | cut -f3	| sort | uniq -c \
				| column -t | sed 's/  */\t/g' > out_mapped_${mapper}.txt
  			samtools view -f 4 ${mapper}_output.sam | cut -f3 |	sort | uniq -c \
				| column -t | sed 's/  */\t/g' > out_unmapped_${mapper}.txt
  			echo -e "counts\tsequence_name" > merge_input_mapped_${mapper}.txt \
				&& cat out_mapped_${mapper}.txt >> merge_input_mapped_${mapper}.txt
  			echo -e "counts\tsequence_name" > merge_input_unmapped_${mapper}.txt \
				&& cat out_unmapped_${mapper}.txt >> merge_input_unmapped_${mapper}.txt

  			rm *_index* out_*mapped_${mapper}.txt
        cd ..
			# And we close the mapper loop here because the next step is independent
			# from the mappers and saved under a separate folder:
			done

			# We cd back into the step_3 directory using the base_directory variable
			# we created at the beginning of the script and the nested for loop
			# variables we generated during the script:
      cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)

			echo -e "++++++++ FINISHED STEP 4: MAPPING FOR TRIMMED READS IN FOLDER $trimming_results/ AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ AND ASSEMBLY IN FOLDER $assembly_results/ ++++++++\n"


			######################### Steps 5 and 6.1: Picking a referencd DB and taxonomic classification ################################

			mkdir $assembly_results/step_5_reference_DB/

			# For loop 4: loop over the reference DBs for taxonomic classification:
			for DB in SILVA NCBI_NT; do
				if [[ $DB == "SILVA" ]]; then
					krakenDB="/hdd1/databases/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020/"
					blastDB="/hdd1/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta"
				else # NCBI_NT
					krakenDB="/hdd1/databases/kraken2_nt_DB"
					blastDB="/hdd1/databases/nt_database_feb_2020_indexed/nt"
				fi

				echo -e "++++++++ START STEP 5 AND 6.1: CLASSIFICATION OF ASSEMBLED SCAFFOLDS FROM TRIMMED READS IN FOLDER $trimming_results/ AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ AND ASSEMBLY IN FOLDER $assembly_results/ USING DATABASE ${DB} +++++++\n"
				mkdir $assembly_results/step_5_reference_DB/${DB}/
				mkdir $assembly_results/step_5_reference_DB/${DB}/step_6_classification/
				cd $assembly_results/step_5_reference_DB/${DB}/step_6_classification/

				echo -e "\n======== RUNNING JUSTBLAST WITH DATABASE $DB ========\n"
				# Run BLAST via justblast
				justblast ../../../../$scaffolds $blastDB --cpus $threads --evalue 1e-05 \
				--outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
				--out_filename blast_output.txt
				echo -e "\n======== JUSTBLAST WITH DATABASE $DB DONE ========\n"

				echo -e "\n======== RUNNING BLAST FIRST HIT ========\n"
				# We run a separate script to filter the BLAST results:
				blast_filtering.bash -i blast_output.txt -f blast -t soft -T $threads
				cp blast_output.txt blast_filtering_results/
				mv blast_filtering_results/ BLAST_FIRST_HIT/
        echo -e "\n======== BLAST FIRST HIT DONE ========\n"

				echo -e "\n======== RUNNING BLAST FILTERED ========\n"
				# We run a separate script to filter the BLAST results:
				blast_filtering.bash -i blast_output.txt -f blast -t strict -T $threads
				mv blast_output.txt blast_filtering_results/
				mv blast_filtering_results/ BLAST_FILTERED/
        echo -e "\n======== BLAST FILTERED DONE========\n"

				echo -e "\n======== RUNNING KRAKEN2 WITH DATABASE $DB ========\n"
				mkdir KRAKEN2/
				cd KRAKEN2/
	     	# Run kraken2
				kraken2 --db $krakenDB --threads $threads ../../../../../$scaffolds \
				> kraken2_output.txt

				if [[ $DB == 'SILVA' ]]; then
		      # Now we're gonna edit the output so that is has the same format as
					# CREST output, since we already have a script to deal with
					# SILVA CREST output/taxonomy
		      # Extract the taxids column of the standard kraken output:
					cut -f3 kraken2_output.txt > kraken2_taxids.txt
		      # Access the SILVA taxonomy file (downloaded taxmap files as in
					# subscript SILVA_SSU_LSU_kraken2_preparation and concatenated them)
					# and generate a file containing one column for each SILVA taxid and
					# one column for the respective SILVA taxonomy path:
		     	tail -n +2 /hdd1/databases/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020/taxmap_slv_ssu_lsu_ref_nr_138.1.txt \
		      | cut -f 4,6 | sort -u > SILVA_paths_and_taxids.txt
		      # Kraken2 spits out the taxid 0 when no hit is found, but 0 doesn't
					# exist in the SILVA taxonomy, so manually add taxid 0 with path
					# “No hits” to the SILVA path file:
		      echo -e "No hits;\t0" > tmp && cat SILVA_paths_and_taxids.txt >> tmp \
		      && mv tmp SILVA_paths_and_taxids.txt
		      # Merge your kraken2 taxids with the SILVA path file to assign a SILVA
					# taxonomy path to every kraken2 hit:
		      mergeFilesOnColumn.pl SILVA_paths_and_taxids.txt kraken2_taxids.txt 2 1 > merged.txt
		      cut -f -2 merged.txt | sed 's/;\t/\t/g' > merged_edit.txt # Edit the output
		      # Extract the sequence names from the kraken2 output and generate a
					# final file with sequence name, taxid, and SILVA path:
		     	cut -f 3 kraken2_output.txt > names.txt
		      paste names.txt merged_edit.txt | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2}' \
				 	> kraken2_SILVA_formatted.txt
		      # This file has now the same format as the output of CREST and can be
					# translated into NCBI taxonomy the same way as CREST output

					# We run a separate script that was initially made to deal with the
					# SILVA taxonomy of CREST output, by translating SILVA taxonomic paths
					# into NCBI taxids, and use that script on the formatted kraken2 output.
					# The files NCBI_staxids_(non_)scientific were generated by the script
					# SILVA_SSU_LSU_makeblastdb_preparation:
		      assign_NCBI_staxids_to_CREST_v4.py /hdd1/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/NCBI_staxids_scientific.txt \
		      /hdd1/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/NCBI_staxids_non_scientific.txt \
		      kraken2_SILVA_formatted.txt kraken2_SILVA_formatted_with_NCBI_taxids.txt
		      mergeFilesOnColumn.pl kraken2_SILVA_formatted_with_NCBI_taxids.txt \
		      kraken2_SILVA_formatted.txt 1 1 | cut -f3 > NCBItaxids.txt # Merge SILVA output with taxids and extract taxids
					# We use a separate script to assign taxonomy to NCBI taxids:
		      assign_taxonomy_to_NCBI_staxids.sh -b NCBItaxids.txt -c 1 \
					-e ~/.etetoolkit/taxa.sqlite
		      sed -i '1d' NCBItaxids_with_taxonomy.txt # Remove header
		      cut -f2 kraken2_output.txt > contig_names.txt # Get contig names from original kraken2 output
	       	paste contig_names.txt NCBItaxids_with_taxonomy.txt \
				 	> contigs_with_NCBItaxids_and_taxonomy.txt # Add contig names to taxonomy file
		     	echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
		     	> kraken2_final.txt && cat contigs_with_NCBItaxids_and_taxonomy.txt | sed 's/Unknown/NA/g' \
		    	>> kraken2_final.txt # Add header

		      # Sort files
			   	mkdir intermediate_files
		      mv kraken2_output.txt kraken2_taxids.txt SILVA_paths_and_taxids.txt merged* \
		      names.txt kraken2_SILVA_formatted* NCBItaxids* contig* intermediate_files/

	      else # NCBI_NT
	      	cut -f 2-3 kraken2_output.txt > kraken2_output_contig_taxid.txt # Isolate contig names and taxids
					# We use a separate script to assign taxonomy to NCBI taxids:
	        assign_taxonomy_to_NCBI_staxids.sh -b kraken2_output_contig_taxid.txt \
					-c 2 -e ~/.etetoolkit/taxa.sqlite
	        sed -i '1d' kraken2_output_contig_taxid_with_taxonomy.txt # Remove header
	        echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
	        > kraken2_final.txt && cat kraken2_output_contig_taxid_with_taxonomy.txt | sed 's/Unknown/NA/g' \
	        >> kraken2_final.txt # Add header

	       	# Sort files
	        mkdir intermediate_files
	        mv kraken2_output* intermediate_files/
	      fi

				cd ..
	      echo -e "\n======== KRAKEN2 WITH DATABASE $DB DONE ========\n"

				######################### Step 6.2: Generating final putput files ################################
				# Each assembler/classification tool output has a different format. We
				# make that format universal with the following code.

				# For loop 5: loop over the classification tools for universal edits:
				for classification_tool in KRAKEN2 BLAST_FILTERED BLAST_FIRST_HIT; do

					echo -e "++++++++ START STEP 6.2: GENERATING FINAL OUTPUT FILES OF READ IN $trimming_results/ AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ AND ASSEMBLY IN FOLDER $assembly_results/ AND CLASSIFICATION IN FOLDER $classification_tool/ ++++++++\n"

					mkdir ${classification_tool}/FINAL_FILES/
					mkdir ${classification_tool}/FINAL_FILES/intermediate_files/

					# We need to open a loop for the mappers again to access their files,
					# since we closed that loop earlier:
					for mapper in BWA BOWTIE2; do

						# We use full paths here to access the separate files.
 					  # Add assembly sequences to BWA and BOWTIE file:
						if [[ $assembly_results == 'SPADES' || $assembly_results == 'METASPADES' \
						|| $assembly_results == 'IDBA_UD' || $assembly_results == 'RNASPADES' \
						|| $assembly_results == 'TRANSABYSS' ]]; then
							# Change the sequences' fasta format to tab delimited:
							fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} \
							> tmp
							# Add header name so that we can later merge on outer with a python script:
							echo -e "sequence_name\tsequence" > ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt \
							&& cat tmp >> ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt
							# And remove intermediate file:
							rm tmp
							# Now the scaffold names can be merged with the mapper output:
							merge_on_outer.py ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_final_${mapper}_merge_ready.txt

						elif [[ $assembly_results == 'MEGAHIT' ]]; then
							# Change the sequences' fasta format to tab delimited:
							fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} \
							| sed 's/ /\t/g' | cut -f1,4,5 | sed 's/len=//g' \
							> tmp
							# Add header name so that we can later merge on outer with a python script:
							echo -e "sequence_name\tsequence_length\tsequence" > ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt \
							&& cat tmp >> ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt
							# And remove intermediate file:
							rm tmp
							# Now the scaffold names can be merged with the mapper output:
							merge_on_outer.py ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_final_${mapper}_merge_ready.txt

						elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
							# Change the sequences' fasta format to tab delimited:
							fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} \
							| sed 's/_/\t/2'  | sed 's/_/\t/3' | sed 's/ /\t/g' | cut -f1,3,5,6 \
							> tmp
							# Add header name so that we can later merge on outer with a python script:
							echo -e "sequence_name\tsequence_length\tcontig_kmer_count\tsequence" > ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt \
							&& cat tmp >> ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt
							# And remove intermediate file:
							rm tmp
							# Now the scaffold names can be merged with the mapper output:
							merge_on_outer.py ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_final_${mapper}_merge_ready.txt

						else # TRANSABYSS
							# Change the sequences' fasta format to tab delimited:
							fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} \
							| sed 's/ /\t/1' | cut -f1,3 > tmp
							# Add header name so that we can later merge on outer with a python script:
							echo -e "sequence_name\tsequence" > ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt \
							&& cat tmp >> ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt
							# And remove intermediate file:
							rm tmp
							# Now the scaffold names can be merged with the mapper output:
							merge_on_outer.py ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_tab_to_merge.txt ${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_final_${mapper}_merge_ready.txt

						fi

						# Add classification data to the mapper output and do
						# final individual edits:

						if [[ $classification_tool == 'BLAST_FILTERED' ]]; then

							# We run a simple python script because we need to merge on "outer":
							merge_on_outer.py \
							${classification_tool}/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_final_${mapper}_merge_ready.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt

							if [ $assembly_results == 'SPADES' ] \
							|| [ $assembly_results == 'METASPADES' ]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
								| sed 's/length_//g' | sed 's/cov_//g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							elif [[ $assembly_results == 'IDBA_UD' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/scaffold_//g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							elif [[ $assembly_results == 'MEGAHIT' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/1' | cut -f2-17 > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							elif [[ $assembly_results == 'RNASPADES' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
								| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
								| cut -f1,2,3,5-21 > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								|	sed 's/contig-[0-9]*_//g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $16, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $17}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							elif [[ $assembly_results == 'TRINITY' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' \
								> ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							else #TRANSABYSS
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							fi

							echo -e "\n======== DONE FINALIZING BLAST_FILTERED FILES FOR MAPPER $mapper =======\n"

						elif [[ $classification_tool == 'BLAST_FIRST_HIT' ]]; then

							# We run a simple python script because we need to merge on "outer":
							merge_on_outer.py ${classification_tool}/blast_output_with_taxonomy_and_best_hit.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_final_${mapper}_merge_ready.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt

							if [ $assembly_results == 'SPADES' ] \
							|| [ $assembly_results == 'METASPADES' ]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
								| sed 's/length_//g' | sed 's/cov_//g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							elif [[ $assembly_results == 'IDBA_UD' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/scaffold_//g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							elif [[ $assembly_results == 'MEGAHIT' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/1' | cut -f2-17 > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							elif [[ $assembly_results == 'RNASPADES' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
								| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
								| cut -f1,2,3,5-20 > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								|	sed 's/contig-[0-9]*_//g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $16, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $17}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							elif [[ $assembly_results == 'TRINITY' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							else #TRANSABYSS
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

							fi

							echo -e "\n======== DONE FINALIZING BLAST_FIRST_HIT FILES FOR MAPPER $mapper ========\n"

						else # kraken2

							# We run a simple python script because we need to merge on "outer":
							merge_on_outer.py ${classification_tool}/kraken2_final.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${assembly_results}_final_${mapper}_merge_ready.txt \
							${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt

							if [ $assembly_results == 'SPADES' ] \
							|| [ $assembly_results == 'METASPADES' ]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
								| sed 's/length_//g' | sed 's/cov_//g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							elif [[ $assembly_results == 'IDBA_UD' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/1' | cut -f2-20 > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $17}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							elif [[ $assembly_results == 'MEGAHIT' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/1' | cut -f2-19 > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tsequence_length\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $16, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $17, $18}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							elif [[ $assembly_results == 'RNASPADES' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/2' |sed 's/_/\t/3' | sed 's/NODE_//g' \
								| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
								| cut -f1,2,3,5-20 > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/1' | cut -f2-20 > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $17, $18, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $19}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							elif [[ $assembly_results == 'TRINITY' ]]; then
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tsequence_length\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $5, $4, $17, $18}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							else #TRANSABYSS
								sed '1d' ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt \
								| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt
								echo -e "sequence_name\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \								> tmp \
								&& cat ${classification_tool}/FINAL_FILES/intermediate_files/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_no_header.txt \
								>> tmp
								awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
								> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt \
								&& rm tmp

							fi

							echo "PWD: $(pwd)"

							echo -e "\n======== DONE FINALIZING KRAKEN2 FILES FOR MAPPER $mapper ========\n"
						fi
					done
				done
				# After each loop, we have to go back to the directory we started the
				# loop in, and we're using realtive paths for that, making use of the
				# variables generated in the previous loops
				cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)
			done
			cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)
		done
		cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/)
	done
	cd $(realpath --relative-to=$(pwd) $base_directory)
done

# We make a final directory and extract all final files into that directory to
# have everything in one place
#mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/
#find . -type f -name "*_pipeline_final.txt" -print0 | xargs -0 -Ifile cp file \
#METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"

# Write output to console and log file
) 2>&1 | tee METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_LOG.txt
