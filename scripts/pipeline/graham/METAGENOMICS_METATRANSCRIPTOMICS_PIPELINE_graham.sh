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

# NOTE FOR SERGIO: Trinity uses the option --max_memory $(echo $(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024)-5)))G
# I assume that doesn't work so needs to be replaced by global variable for job's max memory

cmd="$0 $@" # Make variable containing the entire entered command to print command to logfile
usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq> -p <line from external file with pipeline tools to use> -N <NCBI_NT BLAST database> -S <SILVA BLAST databse> -n <NCBI_NT kraken2 database> -s <SILVA kraken2 database> -e <PATH/TO/.etetoolkit/taxa.sqlite> [-t <n>]
Usage:
	-1 Full path to forward reads in .fastq/.fq format
	-2 Full path to reverse reads in .fastq/.fq format
	-p Pipeline tools
	-N Path to NCBI_NT database for BLAST
	-S Path to SILVA database for BLAST
	-n Path to NCBI_NT database for kraken2
	-s Path to SILVA database for kraken23.
	-e Path to .etetoolkit/taxa.sqlite
	-t Number of threads (default:16)
	-h Display this help and exit"

# Set default options:
threads='16'

# Set specified options:
while getopts ':1:2:p:N:S:n:s:e:t:h' opt; do
 	case "${opt}" in
		1) forward_reads="${OPTARG}" ;;
		2) reverse_reads="${OPTARG}" ;;
		p) pipeline="${OPTARG}" ;;
		N) ncbi_nt_blast_db="${OPTARG}" ;;
		S) silva_blast_db="${OPTARG}" ;;
		n) ncbi_nt_kraken2_db="${OPTARG}" ;;
		s) silva_kraken2_db="${OPTARG}" ;;
		e) etetoolkit="${OPTARG}" ;;
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
if [[ -z $forward_reads || -z $reverse_reads || -z $pipeline || -z $ncbi_nt_blast_db || -z $silva_blast_db || -z $ncbi_nt_kraken2_db || -z $silva_kraken2_db || -z $etetoolkit ]]
then
   echo -e "-1, -2, -p, -N, -S, -n, -s, and -e must be set.\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script.\n"
   exit
fi

# Set pipeline tools to use
trimming=$(echo $pipeline | cut -f1 -d,)
sorting=$(echo $pipeline | cut -f2 -d,)
assembly=$(echo $pipeline | cut -f3 -d,)
mapping=$(echo $pipeline | cut -f4 -d,)
db=$(echo $pipeline | cut -f5 -d,)
classification=$(echo $pipeline | cut -f6 -d,)



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
echo -e "Tools for the pipeline were set to $pipeline.\n"
echo -e "Number of threads was set to $threads.\n"
echo -e "Script started with full command: $cmd\n"



######################### Start of the actual script ################################
echo -e "++++++++ START RUNNING SCRIPT ++++++++\n"

# Make output directory and directory for final files:
mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/
mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/
cd METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/

# Save full current path in variable to make navigation between directories easier:
base_directory=$(pwd)

######################### Step 1: trimming ################################
echo -e "++++++++ START STEP 1: TRIMMING AND ERROR CORRECTION ++++++++\n"
pwd
# Trimming is done with a separate subscript:
fastqc_on_R1_R2_and_optional_trimming.sh \
-T /hdd1/programs_for_pilot/Trimmomatic-0.39/trimmomatic-0.39.jar \
-1 $forward_reads -2 $reverse_reads -t yes -p $threads -P $trimming
mv trimming_with_phred_scores_and_fastqc_report_output/ step_1_trimming/

# Use lines from the script "fastqc_on_R1_R2_and_optional_trimming.sh" to be\
# able to change into the generated directory:
baseout=${forward_reads%_*} # Make basename
pwd
cd step_1_trimming/trimmomatic/trimmed_at_phred_${trimming}_${baseout##*/}

# Running error correction module of SPAdes on all trimmed reads
echo -e "\n======== ERROR-CORRECTING READS ========\n"
spades.py -1 *1P.fastq -2 *2P.fastq \
--only-error-correction --disable-gzip-output -o error_correction \
-t $threads
mv error_correction/corrected/*1P*.fastq \
error_correction/corrected/*2P*.fastq .
# Rename weird name of error-corrected reads:
R1=$(echo *1P.00.0_0.cor.fastq) \
&& 	mv *1P.00.0_0.cor.fastq ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
sed -r -i 's/ BH:.{2,6}//g' ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
R2=$(echo *2P.00.0_0.cor.fastq) \
&& mv *2P.00.0_0.cor.fastq ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
sed -r -i 's/ BH:.{2,6}//g' ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
rm -r error_correction/
echo -e "\n======== FINISHED ERROR-CORRECTING READS ========\n"

echo -e "++++++++ FINISHED STEP 1: TRIMMING AND ERROR CORRECTION ++++++++\n"
pwd
######################### Step 2: rRNA sorting ################################

echo -e "++++++++ START STEP 2: rRNA SORTING OF TRIMMED READS ++++++++\n"

mkdir step_2_rrna_sorting/
cd step_2_rrna_sorting/

if [[ ${sorting} == "barrnap" || ${sorting} == "rrnafilter" ]]; then
	echo -e "\n======== CONVERT READS IN FASTA FORMAT ========\n"
	mkdir reads_in_fasta_format/
	fq2fa ../*1P_error_corrected.fastq reads_in_fasta_format/R1.fa
	fq2fa ../*2P_error_corrected.fastq reads_in_fasta_format/R2.fa
	echo -e "\n======== READS TO FASTA CONVERSION DONE ========\n"
fi

if [[ ${sorting} == "sortmerna" ]]; then
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
	cd SORTMERNA/
	echo -e "\n======== SORTMERNA DONE ========\n"
fi

if [[ ${sorting} == "rrnafilter" ]]; then
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
	echo -e "\n======== rRNAFILTER DONE ========\n"
fi

if [[ ${sorting} == "barrnap" ]]; then
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
	cd BARRNAP/
	echo -e "\n======== BARRNAP DONE ========\n"
fi

if [[ ${sorting} == "unsorted" ]]; then
	echo -e "\n======== MAKING FOLDER UNSORTED/ AND COPYING UNSORTED READS IN THERE TO KEEP THE FOLDER STRUCTURE CONSTANT ========\n"
	mkdir UNSORTED/
	cp ../*1P_error_corrected.fastq ../*2P_error_corrected.fastq UNSORTED/
	cd UNSORTED/
fi

echo -e "\n++++++++ FINISHED STEP 2: rRNA SORTING ++++++++\n"
pwd
######################### Step 3: Assembly ################################

echo -e "++++++++ START STEP 3: ASSEMBLY ++++++++\n"

mkdir step_3_assembly/
cd step_3_assembly/

if [[ ${sorting} == 'rrnafilter' ]]; then
	R1_sorted='rRNAFilter_paired_R1.fa'
	R2_sorted='rRNAFilter_paired_R2.fa'
elif [[ ${sorting} == 'sortmerna' ]]; then
	R1_sorted='out/aligned_R1.fq'
	R2_sorted='out/aligned_R2.fq'
elif [[ ${sorting} == 'barrnap' ]]; then
	R1_sorted='barrnap_paired_R1.fa'
	R2_sorted='barrnap_paired_R2.fa'
else # UNSORTED
	R1_sorted='*1P_error_corrected.fastq'
	R2_sorted='*2P_error_corrected.fastq'
fi

if [[ $assembly == "spades" ]]; then
	echo -e "\n======== RUNNING SPADES ========\n"
	mkdir SPADES/
	spades.py -1 ../$R1_sorted -2 ../$R2_sorted --only-assembler \
	-o SPADES/ -t $threads
	cd SPADES/
	echo -e "\n======== SPADES DONE ========\n"
fi

if [[ $assembly == "metaspades" ]]; then
	echo -e "\n======== RUNNING METASPADES ========\n"
	mkdir METASPADES/
	metaspades.py -1 ../$R1_sorted -2 ../$R2_sorted --only-assembler \
	-o METASPADES/ -t $threads
	cd METASPADES/
	echo -e "\n======== METASPADES DONE ========\n"
fi

if [[ $assembly == "megahit" ]]; then
	echo -e "\n======== RUNNING MEGAHIT ========\n"
	megahit --presets meta-large -t $threads -1 ../$R1_sorted \
	-2 ../$R2_sorted -o MEGAHIT/
	cd MEGAHIT/
	echo -e "\n======== MEGAHIT DONE ========\n"
fi

if [[ $assembly == "idba-ud" ]]; then
	echo -e "\n======== RUNNING IDBA_UD ========\n"
	# Note: we had to edit IDBA prior to compiling it because it didn't work
	# using long reads and the -l option. This seems to be a common problem and
	# can be circumvented following for example the instructions in
	# http://seqanswers.com/forums/showthread.php?t=29109, and see also
	# https://github.com/loneknightpy/idba/issues/26
	# IDBA_UD only takes interleaved fasta files
	fq2fa --merge --filter ../$R1_sorted ../$R2_sorted idba_ud_input.fa
	idba_ud --num_threads $threads --pre_correction -r idba_ud_input.fa \
  -o IDBA_UD/
	mv idba_ud_input.fa IDBA_UD/
	cd IDBA_UD/
	echo -e "\n======== IDBA_UD DONE ========\n"
fi

if [[ $assembly == "rnaspades" ]]; then
	echo -e "\n======== RUNNING RNASPADES ========\n"
	mkdir RNASPADES/
	rnaspades.py -1 ../$R1_sorted -2 ../$R2_sorted --only-assembler \
	-o RNASPADES/ -t $threads
	cd RNASPADES/
	echo -e "\n======== RNASPADES DONE ========\n"
fi

if [[ $assembly == "idba-tran" ]]; then
	echo -e "\n======== RUNNING IDBA_TRAN ========\n"
	# IDBA_TRAN only takes interleaved fasta files
	fq2fa --merge ../$R1_sorted ../$R2_sorted idba_tran_input.fa
	idba_tran --num_threads $threads --pre_correction -l idba_tran_input.fa \
  -o IDBA_TRAN/
	mv idba_tran_input.fa IDBA_TRAN/
	cd IDBA_TRAN
	echo -e "\n======== IDBA_TRAN DONE ========\n"
fi

if [[ $assembly == "trinity" ]]; then
	echo -e "\n======== RUNNING TRINITY ========\n"
	# Barrnap and rRNAFilter output fasta files which has to be indicated to Trinity:
  if [[ ${sorting} == "rrnafilter" || ${sorting} == "barrnap" ]]; then
		Trinity --seqType fa --max_memory $(echo $(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024)-5)))G --left ../$R1_sorted --right \
    ../$R2_sorted --CPU $threads --output TRINITY/ # The max_memory command simply takes the maximum RAM size in GB and subtracts 5GB
  else
    Trinity --seqType fq --max_memory $(echo $(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024)-5)))G --left ../$R1_sorted --right \
    ../$R2_sorted --CPU $threads --output TRINITY/ # The max_memory command simply takes the maximum RAM size in GB and subtracts 5GB
  fi
  cat TRINITY/Trinity.fasta | sed 's/ len/_len/g' \
	> TRINITY/Trinity_with_length.fasta  # Edit for universal format
	cd TRINITY/
	echo -e "\n======== TRINITY DONE ========\n"
fi

if [[ $assembly == "transabyss" ]]; then
	echo -e "\n======== RUNNING TRANSABYSS ========\n"
  transabyss --pe ../$R1_sorted ../$R2_sorted --threads $threads \
  --outdir TRANSABYSS/
  sed 's/ /_/g' TRANSABYSS/transabyss-final.fa \
	> TRANSABYSS/transabyss-final_edited.fa # Edit for universal format
	cd TRANSABYSS/
	echo -e "\n======== TRANSABYSS DONE ========\n"
fi

echo -e "\n++++++++ FINISHED STEP 3: ASSEMBLY ++++++++\n"
pwd
######################### Step 4: Mapping ################################

echo -e "++++++++ START STEP 4: MAPPING ++++++++\n"

mkdir step_4_mapping/
cd step_4_mapping/

if [[ $assembly == 'spades' ]]; then
	scaffolds='scaffolds.fasta'
elif [[ $assembly == 'metaspades' ]]; then
	scaffolds='scaffolds.fasta'
elif [[ $assembly == 'megahit' ]]; then
	scaffolds='final.contigs.fa'
elif [[ $assembly == 'idba_ud' ]]; then
	scaffolds='scaffold.fa'
elif [[ $assembly == 'rnaspades' ]]; then
	scaffolds='transcripts.fasta'
elif [[ $assembly == 'idba_tran' ]]; then
	scaffolds='contig.fa'
elif [[ $assembly == 'trinity' ]]; then
	scaffolds='Trinity_with_length.fasta'
else # TRANSABYSS
	scaffolds='transabyss-final_edited.fa'
fi

if [[ ${mapping} == 'bwa' ]]; then
  echo -e "\n======== Starting bwa index ========\n"
  bwa index -p bwa_index ../$scaffolds
  echo -e "\n======== bwa index complete. Starting bwa mem ========\n"
  bwa mem -t $threads bwa_index ../../../../../*1P_error_corrected.fastq \
	../../../../../*2P_error_corrected.fastq > ${mapping}_output.sam
  rm bwa_index*
	echo -e "\n======== bwa mem complete ========\n"
else # bowtie2
  echo -e "\n======== Starting bowtie2 index ========\n"
	bowtie2-build -f ../$scaffolds bowtie_index
	echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"
	bowtie2 -q -x bowtie_index -1 ../../../../../*1P_error_corrected.fastq \
	-2 ../../../../../*2P_error_corrected.fastq -S ${mapping}_output.sam \
	-p $threads
	rm bowtie_index*
  echo -e "\n======== bowtie2 complete ========\n"
fi

# Editing the mapper outputs:
samtools view -F 4 ${mapping}_output.sam | cut -f3	| sort | uniq -c \
| column -t | sed 's/  */\t/g' > out_mapped_${mapping}.txt
samtools view -f 4 ${mapping}_output.sam | cut -f3 |	sort | uniq -c \
| column -t | sed 's/  */\t/g' > out_unmapped_${mapping}.txt
echo -e "counts\tsequence_name" > merge_input_mapped_${mapping}.txt \
&& cat out_mapped_${mapping}.txt >> merge_input_mapped_${mapping}.txt
echo -e "counts\tsequence_name" > merge_input_unmapped_${mapping}.txt \
&& cat out_unmapped_${mapping}.txt >> merge_input_unmapped_${mapping}.txt

rm out_*mapped_${mapping}.txt

cd ../
# We cd back into the step_3 directory using the base_directory variable
# we created at the beginning of the script and the nested for loop
# variables we generated during the script:
#cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)

echo -e "++++++++ FINISHED STEP 4: MAPPING ++++++++\n"
pwd
######################### Steps 5 and 6.1: Picking a referencd DB and taxonomic classification ################################

echo -e "++++++++ START STEP 5 AND 6.1: CLASSIFICATION OF ASSEMBLED SCAFFOLDS WITH $(echo $db | tr '[:lower:]' '[:upper:]') DATABASE +++++++\n"

mkdir step_5_reference_DB/
mkdir step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/
mkdir step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/step_6_classification/
mkdir step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/step_6_classification/$(echo $classification | tr '[:lower:]' '[:upper:]')/
cd step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/step_6_classification/$(echo $classification | tr '[:lower:]' '[:upper:]')/


if [[ $db == "silva" ]]; then
	krakenDB=$silva_kraken2_db
	blastDB=$silva_blast_db
else # ncbi_nt
	krakenDB=$ncbi_nt_kraken2_db
	blastDB=$ncbi_nt_blast_db
fi

if [[ $classification == "blast_first_hit" || $classification == "blast_filtered" ]]; then
	echo -e "\n======== RUNNING JUSTBLAST WITH DATABASE $blastDB ========\n"
	# Run BLAST via justblast
	justblast ../../../../$scaffolds $blastDB --cpus $threads --evalue 1e-05 \
	--outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
	--out_filename blast_output.txt
	rm -r dask-worker-space/
	echo -e "\n======== JUSTBLAST WITH DATABASE $blastDB DONE ========\n"
fi

if [[ $classification == "blast_first_hit" ]]; then
	echo -e "\n======== RUNNING BLAST FIRST HIT ========\n"
	# We run a separate script to filter the BLAST results:
	blast_filtering.bash -i blast_output.txt -f blast -t soft -T $threads
	mv blast_output.txt blast_filtering_results/
	echo -e "\n======== BLAST FIRST HIT DONE ========\n"
fi

if [[ $classification == "blast_filtered" ]]; then
	echo -e "\n======== RUNNING BLAST FILTERED ========\n"
	# We run a separate script to filter the BLAST results:
	blast_filtering.bash -i blast_output.txt -f blast -t strict -T $threads
	mv blast_output.txt blast_filtering_results/
	echo -e "\n======== BLAST FILTERED DONE========\n"
fi

if [[ $classification == "kraken2" ]]; then
	echo -e "\n======== RUNNING KRAKEN2 WITH DATABASE $krakenDB ========\n"
	# Run kraken2
	kraken2 --db $krakenDB --threads $threads ../../../../$scaffolds \
	> kraken2_output.txt
	if [[ $db == 'silva' ]]; then
    # Now we're gonna edit the output so that is has the same format as
		# CREST output, since we already have a script to deal with
		# SILVA CREST output/taxonomy
    # Extract the taxids column of the standard kraken output:
		cut -f3 kraken2_output.txt > kraken2_taxids.txt

    # Access the SILVA taxonomy file (downloaded taxmap files as in
		# subscript SILVA_SSU_LSU_kraken2_preparation.sh and concatenated them)
		# and generate a file containing one column for each SILVA taxid and
		# one column for the respective SILVA taxonomy path:
			# Note: since we work on the graham cluster, rather than accessing the files generated from the subscript SILVA_SSU_LSU_kraken2_preparation.sh, we perform the respecitve code here:
				# As of 04 Sep 2020, the available SILVA LSU and SSU taxmap files contain
				# duplicate accession IDs with different taxids in the SSU and LSU files. We're
				# going to use all accession IDs from the SSU file, check which additional ones
				# are in the LSU file (about 23,000 are not in the SSU file) and just take these
				# extra ones from the LSU taxmap file to not overwrite SSU taxids with LSU taxids:
				wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz
				wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz
				gunzip taxmap_slv_*su_ref_nr_138.1.txt.gz
				cat taxmap_slv_ssu_ref_nr_138.1.txt > taxmap_slv_ssu_lsu_ref_nr_138.1.txt
				tail -n +2 taxmap_slv_ssu_ref_nr_138.1.txt | cut -f1 > grep_list.txt
				tail -n +2 taxmap_slv_lsu_ref_nr_138.1.txt | grep -v -f grep_list.txt >> taxmap_slv_ssu_lsu_ref_nr_138.1.txt
		tail -n +2 taxmap_slv_ssu_lsu_ref_nr_138.1.txt | cut -f 4,6 | sort -u \
		> SILVA_paths_and_taxids.txt
		rm taxmap_slv_*_ref_nr_138.1.txt grep_list.txt

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
		# SILVA_SSU_LSU_makeblastdb_preparation.sh:
			# Note: since we work on the graham cluster, rather than accessing the files generated from the subscript SILVA_SSU_LSU_makeblastdb_preparation.sh, we perform the respecitve code here:
				# Get NCBI taxonomy files
				wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
				unzip taxdmp.zip
				# Edit names.dmp file into a file only containing scientific names
				sed "s/ <[a-zA-Z -,.&:'0-9]*>//g" names.dmp | grep 'scientific name' \
				| cut -f 1,3 | awk -F $'\t' ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' \
				| grep -v 'environmental\|uncultured\|unidentified\|metagenome' \
				> NCBI_staxids_scientific.txt
				  # 1. Remove <genus>, <family>, and other strings in <> brackets from taxonomy
				  #    (otherwise the NCBI taxonomy won't match the SILVA taxonomy)
				  # 2. Extract only scientific names and staxids
				  # 3. Cut out columns we need
				  # 4. Invert columns for script to work
				  # 5. Removes lines containing "environmental", "uncultured", "unidentified",
				  #    and "metagenome"
				    # Needed because
				      # SILVA taxonomy can have the same taxonomic ranks (e.g., "environmental
				      # sample") for different higher ranks (e.g., "nematode; environmental
				      # sample" and "bacteria;environmental sample"), which would, however, be
				      # assigned to the same staxid because the lower rank "environmental sample"
				      # is similar
				      # NCBI taxonomy can have different staxids for the same taxonomic name,
				      # which will cause issues when matching
				# We match SILVA taxonomy against this file (against scientific names) first

				# Edit names.dmp file into a file only containing non-scientific names
				sed "s/ <[a-zA-Z -,.&:'0-9]*>//g" names.dmp | grep -v 'scientific name' \
				| cut -f 1,3 | awk -F $'\t' ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' \
				| grep -v 'environmental\|uncultured\|unidentified\|metagenome' \
				> NCBI_staxids_non_scientific.txt
					# 1. Remove <genus>, <family>, and other strings in <> brackets from taxonomy
				  #    (otherwise the NCBI taxonomy won't match the SILVA taxonomy)
					# 2. Extract only non-scientific names and staxids
					# 3. Cut out columns we need
					# 4. Invert columns for script to work
					# 5. Removes lines containing "environmental", "uncultured", "unidentified",
				  #    and "metagenome"
						# Needed because
							# SILVA taxonomy can have the same taxonomic ranks (e.g., "environmental
				      # sample") for different higher ranks (e.g., "nematode; environmental
				      # sample" and "bacteria;environmental sample"), wich would, however, be
				      # assigned to the same staxid because the lower rank "environmental sample"
				      # is similar
							# NCBI taxonomy can different staxids for the same taxonomic name, which
				      # will cause issue when matching
				# If SILVA taxonomy is not in exact scientific names, then we match against these
				# to check for synonyms etc.
				rm *.dmp readme.txt taxdmp.zip gc.prt
    assign_NCBI_staxids_to_CREST_v4.py NCBI_staxids_scientific.txt \
  	NCBI_staxids_non_scientific.txt \
    kraken2_SILVA_formatted.txt kraken2_SILVA_formatted_with_NCBI_taxids.txt
		rm NCBI_staxids_*scientific.txt
    mergeFilesOnColumn.pl kraken2_SILVA_formatted_with_NCBI_taxids.txt \
    kraken2_SILVA_formatted.txt 1 1 | cut -f3 > NCBItaxids.txt # Merge SILVA output with taxids and extract taxids
		# We use a separate script to assign taxonomy to NCBI taxids:
    assign_taxonomy_to_NCBI_staxids.sh -b NCBItaxids.txt -c 1 \
		-e $etetoolkit
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

  else # ncbi_nt
  	cut -f 2-3 kraken2_output.txt > kraken2_output_contig_taxid.txt # Isolate contig names and taxids
		# We use a separate script to assign taxonomy to NCBI taxids:
    assign_taxonomy_to_NCBI_staxids.sh -b kraken2_output_contig_taxid.txt \
		-c 2 -e $etetoolkit
    sed -i '1d' kraken2_output_contig_taxid_with_taxonomy.txt # Remove header
    echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
    > kraken2_final.txt && cat kraken2_output_contig_taxid_with_taxonomy.txt | sed 's/Unknown/NA/g' \
    >> kraken2_final.txt # Add header

   	# Sort files
    mkdir intermediate_files
    mv kraken2_output* intermediate_files/
  fi
echo -e "\n======== KRAKEN2 WITH DATABASE $krakenDB DONE ========\n"
fi

######################### Step 6.2: Generating final putput files ################################
pwd
echo -e "++++++++ START STEP 6.2: GENERATING FINAL OUTPUT FILES ++++++++\n"

# Each assembler/classification tool output has a different format. We
# make that format universal with the following code.

mkdir FINAL_FILES/
mkdir FINAL_FILES/intermediate_files/
cd FINAL_FILES/intermediate_files/
pwd
if [[ $assembly == 'spades' || $assembly == 'metaspades' \
|| $assembly == 'idba_ud' || $assembly == 'rnaspades' \
|| $assembly == 'transabyss' ]]; then
	# Change the sequences' fasta format to tab delimited:
	fasta_to_tab ../../../../../../${scaffolds} > tmp
	# Add header name so that we can later merge on outer with a python script:
	echo -e "sequence_name\tsequence" > ${assembly}_tab_to_merge.txt \
	&& cat tmp >> ${assembly}_tab_to_merge.txt
	# And remove intermediate file:
	rm tmp
	# Now the scaffold names can be merged with the mapper output:
	merge_on_outer.py ../../../../../../step_4_mapping/${mapping}/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt

elif [[ $assembly == 'megahit' ]]; then
	# Change the sequences' fasta format to tab delimited:
	fasta_to_tab ../../../../../../${scaffolds} \
	| sed 's/ /\t/g' | cut -f1,4,5 | sed 's/len=//g' \
	> tmp
	# Add header name so that we can later merge on outer with a python script:
	echo -e "sequence_name\tsequence_length\tsequence" > ${assembly}_tab_to_merge.txt \
	&& cat tmp >> ${assembly}_tab_to_merge.txt
	# And remove intermediate file:
	rm tmp
	# Now the scaffold names can be merged with the mapper output:
	merge_on_outer.py ../../../../../../step_4_mapping/${mapping}/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt

elif [[ $assembly == 'idba_tran' ]]; then
	# Change the sequences' fasta format to tab delimited:
	fasta_to_tab ../../../../../../${scaffolds} \
	| sed 's/_/\t/2'  | sed 's/_/\t/3' | sed 's/ /\t/g' | cut -f1,3,5,6 \
	> tmp
	# Add header name so that we can later merge on outer with a python script:
	echo -e "sequence_name\tsequence_length\tcontig_kmer_count\tsequence" > ${assembly}_tab_to_merge.txt \
	&& cat tmp >> ${assembly}_tab_to_merge.txt
	# And remove intermediate file:
	rm tmp
	# Now the scaffold names can be merged with the mapper output:
	merge_on_outer.py ../../../../../../step_4_mapping/${mapping}/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt

else # trinity
	# Change the sequences' fasta format to tab delimited:
	fasta_to_tab ../../../../../../${scaffolds} \
	| sed 's/ /\t/1' | cut -f1,3 > tmp
	# Add header name so that we can later merge on outer with a python script:
	echo -e "sequence_name\tsequence" > ${assembly}_tab_to_merge.txt \
	&& cat tmp >> ${assembly}_tab_to_merge.txt
	# And remove intermediate file:
	rm tmp
	# Now the scaffold names can be merged with the mapper output:
	merge_on_outer.py ../../../../../../step_4_mapping/${mapping}/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt

fi

# Add classification data to the mapper output and do
# final individual edits:

if [[ $classification == 'blast_filtered' ]]; then

	# We run a simple python script because we need to merge on "outer":
	merge_on_outer.py \
	../../blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt \
	${assembly}_final_${mapping}_merge_ready.txt \
	trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	if [[ $assembly == 'spades' || $assembly == 'metaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/scaffold_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-17 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	elif [[ $assembly == 'rnaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
		| cut -f1,2,3,5-21 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		|	sed 's/contig-[0-9]*_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $16, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $17}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' \
		> trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	else #transabyss
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt
	fi

	echo -e "\n======== DONE FINALIZING BLAST_FILTERED FILES =======\n"

elif [[ $classification_tool == 'blast_first_hit' ]]; then

	# We run a simple python script because we need to merge on "outer":
	merge_on_outer.py ../../blast_output_with_taxonomy_and_best_hit.txt \
	${assembly}_final_${mapping}_merge_ready.txt \
	trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	if [[ $assembly == 'spades' || $assembly == 'metaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/scaffold_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-17 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	elif [[ $assembly == 'rnaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
		| cut -f1,2,3,5-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		|	sed 's/contig-[0-9]*_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $16, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $17}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	else #transabyss
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt
	fi

	echo -e "\n======== DONE FINALIZING BLAST_FIRST_HIT FILES ========\n"

else # kraken2

	# We run a simple python script because we need to merge on "outer":
	merge_on_outer.py ../../kraken2_final.txt \
	${assembly}_final_${mapping}_merge_ready.txt \
	trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	if [[ $assembly == 'spades' || $assembly == 'metaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $17}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-19 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tsequence_length\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $16, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $17, $18}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	elif [[ $assembly == 'rnaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' |sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
		| cut -f1,2,3,5-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $17, $18, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $19}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $5, $4, $17, $18}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp

	else #transabyss
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \								> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		&& rm tmp
	fi

	echo -e "\n======== DONE FINALIZING KRAKEN2 FILES ========\n"
fi

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"

# Write output to console and log file
) 2>&1 | tee METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_LOG.txt
