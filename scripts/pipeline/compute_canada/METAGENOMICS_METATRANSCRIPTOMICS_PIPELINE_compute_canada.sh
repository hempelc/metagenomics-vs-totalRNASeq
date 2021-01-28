#!/bin/bash

# Version 0.1
# Written by Natalie Wright (nwrigh06@uoguelph.ca) and Chris Hempel (hempelc@uoguelph.ca)

# This is a pipeline for Chris Hempel's first PhD chapter to be run on the Compute Canada servers

# It takes in lines of a separate file containing combiantions of the pipeline steps that are to be run

# The output is a folder called METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/
# that contains tab-separated, taxonomically annotated scaffolds and read counts
# for the specified pipeline.

# To run every possible combination of tools, the pipeline requires the following
# subscripts, which are all located in the subscripts/ directory:
	# assign_NCBI_staxids_to_CREST_v4.py, fasta_to_tab, mergeFilesOnColumn.pl,
	# assign_taxonomy_to_NCBI_staxids.sh  fastqc_on_R1_R2_and_optional_trimming.sh,
	# merge_on_outer.py, blast_filtering.bash, filter-fasta.awk,
	# deinterleave_fastq_reads.sh, LookupTaxonDetails3.py

# To run every possible combination of tools, the pipeline requires the following
# programs/python packages (versions we used when writing this script are
# indicated in brackets):
	# FastQC (0.11.5), Trimmomatic (0.33), sortmeRNA (4.0.0), barrnap (0.9),
	# rRNAFILTER (1.1), SPADES (3.14.0)[note: runs with the --meta and --rna
	# options for METASPADES and RNASPADES], MEGAHIT (1.2.9), IDBA-UD (1.1.1),
	# IDBA_tran (1.1.1), Trinity (2.10.0),	Trans-ABySS (2.0.1), bowtie2 (2.3.3.1),
	# bwa (0.7.17), blast+ (2.10.0+), kraken2 (2.1.1), seqtk (1.2-r94),
	# samtools (1.10), python module ete3 (3.1.2)

	# Note: we had to edit IDBA prior to compiling it because it didn't work
	# using long reads and the -l option. This seems to be a common problem and
	# can be circumvented following for example the instructions in
	# http://seqanswers.com/forums/showthread.php?t=29109, and see also
	# https://github.com/loneknightpy/idba/issues/26

cmd="$0 $@" # Make variable containing the entire entered command to print command to logfile
usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq> \
-P <line from external file with pipeline tools to use> -N <NCBI_NT BLAST database> \
-S <SILVA BLAST databse> -n <NCBI_NT kraken2 database> -s <SILVA kraken2 database> \
-B <SILVA SortMeRNA LSU bacteria database> -b <SILVA SortMeRNA SSU bacteria database> \
-A <SILVA SortMeRNA LSU archaea database> -a <SILVA SortMeRNA SSU archaea database> \
-E <SILVA SortMeRNA LSU eukaryota database> -e <SILVA SortMeRNA SSU eukaryota database> \
-R <SILVA SortMeRNA rfam 5.8S database> -r <SILVA SortMeRNA rfam 5S database> \
-x <SILVA path taxid file> -F <NCBI staxid scientific file> \
-f <NCBI staxid non-scientific file> -T <PATH/TO/trimmomatic-<version>.jar)> \
-t <PATH/TO/.etetoolkit/taxa.sqlite> -m <nnnG>[-p <n>]

Usage:
	-1 Full path to forward reads in .fastq/.fq format
	-2 Full path to reverse reads in .fastq/.fq format
	-P Pipeline tools
	-N Path to NCBI_NT database for BLAST
	-S Path to SILVA database for BLAST
	-n Path to NCBI_NT database for kraken2
	-s Path to SILVA database for kraken2
	-B Path to bacteria LSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-bac-23s-id98.fasta)
	-b Path to bacteria SSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-bac-16s-id90.fasta)
	-A Path to archaea LSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-arc-23s-id98.fasta)
	-a Path to archaea SSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-arc-16s-id95.fasta)
  -E Path to eukaryota LSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-euk-28s-id98.fasta)
	-e Path to eukaryota SSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-euk-18s-id95.fasta)
	-R Path to rfam 5.8S database for SortMeRNA (comes with SortMeRNA, rfam-5.8s-database-id98.fasta)
	-r Path to rfam 5S database for SortMeRNA (comes with SortMeRNA, rfam-5s-database-id98.fasta)
	-x Path to SILVA_paths_and_taxids.txt
	-F Path to NCBI_staxids_scientific.txt
	-f Path to NCBI_staxids_non_scientific.txt
	-T Path to trimmomatic application (trimmomatic-<version>.jar)
	-t Path to .etetoolkit/taxa.sqlite
	-m Maximum memory (format: XXXG, where XXX is a numerical value for teh emmory in Gigabyte)
	-p Number of threads (default:16)
	-h Display this help and exit"

# Set default options:
threads='16'

# Set specified options:
while getopts ':1:2:P:N:S:n:s:B:b:A:a:E:e:R:r:x:F:f:T:t:m:p:h' opt; do
 	case "${opt}" in
		1) forward_reads="${OPTARG}" ;;
		2) reverse_reads="${OPTARG}" ;;
		P) pipeline="${OPTARG}" ;;
		N) ncbi_nt_blast_db="${OPTARG}" ;;
		S) silva_blast_db="${OPTARG}" ;;
		n) ncbi_nt_kraken2_db="${OPTARG}" ;;
		s) silva_kraken2_db="${OPTARG}" ;;
		B) silva_sortmerna_bac_lsu="${OPTARG}" ;;
		b) silva_sortmerna_bac_ssu="${OPTARG}" ;;
		A) silva_sortmerna_arc_lsu="${OPTARG}" ;;
		a) silva_sortmerna_arc_ssu="${OPTARG}" ;;
		E) silva_sortmerna_euk_lsu="${OPTARG}" ;;
		e) silva_sortmerna_euk_ssu="${OPTARG}" ;;
		R) silva_sortmerna_rfam_5="${OPTARG}" ;;
		r) silva_sortmerna_rfam_5_8="${OPTARG}" ;;
		x) silva_path_taxid="${OPTARG}" ;;
		F) ncbi_scientific="${OPTARG}" ;;
		f) ncbi_non_scientific="${OPTARG}" ;;
		T) trimmomatic="${OPTARG}" ;;
		t) etetoolkit="${OPTARG}" ;;
		m) memory="${OPTARG}" ;;
		p) threads="${OPTARG}" ;;
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
if [[ -z $forward_reads || -z $reverse_reads || -z $pipeline \
|| -z $ncbi_nt_blast_db || -z $silva_blast_db || -z $ncbi_nt_kraken2_db \
|| -z $silva_kraken2_db || -z $etetoolkit || -z $silva_sortmerna_bac_lsu \
|| -z $silva_sortmerna_bac_ssu || -z $silva_sortmerna_arc_lsu \
|| -z $silva_sortmerna_arc_ssu || -z $silva_sortmerna_euk_lsu \
|| -z $silva_sortmerna_euk_ssu || -z $silva_sortmerna_rfam_5 \
|| -z $silva_sortmerna_rfam_5_8 || -z $silva_path_taxid || -z $ncbi_scientific \
|| -z $ncbi_non_scientific || -z $trimmomatic || -z $memory ]]; then
   echo -e "-1, -2, -P, -N, -S, -n, -s, -B, -b, -A, -a, -E, -e, -R, -r, -x, -F, -f, -T, -m, and -p must be set.\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script.\n"
   exit
fi

# Set pipeline tools to use
trimming=$(echo $pipeline | cut -f1 -d-)
sorting=$(echo $pipeline | cut -f2 -d-)
assembly=$(echo $pipeline | cut -f3 -d-)
mapping=$(echo $pipeline | cut -f4 -d-)
db=$(echo $pipeline | cut -f5 -d-)
classification=$(echo $pipeline | cut -f6 -d-)

##################### Write start time and options to output ######################

# Define starting time of script for total runtime calculation:
start=$(date +%s)
echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
echo -e "=================================================================\n\n"


# Output specified options:
echo -e "======== [$(date +%H:%M:%S)] OPTIONS [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

echo -e "Forward reads were defined as $forward_reads.\n"
echo -e "Reverse reads were defined as $reverse_reads.\n"
echo -e "Tools for the pipeline were set to $pipeline.\n"
echo -e "Number of threads was set to $threads.\n"
echo -e "Script started with full command: $cmd\n"



######################### Start of the actual script ################################
echo -e "++++++++ [$(date +%H:%M:%S)] START RUNNING SCRIPT [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"

# Activate the conda ete3 environment within this script to be able to run ete3.
# I found this solution # to activate conda environments in scripts here:
# https://github.com/conda/conda/issues/7980.
#val "$(conda shell.bash hook)" # Without this, the conda environment cannot be
# activated within the script
#conda activate ete3 # ete3 is our conda environemnt in which we installed ete3
# NOTE: outcommented to be run on graham, not needed

# Make output directory and directory for final files:
mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/
mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/
cd METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/

# Save full current path in variable to make navigation between directories easier:
base_directory=$(pwd)

######################### Step 1: trimming ################################

echo -e "++++++++ [$(date +%H:%M:%S)] START STEP 1: TRIMMING AND ERROR CORRECTION [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"
# Trimming is done with a separate subscript:
fastqc_on_R1_R2_and_optional_trimming.sh \
-T $trimmomatic -1 $forward_reads -2 $reverse_reads -t yes -p $threads -P $trimming
mv fastqc_on_R1_R2_and_optional_trimming_output/ step_1_trimming/

# Use a line from the script "fastqc_on_R1_R2_and_optional_trimming.sh" to
# generate the variable baseout and change into the generated directory:
baseout=${forward_reads%_*} # Make basename
cd step_1_trimming/trimmomatic/trimmed_at_phred_${trimming}_${baseout##*/}

# Running error correction module of SPAdes on all trimmed reads
echo -e "\n======== [$(date +%H:%M:%S)] ERROR-CORRECTING READS [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
spades.py -1 *1P.fastq -2 *2P.fastq \
--only-error-correction --disable-gzip-output -o error_correction \
-t $threads
mv error_correction/corrected/*1P*.fastq \
error_correction/corrected/*2P*.fastq .
# Rename cryptic name of error-corrected reads:
R1=$(echo *1P.00.0_0.cor.fastq) \
&& 	mv *1P.00.0_0.cor.fastq ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
sed -r -i 's/ BH:.{2,6}//g' ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
R2=$(echo *2P.00.0_0.cor.fastq) \
&& mv *2P.00.0_0.cor.fastq ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
sed -r -i 's/ BH:.{2,6}//g' ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
rm -r error_correction/
echo -e "\n======== [$(date +%H:%M:%S)] FINISHED ERROR-CORRECTING READS [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

echo -e "++++++++ [$(date +%H:%M:%S)] FINISHED STEP 1: TRIMMING AND ERROR CORRECTION [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"

######################### Step 2: rRNA sorting ################################

echo -e "++++++++ [$(date +%H:%M:%S)] START STEP 2: rRNA SORTING OF TRIMMED READS [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"

mkdir step_2_rrna_sorting/
cd step_2_rrna_sorting/

if [[ ${sorting} == "barrnap" || ${sorting} == "rrnafilter" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] CONVERT READS IN FASTA FORMAT [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	mkdir reads_in_fasta_format/
	fq2fa ../*1P_error_corrected.fastq reads_in_fasta_format/R1.fa
	fq2fa ../*2P_error_corrected.fastq reads_in_fasta_format/R2.fa
	echo -e "\n======== [$(date +%H:%M:%S)] READS TO FASTA CONVERSION DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
fi

if [[ ${sorting} == "sortmerna" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING SORTMERNA [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	mkdir SORTMERNA/
	sortmerna --ref $silva_sortmerna_bac_lsu \
	--ref $silva_sortmerna_bac_ssu \
	--ref $silva_sortmerna_arc_lsu \
	--ref $silva_sortmerna_arc_ssu \
	--ref $silva_sortmerna_euk_lsu \
	--ref $silva_sortmerna_euk_ssu \
	--ref $silva_sortmerna_rfam_5 \
	--ref $silva_sortmerna_rfam_5_8 \
  --reads ../*1P_error_corrected.fastq --reads ../*2P_error_corrected.fastq \
	--paired_in	--out2 -other -fastx 1 -num_alignments 1 -v -workdir SORTMERNA/ \
	--threads 1:1:$threads
	cd SORTMERNA/
	echo -e "\n======== [$(date +%H:%M:%S)] SORTMERNA DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ ${sorting} == "rrnafilter" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING rRNAFILTER [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	mkdir rRNAFILTER/
	cd rRNAFILTER/
	# rRNAFilter only worked for us when we started it within the directory
	# containing the .jar file. To simplify switching to that directory, we copy
	# it from its location to the pwd:
	cp -r ~/scratch/chris_pilot_project/programs/rRNAFilter .
	cd rRNAFilter/
	# We use 7GB for the rRNAFilter .jar, as shown in the rRNAFilter manual:
	java -jar -Xmx7g rRNAFilter_commandline.jar \
	-i ../../reads_in_fasta_format/R1.fa -r 0
	java -jar -Xmx7g rRNAFilter_commandline.jar \
	-i ../../reads_in_fasta_format/R2.fa -r 0
	mv ../../reads_in_fasta_format/R*.fa_rRNA ..
	cd ..
	rm -r rRNAFilter
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
	echo -e "\n======== [$(date +%H:%M:%S)] rRNAFILTER DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ ${sorting} == "barrnap" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING BARRNAP [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	mkdir BARRNAP/
	for kingdom in euk bac arc; do # barrnap needs to be run on kingdoms separately
		echo -e "\n======== [$(date +%H:%M:%S)] RUNNING BARRNAP ON KINGDOM $kingdom AND R1 READS [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
		~/scratch/chris_pilot_project/programs/barrnap/bin/barrnap \
		--quiet --lencutoff 0.000001 --reject 0.000001 --kingdom $kingdom \
		--threads $threads --outseq BARRNAP/${kingdom}_reads1.fa \
		reads_in_fasta_format/R1.fa
		echo -e "\n======== [$(date +%H:%M:%S)] RUNNING BARRNAP ON KINGDOM $kingdom AND R2 READS [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
		~/scratch/chris_pilot_project/programs/barrnap/bin/barrnap \
		--quiet --lencutoff 0.000001 --reject 0.000001 --kingdom $kingdom \
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
	echo -e "\n======== [$(date +%H:%M:%S)] BARRNAP DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ ${sorting} == "unsorted" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] MAKING FOLDER UNSORTED/ AND COPYING UNSORTED READS IN THERE TO KEEP THE FOLDER STRUCTURE CONSTANT [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	mkdir UNSORTED/
	cp ../*1P_error_corrected.fastq ../*2P_error_corrected.fastq UNSORTED/
	cd UNSORTED/
fi

echo -e "\n++++++++ [$(date +%H:%M:%S)] FINISHED STEP 2: rRNA SORTING [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"

######################### Step 3: Assembly ################################

echo -e "++++++++ [$(date +%H:%M:%S)] START STEP 3: ASSEMBLY [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"

mkdir step_3_assembly/
cd step_3_assembly/

if [[ $sorting == 'rrnafilter' ]]; then
	R1_sorted='rRNAFilter_paired_R1.fa'
	R2_sorted='rRNAFilter_paired_R2.fa'
elif [[ $sorting == 'sortmerna' ]]; then
	R1_sorted='out/aligned_fwd.fastq'
	R2_sorted='out/aligned_rev.fastq'
elif [[ $sorting == 'barrnap' ]]; then
	R1_sorted='barrnap_paired_R1.fa'
	R2_sorted='barrnap_paired_R2.fa'
elif [[ $sorting == 'unsorted' ]]; then
	R1_sorted='*1P_error_corrected.fastq'
	R2_sorted='*2P_error_corrected.fastq'
fi

if [[ $assembly == "spades" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING SPADES [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	mkdir SPADES/
	spades.py -1 ../$R1_sorted -2 ../$R2_sorted --only-assembler \
	-o SPADES/ -t $threads
	cd SPADES/
	echo -e "\n======== [$(date +%H:%M:%S)] SPADES DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ $assembly == "metaspades" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING METASPADES [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	mkdir METASPADES/
	spades.py --meta -1 ../$R1_sorted -2 ../$R2_sorted --only-assembler \
	-o METASPADES/ -t $threads
	cd METASPADES/
	echo -e "\n======== [$(date +%H:%M:%S)] METASPADES DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ $assembly == "megahit" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING MEGAHIT [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	megahit --presets meta-large -t $threads -1 ../$R1_sorted \
	-2 ../$R2_sorted -o MEGAHIT/
	cd MEGAHIT/
	echo -e "\n======== [$(date +%H:%M:%S)] MEGAHIT DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ $assembly == "idba_ud" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING IDBA_UD [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	# Note: we had to edit IDBA prior to compiling it because it didn't work
	# using long reads and the -l option. This seems to be a common problem and
	# can be circumvented following for example the instructions in
	# http://seqanswers.com/forums/showthread.php?t=29109, and see also
	# https://github.com/loneknightpy/idba/issues/26
	# IDBA_UD only takes interleaved fasta files
	fq2fa --merge --filter ../$R1_sorted ../$R2_sorted idba_ud_input.fa
	~/scratch/chris_pilot_project/programs/idba/bin/idba_ud \
	--num_threads $threads -r idba_ud_input.fa -o IDBA_UD/
	mv idba_ud_input.fa IDBA_UD/
	cd IDBA_UD/
	echo -e "\n======== [$(date +%H:%M:%S)] IDBA_UD DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ $assembly == "rnaspades" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING RNASPADES [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	mkdir RNASPADES/
	spades.py --rna -1 ../$R1_sorted -2 ../$R2_sorted --only-assembler \
	-o RNASPADES/ -t $threads
	cd RNASPADES/
	echo -e "\n======== [$(date +%H:%M:%S)] RNASPADES DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ $assembly == "idba_tran" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING IDBA_TRAN [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	# IDBA_TRAN only takes interleaved fasta files
	fq2fa --merge ../$R1_sorted ../$R2_sorted idba_tran_input.fa
	~/scratch/chris_pilot_project/programs/idba/bin/idba_tran \
	--num_threads $threads -l idba_tran_input.fa -o IDBA_TRAN/
	mv idba_tran_input.fa IDBA_TRAN/
	cd IDBA_TRAN
	echo -e "\n======== [$(date +%H:%M:%S)] IDBA_TRAN DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ $assembly == "trinity" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING TRINITY [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	# Barrnap and rRNAFilter output fasta files which has to be indicated to Trinity:
  if [[ $sorting == "rrnafilter" || $sorting == "barrnap" ]]; then
		Trinity --seqType fa --max_memory $memory --left ../$R1_sorted --right \
    ../$R2_sorted --CPU $threads --output TRINITY/ --NO_SEQTK
  else
    Trinity --seqType fq --max_memory $memory --left ../$R1_sorted --right \
    ../$R2_sorted --CPU $threads --output TRINITY/ --NO_SEQTK
  fi
  cat TRINITY/Trinity.fasta | sed 's/ len/_len/g' \
	> TRINITY/Trinity_with_length.fasta  # Edit for universal format
	cd TRINITY/
	echo -e "\n======== [$(date +%H:%M:%S)] TRINITY DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ $assembly == "transabyss" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING TRANSABYSS [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
  transabyss --pe ../$R1_sorted ../$R2_sorted --threads $threads \
  --outdir TRANSABYSS/
  sed 's/ /_/g' TRANSABYSS/transabyss-final.fa \
	> TRANSABYSS/transabyss-final_edited.fa # Edit for universal format
	cd TRANSABYSS/
	echo -e "\n======== [$(date +%H:%M:%S)] TRANSABYSS DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
fi

echo -e "\n++++++++ [$(date +%H:%M:%S)] FINISHED STEP 3: ASSEMBLY [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"

######################## Step 4: Mapping ################################

echo -e "++++++++ [$(date +%H:%M:%S)] START STEP 4: MAPPING [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"

mkdir step_4_mapping/
mkdir step_4_mapping/$(echo $mapping | tr '[:lower:]' '[:upper:]')/
cd step_4_mapping/$(echo $mapping | tr '[:lower:]' '[:upper:]')/

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
elif [[ $assembly == 'transabyss' ]]; then
	scaffolds='transabyss-final_edited.fa'
fi

if [[ $mapping == 'bwa' ]]; then
  echo -e "\n======== [$(date +%H:%M:%S)] Starting bwa index [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
  bwa index -p bwa_index ../../$scaffolds
  echo -e "\n======== [$(date +%H:%M:%S)] bwa index complete. Starting bwa mem [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
  bwa mem -t $threads bwa_index ../../../../../../*1P_error_corrected.fastq \
	../../../../../../*2P_error_corrected.fastq > ${mapping}_output.sam
  rm bwa_index*
	echo -e "\n======== [$(date +%H:%M:%S)] bwa mem complete [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
elif [[ $mapping == 'bowtie2' ]]; then
  echo -e "\n======== [$(date +%H:%M:%S)] Starting bowtie2 index [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	bowtie2-build -f ../../$scaffolds bowtie_index
	echo -e "\n======== [$(date +%H:%M:%S)] bowtie2 index complete. Starting bowtie2 [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	bowtie2 -q -x bowtie_index -1 ../../../../../../*1P_error_corrected.fastq \
	-2 ../../../../../../*2P_error_corrected.fastq -S ${mapping}_output.sam \
	-p $threads
	rm bowtie_index*
  echo -e "\n======== [$(date +%H:%M:%S)] bowtie2 complete [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
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

# Moving back to assembler directory:
cd ../../

echo -e "++++++++ [$(date +%H:%M:%S)] FINISHED STEP 4: MAPPING [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"

######################### Steps 5 and 6.1: Picking a referencd DB and taxonomic classification ################################

echo -e "++++++++ [$(date +%H:%M:%S)] START STEP 5 AND 6.1: CLASSIFICATION OF ASSEMBLED SCAFFOLDS WITH $(echo $db | tr '[:lower:]' '[:upper:]') DATABASE +++++++\n"

mkdir step_5_reference_DB/
mkdir step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/
mkdir step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/step_6_classification/
mkdir step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/step_6_classification/$(echo $classification | tr '[:lower:]' '[:upper:]')/
cd step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/step_6_classification/$(echo $classification | tr '[:lower:]' '[:upper:]')/


if [[ $db == "silva" ]]; then
	krakenDB=$silva_kraken2_db
	blastDB=$silva_blast_db
elif [[ $db == 'ncbi_nt' ]]; then
	krakenDB=$ncbi_nt_kraken2_db
	blastDB=$ncbi_nt_blast_db
fi

if [[ $classification == "blast_first_hit" || $classification == "blast_filtered" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING BLAST WITH DATABASE $blastDB [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	blastn -query ../../../../$scaffolds -db $blastDB -out blast_output.txt \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
	-evalue 1e-05 -num_threads $threads
	echo -e "\n======== [$(date +%H:%M:%S)] BLAST WITH DATABASE $blastDB DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
fi

if [[ $classification == "blast_first_hit" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING BLAST FIRST HIT [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	# We run a separate script to filter the BLAST results:
	blast_filtering.bash -i blast_output.txt -f blast -t soft -T $threads -e $etetoolkit
	mv blast_output.txt blast_filtering_results/
	echo -e "\n======== [$(date +%H:%M:%S)] BLAST FIRST HIT DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ $classification == "blast_filtered" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING BLAST FILTERED [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	# We run a separate script to filter the BLAST results:
	blast_filtering.bash -i blast_output.txt -f blast -t strict -T $threads -e $etetoolkit
	mv blast_output.txt blast_filtering_results/
	echo -e "\n======== [$(date +%H:%M:%S)] BLAST FILTERED DONE========\n"

elif [[ $classification == "kraken2" ]]; then
	echo -e "\n======== [$(date +%H:%M:%S)] RUNNING KRAKEN2 WITH DATABASE $krakenDB [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
	# Run kraken2
	echo $krakenDB
	kraken2 --db $krakenDB --threads $threads ../../../../$scaffolds \
	> kraken2_output.txt
	if [[ $db == 'silva' ]]; then
    # Now we're gonna edit the output so that it has the same format as
		# CREST output, since we already have a script to deal with
		# SILVA CREST output/taxonomy
    # Extract the taxids column of the standard kraken output:
		cut -f3 kraken2_output.txt > kraken2_taxids.txt

    # Merge your kraken2 taxids with the SILVA path file to assign a SILVA
		# taxonomy path to every kraken2 hit:
    mergeFilesOnColumn.pl $silva_path_taxid kraken2_taxids.txt 2 1 > merged.txt
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
		# into NCBI taxids, and use that script on the formatted kraken2 output:
    assign_NCBI_staxids_to_CREST_v4.py $ncbi_scientific $ncbi_non_scientific \
    kraken2_SILVA_formatted.txt kraken2_SILVA_formatted_with_NCBI_taxids.txt
    mergeFilesOnColumn.pl kraken2_SILVA_formatted_with_NCBI_taxids.txt \
    kraken2_SILVA_formatted.txt 1 1 | cut -f3 > NCBItaxids.txt # Merge SILVA output with taxids and extract taxids
		# We use a separate script to assign taxonomy to NCBI taxids:
    assign_taxonomy_to_NCBI_staxids.sh -b NCBItaxids.txt -c 1	-e $etetoolkit
    sed -i '1d' NCBItaxids_with_taxonomy.txt # Remove header
    cut -f2 kraken2_output.txt > contig_names.txt # Get contig names from original kraken2 output
   	paste contig_names.txt NCBItaxids_with_taxonomy.txt \
	 	> contigs_with_NCBItaxids_and_taxonomy.txt # Add contig names to taxonomy file
   	echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
   	> kraken2_final.txt && cat contigs_with_NCBItaxids_and_taxonomy.txt | sed 's/Unknown/NA/g' \
  	>> kraken2_final.txt # Add header

    # Sort files
   	mkdir intermediate_files
    mv kraken2_output.txt kraken2_taxids.txt merged* names.txt \
		kraken2_SILVA_formatted* NCBItaxids* contig* intermediate_files/

  elif [[ $db == 'ncbi_nt' ]]; then
  	cut -f 2-3 kraken2_output.txt > kraken2_output_contig_taxid.txt # Isolate contig names and taxids
		# We use a separate script to assign taxonomy to NCBI taxids:
    assign_taxonomy_to_NCBI_staxids.sh -b kraken2_output_contig_taxid.txt \
		-c 2 -e $etetoolkit
    sed -i '1d' kraken2_output_contig_taxid_with_taxonomy.txt # Remove header
    echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
    > kraken2_final.txt && cat kraken2_output_contig_taxid_with_taxonomy.txt \
		| sed 's/Unknown/NA/g' >> kraken2_final.txt # Add header

   	# Sort files
    mkdir intermediate_files
    mv kraken2_output* intermediate_files/
  fi
echo -e "\n======== [$(date +%H:%M:%S)] KRAKEN2 WITH DATABASE $krakenDB DONE [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
fi

######################### Step 6.2: Generating final putput files ################################

echo -e "++++++++ [$(date +%H:%M:%S)] START STEP 6.2: GENERATING FINAL OUTPUT FILES [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"

# Each assembler/classification tool output has a different format. We
# make that format universal with the following code.

mkdir FINAL_FILES/
mkdir FINAL_FILES/intermediate_files/
cd FINAL_FILES/intermediate_files/

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
	merge_on_outer.py ../../../../../../step_4_mapping/$(echo $mapping | tr '[:lower:]' '[:upper:]')/merge_input_mapped_${mapping}.txt \
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
	merge_on_outer.py ../../../../../../step_4_mapping/$(echo $mapping | tr '[:lower:]' '[:upper:]')/merge_input_mapped_${mapping}.txt \
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
	merge_on_outer.py ../../../../../../step_4_mapping/$(echo $mapping | tr '[:lower:]' '[:upper:]')/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt

elif [[ $assembly == 'trinity' ]]; then
	# Change the sequences' fasta format to tab delimited:
	fasta_to_tab ../../../../../../${scaffolds} \
	| sed 's/ /\t/1' | cut -f1,3 > tmp
	# Add header name so that we can later merge on outer with a python script:
	echo -e "sequence_name\tsequence" > ${assembly}_tab_to_merge.txt \
	&& cat tmp >> ${assembly}_tab_to_merge.txt
	# And remove intermediate file:
	rm tmp
	# Now the scaffold names can be merged with the mapper output:
	merge_on_outer.py ../../../../../../step_4_mapping/$(echo $mapping | tr '[:lower:]' '[:upper:]')/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt
fi

# Add classification data to the mapper output and do
# final individual edits:

if [[ $classification == 'blast_filtered' ]]; then

	# We run a simple python script because we need to merge on "outer":
	merge_on_outer.py \
	../../blast_filtering_results/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt \
	${assembly}_final_${mapping}_merge_ready.txt \
	trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	if [[ $assembly == 'spades' || $assembly == 'metaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/scaffold_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-17 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'rnaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
		| cut -f1,2,3,5-21 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		|	sed 's/contig-[0-9]*_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $16, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $17}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' \
		> trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'transabyss' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt
	fi

	echo -e "\n======== [$(date +%H:%M:%S)] DONE FINALIZING BLAST_FILTERED FILES =======\n"

elif [[ $classification == 'blast_first_hit' ]]; then

	# We run a simple python script because we need to merge on "outer":
	merge_on_outer.py ../../blast_filtering_results/blast_output_with_taxonomy_and_best_hit.txt \
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
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/scaffold_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-17 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
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
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		|	sed 's/contig-[0-9]*_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $16, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $17}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'transabyss' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp
	fi

	echo -e "\n======== [$(date +%H:%M:%S)] DONE FINALIZING BLAST_FIRST_HIT FILES [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"

elif [[ $classification == 'kraken2' ]]; then

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
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $17}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-19 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tsequence_length\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $17, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $18}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
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
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $17, $18, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $19}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $5, $4, $17, $18}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'transabyss' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \								> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp
	fi

	echo -e "\n======== [$(date +%H:%M:%S)] DONE FINALIZING KRAKEN2 FILES [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n"
fi

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"
