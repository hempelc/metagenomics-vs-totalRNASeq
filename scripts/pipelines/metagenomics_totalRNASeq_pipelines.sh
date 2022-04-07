#!/bin/bash

# Version 0.1
# Written by Natalie Wright (nwrigh06@uoguelph.ca) and Chris Hempel (hempelc@uoguelph.ca)

# This script was designed for Chris Hempel's first PhD chapter to run any combination
# of specified data-processing tools on metagenomcis and total RNA-Seq data

# The output is a folder named after the specified pipeline
# that contains tab-separated, taxonomically annotated scaffolds and scaffold coverage
# for the specified pipeline.

# The pipleine that is to be run should be given in the following format:
# "PHREDScore-rRNASortingTool-AssemblyTool-Mapper-ClassificationTool"
# You can choose from the following options:
	# PHREDScore: 5, 10, 15, 20
	# rRNASortingTool: barrnap, rrnafilter, sortmerna, unsorted
	# AssemblyTool: megahit, metaspades, spades, rnaspades, trinity, idba_ud, idba_tran, transabyss
	# Mapper: bwa, bowtie2
	# ClassificationTool: blast_filtered, blast_first_hit, kraken2


# The pipelines require the following subscripts, which are all located in the
# subscripts/ directory:
	# assign_taxonomy_to_NCBI_staxids.sh, blast_filter.py, fasta_to_tab,
	# fastqc_on_R1_R2_and_optional_trimming.sh, filter-fasta.awk,
	# deinterleave_fastq_reads.sh, LookupTaxonDetails3.py, merge_on_outer.py,
	# mergeFilesOnColumn.pl


# To run every possible combination of tools, the pipelines require the following
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

# The pipelines also require the file <PATH/TO/.etetoolkit/taxa.sqlite>,
# which is generated in the home directory when ete3 is run for the first time.
# Note that we installed ete3 in a conda environment and activate the conda
# environment within this script.

# rRNAFilter only worked for us when we started it within the directory
# containing the rRNAFilter.jar file. Therefore, we cope it from it's current
# location within the script.

# Adapt default options below as necessary.


cmd="$0 $@" # Make variable containing the entire entered command to print command to logfile
usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq> \
-P <pipeline> -S <SILVA BLAST database> -s <SILVA kraken2 database> \
-B <SILVA SortMeRNA LSU bacteria database> -b <SILVA SortMeRNA SSU bacteria database> \
-A <SILVA SortMeRNA LSU archaea database> -a <SILVA SortMeRNA SSU archaea database> \
-E <SILVA SortMeRNA LSU eukaryota database> -e <SILVA SortMeRNA SSU eukaryota database> \
-R <SILVA SortMeRNA rfam 5.8S database> -r <SILVA SortMeRNA rfam 5S database> \
-T <PATH/TO/trimmomatic-<version>.jar)> -i <PATH/TO/rRNAFilter directory)>\
-t <PATH/TO/.etetoolkit/taxa.sqlite> -m <nnnG> [-p <n>]

Usage:
	-1 Full path to forward reads in .fastq/.fq format
	-2 Full path to reverse reads in .fastq/.fq format
	-P Pipeline tools
	-S Path to SILVA database for BLAST
	-s Path to SILVA database for kraken2
	-B Path to bacteria LSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-bac-23s-id98.fasta)
	-b Path to bacteria SSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-bac-16s-id90.fasta)
	-A Path to archaea LSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-arc-23s-id98.fasta)
	-a Path to archaea SSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-arc-16s-id95.fasta)
  -E Path to eukaryota LSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-euk-28s-id98.fasta)
	-e Path to eukaryota SSU SILVA database for SortMeRNA (comes with SortMeRNA, silva-euk-18s-id95.fasta)
	-R Path to rfam 5.8S database for SortMeRNA (comes with SortMeRNA, rfam-5.8s-database-id98.fasta)
	-r Path to rfam 5S database for SortMeRNA (comes with SortMeRNA, rfam-5s-database-id98.fasta)
	-T Path to trimmomatic application (trimmomatic-<version>.jar)
	-t Path to .etetoolkit/taxa.sqlite
	-i Path to rRNAFilter directory
	-m Maximum memory (format: XXXG, where XXX is a numerical value for teh emmory in Gigabyte)
	-p Number of threads (default:16)
	-h Display this help and exit"

# Set default options:
threads='16'
pipeline="20-unsorted-spades-bwa-silva-kraken2"
silva_blast_db="/hdd2/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Jul_2021/blastdb"
silva_kraken2_db="/hdd2/databases/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Jul_2021/"
silva_sortmerna_bac_lsu="/hdd2/databases/sortmerna_silva_databases/silva-bac-23s-id98.fasta"
silva_sortmerna_bac_ssu="/hdd2/databases/sortmerna_silva_databases/silva-bac-16s-id90.fasta"
silva_sortmerna_arc_lsu="/hdd2/databases/sortmerna_silva_databases/silva-arc-23s-id98.fasta"
silva_sortmerna_arc_ssu="/hdd2/databases/sortmerna_silva_databases/silva-arc-16s-id95.fasta"
silva_sortmerna_euk_lsu="/hdd2/databases/sortmerna_silva_databases/silva-euk-28s-id98.fasta"
silva_sortmerna_euk_ssu="/hdd2/databases/sortmerna_silva_databases/silva-euk-18s-id95.fasta"
silva_sortmerna_rfam_5="/hdd2/databases/sortmerna_silva_databases/rfam-5s-database-id98.fasta"
silva_sortmerna_rfam_5_8="/hdd2/databases/sortmerna_silva_databases/rfam-5.8s-database-id98.fasta"
trimmomatic="/hdd1/programs_for_pilot/Trimmomatic-0.39/trimmomatic-0.39.jar"
etetoolkit="/hdd1/chempel/.etetoolkit/taxa.sqlite"
rrnafil="/hdd1/programs_for_pilot/rRNAFilter"
memory="120G"


# Set specified options:
while getopts ':1:2:P:N:S:n:s:B:b:A:a:E:e:R:r:x:F:f:T:t:i:m:p:h' opt; do
 	case "${opt}" in
		1) forward_reads="${OPTARG}" ;;
		2) reverse_reads="${OPTARG}" ;;
		P) pipeline="${OPTARG}" ;;
		S) silva_blast_db="${OPTARG}" ;;
		s) silva_kraken2_db="${OPTARG}" ;;
		B) silva_sortmerna_bac_lsu="${OPTARG}" ;;
		b) silva_sortmerna_bac_ssu="${OPTARG}" ;;
		A) silva_sortmerna_arc_lsu="${OPTARG}" ;;
		a) silva_sortmerna_arc_ssu="${OPTARG}" ;;
		E) silva_sortmerna_euk_lsu="${OPTARG}" ;;
		e) silva_sortmerna_euk_ssu="${OPTARG}" ;;
		R) silva_sortmerna_rfam_5="${OPTARG}" ;;
		r) silva_sortmerna_rfam_5_8="${OPTARG}" ;;
		T) trimmomatic="${OPTARG}" ;;
		t) etetoolkit="${OPTARG}" ;;
		i) rrnafil="${OPTARG}" ;;
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
|| -z $silva_blast_db || -z $silva_kraken2_db || -z $silva_sortmerna_bac_lsu \
|| -z $silva_sortmerna_bac_ssu || -z $silva_sortmerna_arc_lsu \
|| -z $silva_sortmerna_arc_ssu || -z $silva_sortmerna_euk_lsu \
|| -z $silva_sortmerna_euk_ssu || -z $silva_sortmerna_rfam_5 \
|| -z $silva_sortmerna_rfam_5_8 || -z $trimmomatic || -z $etetoolkit || -z $rrnafil \
|| -z $memory ]]; then
   echo -e "-1, -2, -P, -S, -s, -B, -b, -A, -a, -E, -e, -R, -r, -T, -t, -i, -m, and -p must be set.\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script.\n"
   exit
fi

# Set pipeline tools to use
trimming=$(echo $pipeline | cut -f1 -d-)
sorting=$(echo $pipeline | cut -f2 -d-)
assembly=$(echo $pipeline | cut -f3 -d-)
mapping=$(echo $pipeline | cut -f4 -d-)
database=$(echo $pipeline | cut -f5 -d-)
classification=$(echo $pipeline | cut -f6 -d-)

# Define functions to print steps with time
start=$(date +%s)

step_description_and_time_first () {
	echo -e "\n++++++++ [$(date +%H:%M:%S)] ${1} [Runtime: $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"
}

step_description_and_time_second () {
	echo -e "\n======== [$(date +%H:%M:%S)] ${1} [Runtime: $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n" #" adding outcommented quote here to fix bug in colouring scheme of personal text editor
}

# Make output directory and directory for final files:
mkdir "${pipeline}"/
cd "${pipeline}"/

# Make open bracket to later tell script to write everything that follows into a logfile:
(

##################### Write start time and options to output ######################

# Define starting time of script for total runtime calculation:
echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
echo -e "=================================================================\n\n"


# Output specified options:
step_description_and_time_first "OPTIONS"

echo -e "Forward reads were defined as $forward_reads.\n"
echo -e "Reverse reads were defined as $reverse_reads.\n"
echo -e "Tools for the pipeline were set to $pipeline.\n"
echo -e "Number of threads was set to $threads.\n"
echo -e "Script started with full command: $cmd\n"



######################### Start of the actual script ################################
step_description_and_time_first "START RUNNING SCRIPT"

# Activate the conda ete3 environment within this script to be able to run ete3.
# Based on https://github.com/conda/conda/issues/7980
source /hdd1/chempel/programs/miniconda/etc/profile.d/conda.sh
conda activate ete3 # ete3 is our conda environemnt in which we installed ete3

# Save full current path in variable to make navigation between directories easier:
base_directory=$(pwd)

######################### Step 1: trimming ################################

step_description_and_time_first "START STEP 1: TRIMMING AND ERROR CORRECTION"
# Trimming is done with a separate subscript:
fastqc_on_R1_R2_and_optional_trimming.sh \
-T $trimmomatic -1 $forward_reads -2 $reverse_reads -t yes -p $threads -P $trimming
mv fastqc_on_R1_R2_and_optional_trimming_output/ step_1_trimming/

# Use a line from the script "fastqc_on_R1_R2_and_optional_trimming.sh" to
# generate the variable baseout and change into the generated directory:
baseout=${forward_reads%_*} # Make basename
cd step_1_trimming/trimmomatic/trimmed_at_phred_${trimming}_${baseout##*/}

# Running error correction module of SPAdes on all trimmed reads
step_description_and_time_first "ERROR-CORRECTING READS"
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
step_description_and_time_first "FINISHED ERROR-CORRECTING READS"

step_description_and_time_first "FINISHED STEP 1: TRIMMING AND ERROR CORRECTION"

######################### Step 2: rRNA sorting ################################

step_description_and_time_first "START STEP 2: rRNA SORTING OF TRIMMED READS"

mkdir step_2_rrna_sorting/
cd step_2_rrna_sorting/

if [[ ${sorting} == "barrnap" || ${sorting} == "rrnafilter" ]]; then
	step_description_and_time_first "CONVERT READS IN FASTA FORMAT"
	mkdir reads_in_fasta_format/
	fq2fa ../*1P_error_corrected.fastq reads_in_fasta_format/R1.fa
	fq2fa ../*2P_error_corrected.fastq reads_in_fasta_format/R2.fa
	step_description_and_time_first "READS TO FASTA CONVERSION DONE"
fi

if [[ ${sorting} == "sortmerna" ]]; then
	step_description_and_time_first "RUNNING SORTMERNA"
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
	step_description_and_time_first "SORTMERNA DONE"

elif [[ ${sorting} == "rrnafilter" ]]; then
	step_description_and_time_first "RUNNING rRNAFILTER"
	# rRNAFilter only worked for us when we started it within the directory
	# containing the .jar file. To simplify switching to that directory, we copy
	# it from its location to the pwd:
	cp -r $rrnafil .
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
	step_description_and_time_first "rRNAFILTER DONE"

elif [[ ${sorting} == "barrnap" ]]; then
	step_description_and_time_first "RUNNING BARRNAP"
	mkdir BARRNAP/
	for kingdom in euk bac arc; do # barrnap needs to be run on kingdoms separately
		step_description_and_time_first "RUNNING BARRNAP ON KINGDOM $kingdom AND R1 READS"
		barrnap \
		--quiet --lencutoff 0.000001 --reject 0.000001 --kingdom $kingdom \
		--threads $threads --outseq BARRNAP/${kingdom}_reads1.fa \
		reads_in_fasta_format/R1.fa
		step_description_and_time_first "RUNNING BARRNAP ON KINGDOM $kingdom AND R2 READS"
		barrnap \
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
	step_description_and_time_first "BARRNAP DONE"

elif [[ ${sorting} == "unsorted" ]]; then
	step_description_and_time_first "MAKING FOLDER UNSORTED/ AND COPYING UNSORTED READS IN THERE TO KEEP THE FOLDER STRUCTURE CONSTANT"
	mkdir UNSORTED/
	cp ../*1P_error_corrected.fastq ../*2P_error_corrected.fastq UNSORTED/
	cd UNSORTED/
fi

step_description_and_time_first "FINISHED STEP 2: rRNA SORTING"

######################### Step 3: Assembly ################################

step_description_and_time_first "START STEP 3: ASSEMBLY"

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
	step_description_and_time_first "RUNNING SPADES"
	mkdir SPADES/
	spades.py -1 ../$R1_sorted -2 ../$R2_sorted --only-assembler \
	-o SPADES/ -t $threads
	cd SPADES/
	step_description_and_time_first "SPADES DONE"

elif [[ $assembly == "metaspades" ]]; then
	step_description_and_time_first "RUNNING METASPADES"
	mkdir METASPADES/
	spades.py --meta -1 ../$R1_sorted -2 ../$R2_sorted --only-assembler \
	-o METASPADES/ -t $threads
	cd METASPADES/
	step_description_and_time_first "METASPADES DONE"

elif [[ $assembly == "megahit" ]]; then
	step_description_and_time_first "RUNNING MEGAHIT"
	megahit --presets meta-large -t $threads -1 ../$R1_sorted \
	-2 ../$R2_sorted -o MEGAHIT/
	cd MEGAHIT/
	step_description_and_time_first "MEGAHIT DONE"

elif [[ $assembly == "idba_ud" ]]; then
	step_description_and_time_first "RUNNING IDBA_UD"
	# Note: we had to edit IDBA prior to compiling it because it didn't work
	# using long reads and the -l option. This seems to be a common problem and
	# can be circumvented following for example the instructions in
	# http://seqanswers.com/forums/showthread.php?t=29109, and see also
	# https://github.com/loneknightpy/idba/issues/26
	# IDBA_UD only takes interleaved fasta files
	fq2fa --merge --filter ../$R1_sorted ../$R2_sorted idba_ud_input.fa
	idba_ud --num_threads $threads -r idba_ud_input.fa -o IDBA_UD/
	mv idba_ud_input.fa IDBA_UD/
	cd IDBA_UD/
	step_description_and_time_first "IDBA_UD DONE"

elif [[ $assembly == "rnaspades" ]]; then
	step_description_and_time_first "RUNNING RNASPADES"
	mkdir RNASPADES/
	spades.py --rna -1 ../$R1_sorted -2 ../$R2_sorted --only-assembler \
	-o RNASPADES/ -t $threads
	cd RNASPADES/
	step_description_and_time_first "RNASPADES DONE"

elif [[ $assembly == "idba_tran" ]]; then
	step_description_and_time_first "RUNNING IDBA_TRAN"
	# IDBA_TRAN only takes interleaved fasta files
	fq2fa --merge ../$R1_sorted ../$R2_sorted idba_tran_input.fa
	idba_tran --num_threads $threads -l idba_tran_input.fa -o IDBA_TRAN/
	mv idba_tran_input.fa IDBA_TRAN/
	cd IDBA_TRAN
	step_description_and_time_first "IDBA_TRAN DONE"

elif [[ $assembly == "trinity" ]]; then
	step_description_and_time_first "RUNNING TRINITY"
	# Barrnap and rRNAFilter output fasta files which has to be indicated to Trinity:
  if [[ $sorting == "rrnafilter" || $sorting == "barrnap" ]]; then
		Trinity --seqType fa --max_memory 20G --left ../$R1_sorted --right \
    ../$R2_sorted --CPU $threads --output TRINITY/ --NO_SEQTK
  else
    Trinity --seqType fq --max_memory 20G --left ../$R1_sorted --right \
    ../$R2_sorted --CPU $threads --output TRINITY/ --NO_SEQTK
  fi
  cat TRINITY/Trinity.fasta | sed 's/ len/_len/g' \
	> TRINITY/Trinity_with_length.fasta  # Edit for universal format
	cd TRINITY/
	step_description_and_time_first "TRINITY DONE"

elif [[ $assembly == "transabyss" ]]; then
	step_description_and_time_first "RUNNING TRANSABYSS"
  transabyss --pe ../$R1_sorted ../$R2_sorted --threads $threads \
  --outdir TRANSABYSS/
  sed 's/ /_/g' TRANSABYSS/transabyss-final.fa \
	> TRANSABYSS/transabyss-final_edited.fa # Edit for universal format
	cd TRANSABYSS/
	step_description_and_time_first "TRANSABYSS DONE"
fi

step_description_and_time_first "FINISHED STEP 3: ASSEMBLY"

######################## Step 4: Mapping ################################

step_description_and_time_first "START STEP 4: MAPPING"

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
  step_description_and_time_first "Starting bwa index"
  bwa index -p bwa_index ../../$scaffolds
  step_description_and_time_first "bwa index complete. Starting bwa mem"
  bwa mem -t $threads bwa_index ../../../../../../*1P_error_corrected.fastq \
	../../../../../../*2P_error_corrected.fastq > ${mapping}_output.sam
	step_description_and_time_first "bwa mem complete"
elif [[ $mapping == 'bowtie2' ]]; then
  step_description_and_time_first "Starting bowtie2 index"
	bowtie2-build -f ../../$scaffolds bowtie_index
	step_description_and_time_first "bowtie2 index complete. Starting bowtie2"
	bowtie2 -q -x bowtie_index -1 ../../../../../../*1P_error_corrected.fastq \
	-2 ../../../../../../*2P_error_corrected.fastq -S ${mapping}_output.sam \
	-p $threads
  step_description_and_time_first "bowtie2 complete"
fi

# Editing the mapper outputs to detemine coverge of each contig:
samtools sort ${mapping}_output.sam > ${mapping}_output_sorted.sam
samtools coverage ${mapping}_output_sorted.sam | cut -f 1,7 | tail -n +2 \
| head -n -1 | awk '{ print $2, $1}' OFS=$'\t' > coverage_${mapping}.tsv
echo -e "coverage\tsequence_name" > merge_input_mapped_${mapping}.txt \
&& cat coverage_${mapping}.tsv >> merge_input_mapped_${mapping}.txt

rm ${mapping}_output* *index* coverage_${mapping}.tsv

# Moving back to assembler directory:
cd ../../

step_description_and_time_first "FINISHED STEP 4: MAPPING"

######################### Steps 5 and 6.1: Picking a referencd DB and taxonomic classification ################################

step_description_and_time_first "START STEP 5.1: CLASSIFICATION OF ASSEMBLED SCAFFOLDS"


mkdir step_5_classification/
mkdir step_5_classification/$(echo $classification | tr '[:lower:]' '[:upper:]')/
cd step_5_classification/$(echo $classification | tr '[:lower:]' '[:upper:]')/


if [[ $classification == "blast_first_hit" || $classification == "blast_filtered" ]]; then
	step_description_and_time_first "RUNNING BLAST WITH DATABASE $silva_blast_db"
	blastn -query ../../../../$scaffolds -db $silva_blast_db -out blast_output.txt \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
	-evalue 1e-05 -num_threads $threads
	step_description_and_time_first "BLAST WITH DATABASE $silva_blast_db DONE"
fi

if [[ $classification == "blast_first_hit" ]]; then
	step_description_and_time_first "RUNNING BLAST FIRST HIT"
	# We run a separate script to add taxonomy and filter the BLAST results:
	assign_taxonomy_to_NCBI_staxids.sh -b blast_output.txt -c 13 \
	-e $etetoolkit
	sed -i 's/Unknown/NA/g' blast_output_with_taxonomy.txt
	conda deactivate # ete3 env incompatible with blast filter script
	blast_filter.py blast_output_with_taxonomy.txt soft
	step_description_and_time_first "BLAST FIRST HIT DONE"

elif [[ $classification == "blast_filtered" ]]; then
	step_description_and_time_first "RUNNING BLAST FILTERED"
	# We run a separate script to add taxonomy and filter the BLAST results:
	assign_taxonomy_to_NCBI_staxids.sh -b blast_output.txt -c 13 \
	-e $etetoolkit
	sed -i 's/Unknown/NA/g' blast_output_with_taxonomy.txt
	conda deactivate # ete3 env incompatible with blast filter script
  blast_filter.py blast_output_with_taxonomy.txt strict
	step_description_and_time_first "BLAST FILTERED DONE"

elif [[ $classification == "kraken2" ]]; then
	step_description_and_time_first "RUNNING KRAKEN2 WITH DATABASE $silva_kraken2_db"
	# Run kraken2
	kraken2 --db $silva_kraken2_db --threads $threads ../../../../$scaffolds \
	> kraken2_output.txt

	cut -f 2-3 kraken2_output.txt > kraken2_output_contig_taxid.txt # Isolate contig names and taxids
	# We use a separate script to assign taxonomy to NCBI taxids:
	assign_taxonomy_to_NCBI_staxids.sh -b kraken2_output_contig_taxid.txt \
	-c 2 -e $etetoolkit
	sed -i '1d' kraken2_output_contig_taxid_with_taxonomy.txt # Remove header
	sed -i 's/Unknown/NA/g' kraken2_output_contig_taxid_with_taxonomy.txt # Change unknown to NA

	# Turn lowest_hit into species
	echo -e "sequence_name\tstaxid\tlowest_rank\tspecies\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
	> kraken2_final.txt
	while read line; do
		pre=$(cut -f 1-3 <<< $line)
		spec=$(cut -f 4 <<< $line | cut -f 1-2 -d ' ') # Cut down to first two words
		post=$(cut -f 5- <<< $line)
		if [[ "${spec:0:1}" =~ [a-z] ]]; then # if first letter is not capitalized (not in format "Genus species")
			spec="NA"
		elif [[ $(wc -w <<< "${spec}") != 2 ]]; then # if only one word (not in format "Genus species")
			spec="NA"
		fi
		echo -e "${pre}\t${spec}\t${post}" >> kraken2_final.txt
	done <kraken2_output_contig_taxid_with_taxonomy.txt

	# Sort files
	mkdir intermediate_files
	mv kraken2_output* intermediate_files/

step_description_and_time_first "KRAKEN2 WITH DATABASE $silva_kraken2_db DONE"
fi

######################### Step 6.2: Generating final putput files ################################

step_description_and_time_first "START STEP 6.2: GENERATING FINAL OUTPUT FILES"

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
	merge_on_outer.py ../../blast_filtered.txt \
	${assembly}_final_${mapping}_merge_ready.txt \
	trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt

	if [[ $assembly == 'spades' || $assembly == 'metaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/scaffold_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-17 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tsequence_length\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'rnaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
		| cut -f1,2,3,5-21 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		|	sed 's/contig-[0-9]*_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $16, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $17}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' \
		> trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt

	elif [[ $assembly == 'transabyss' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt
	fi

	step_description_and_time_first "DONE FINALIZING BLAST_FILTERED FILES"

elif [[ $classification == 'blast_first_hit' ]]; then

	# We run a simple python script because we need to merge on "outer":
	merge_on_outer.py ../../blast_filtered.txt \
	${assembly}_final_${mapping}_merge_ready.txt \
	trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt

	if [[ $assembly == 'spades' || $assembly == 'metaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/scaffold_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-17 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tsequence_length\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'rnaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
		| cut -f1,2,3,5-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		|	sed 's/contig-[0-9]*_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $16, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $17}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'transabyss' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcoverage\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp
	fi

	step_description_and_time_first "DONE FINALIZING BLAST_FIRST_HIT FILES"

elif [[ $classification == 'kraken2' ]]; then

	# We run a simple python script because we need to merge on "outer":
	merge_on_outer.py ../../kraken2_final.txt \
	${assembly}_final_${mapping}_merge_ready.txt \
	trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt

	if [[ $assembly == 'spades' || $assembly == 'metaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tstaxid\tlowest_rank\tspecies\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcoverage\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tstaxid\tlowest_rank\tspecies\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcoverage\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $17}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-19 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tstaxid\tlowest_rank\tspecies\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcoverage\tsequence_length\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $17, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $18}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'rnaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/2' |sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' \
		| cut -f1,2,3,5-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tstaxid\tlowest_rank\tspecies\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcoverage\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-20 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tstaxid\tlowest_rank\tspecies\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcoverage\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $17, $18, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $4, $2, $16, $19}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tstaxid\tlowest_rank\tspecies\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcoverage\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $5, $4, $17, $18}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'transabyss' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt
		echo -e "sequence_name\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tspecies\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcoverage\tassembly_sequence" \								> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $6, $4, $18, $19}' tmp \
		> ${base_directory}/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${classification}_final.txt \
		&& rm tmp
	fi

	step_description_and_time_first "DONE FINALIZING KRAKEN2 FILES"
fi

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"

# Write output to console and log file
) 2>&1 | tee "${pipeline}"/log.txt
