#!/bin/bash

# Version 1.0, made on 31 Mar 2020 by Chris Hempel (hempelc@uoguelph.ca)

# Takes R1 and R2 reads, performs FastQC on the raw reads, then trims them using Trimmomatic
# with a list of PHRED scores and performs FastQC on all trimmed read sets

# FastQC must be installed and in path

# A folder with adapters to trim must be located in the same folder as the trimmomatic .jar application
# and called "adapters" (that's usually the case when you install trimmomatic)

# For some reason, the limit of the PHRED score can only be set to 38, FastQC is not able to deal with data
# generated with higher PHRED scores (in my test dataset)

usage="$(basename "$0") -T <path/to/trimmomatic.jar> -1 <R1.fastq> -2 <R2.fastq> [-P <'score score score ...'>] [-l <length>] [-t <threads>]

Usage:
	-T  Path to trimmomatic java application (.jar)
	-1  Reads1 (can also be .gz compressed)
	-2  Reads2 (can also be .gz compressed)
	-P  PHRED scores to trim on (default: '5 10 15 20' (numbers must be surrounded by ' ' and separated only by space); based on recommendations of MacMarnes 2014 and a literature review of papers that cite them)
	-l  Minimum read length to keep (default: 25 , based on recommendations of MacMarnes 2014 and a literature review of papers that cite them)
	-t  Number of threads used	(default: 16)
	-h  Display this help and exit"

# Set default options
PHRED='5 10 15 20'
min_length='25'
threads='16'


# Set specified options
while getopts ':T:1:2:P:l:t:h' opt; do
  case "${opt}" in
  	T) trimmomatic="${OPTARG}" ;;
    1) R1="${OPTARG}" ;;
    2) R2="${OPTARG}" ;;
	P) PHRED="${OPTARG}" ;;
	l) min_length="${OPTARG}" ;;
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


# Check if required options are set
if [[ -z "$trimmomatic" || -z "$R1" || -z "$R2" ]]
then
   echo -e "-T, -1, and -2 must be set.\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script.\n"
   exit
fi


######################## Write time, options etc. to output #########################

# Making output directory for this script and open bracket to later tell script to write logfile
mkdir trimming_with_phred_scores_and_fastqc_report_output/
(


# Define starting time of script for total runtime calculation
start=$(date +%s)
echo -e "\nSTART RUNNING SCRIPT AT $(date)\n"


# Output specified options
echo -e "~~~~~~~~~~ OPTIONS ~~~~~~~~~~\n"

echo -e "Trimmomatic is located in $trimmomatic."
echo -e "R1 was defined as $R1."
echo -e "R2 was defined as $R2."
echo -e "PHRED scores to trim on are $PHRED."
echo -e "Minimum length of reads is $min_length."
echo -e "Number of threads was set to $threads."

# Running Trimmomatic
echo -e "\n\n~~~~~~~~~~ RUNNING TRIMMOMATIC ~~~~~~~~~~\n"
mkdir trimming_with_phred_scores_and_fastqc_report_output/trimmomatic/
baseout=$(echo ${R1%%_*})
for score in $PHRED; do
	mkdir trimming_with_phred_scores_and_fastqc_report_output/trimmomatic/$(echo ${baseout##*/})_trimmed_at_phred_$(echo $score)
	java -jar $trimmomatic PE $R1 $R2 ILLUMINACLIP:$(echo ${trimmomatic%/*})/adapters/TruSeq3-PE.fa:2:30:10 LEADING:$score TRAILING:$score SLIDINGWINDOW:4:$score MINLEN:$min_length -baseout $(echo ${baseout##*/})_trimmed_at_phred_$(echo $score).fastq -threads $threads
	echo -e '\n'
	mv $(echo ${baseout##*/})_trimmed_at_phred_$(echo $score)*.fastq trimming_with_phred_scores_and_fastqc_report_output/trimmomatic/$(echo ${baseout##*/})_trimmed_at_phred_$(echo $score)
done

# Running FastQC
echo -e "~~~~~~~~~~ RUNNING FASTQC ~~~~~~~~~~\n"
mkdir trimming_with_phred_scores_and_fastqc_report_output/fastqc_reports
mkdir trimming_with_phred_scores_and_fastqc_report_output/fastqc_reports/${baseout##*/}_untrimmed
fastqc $R1 -o trimming_with_phred_scores_and_fastqc_report_output/fastqc_reports/${baseout##*/}_untrimmed
echo -e "\n"
fastqc $R2 -o trimming_with_phred_scores_and_fastqc_report_output/fastqc_reports/${baseout##*/}_untrimmed
echo -e "\n"

for reads in trimming_with_phred_scores_and_fastqc_report_output/trimmomatic/$(echo ${baseout##*/})_trimmed_at_phred_*/*_1P.fastq; do
	dir_name=${reads%_*} # making variable so that read names can be used to make directory
	mkdir trimming_with_phred_scores_and_fastqc_report_output/fastqc_reports/$(echo ${dir_name##*/})
	fastqc $reads -o trimming_with_phred_scores_and_fastqc_report_output/fastqc_reports/$(echo ${dir_name##*/}) -t $threads
	echo -e "\n"
done

for reads in trimming_with_phred_scores_and_fastqc_report_output/trimmomatic/$(echo ${baseout##*/})_trimmed_at_phred_*/*_2P.fastq; do
	dir_name=${reads%_*} # making variable so that read names can be used to sort into directory
	fastqc $reads -o trimming_with_phred_scores_and_fastqc_report_output/fastqc_reports/$(echo ${dir_name##*/}) -t $threads
	echo -e "\n"
done


# Display runtime
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nSCRIPT RUNTIME: $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"


# Write output to console and log file
) 2>&1 | tee trimming_with_phred_scores_and_fastqc_report_output/trimming_with_phred_scores_and_fastqc_report_log.txt
