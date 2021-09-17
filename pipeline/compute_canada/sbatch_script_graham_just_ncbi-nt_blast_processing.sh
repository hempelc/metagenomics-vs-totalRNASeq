#!/bin/bash

#SBATCH --account=def-dsteinke
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --time=4:00:00
#SBATCH --array=1-48

# A script to run Chris Hempel's METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE in
# parallel on graham

# Uses full graham nodes with 32 cores and all available memory (903 nodes available)

# --array=1-48 is set to 48 to generate 48 jobs for 48 bundled pipeline combinations


# Load modules
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 trimmomatic/0.39 fastqc/0.11.9 \
spades/3.14.1 bowtie/1.3.0 trinity/2.11.0 jellyfish/2.3.0 salmon/1.3.0 \
bowtie2/2.4.1 bwa/0.7.17 boost/1.72.0 kraken2/2.1.1 blast+/2.11.0 \
seqtk/1.3 samtools/1.10 sortmerna/4.2.0 qt/5.12.8 scipy-stack/2020b \
leveldb/1.22 trans-abyss/2.0.1 megahit/1.2.9 bedtools/2.29.2

# Set some general variables
BASE="/home/hempelc/projects/def-dsteinke/hempelc/pilot_project"
start=$(date +%s)

# Echo array ID
echo -e "Job array ID is ${SLURM_ARRAY_TASK_ID}"

# Activate copied environment
source ${BASE}/programs/ete3_env/bin/activate

# Assign each job in array to bundle of pipelines
pipeline=$(cat file_chunk${SLURM_ARRAY_TASK_ID}.txt)
etetoolkit=${HOME}/.etetoolkit/taxa.sqlite

# Save path for directory in which pipeline was started
cwd1=${PWD}

mkdir -p ${pipeline}
cd ${pipeline}

echo -e "\n\n ================ [$(date +%H:%M:%S)] START PIPELINE ${pipeline} [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m] ==============\n\n"


# Set pipeline tools to use
trimming=$(echo $pipeline | cut -f1 -d-)
sorting=$(echo $pipeline | cut -f2 -d-)
assembly=$(echo $pipeline | cut -f3 -d-)
mapping=$(echo $pipeline | cut -f4 -d-)
db=$(echo $pipeline | cut -f5 -d-)
classification=$(echo $pipeline | cut -f6 -d-)


step_description_and_time_first () {
	echo -e "\n++++++++ [$(date +%H:%M:%S)] ${1} [Runtime: $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"
}

step_description_and_time_second () {
	echo -e "\n======== [$(date +%H:%M:%S)] ${1} [Runtime: $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n" #" adding outcommented quote here to fix bug in colouring scheme of personal text editor
}

##################### Write start time and options to output ######################

# Define starting time of script for total runtime calculation:
echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
echo -e "=================================================================\n\n"

######################### Start of the actual script ################################
step_description_and_time_first "START RUNNING SCRIPT"

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

# Copy scaffolds and mapped files into base dir
cp ${scaffolds} merge_input_mapped_${mapping}.txt blast_output.txt METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/
cd METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/

# Save full current path in variable to make navigation between directories easier:
base_directory=$(pwd)

step_description_and_time_first "START STEP 5 AND 6.1: CLASSIFICATION OF ASSEMBLED SCAFFOLDS WITH $(echo $db | tr '[:lower:]' '[:upper:]') DATABASE"

mkdir step_5_reference_DB/
mkdir step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/
mkdir step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/step_6_classification/
mkdir step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/step_6_classification/$(echo $classification | tr '[:lower:]' '[:upper:]')/
cd step_5_reference_DB/$(echo $db | tr '[:lower:]' '[:upper:]')/step_6_classification/$(echo $classification | tr '[:lower:]' '[:upper:]')/
cp ${base_directory}/blast_output.txt $(pwd)


if [[ $classification == "blast_first_hit" ]]; then
	step_description_and_time_first "RUNNING BLAST FIRST HIT"
	# We run a separate script to add taxonomy and filter the BLAST results:
	assign_taxonomy_to_NCBI_staxids.sh -b blast_output.txt -c 13 \
	-e $etetoolkit
	sed -i 's/Unknown/NA/g' blast_output_with_taxonomy.txt
	blast_filter.py blast_output_with_taxonomy.txt soft
	if [ ! -f "blast_filtered.txt" ]; then
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies" \
		> blast_filtered.txt
	fi
	step_description_and_time_first "BLAST FIRST HIT DONE"

elif [[ $classification == "blast_filtered" ]]; then
	step_description_and_time_first "RUNNING BLAST FILTERED"
	# We run a separate script to add taxonomy and filter the BLAST results:
	assign_taxonomy_to_NCBI_staxids.sh -b blast_output.txt -c 13 \
	-e $etetoolkit
	sed -i 's/Unknown/NA/g' blast_output_with_taxonomy.txt
	blast_filter.py blast_output_with_taxonomy.txt strict
	if [ ! -f "blast_filtered.txt" ]; then
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies" \
		> blast_filtered.txt
	fi
fi
	step_description_and_time_first "BLAST FILTERED DONE"

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
	fasta_to_tab ${base_directory}/${scaffolds} > tmp
	# Add header name so that we can later merge on outer with a python script:
	echo -e "sequence_name\tsequence" > ${assembly}_tab_to_merge.txt \
	&& cat tmp >> ${assembly}_tab_to_merge.txt
	# And remove intermediate file:
	rm tmp
	# Now the scaffold names can be merged with the mapper output:
	merge_on_outer.py ${base_directory}/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt

elif [[ $assembly == 'megahit' ]]; then
	# Change the sequences' fasta format to tab delimited:
	fasta_to_tab ${base_directory}/${scaffolds} \
	| sed 's/ /\t/g' | cut -f1,4,5 | sed 's/len=//g' \
	> tmp
	# Add header name so that we can later merge on outer with a python script:
	echo -e "sequence_name\tsequence_length\tsequence" > ${assembly}_tab_to_merge.txt \
	&& cat tmp >> ${assembly}_tab_to_merge.txt
	# And remove intermediate file:
	rm tmp
	# Now the scaffold names can be merged with the mapper output:
	merge_on_outer.py ${base_directory}/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt

elif [[ $assembly == 'idba_tran' ]]; then
	# Change the sequences' fasta format to tab delimited:
	fasta_to_tab ${base_directory}/${scaffolds} \
	| sed 's/_/\t/2'  | sed 's/_/\t/3' | sed 's/ /\t/g' | cut -f1,3,5,6 \
	> tmp
	# Add header name so that we can later merge on outer with a python script:
	echo -e "sequence_name\tsequence_length\tcontig_kmer_count\tsequence" > ${assembly}_tab_to_merge.txt \
	&& cat tmp >> ${assembly}_tab_to_merge.txt
	# And remove intermediate file:
	rm tmp
	# Now the scaffold names can be merged with the mapper output:
	merge_on_outer.py ${base_directory}/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt

elif [[ $assembly == 'trinity' ]]; then
	# Change the sequences' fasta format to tab delimited:
	fasta_to_tab ${base_directory}/${scaffolds} \
	| sed 's/ /\t/1' | cut -f1,3 > tmp
	# Add header name so that we can later merge on outer with a python script:
	echo -e "sequence_name\tsequence" > ${assembly}_tab_to_merge.txt \
	&& cat tmp >> ${assembly}_tab_to_merge.txt
	# And remove intermediate file:
	rm tmp
	# Now the scaffold names can be merged with the mapper output:
	merge_on_outer.py ${base_directory}/merge_input_mapped_${mapping}.txt \
	${assembly}_tab_to_merge.txt ${assembly}_final_${mapping}_merge_ready.txt
fi

# Add classification data to the mapper output and do
# final individual edits:

if [[ $classification == 'blast_filtered' ]]; then

	# We run a simple python script because we need to merge on "outer":
	merge_on_outer.py ../../blast_filtered.txt \
	${assembly}_final_${mapping}_merge_ready.txt \
	trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	if [[ $assembly == 'spades' || $assembly == 'metaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/scaffold_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-17 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tsequence_length\tassembly_sequence" \
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
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		|	sed 's/contig-[0-9]*_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
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
		echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'transabyss' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt
	fi

	step_description_and_time_first "DONE FINALIZING BLAST_FILTERED FILES"

elif [[ $classification == 'blast_first_hit' ]]; then

	# We run a simple python script because we need to merge on "outer":
	merge_on_outer.py ../../blast_filtered.txt \
	${assembly}_final_${mapping}_merge_ready.txt \
	trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt

	if [[ $assembly == 'spades' || $assembly == 'metaspades' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' \
		| sed 's/length_//g' | sed 's/cov_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \

	elif [[ $assembly == 'idba_ud' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/scaffold_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt

	elif [[ $assembly == 'megahit' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | cut -f2-17 > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tsequence_length\tassembly_sequence" \
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
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'idba_tran' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		|	sed 's/contig-[0-9]*_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tsequence_length\tcontig_kmer_count\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $15, $16, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $17}' tmp \
		> ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'trinity' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp

	elif [[ $assembly == 'transabyss' ]]; then
		sed '1d' trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_merged.txt \
		| sed 's/_/\t/1' | sed 's/_/\t/1' | sed -r 's/_[0-9,.+-]*\t/\t/g' > trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt
		echo -e "sequence_name\tsequence_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tspecies\tcounts\tassembly_sequence" \
		> tmp \
		&& cat trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_no_header.txt \
		>> tmp
		cat tmp > ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/trimmed_at_phred_${trimming}_${sorting}_${assembly}_${mapping}_${db}_${classification}_final.txt \
		&& rm tmp
	fi

	step_description_and_time_first "DONE FINALIZING BLAST_FIRST_HIT FILES"

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"




cd ${cwd1}
cp ${pipeline}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/* ${pipeline}
rm -r ${pipeline}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/

echo -e "\n\n ================ [$(date +%H:%M:%S)] END PIPELINE ${pipeline} [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m] ==============\n\n"
