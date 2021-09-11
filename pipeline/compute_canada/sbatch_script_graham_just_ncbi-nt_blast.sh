#!/bin/bash

#SBATCH --account=def-dsteinke
#SBATCH --cpus-per-task=2
#SBATCH --mem=124G
#SBATCH --time=4:00:00
#SBATCH --array=1-56

# A script to run Chris Hempel's METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE in
# parallel on graham

# Uses full graham nodes with 32 cores and all available memory (903 nodes available)

# --array=1-64 is set to 64 to generate 64 jobs for 64 bundled pipeline combinations

# Load modules
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 trimmomatic/0.39 fastqc/0.11.9 \
spades/3.14.1 bowtie/1.3.0 trinity/2.11.0 jellyfish/2.3.0 salmon/1.3.0 \
bowtie2/2.4.1 bwa/0.7.17 boost/1.72.0 kraken2/2.1.1 blast+/2.11.0 \
seqtk/1.3 samtools/1.10 sortmerna/4.2.0 qt/5.12.8 scipy-stack/2020b \
leveldb/1.22 trans-abyss/2.0.1 megahit/1.2.9 bedtools/2.29.2

# Set some general variables
BASE="/home/hempelc/projects/def-dsteinke/hempelc/pilot_project"
start=$(date +%s)
DBS=${BASE}/databases/nt_database_14_Feb_2021_indexed/nt
jobfile=${BASE}/split_files_8_pips_just_ncbi_nt/file_chunk${SLURM_ARRAY_TASK_ID}.txt
threads=${SLURM_CPUS_PER_TASK}

# Run pipeline for each line in chunk file, i.e., each bundled pipeline
echo "[$(date +%H:%M:%S)] Bundled pipelines started [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m]"
echo "Pipelines that are to be run are:"
cat ${jobfile}

cwd=${PWD}
mkdir -p ${pipeline}
cd ${pipeline}

echo -e "\n\n ================ [$(date +%H:%M:%S)] START PIPELINE ${pipeline} [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m] ==============\n\n"

blastn -query *.fa* -db $DBS -out blast_output.txt \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
-evalue 1e-05 -num_threads $threads

cd $cwd
echo -e "\n\n ================ [$(date +%H:%M:%S)] END PIPELINE ${pipeline} [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m] ==============\n\n"

echo "[$(date +%H:%M:%S)] Bundled pipelines ended in $((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m"
