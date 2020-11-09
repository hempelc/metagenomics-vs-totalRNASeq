#!/bin/bash

#SBATCH --account=def-dsteinke
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=120G
#SBATCH --array=1-512
#SBATCH --time=10:00:00

# A script to run Chris Hempel's METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE in
# parallel on graham

# Uses full graham nodes with 32 cores and 125G memory (- 5G buffer) (903 nodes available)

# --array=1-512 is set to 512 to generate 512 jobs for 512 pipeline combinations

# Set slurm options:
memory="$((${SLURM_MEM_PER_NODE} / 1024))G" # $SLURM_MEM_PER_NODE is in Megabyte,
# so we transform it in Gigabyte and add "G" as required by Trinity
threads=$SLURM_NTASKS_PER_NODE

# Load required modules:
module load StdEnv/2018.3 gcc/7.3.0 nixpkgs/16.09 openmpi/3.1.4 \
trimmomatic/0.39 fastqc/0.11.9 spades/3.13.1 bowtie/1.1.2 trinity/2.9.0 \
bowtie2/2.3.4.3 bwa/0.7.17 boost/1.68.0 kraken2/2.0.8-beta blast+/2.10.1 \
seqtk/1.3 samtools/1.10 sortmerna/2.1 qt/5.12.3 scipy-stack/2018b leveldb/1.20 \
cuda/10.0.130 megahit/1.2.7 trans-abyss/2.0.1 bedtools/2.29.2

# This is the code I used to generate a virtual environment for ete3 and justblast:
#virtualenv --no-download ~/scratch/chris_pilot_project/programs/pipeline_environment
#source ~/scratch/chris_pilot_project/programs/pipeline_environment/bin/activate
#pip install --no-index ete3 --upgrade
#pip install git+https://github.com/jshleap/justblast.git

# activate the environment for justblast and ete3
source ~/scratch/chris_pilot_project/programs/pipeline_environment/bin/activate

# Set pipeline based on SLURM_ARRAY_TASK_ID (1-512 for each job created in the job array):
pipeline=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $file)

# Each job in the array will create the same output directory by default. To not
# let the jobs overwrite each other, each job will create a directory named after
# the pipeline and change into that directory. Once all jobs are done, we can
# copy the final files out of these:
mkdir ${pipeline}
cd ${pipeline}

# Run pipeline:
~/chrisnatjulia/scripts/pipeline/graham/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_graham.sh \
-1 $R1 -2 $R2 -P $pipeline \
-N /home/hempelc/scratch/chris_pilot_project/databases/nt_database_feb_2020_indexed/nt \
-S /home/hempelc/scratch/chris_pilot_project/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta \
-s /home/hempelc/scratch/chris_pilot_project/databases/databases/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020 \
-n /home/hempelc/scratch/chris_pilot_project/databases/kraken2_nt_DB \
-a /home/hempelc/scratch/chris_pilot_project/databases/sortmerna_silva_databases/silva-arc-16s-id95.fasta \
-b /home/hempelc/scratch/chris_pilot_project/databases/sortmerna_silva_databases/silva-bac-16s-id90.fasta \
-e /home/hempelc/scratch/chris_pilot_project/databases/sortmerna_silva_databases/silva-euk-18s-id95.fasta \
-E /home/hempelc/scratch/chris_pilot_project/databases/sortmerna_silva_databases/silva-euk-28s-id98.fasta \
-A /home/hempelc/scratch/chris_pilot_project/databases/sortmerna_silva_databases/silva-arc-23s-id98.fasta \
-B /home/hempelc/scratch/chris_pilot_project/databases/sortmerna_silva_databases/silva-bac-23s-id98.fasta \
-R /home/hempelc/scratch/chris_pilot_project/databases/sortmerna_silva_databases/rfam-5.8s-database-id98.fasta \
-r /home/hempelc/scratch/chris_pilot_project/databases/sortmerna_silva_databases/rfam-5s-database-id98.fasta \
-t ~/.etetoolkit/taxa.sqlite -T $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar -m $memory -p $threads
