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

# Record env
parallel --record-env

# Load modules
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 trimmomatic/0.39 fastqc/0.11.9 \
spades/3.14.1 bowtie/1.3.0 trinity/2.11.0 jellyfish/2.3.0 salmon/1.3.0 \
bowtie2/2.4.1 bwa/0.7.17 boost/1.72.0 kraken2/2.1.1 blast+/2.11.0 \
seqtk/1.3 samtools/1.10 sortmerna/4.2.0 qt/5.12.8 scipy-stack/2020b \
leveldb/1.22 trans-abyss/2.0.1 megahit/1.2.9 bedtools/2.29.2

# Set some general variables
memory="$((${SLURM_MEM_PER_NODE} / 1024))G" # $SLURM_MEM_PER_NODE is in Megabyte
BASE="/home/hempelc/projects/def-dsteinke/hempelc/pilot_project"
start=$(date +%s)
pip_num="8"

# Echo array ID
echo -e "Job array ID is ${SLURM_ARRAY_TASK_ID}"


# Copy all necessary DBs, reads, and files to temporary dir on server (SLURM_TMPDIR)
echo "[$(date +%H:%M:%S)] Copying started [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m]"
cp -r ${BASE}/databases/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Jul_2021 \
${BASE}/databases/NCBI_staxids*_scientific.txt ${BASE}/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Jul_2021 \
${BASE}/databases/SILVA_paths_and_taxids.txt ${BASE}/databases/sortmerna_silva_databases \
${BASE}/programs/ete3_env ${HOME}/.etetoolkit \
${R1} ${R2} ${BASE}/split_files_${pip_num}_pips/file_chunk${SLURM_ARRAY_TASK_ID}.txt \
${BASE}/programs/rRNAFilter ${SLURM_TMPDIR}
echo "[$(date +%H:%M:%S)] Copying finished [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m]"

# Set some directory-specific variables
R1=${SLURM_TMPDIR}/$(basename ${R1})
R2=${SLURM_TMPDIR}/$(basename ${R2})
DBS=${SLURM_TMPDIR}/databases

# Activate copied environment
source ${SLURM_TMPDIR}/ete3_env/bin/activate

# Assign each job in array to bundle of pipelines
jobfile=${SLURM_TMPDIR}/file_chunk${SLURM_ARRAY_TASK_ID}.txt

# Run pipeline for each line in chunk file, i.e., each bundled pipeline
echo "[$(date +%H:%M:%S)] Bundled pipelines started [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m]"
echo "Pipelines that are to be run are:"
cat ${jobfile}

# Define function to run pipelines in parallel
run_it(){
  pipeline=$1
  R1=$2
  R2=$3
  memory=$4
  threads=$5
  start=$6

  # Save path for directory in which pipeline was started
  cwd1=${PWD}

  mkdir -p ${pipeline}
  cd ${SLURM_TMPDIR}
  mkdir -p ${pipeline}
  cd ${pipeline}

  echo -e "\n\n ================ [$(date +%H:%M:%S)] START PIPELINE ${pipeline} [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m] ==============\n\n"

  METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_compute_canada_no_ncbi-nt.sh \
  -1 ${R1} -2 ${R2} -P ${pipeline} \
  -N ${SLURM_TMPDIR}/nt_database_feb_2020_indexed/nt \
  -S ${SLURM_TMPDIR}/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Jul_2021/blastdb \
  -s ${SLURM_TMPDIR}/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Jul_2021 \
  -n ${SLURM_TMPDIR}/kraken2_nt_DB \
  -a ${SLURM_TMPDIR}/sortmerna_silva_databases/silva-arc-16s-id95.fasta \
  -b ${SLURM_TMPDIR}/sortmerna_silva_databases/silva-bac-16s-id90.fasta \
  -e ${SLURM_TMPDIR}/sortmerna_silva_databases/silva-euk-18s-id95.fasta \
  -E ${SLURM_TMPDIR}/sortmerna_silva_databases/silva-euk-28s-id98.fasta \
  -A ${SLURM_TMPDIR}/sortmerna_silva_databases/silva-arc-23s-id98.fasta \
  -B ${SLURM_TMPDIR}/sortmerna_silva_databases/silva-bac-23s-id98.fasta \
  -R ${SLURM_TMPDIR}/sortmerna_silva_databases/rfam-5.8s-database-id98.fasta \
  -r ${SLURM_TMPDIR}/sortmerna_silva_databases/rfam-5s-database-id98.fasta \
  -x ${SLURM_TMPDIR}/SILVA_paths_and_taxids.txt \
  -F ${SLURM_TMPDIR}/NCBI_staxids_scientific.txt \
  -f ${SLURM_TMPDIR}/NCBI_staxids_non_scientific.txt \
  -t ${SLURM_TMPDIR}/.etetoolkit/taxa.sqlite \
  -T ${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar \
  -i ${SLURM_TMPDIR}/rRNAFilter \
  -m ${memory} \
  -p ${threads}

  cp ${SLURM_TMPDIR}/${pipeline}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/* ${cwd1}/${pipeline}
  rm -r ${SLURM_TMPDIR}/${pipeline}/

  echo -e "\n\n ================ [$(date +%H:%M:%S)] END PIPELINE ${pipeline} [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m] ==============\n\n"
}

# Set threads specifically before exporting the function
threads=${SLURM_CPUS_PER_TASK}

# Exporting the function
export -f run_it

# And dividing the numbers of threads by 8 for 8 parallel processes
in_threads=$(( ${threads} / ${pip_num} ))

parallel --env _ -j ${pip_num} run_it {} ${R1} ${R2} ${memory} ${in_threads} ${start} :::: ${jobfile}

echo "[$(date +%H:%M:%S)] Bundled pipelines ended in $((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m"
