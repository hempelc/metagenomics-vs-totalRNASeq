#!/bin/bash

#SBATCH --account=def-dsteinke
#SBATCH --ntasks-per-node=32
#SBATCH --mem=0
#SBATCH --array=1-140
#SBATCH --time=4:00:00

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
memory="$((${SLURM_MEM_PER_NODE} / 1024))G" # $SLURM_MEM_PER_NODE is in Megabyte
threads=${SLURM_NTASKS_PER_NODE}
BASE="/home/hempelc/projects/def-dsteinke/hempelc/pilot_project"
start=$(date +%s)

# Echo array ID
echo -e "Job array ID is ${SLURM_ARRAY_TASK_ID}"

# Copy all necessary DBs and reads to temporary dir on server (SLURM_TMPDIR)
echo "[$(date +%H:%M:%S)] Copying started [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m]"
cp -r ${BASE}/databases/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020 \
${BASE}/databases/NCBI_staxids*_scientific.txt ${BASE}/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020 \
${BASE}/databases/SILVA_paths_and_taxids.txt ${BASE}/databases/sortmerna_silva_databases \
${BASE}/programs/ete3_env ${HOME}/.etetoolkit \
${R1} ${R2} ${BASE}/split_files_low/file_chunk_${SLURM_ARRAY_TASK_ID} \
${BASE}/programs/rRNAFilter ${SLURM_TMPDIR}
echo "[$(date +%H:%M:%S)] Copying finished [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m]"

# Set some directory-specific variables
R1=${SLURM_TMPDIR}/$(basename ${R1})
R2=${SLURM_TMPDIR}/$(basename ${R2})
DBS=${SLURM_TMPDIR}/databases

# Activate copies environment
source ${SLURM_TMPDIR}/ete3_env/bin/activate

# Assign each job in array to bundle of pipelines
jobfile=${SLURM_TMPDIR}/file_chunk_${SLURM_ARRAY_TASK_ID}

# Save path for directory in which pipeline was started
cwd1=${PWD}

# Run pipeline for each line in chunk file, i.e., each bundled pipeline
echo "[$(date +%H:%M:%S)] Bundled pipelines started [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m]"
echo "Pipelines that are to be run are:"
cat ${jobfile}

run_it(){
  jobfile=$1
  line=$2
  R1=$3
  R2=$4
  memory=$5
  threads=$6
  pipeline=$(sed -n ${line}p ${jobfile})
  mkdir -p ${pipeline}
  cwd2=${PWD}/${pipeline}
  cd ${SLURM_TMPDIR}

  echo -e "\n\n ================ [$(date +%H:%M:%S)] START PIPELINE ${pipeline} [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ==============\n\n"

  METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_compute_canada.sh \
  -1 ${R1} -2 ${R2} -P ${pipeline} \
  -N ${SLURM_TMPDIR}/nt_database_feb_2020_indexed/nt \
  -S ${SLURM_TMPDIR}/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta \
  -s ${SLURM_TMPDIR}/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020 \
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
  -i ${SLURM_TMPDIR}/rRNAFilter
  -m ${memory} \
  -p ${threads}

  mv METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/* ${cwd2}
  rm -r METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/

  echo -e "\n\n ================ [$(date +%H:%M:%S)] END PIPELINE ${pipeline} [$((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m] ==============\n\n"
}

export -f run_it
in_threads=$(( ${threads} / 8 )

parallel -j 8 run_it ${jobfile} {} ${R1} ${R2} ${memory} ${in_threads} :::: ${jobfile}

echo "[$(date +%H:%M:%S)] Bundled pipelines ended in $((($(date +%s)-${start})/3600))h $(((($(date +%s)-${start})%3600)/60))m"