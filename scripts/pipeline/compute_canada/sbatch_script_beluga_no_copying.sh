#!/bin/bash

#SBATCH --account=def-dsteinke
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=186G
#SBATCH --array=1-64
#SBATCH --time=10:00:00

# A script to run Chris Hempel's METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE in
# parallel on beluga

# Uses full beluga nodes with 40 cores and all available memory (516 nodes available)

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
echo "[$(date +%H:%M:%S)] Copying started [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m]"
cp -r ${BASE}/databases/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020 \
${BASE}/databases/kraken2_nt_DB \
${BASE}/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020 \
${BASE}/databases/SILVA_paths_and_taxids.txt ${BASE}/databases/sortmerna_silva_databases \
${HOME}/ete3_env ${HOME}/.etetoolkit /home/hempelc/scratch/nt/ \
${R1} ${R2} ${BASE}/split_files/file_chunk_${SLURM_ARRAY_TASK_ID} ${SLURM_TMPDIR}
echo "[$(date +%H:%M:%S)] Copying finished [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m]"

# # Copy all necessary DBs and reads to temporary dir on server (SLURM_TMPDIR)
# echo "[$(date +%H:%M:%S)] Copying started [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m]"
# cp -r ${BASE}/databases ${BASE}/programs/pipeline_environment ${HOME}/.etetoolkit \
# ${R1} ${R2} ${BASE}/split_files/file_chunk_${SLURM_ARRAY_TASK_ID} ${SLURM_TMPDIR}
# echo "[$(date +%H:%M:%S)] Copying finished [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m]"

# Set some directory-specific variables
R1=${SLURM_TMPDIR}/$(basename ${R1})
R2=${SLURM_TMPDIR}/$(basename ${R2})
#DBS=${SLURM_TMPDIR}/databases

# Activate copies environment
#source ${SLURM_TMPDIR}/pipeline_environment/bin/activate
source ${SLURM_TMPDIR}/ete3_env/bin/activate

# Assign each job in array to bundle of pipelines
jobfile=${SLURM_TMPDIR}/file_chunk_${SLURM_ARRAY_TASK_ID}

# Save path for directory in which pipeline was started
cwd1=${PWD}

# Run pipeline for each line in chunk file, i.e., each bundled pipeline
echo "[$(date +%H:%M:%S)] Bundled pipelines started [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m]"
echo "Pipelines that are to be run are:"
cat ${jobfile}

for line in {1..8}; do
  pipeline=$(sed -n ${line}p ${jobfile})
  echo -e "\n\n ================ [$(date +%H:%M:%S)] START PIPELINE ${pipeline} [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ==============\n\n"
  mkdir -p ${pipeline}
  cd ${pipeline}
  cwd2=${PWD}
  cd ${SLURM_TMPDIR}

  # manually changed ${DBS} to ${SLURM_TMPDIR}, needs to be changed once once the directory structure is adapted to  include new BLAST nt DB in database/
    METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_compute_canada.sh \
    -1 $R1 -2 $R2 -P $pipeline \
    -N ${SLURM_TMPDIR}/nt/nt \
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
    -F ${SLURM_TMPDIR}/nt/NCBI_staxids_scientific.txt \
    -f ${SLURM_TMPDIR}/nt/NCBI_staxids_non_scientific.txt \
    -t ${SLURM_TMPDIR}/.etetoolkit/taxa.sqlite \
    -T ${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar \
    -m ${memory} \
    -p ${threads}
  # -N ${SLURM_TMPDIR}/nt_database_feb_2020_indexed/nt


  # METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_compute_canada.sh \
  # -1 $R1 -2 $R2 -P $pipeline \
  # -N ${DBS}/nt_database_feb_2020_indexed/nt \
  # -S ${DBS}/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta \
  # -s ${DBS}/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020 \
  # -n ${DBS}/kraken2_nt_DB \
  # -a ${DBS}/sortmerna_silva_databases/silva-arc-16s-id95.fasta \
  # -b ${DBS}/sortmerna_silva_databases/silva-bac-16s-id90.fasta \
  # -e ${DBS}/sortmerna_silva_databases/silva-euk-18s-id95.fasta \
  # -E ${DBS}/sortmerna_silva_databases/silva-euk-28s-id98.fasta \
  # -A ${DBS}/sortmerna_silva_databases/silva-arc-23s-id98.fasta \
  # -B ${DBS}/sortmerna_silva_databases/silva-bac-23s-id98.fasta \
  # -R ${DBS}/sortmerna_silva_databases/rfam-5.8s-database-id98.fasta \
  # -r ${DBS}/sortmerna_silva_databases/rfam-5s-database-id98.fasta \
  # -x ${DBS}/SILVA_paths_and_taxids.txt \
  # -F ${DBS}/NCBI_staxids_scientific.txt \
  # -f ${DBS}/NCBI_staxids_non_scientific.txt \
  # -t ${SLURM_TMPDIR}/.etetoolkit/taxa.sqlite \
  # -T ${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar \
  # -m ${memory} \
  # -p ${threads}

  cp METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/* ${cwd2}
  rm -r METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/
  cd ${cwd1}
  echo -e "\n\n ================ [$(date +%H:%M:%S)] END PIPELINE ${pipeline} [$((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ==============\n\n"
done
echo "[$(date +%H:%M:%S)] Bundled pipelines ended in $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"
