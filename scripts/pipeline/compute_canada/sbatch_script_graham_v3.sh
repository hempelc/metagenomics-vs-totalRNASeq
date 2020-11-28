#!/bin/bash

#SBATCH --account=def-dsteinke
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=120G
#SBATCH --array=1-64
#SBATCH --time=10:00:00

copy_structure(){
# USAGE: copy_structure parent2retain line dest
prnt=${1}
line=${2}
dest=${3}
file=$(basename ${line})
newline=$(echo "${line}" | sed -e "s/.*${prnt}\(.*\)${file}/\1/" -e "s|^/||g" -e "s|/$||g")
mkdir -p ${dest}/${prnt}
IFS=/ read -r -a array <<< "${newline}"
for subfldr in "${array[@]}"
do
  mkdir -p ${dest}/${prnt}/${subfldr}
done
rsync -Larv ${line} ${dest}/${prnt}/${subfldr}
}

source ${HOME}/.bashrc

export -f copy_structure


memory="$((${SLURM_MEM_PER_NODE} / 1024))G" # $SLURM_MEM_PER_NODE is in Megabyte,
threads=${SLURM_NTASKS_PER_NODE}
BASE="/scratch/hempelc/chris_pilot_project"

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 trimmomatic/0.39 fastqc/0.11.9 \
spades/3.14.1 bowtie/1.3.0 trinity/2.11.0 jellyfish/2.3.0 salmon/1.3.0 \
bowtie2/2.4.1 bwa/0.7.17 boost/1.72.0 kraken2/2.1.1 blast+/2.11.0 \
seqtk/1.3 samtools/1.10 sortmerna/4.2.0 qt/5.12.8 scipy-stack/2020b \
leveldb/1.22 trans-abyss/2.0.1 megahit/1.2.9 bedtools/2.29.2

#/usr/bin/time -v -o copy_databases.log rsync -Lar ${BASE}/databases ${BASE}/programs/pipeline_environment ${HOME}/.etetoolkit/taxa.sqlite ${SLURM_TMPDIR}
#/usr/bin/time -v -o copy_input.log rsync -Lar ${R1} ${R2} ${SLURM_TMPDIR}

for prnt in databases programs/pipeline_environment
do
   parallel -j 32 --will-cite copy_structure $(basename ${prnt}) {} ${SLURM_TMPDIR} ::: $(find ${BASE}/${prnt} -type f)
done

rsync -Lar ${HOME}/.etetoolkit/taxa.sqlite ${R1} ${R2} ${BASE}/split_files/file_chunk_${SLURM_ARRAY_TASK_ID} ${SLURM_TMPDIR}

R1=${SLURM_TMPDIR}/$(basename ${R1})
R2=${SLURM_TMPDIR}/$(basename ${R2})


# echo -e "\n----------\nContents of SLURM_TMPDIR" >> copy_databases.log
#ls ${SLURM_TMPDIR} >> copy_databases.log


# Assign each job in array to nubdle of pipelines
jobfile=${SLURM_TMPDIR}/file_chunk_${SLURM_ARRAY_TASK_ID}

# Run pipeline for each line in chunk file, i.e., each bundled pipeline
while read pipeline; do
  mkdir -p ${pipeline}
  cd ${pipeline}
  cwd=${PWD}

  cd ${SLURM_TMPDIR}
  DBS=${SLURM_TMPDIR}/databases
  source ${SLURM_TMPDIR}/pipeline_environment/bin/activate

  /usr/bin/time -v -o time_pipeline_external.log METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_compute_canada.sh \
  -1 $R1 -2 $R2 -P $pipeline \
  -N ${DBS}/nt_database_feb_2020_indexed/nt \
  -S ${DBS}/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta \
  -s ${DBS}/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020 \
  -n ${DBS}/kraken2_nt_DB \
  -a ${DBS}/sortmerna_silva_databases/silva-arc-16s-id95.fasta \
  -b ${DBS}/sortmerna_silva_databases/silva-bac-16s-id90.fasta \
  -e ${DBS}/sortmerna_silva_databases/silva-euk-18s-id95.fasta \
  -E ${DBS}/sortmerna_silva_databases/silva-euk-28s-id98.fasta \
  -A ${DBS}/sortmerna_silva_databases/silva-arc-23s-id98.fasta \
  -B ${DBS}/sortmerna_silva_databases/silva-bac-23s-id98.fasta \
  -R ${DBS}/sortmerna_silva_databases/rfam-5.8s-database-id98.fasta \
  -r ${DBS}/sortmerna_silva_databases/rfam-5s-database-id98.fasta \
  -x ${DBS}/SILVA_paths_and_taxids.txt \
  -F ${DBS}/NCBI_staxids_scientific.txt \
  -f ${DBS}/NCBI_staxids_non_scientific.txt \
  -t ${SLURM_TMPDIR}/taxa.sqlite \
  -T ${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar \
  -m ${memory} \
  -p ${threads}

  /usr/bin/time -v -o copy_cp_results.log cp METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/* ${cwd}
  echo -e "\n----------\nContents of ${cwd}" >> copy_cp_results.log
  ls ${SLURM_TMPDIR} >> copy_cp_results.log
  rsync -a copy_cp_results.log time_pipeline_external.log ${cwd}
done < ${jobfile}

