#!/bin/bash

# Version 0.1
# Written by Natalie Wright (nwrigh06@uoguelph.ca) and Chris Hempel (hempelc@uoguelph.ca)

# This is a pipeline for Chris Hempel's first PhD chapter

# It trims raw paired-end input sequences at 4 PHRED scores, filters rRNA
# with 4 approaches, uses 8 assemblers, maps trimmed reads back to scaffolds
# using 2 mappers, and assigns taxonomy to scaffolds using 2 databases with
# 3 classification approaches.

# The output is a folder called METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/
# that contains tab-separated, taxonomically annotated scaffolds and read counts
# for all combinations of the above described steps.

# The pipeline requires the following subscripts, which are all located in the
# subscripts/ directory:
# assign_NCBI_staxids_to_CREST_v4.py, fasta_to_tab, mergeFilesOnColumn.pl,
# assign_taxonomy_to_NCBI_staxids.sh  fastqc_on_R1_R2_and_optional_trimming.sh,
# merge_mapped_reads_and_contigs.py, blast_filtering.bash, filter-fasta.awk,
# SILVA_SSU_LSU_kraken2_preparation.sh, deinterleave_fastq_reads.sh,
# LookupTaxonDetails3.py, SILVA_SSU_LSU_makeblastdb_preparation.sh

# The pipeline requires the following programs/python packages/commands (versions
# we used when writing this script are indicated in brackets):
# FastQC (0.11.5), Trimmomatic (0.33), sortmeRNA (4.0.0), barrnap (0.9),
# rRNAFILTER (1.1), SPADES (3.14.0), METASPADES (3.14.0), RNASPADES (3.14.0),
# MEGAHIT (1.2.9), IDBA-UD (1.1.1), IDBA-TRAN (1.1.1), Trinity (2.10.0),
# bowtie2 (2.3.3.1), bwa (0.7.17), blast (2.10.0+), justblast (2020.0.3),
# seqtk (1.2-r94)

# Note: we had to edit IDBA prior to compiling it because it didn't work
# using long reads and the -l option. This seems to be a common problem and
# can be circumvented following for example the instructions in
# http://seqanswers.com/forums/showthread.php?t=29109, and see also
# https://github.com/loneknightpy/idba/issues/26

cmd="$0 $@" # Make variable containing the entire entered command to print command to logfile
usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq> [-t <n>]
Usage:
	-1 Full path to forward reads in .fastq/.fq format
	-2 Full path to reverse reads in .fastq/.fq format
	-t Number of threads (default:16)
	-h Display this help and exit"

# Set default options:
threads='16'

# Set specified options:
while getopts ':1:2:t:h' opt; do
  case "${opt}" in
  1) forward_reads="${OPTARG}" ;;
  2) reverse_reads="${OPTARG}" ;;
  t) threads="${OPTARG}" ;;
  h)
    echo "$usage"
    exit
    ;;
  :)
    printf "Option -$OPTARG requires an argument."
    echo -e "\n$usage"
    exit
    ;;
  \?)
    printf "Invalid option: -$OPTARG"
    echo -e "\n$usage"
    exit
    ;;
  esac
done
shift $((OPTIND - 1))

# Check if required options are set:
if [[ ! -s "${forward_reads}" || ! -s "${reverse_reads}" ]]; then
  echo -e "-1 and -2 must be set, present and not empty\n"
  echo -e "${usage}\n\n"
  echo -e "Exiting script.\n"
  exit
fi

##################### Write start time and options to output ######################

# Make open bracket to later tell script to write everything that follows into
# a logfile:
# SERGIO: this creates more problems that it solves, just redirect the stdout
# to file when calling the script
#(

# initialize the especial variable SECONDS so convertsecs will measure
SECONDS=0
# Import utility functions
SCRIPT=$(realpath "${0}")
SCRIPTPATH=$(dirname "${SCRIPT}")
source "${SCRIPTPATH}"/utils.sh

# Define starting time of script for total runtime calculation:
echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
echo -e "=================================================================\n\n"

# Output specified options:
echo -e "======== OPTIONS ========\n"

echo -e "Forward reads were defined as ${forward_reads}.\n"
echo -e "Reverse reads were defined as ${reverse_reads}.\n"
echo -e "Number of threads was set to ${threads}.\n"
echo -e "Script started with full command: ${cmd}\n"

######################### Start of the actual script ################################
echo_section START RUNNING SCRIPT

# Make output directory:
mkdir -p METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE &&
  cd METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE || exit

# Save full current path in variable to make navigation between directories easier:
base_directory=${PWD}

######################### Step 1: trimming ################################
echo_section "START STEP 1: TRIMMING AND ERROR CORRECTION"

# Trimming is done with a separate subscript:
fastqc_on_R1_R2_and_optional_trimming.sh \
  -T /hdd1/programs_for_pilot/Trimmomatic-0.39/trimmomatic-0.39.jar \
  -1 ${forward_reads} -2 ${reverse_reads} -t yes -p $threads

# SERGIO: WHY? just pass the name of the folder to the fastqc script (first
# make it an option)
mv trimming_with_phred_scores_and_fastqc_report_output/ step_1_trimming/

# Running error correction module of SPAdes on all trimmed reads
for trimming_results in step_1_trimming/trimmomatic/*; do
  echo_subsection "ERROR-CORRECTING READS IN FOLDER ${trimming_results}"
  input1="${trimming_results}/*1P.fastq"
  input2="${trimming_results}/*2P.fastq"
  output="${trimming_results}"/error_correction
  spades.py -1 "${input1}" -2 "${input2}" --only-error-correction \
    --disable-gzip-output -o "${output}" -t "${threads}"

  # SERGIO: Again, I do not understand why are you moving stuff around. Also
  # this overwrites the 2P
  # mv "${trimming_results}"/error_correction/corrected/*1P*.fastq \
  #    "${trimming_results}"/error_correction/corrected/*2P*.fastq "${trimming_results}"
  mv ${output}/corrected/*[12]P*.fastq
  # Rename weird name of error-corrected reads:
  rename 's/.00.0_0.cor.fastq/_error_corrected.fastq/' "${trimming_results}"/*.00.0_0.cor.fastq
  for f in "${trimming_results}"/*_error_corrected.fastq; do
    sed -r -i 's/ BH:.{2,6}//g' "${f}"
  done
  rm -rf "${trimming_results}"/error_correction/
  echo_subsection FINISHED ERROR-CORRECTING READS IN FOLDER "${trimming_results}"
done

echo_section FINISHED STEP 1: TRIMMING AND ERROR CORRECTION

######################### Step 2: rRNA sorting ################################

# NOTE: To enable each combination of programs in each step, we run nested for
# loops and close them all at the very end of the script.

# For loop 1: loop over the folders for trimmed data at different PHRED scores
# for rRNA filtration:
for trimming_results in step_1_trimming/trimmomatic/*; do
  filter_dir="${trimming_results}/step_2_rrna_sorting/"
  fastas_dir="${filter_dir}/reads_in_fasta_format"
  mkdir -p "${fastas_dir}"
  #  cd "${filter_dir}" || exit

  echo_section START STEP 2: rRNA SORTING OF TRIMMED READS IN FOLDER \
    "${trimming_results}"

  echo_subsection CONVERT READS IN FASTA FORMAT FOR rRNAFILTER AND BARRNAP
  #  mkdir reads_in_fasta_format/
  fq2fa "${trimming_results}"/*1P_error_corrected.fastq "${fastas_dir}"/R1.fa
  fq2fa "${trimming_results}"/*2P_error_corrected.fastq "${fastas_dir}"/R2.fa
  echo_subsection READS TO FASTA CONVERSION DONE

  echo_subsection RUNNING SORTMERNA
  sortmerna_dir="${filter_dir}/SORTMERNA"
  mkdir -p "${sortmerna_dir}"
  # SERGIO: Hardcoding the paths will make it unusable for others and on the clusters
  # SERGIO: Are you sure that reads takes a
  sortmerna --ref /hdd1/databases/sortmerna_silva_databases/silva-bac-16s-id90.fasta \
    --ref /hdd1/databases/sortmerna_silva_databases/silva-arc-16s-id95.fasta \
    --ref /hdd1/databases/sortmerna_silva_databases/silva-euk-18s-id95.fasta \
    --reads "${trimming_results}"/*1P_error_corrected.fastq \
    --reads "${trimming_results}"/*2P_error_corrected.fastq \
    --paired_in -other -fastx 1 -num_alignments 1 -v \
    -workdir "${sortmerna_dir}" --threads 1:1:"${threads}"
  # SortMeRNA interleaves reads, which we don't want, so we deinterleave them:
  deinterleave_fastq_reads.sh <"${sortmerna_dir}"/out/aligned.fastq \
    "${sortmerna_dir}"/out/aligned_R1.fq \
    "${sortmerna_dir}"/out/aligned_R2.fq
  echo_subsection SORTMERNA DONE

  echo_subsection RUNNING rRNAFILTER
  rrnafilter_dir="${filter_dir}/rRNAFILTER"
  mkdir -p "${rrnafilter_dir}"
  #  cd rRNAFILTER/
  # rRNAFilter only worked for us when we started it within the directory
  # containing the .jar file. To simplify switching to that directory, we simply
  # download the small program within the script and delete it after usage:
  # SERGIO: I WILL AVOID DOING THIS. If is truly necessary to execute in the
  # same folder, you are better off having the zip file on path (or the folder)
  # and either link it or copy it.
  wget http://hulab.ucf.edu/research/projects/rRNAFilter/software/rRNAFilter.zip
  unzip rRNAFilter.zip
  #  cd rRNAFilter/
  # We use 7GB for the rRNAFilter .jar, as shown in the rRNAFilter manual:
  # SERGIO: The rRNAFilter_commandline does not take regex?
  java -jar -Xmx7g rRNAFilter_commandline.jar -i "${fastas_dir}/R1.fa" -r 0
  java -jar -Xmx7g rRNAFilter_commandline.jar -i "${fastas_dir}/R2.fa" -r 0
  mv "${fastas_dir}"/R*.fa_rRNA "${filter_dir}"

  rm -r rRNAFilter rRNAFilter.zip
  # We want to keep paired reads, so we extract all rRNA read names that were
  # found in R1 and R2, save them as one list, and extract all reads from both
  # R1 and R2 reads. That way, even if only one read from a pair was identified
  # as rRNA, we keep the pair of reads:
  fasta_to_tab "${filter_dir}"/R1.fa_rRNA | cut -f 1 | cut -f1 -d " " >names.txt
  fasta_to_tab "${filter_dir}"/R2.fa_rRNA | cut -f 1 | cut -f1 -d " " >>names.txt
  sort -u names.txt >names_sorted.txt
  seqtk subseq "${fastas_dir}"/R1.fa names_sorted.txt >rRNAFilter_paired_R1.fa
  seqtk subseq "${fastas_dir}"/R2.fa names_sorted.txt >rRNAFilter_paired_R2.fa
  rm names_sorted.txt names.txt
  #  cd ..
  echo_subsection rRNAFILTER DONE
  echo_subsection RUNNING BARRNAP
  barnap_dir="${filter_dir}/BARRNAP"
  mkdir -p "${barnap_dir}"
  for kingdom in euk bac arc; do # barrnap needs to be run on kingdoms separately
    echo_subsection RUNNING BARRNAP ON KINGDOM "${kingdom}" AND R1 READS
    barrnap --lencutoff 0.000001 --reject 0.000001 --kingdom "${kingdom}" \
      --threads "${threads}" --outseq "${barnap_dir}/${kingdom}_reads1.fa" \
      "${fastas_dir}/R1.fa"
    echo_subsection RUNNING BARRNAP ON KINGDOM ${kingdom} AND R2 READS
    barrnap --lencutoff 0.000001 --reject 0.000001 --kingdom "${kingdom}" \
      --threads "${threads}" --outseq "${barnap_dir}"/${kingdom}_reads2.fa \
      reads_in_fasta_format/R2.fa
    rm reads_in_fasta_format/*.fai
    sed 's/.*::/>/g' "${barnap_dir}"/${kingdom}_reads1.fa | sed 's/:[^:]*$//g' \
      >"${barnap_dir}"/${kingdom}_reads1_edited.fa
    sed 's/.*::/>/g' "${barnap_dir}"/${kingdom}_reads2.fa | sed 's/:[^:]*$//g' \
      >"${barnap_dir}"/${kingdom}_reads2_edited.fa
  done
  # Concatenating results from the three kingdoms and R1 and R2 files
  cat "${barnap_dir}"/*edited.fa >"${barnap_dir}"/all_reads.fa
  # We want to keep paired reads, so we extract all rRNA read names that were
  # found in R1 and R2 for the three kingdoms (in all_reads.fa), save them as
  # one list, and extract all reads from both R1 and R2 reads. That way, even if
  # only one read from a pair was identified as rRNA, we keep the pair of reads:
  fasta_to_tab "${barnap_dir}"/all_reads.fa | cut -f 1 | cut -f1 -d " " | sort -u \
    >"${barnap_dir}"/names_sorted.txt
  seqtk subseq reads_in_fasta_format/R1.fa "${barnap_dir}"/names_sorted.txt \
    >"${barnap_dir}"/barrnap_paired_R1.fa
  seqtk subseq reads_in_fasta_format/R2.fa "${barnap_dir}"/names_sorted.txt \
    >"${barnap_dir}"/barrnap_paired_R2.fa
  rm "${barnap_dir}"/names_sorted.txt

  echo_subsection BARRNAP DONE
  echo_subsection MAKING FOLDER UNSORTED/ AND COPYING UNSORTED READS IN THERE TO KEEP THE FOLDER STRUCTURE CONSTANT
  unsorted_dir="${barnap_dir}"/UNSORTED
  mkdir -p "${unsorted_dir}"
  # SERGIO: again better absolute path... check if I did not mess that up
  cp "${filter_dir}"/*[12]P_error_corrected.fastq "${unsorted_dir}"

  echo_section FINISHED STEP 2: rRNA SORTING OF TRIMMED READS IN FOLDER "${trimming_results}"

  ######################### Step 3: Assembly ################################

  # For loop 2: loop over the folders for rRNA filtered data for assembly:
  for rrna_filter_results in rRNAFILTER SORTMERNA BARRNAP UNSORTED; do
    if [[ $rrna_filter_results == 'rRNAFILTER' ]]; then
      R1_sorted="${rrnafilter_dir}/rRNAFilter_paired_R1.fa"
      R2_sorted="${rrnafilter_dir}/rRNAFilter_paired_R2.fa"
    elif [[ $rrna_filter_results == 'SORTMERNA' ]]; then
      R1_sorted="${sortmerna_dir}/out/aligned_R1.fq"
      R2_sorted="${sortmerna_dir}/out/aligned_R2.fq"
    elif [[ $rrna_filter_results == 'BARRNAP' ]]; then
      R1_sorted="${barnap_dir}/barrnap_paired_R1.fa"
      R1_sorted="${barnap_dir}//barrnap_paired_R2.fa"
    else # UNSORTED
      R1_sorted="${unsorted_dir}/*1P_error_corrected.fastq"
      R2_sorted="${unsorted_dir}/*2P_error_corrected.fastq"
    fi

    echo_section START STEP 3: ASSEMBLY OF TRIMMED READS IN FOLDER "${trimming_results}" \
      AND rRNA FILTERED READS IN FOLDER "${rrna_filter_results}"

    assembly_folder="${rrna_filter_results}/step_3_assembly"
    mkdir -p "${assembly_folder}"
    #    cd ${assembly_folder} || exit

    echo_subsection RUNNING SPADES
    spades_dir="${assembly_folder}/SPADES"
    mkdir -p "${spades_dir}"
    spades.py -1 "${R1_sorted}" -2 "${R2_sorted}" --only-assembler \
      -o "${spades_dir}" -t "${threads}"
    echo_subsection SPADES DONE

    echo_subsection RUNNING METASPADES
    metaspades_dir="${assembly_folder}/METASPADES"
    mkdir -p "${metaspades_dir}"
    metaspades.py -1 "${R1_sorted}" -2 "${R2_sorted}" --only-assembler \
      -o "${metaspades_dir}" -t "${threads}"
    echo_subsection METASPADES DONE

    echo_subsection RUNNING MEGAHIT
    # SERGIO: Does not need to be created?
    megahit_dir="${assembly_folder}/MEGAHIT"
    megahit --presets meta-large -t "${threads}" -1 "${R1_sorted}" \
      -2 "${R2_sorted}" -o "${megahit_dir}"
    echo_subsectionMEGAHIT DONE

    echo_subsection RUNNING IDBA_UD
    # Note: we had to edit IDBA prior to compiling it because it didn't work
    # using long reads and the -l option. This seems to be a common problem and
    # can be circumvented following for example the instructions in
    # http://seqanswers.com/forums/showthread.php?t=29109, and see also
    # https://github.com/loneknightpy/idba/issues/26
    # IDBA_UD only takes interleaved fasta files

    #SERGIO: does not need to be created?
    idba_dir="${assembly_folder}/IDBA_UD"
    idba_input="${idba_dir}/idba_ud_input.fa"
    fq2fa --merge --filter "${R1_sorted}" "${R2_sorted}" "${idba_input}"
    idba_ud --num_threads "${threads}" --pre_correction -r "${idba_input}" \
      -o "${idba_dir}"
    echo_subsection IDBA_UD DONE

    echo_subsection RUNNING RNASPADES
    rnaspades_dir="${assembly_folder}/RNASPADES"
    mkdir -p "${rnaspades_dir}"
    rnaspades.py -1 "${R1_sorted}" -2 "${R2_sorted}" --only-assembler \
      -o "${rnaspades_dir}" -t "${threads}"
    echo_subsection RNASPADES DONE

    echo_subsection RUNNING IDBA_TRAN
    idbatran_dir="${assembly_folder}/IDBA_TRAN"
    idbat_input="${idbatran_dir}/idba_tran_input.fa"
    # IDBA_TRAN only takes interleaved fasta files
    fq2fa --merge "${R1_sorted}" "${R2_sorted}" "${idbat_input}"
    idba_tran --num_threads $threads --pre_correction -l "${idbat_input}" \
      -o IDBA_TRAN/
    echo_subsection IDBA_TRAN DONE

    echo_subsection RUNNING TRINITY
    trinity_dir="${assembly_folder}/TRINITY"
    # SERGIO: to run in shared systems this cannot be set this way. Have to be
    # user provided
    max_mem=$(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024) - 5))
    # Barrnap and rRNAFilter output fasta files which has to be indicated to Trinity:
    if [[ "${rrna_filter_results}" == "rRNAFILTER" || "${rrna_filter_results}" == "BARRNAP" ]]; then
      t=fa
    else
      t=fq
    fi
    Trinity --seqType ${t} --max_memory "${max_mem}"G --left "${R1_sorted}" \
      --right "${R2_sorted}" --CPU "${threads}" --output "${trinity_dir}" # The max_memory command simply takes the maximum RAM size in GB and subtracts 5GB
    sed 's/ len/_len/g' "${trinity_dir}/Trinity.fasta" >"${trinity_dir}/Trinity_with_length.fasta" # Edit for universal format
    echo_subsection TRINITY DONE

    echo_subsection RUNNING TRANSABYSS
    tabyss_dir="${assembly_folder}/TRANSABYSS"
    transabyss --pe "${R1_sorted}" "${R2_sorted}" --threads "${threads}" \
      --outdir "${tabyss_dir}"
    sed 's/ /_/g' "${tabyss_dir}/transabyss-final.fa" >"${tabyss_dir}/transabyss-final_edited.fa" # Edit for universal format
    echo_subsection TRANSABYSS DONE

    echo_section FINISHED STEP 3: ASSEMBLY OF TRIMMED READS IN FOLDER "${trimming_results}" \
      AND rRNA FILTERED READS IN FOLDER "${rrna_filter_results}"

    ######################### Step 4: Mapping ################################

    # For loop 3: loop over the folders for assembled data for mapping:
    for assembly_results in SPADES METASPADES MEGAHIT IDBA_UD RNASPADES IDBA_TRAN TRINITY TRANSABYSS; do
      if [[ "${assembly_results}" == 'SPADES' ]]; then
        scaffolds="${spades_dir}/scaffolds.fasta"
      elif [[ "${assembly_results}" == 'METASPADES' ]]; then
        scaffolds="${metaspades_dir}/scaffolds.fasta"
      elif [[ "${assembly_results}" == 'MEGAHIT' ]]; then
        scaffolds="${megahit_dir}/final.contigs.fa"
      elif [[ "${assembly_results}" == 'IDBA_UD' ]]; then
        scaffolds="${idba_dir}/scaffold.fa"
      elif [[ "${assembly_results}" == 'RNASPADES' ]]; then
        scaffolds="${rnaspades_dir}/transcripts.fasta"
      elif [[ "${assembly_results}" == 'IDBA_TRAN' ]]; then
        scaffolds="${idbatran_dir}/contig.fa"
      elif [[ "${assembly_results}" == 'TRINITY' ]]; then
        scaffolds="${trinity_dir}/Trinity_with_length.fasta"
      else # TRANSABYSS
        scaffolds="${tabyss_dir}/transabyss-final_edited.fa"
      fi

      echo_section START STEP 4: MAPPING OF TRIMMED READS IN FOLDER \
        "${trimming_results}" AND rRNA FILTERED READS IN FOLDER \
        "${rrna_filter_results}" AND ASSEMBLY IN FOLDER "${assembly_results}"

      mapping_dir="${assembly_results}/step_4_mapping"
      mkdir -p "${mapping_dir}"
      #      cd "${mapping_dir}" || exit

      # We insert a loop here for the mappers because their output is edited
      # with identical commands:
      for mapper in BWA BOWTIE2; do
        mapper="${mapping_dir}/${mapper}"
        mkdir -p "${mapper}"
        #        cd "${mapper}" || exit

        if [[ "${mapper}" == 'BWA' ]]; then
          echo_subsection Starting bwa index
          index_pref="${mapper}/bwa_index"
          bwa index -p "${index_pref}" ${scaffolds}
          echo_subsection bwa index complete. Starting bwa mem
          bwa mem -t "${threads}" ${index_pref} "${trimming_results}"/*1P_error_corrected.fastq \
            "${trimming_results}"/*2P_error_corrected.fastq >${mapper}_output.sam
          rm "${index_pref}*"
          echo_subsection bwa mem complete
        else
          echo_subsection Starting bowtie2 index
          index_pref="${mapper}/bowtie_index"
          bowtie2-build -f "${scaffolds}" "${index_pref}"
          echo_subsection bowtie2 index complete. Starting bowtie2
          bowtie2 -q -x "${index_pref}" -1 "${trimming_results}"/*1P_error_corrected.fastq \
            -2 "${trimming_results}"/*2P_error_corrected.fastq -S ${mapper}_output.sam \
            -p "${threads}"
          rm "${index_pref}"*
          echo_subsection bowtie2 complete
        fi

        # Editing the mapper outputs:
        name=$(basename ${mapper})
        otm="${mapper}/out_mapped_${name}.txt"
        otu="${mapper}/out_unmapped_${name}.txt"
        mim="${mapper}/merge_input_mapped_${name}.txt"
        miu="${mapper}/merge_input_unmapped_${name}.txt"
        samtools view -F 4 ${mapper}_output.sam | cut -f3 | sort | uniq -c |
          column -t | sed 's/  */\t/g' >"${otm}"
        samtools view -f 4 ${mapper}_output.sam | cut -f3 | sort | uniq -c |
          column -t | sed 's/  */\t/g' >"${otu}"
        echo -e "counts\tcontig_number" >"${mim}" && cat ${otm} >>"${mim}"
        echo -e "counts\tcontig_number" >"${miu}" && cat "${otu}" >>"${miu}"

        rm "${mapper}"/*_index* "${mapper}"/out_*mapped_${name}.txt
        #        cd ..
        # And we close the mapper loop here because the next step is independent
        # from the mappers and saved under a separate folder:
      done

      # We cd back into the step_3 directory using the base_directory variable
      # we created at the beginning of the script and the nested for loop
      # variables we generated during the script:
      #      cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)
      cd "${assembly_folder}" || exit

      echo_section FINISHED STEP 4: MAPPING FOR TRIMMED READS IN FOLDER \
        "${trimming_results}"/ AND rRNA FILTERED READS IN FOLDER \
        "${rrna_filter_results}" AND ASSEMBLY IN FOLDER "${assembly_results}"

      ######################### Steps 5 and 6.1: Picking a reference DB and taxonomic classification ################################
      refdb_dir="${assembly_results}/step_5_reference_DB/"
      mkdir - p "${refdb_dir}"

      # For loop 4: loop over the reference DBs for taxonomic classification:
      # SERGIO: Hardcoded paths are not transferable
      for DB in SILVA NCBI_NT; do
        if [[ "${DB}" == "SILVA" ]]; then
          krakenDB="/hdd1/databases/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020/"
          blastDB="/hdd1/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc.fasta"
        else # NCBI_NT
          krakenDB="/hdd1/databases/kraken2_nt_DB"
          blastDB="/hdd1/databases/nt_database_feb_2020_indexed/nt"
        fi

        echo_section START STEP 5 AND 6.1: CLASSIFICATION OF ASSEMBLED SCAFFOLDS \
          FROM TRIMMED READS IN FOLDER "${trimming_results}" AND rRNA FILTERED READS \
          IN FOLDER "${rrna_filter_results}" AND ASSEMBLY IN FOLDER "${assembly_results}"
        class_dir="${refdb_dir}/${DB}/step_6_classification"
        mkdir -p ${class_dir}
        #        cd $assembly_results/step_5_reference_DB/${DB}/step_6_classification/

        echo_subsection RUNNING JUSTBLAST WITH DATABASE "${DB}"
        # Run BLAST via justblast
        # SERGIO: Do you really need all those headers?
        justblast "${scaffolds}" "${blastDB}" --cpus ${threads} --evalue 1e-05 \
          --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
          --out_filename "${class_dir}/blast_output.txt"
        echo_subsection JUSTBLAST WITH DATABASE "${DB}" DONE

        echo -e "\n======== RUNNING BLAST FIRST HIT ========\n"
        # We run a separate script to filter the BLAST results:
        blast_filtering.bash -i blast_output.txt -f blast -t soft -T $threads
        cp blast_output.txt blast_filtering_results/
        mv blast_filtering_results/ BLAST_FIRST_HIT/
        echo -e "\n======== BLAST FIRST HIT DONE ========\n"

        echo -e "\n======== RUNNING BLAST FILTERED ========\n"
        # We run a separate script to filter the BLAST results:
        blast_filtering.bash -i blast_output.txt -f blast -t strict -T $threads
        mv blast_output.txt blast_filtering_results/
        mv blast_filtering_results/ BLAST_FILTERED/
        echo -e "\n======== BLAST FILTERED DONE========\n"

        echo -e "\n======== RUNNING KRAKEN2 WITH DATABASE $DB ========\n"
        mkdir KRAKEN2/
        cd KRAKEN2/
        # Run kraken2
        kraken2 --db $krakenDB --threads $threads ../../../../../$scaffolds \
          >kraken2_output.txt

        if [[ $DB == 'SILVA' ]]; then
          # Now we're gonna edit the output so that is has the same format as
          # CREST output, since we already have a script to deal with
          # SILVA CREST output/taxonomy
          # Extract the taxids column of the standard kraken output:
          cut -f3 kraken2_output.txt >kraken2_taxids.txt
          # Access the SILVA taxonomy file (downloaded taxmap files as in
          # subscript SILVA_SSU_LSU_kraken2_preparation and concatenated them)
          # and generate a file containing one column for each SILVA taxid and
          # one column for the respective SILVA taxonomy path:
          tail -n +2 /hdd1/databases/kraken2_SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_DB_Sep_2020/taxmap_slv_ssu_lsu_ref_nr_138.1.txt |
            cut -f 4,6 | sort -u >SILVA_paths_and_taxids.txt
          # Kraken2 spits out the taxid 0 when no hit is found, but 0 doesn't
          # exist in the SILVA taxonomy, so manually add taxid 0 with path
          # “No hits” to the SILVA path file:
          echo -e "No hits;\t0" >tmp && cat SILVA_paths_and_taxids.txt >>tmp &&
            mv tmp SILVA_paths_and_taxids.txt
          # Merge your kraken2 taxids with the SILVA path file to assign a SILVA
          # taxonomy path to every kraken2 hit:
          mergeFilesOnColumn.pl SILVA_paths_and_taxids.txt kraken2_taxids.txt 2 1 >merged.txt
          cut -f -2 merged.txt | sed 's/;\t/\t/g' >merged_edit.txt # Edit the output
          # Extract the sequence names from the kraken2 output and generate a
          # final file with sequence name, taxid, and SILVA path:
          cut -f 3 kraken2_output.txt >names.txt
          paste names.txt merged_edit.txt | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2}' \
            >kraken2_SILVA_formatted.txt
          # This file has now the same format as the output of CREST and can be
          # translated into NCBI taxonomy the same way as CREST output

          # We run a separate script that was initially made to deal with the
          # SILVA taxonomy of CREST output, by translating SILVA taxonomic paths
          # into NCBI taxids, and use that script on the formatted kraken2 output.
          # The files NCBI_staxids_(non_)scientific were generated by the script
          # SILVA_SSU_LSU_makeblastdb_preparation:
          assign_NCBI_staxids_to_CREST_v4.py /hdd1/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/NCBI_staxids_scientific.txt \
            /hdd1/databases/SILVA_138.1_SSU_LSURef_NR99_tax_silva_trunc_BLAST_DB_Sep_2020/NCBI_staxids_non_scientific.txt \
            kraken2_SILVA_formatted.txt kraken2_SILVA_formatted_with_NCBI_taxids.txt
          mergeFilesOnColumn.pl kraken2_SILVA_formatted_with_NCBI_taxids.txt \
            kraken2_SILVA_formatted.txt 1 1 | cut -f3 >NCBItaxids.txt # Merge SILVA output with taxids and extract taxids
          # We use a separate script to assign taxonomy to NCBI taxids:
          assign_taxonomy_to_NCBI_staxids.sh -b NCBItaxids.txt -c 1 \
            -e ~/.etetoolkit/taxa.sqlite
          sed -i '1d' NCBItaxids_with_taxonomy.txt     # Remove header
          cut -f2 kraken2_output.txt >contig_names.txt # Get contig names from original kraken2 output
          paste contig_names.txt NCBItaxids_with_taxonomy.txt \
            >contigs_with_NCBItaxids_and_taxonomy.txt # Add contig names to taxonomy file
          echo -e "sequence\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
            >kraken2_final.txt && cat contigs_with_NCBItaxids_and_taxonomy.txt \
            >>kraken2_final.txt # Add header

          # Sort files
          mkdir intermediate_files
          mv kraken2_output.txt kraken2_taxids.txt SILVA_paths_and_taxids.txt merged* \
            names.txt kraken2_SILVA_formatted* NCBItaxids* contig* intermediate_files/

        else # NCBI_NT
          cut -f 2-3 kraken2_output.txt >kraken2_output_contig_taxid.txt # Isolate contig names and taxids
          # We use a separate script to assign taxonomy to NCBI taxids:
          assign_taxonomy_to_NCBI_staxids.sh -b kraken2_output_contig_taxid.txt \
            -c 2 -e ~/.etetoolkit/taxa.sqlite
          sed -i '1d' kraken2_output_contig_taxid_with_taxonomy.txt # Remove header
          echo -e "sequence\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
            >kraken2_final.txt && cat kraken2_output_contig_taxid_with_taxonomy.txt \
            >>kraken2_final.txt # Add header

          # Sort files
          mkdir intermediate_files
          mv kraken2_output* intermediate_files/
        fi

        cd ..
        echo -e "\n======== KRAKEN2 WITH DATABASE $DB DONE========\n"

        ######################### Step 6.2: Generating final putput files ################################
        # Each assembler/classification tool output has a different format. We
        # make that format universal with the following code.

        # For loop 5: loop over the classification tools for universal edits:
        for classification_tool in KRAKEN2 BLAST_FILTERED BLAST_FIRST_HIT; do

          echo -e "++++++++ START STEP 6.2: GENERATING FINAL OUTPUT FILES OF READ IN "${trimming_results}"/ AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ AND ASSEMBLY IN FOLDER $assembly_results/ AND CLASSIFICATION IN FOLDER $classification_tool/ ++++++++\n"

          # We need to open a loop for the mappers again to access their files,
          # since we closed that loop earlier:
          for mapper in BWA BOWTIE2; do
            cd $classification_tool

            # We use full paths here to access the separate files.
            # Add assembly sequences to BWA and BOWTIE file:
            if [[ $assembly_results == 'SPADES' || $assembly_results == 'METASPADES' || $assembly_results == 'IDBA_UD' || $assembly_results == 'RNASPADES' || $assembly_results == 'TRANSABYSS' ]]; then
              # Change the sequences' fasta format to tab delimited:
              fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} \
                >./${assembly_results}_tab_to_merge.txt
              # Now the scaffold names can be merged with the mapper output:
              mergeFilesOnColumn.pl ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt \
                ./${assembly_results}_tab_to_merge.txt 2 1 | cut -f1,2,4 |
                awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' \
                  >./final_order_${mapper}_${assembly_results}.txt
              # Final edits:
              echo -e "sequence\tcounts\tassembly_sequence" \
                >./${assembly_results}_final_${mapper}_merge_ready.txt &&
                cat ./final_order_${mapper}_${assembly_results}.txt \
                  >>./${assembly_results}_final_${mapper}_merge_ready.txt
              # Sorting:
              rm ./final_order_${mapper}_${assembly_results}.txt

            elif [[ $assembly_results == 'MEGAHIT' ]]; then
              # Change the sequences' fasta format to tab delimited:
              fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} |
                sed 's/ /\t/g' | cut -f1,4,5 | sed 's/len=//g' \
                  >./${assembly_results}_tab_to_merge.txt
              # Now the scaffold names can be merged with the mapper output:
              mergeFilesOnColumn.pl ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt \
                ./${assembly_results}_tab_to_merge.txt 2 1 | cut -f1,2,4,5 |
                awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $1, $4}' \
                  >./final_order_${mapper}_${assembly_results}.txt
              # Final edits:
              echo -e "sequence\tcontig_length\tcounts\tassembly_sequence" \
                >./${assembly_results}_final_${mapper}_merge_ready.txt &&
                cat ./final_order_${mapper}_${assembly_results}.txt \
                  >>./${assembly_results}_final_${mapper}_merge_ready.txt
              # Sorting:
              rm ./final_order_${mapper}_${assembly_results}.txt

            elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
              # Change the sequences' fasta format to tab delimited:
              fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} |
                sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/ /\t/g' | cut -f1,3,5,6 \
                  >./${assembly_results}_tab_to_merge.txt
              # Now the scaffold names can be merged with the mapper output:
              mergeFilesOnColumn.pl ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt \
                ./${assembly_results}_tab_to_merge.txt 2 1 | cut -f1,2,4,5,6 |
                awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $4, $1, $5}' \
                  >./final_order_${mapper}_${assembly_results}.txt
              # Final edits:
              echo -e "sequence\tcontig_length\tcoverage\tcounts\tassembly_sequence" \
                >./${assembly_results}_final_${mapper}_merge_ready.txt &&
                cat ./final_order_${mapper}_${assembly_results}.txt \
                  >>./${assembly_results}_final_${mapper}_merge_ready.txt
              # Sorting:
              rm ./final_order_${mapper}_${assembly_results}.txt

            else # TRANSABYSS
              # Change the sequences' fasta format to tab delimited:
              fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} |
                sed 's/ /\t/1' | cut -f1,3 >./${assembly_results}_tab_to_merge.txt
              # Now the scaffold names can be merged with the mapper output:
              mergeFilesOnColumn.pl ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt \
                ./${assembly_results}_tab_to_merge.txt 2 1 | cut -f1,2,4 |
                awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' \
                  >./final_order_${mapper}_${assembly_results}.txt
              # Final edits:
              echo -e "sequence\tcounts\tassembly_sequence" \
                >./${assembly_results}_final_${mapper}_merge_ready.txt &&
                cat ./final_order_${mapper}_${assembly_results}.txt \
                  >>./${assembly_results}_final_${mapper}_merge_ready.txt
              # Sorting:
              rm ./final_order_${mapper}_${assembly_results}.txt
            fi

            cd ..

            # Add classification data to the mapper output and do
            # final individual edits:
            if [[ $classification_tool == 'BLAST_FILTERED' ]]; then

              cd $classification_tool
              mkdir FINAL_FILES/

              # We run a simple python script because we need to merge on "outer":
              merge_mapped_reads_and_contigs.py \
                ./blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt \
                ./${assembly_results}_final_${mapper}_merge_ready.txt \
                ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt

              if [ $assembly_results == 'SPADES' ] ||
                [ $assembly_results == 'METASPADES' ]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' |
                  sed 's/length_//g' | sed 's/cov_//g' >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'IDBA_UD' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | cut -f2-16 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'MEGAHIT' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | cut -f2-17 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcontig_length\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'RNASPADES' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' |
                  sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' |
                  cut -f1,2,3,5-21 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | cut -f2-18 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'TRINITY' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' \
                    >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              else #TRANSABYSS
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | sed 's/_/\t/1' >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              fi

              # Final sorting
              mkdir FINAL_FILES/intermediate_files/
              mv *merge* new.txt FINAL_FILES/intermediate_files/

              cd ..
              echo -e "\n======== DONE FINALIZING BLAST_FILTERED FILES FOR MAPPER $mapper =======\n"

            elif [[ $classification_tool == 'BLAST_FIRST_HIT' ]]; then

              cd $classification_tool
              mkdir FINAL_FILES/

              # We run a simple python script because we need to merge on "outer":
              merge_mapped_reads_and_contigs.py ./blast_output_with_taxonomy_and_best_hit.txt \
                ./${assembly_results}_final_${mapper}_merge_ready.txt \
                ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt

              if [ $assembly_results == 'SPADES' ] ||
                [ $assembly_results == 'METASPADES' ]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' |
                  sed 's/length_//g' | sed 's/cov_//g' >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'IDBA_UD' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | cut -f2-16 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'MEGAHIT' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | cut -f2-17 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'RNASPADES' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' |
                  sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' |
                  cut -f1,2,3,5-20 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | cut -f2-18 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'TRINITY' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              else #TRANSABYSS
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | sed 's/_/\t/1' >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              fi

              # Final sorting
              mkdir FINAL_FILES/intermediate_files/
              mv *merge* new.txt FINAL_FILES/intermediate_files/

              cd ..
              echo -e "\n======== DONE FINALIZING BLAST_FIRST_HIT FILES FOR MAPPER $mapper ========\n"

            else # kraken2

              cd $classification_tool
              mkdir FINAL_FILES/

              # We run a simple python script because we need to merge on "outer":
              merge_mapped_reads_and_contigs.py ./kraken2_final.txt \
                ./${assembly_results}_final_${mapper}_merge_ready.txt \
                ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt

              if [ $assembly_results == 'SPADES' ] ||
                [ $assembly_results == 'METASPADES' ]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' |
                  sed 's/length_//g' | sed 's/cov_//g' >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'IDBA_UD' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | cut -f2-20 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'MEGAHIT' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | cut -f2-19 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'RNASPADES' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' |
                  sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' |
                  cut -f1,2,3,5-20 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | cut -f2-20 >./FINAL_FILES/new.txt
                echo -e "sequence_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              elif [[ $assembly_results == 'TRINITY' ]]; then
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' >./FINAL_FILES/new.txt
                echo -e "contig_number\tcontig_length\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              else #TRANSABYSS
                sed '1d' ./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_merged.txt |
                  sed 's/_/\t/1' | sed 's/_/\t/1' >./FINAL_FILES/new.txt
                echo -e "sequence_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" \
                  >./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt &&
                  cat ./FINAL_FILES/new.txt \
                    >>./FINAL_FILES/${trimming_results##*/}_${rrna_filter_results}_${assembly_results}_${mapper}_${DB}_${classification_tool}_pipeline_final.txt

              fi

              # Final sorting
              mkdir FINAL_FILES/intermediate_files/
              mv *merge* new.txt FINAL_FILES/intermediate_files/

              cd ..
              echo -e "\n======== DONE FINALIZING KRAKEN2 FILES FOR MAPPER $mapper ========\n"
            fi
          done
        done
        # After each loop, we have to go back to the directory we started the
        # loop in, and we're using realtive paths for that, making use of the
        # variables generated in the previous loops
        cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)
      done
      cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)
    done
    cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/)
  done
  cd $(realpath --relative-to=$(pwd) $base_directory)
done

# We make a final directory and extract all final files into that directory to
# have everything in one place
mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/
find . -type f -name "*_pipeline_final.txt" -print0 | xargs -0 -Ifile cp file \
  METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_FINAL_FILES/

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $(convertsecs)"
