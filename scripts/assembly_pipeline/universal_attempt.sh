#!/bin/bash

cmd="$0 $@" # Make variable containing full used command to print command in logfile
usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq>
Usage:
	-1 Forward reads trimmed - must state full path from root to the file
	-2 Reverse reads trimmed - must state full path from root to the file
	-h Display this help and exit"

# Set default options
forward_reads=''
reverse_reads=''

# Set specified options
while getopts ':1:2:aDRSMmUrtTBCfsh' opt; do
 	case "${opt}" in
		1) forward_reads="${OPTARG}" ;;
		2) reverse_reads="${OPTARG}" ;;
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
if [[ -z "$forward_reads" || -z "$reverse_reads" ]]
then
   echo -e "-1 and -2 must be set.\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script.\n"
   exit
fi


##################### Write time, options etc. to output ######################
# Define starting time of script for total runtime calculation
start=$(date +%s)
echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
echo -e "=================================================================\n\n"


# Output specified options
echo -e "======== OPTIONS ========\n"

echo -e "Forward reads were defined as $forward_reads.\n"
echo -e "Reverse reads were defined as $reverse_reads.\n"
echo -e "Script started with full command: $cmd\n"

echo -e "++++++++ START RUNNING SCRIPT ========\n"

# Save current path in variable to make navigation between directories easier
base_directory=$(pwd)

######################### Step 1: trimming ################################
echo -e "++++++++ START STEP 1: TRIMMING AND ERROR CORRECTION ========\n"

# Trimming is done with separate script:
trimming_with_phred_scores_and_fastqc_report.sh \
-T /hdd1/programs_for_pilot/Trimmomatic-0.39/trimmomatic-0.39.jar \
-1 $forward_reads -2 $reverse_reads
mv trimming_with_phred_scores_and_fastqc_report_output/ step_1_trimming/

# Running error correction module of SPAdes on trimmed reads
for trimming_results in step_1_trimming/trimmomatic/*; do
	echo -e "\n======== ERROR-CORRECTING READS IN FOLDER $trimming_results ++++++++\n"
	spades.py -1 $trimming_results/*1P.fastq -2 $trimming_results/*2P.fastq \
	--only-error-correction --disable-gzip-output -o $trimming_results/error_correction
	mv $trimming_results/error_correction/corrected/*1P*.fastq \
	$trimming_results/error_correction/corrected/*2P*.fastq $trimming_results
	R1=$(echo $trimming_results/*1P.00.0_0.cor.fastq) \
  && 	mv $trimming_results/*1P.00.0_0.cor.fastq ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
  sed -r -i 's/ BH:.{2,6}//g' ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
	R2=$(echo $trimming_results/*2P.00.0_0.cor.fastq) \
  && mv $trimming_results/*2P.00.0_0.cor.fastq ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
  sed -r -i 's/ BH:.{2,6}//g' ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
	rm -r $trimming_results/error_correction/
	echo -e "\n++++++++ FINISHED ERROR-CORRECTING READS IN FOLDER $trimming_results ++++++++\n"
done

echo -e "++++++++ FINISHED STEP 1: TRIMMING AND ERROR CORRECTION ========\n"

######################### Step 2: rRNA sorting ################################

for trimming_results in step_1_trimming/trimmomatic/*; do
	mkdir $trimming_results/step_2_rrna_sorting/
	cd $trimming_results/step_2_rrna_sorting/

	echo -e "++++++++ START STEP 2: rRNA SORTING OF TRIMMED READS IN FOLDER $trimming_results ++++++++\n"

	echo -e "\n======== RUNNING SORTMERNA ========\n"
	mkdir SORTMERNA/
	sortmerna --ref /hdd1/databases/sortmerna_silva_databases/silva-bac-16s-id90.fasta \
	--ref /hdd1/databases/sortmerna_silva_databases/silva-arc-16s-id95.fasta \
	--ref /hdd1/databases/sortmerna_silva_databases/silva-euk-18s-id95.fasta \
	--reads ../*1P_error_corrected.fastq --reads ../*2P_error_corrected.fastq \
	--paired_in -other -fastx 1 \
	-num_alignments 1 -v \
	-workdir SORTMERNA/
  pwd
	deinterleave_fastq_reads.sh < SORTMERNA/out/aligned.fastq \
	SORTMERNA/out/aligned_R1.fq SORTMERNA/out/aligned_R2.fq
	echo -e "\n======== SORTMERNA DONE ========\n"

	echo -e "\n======== RUNNING rRNAFILTER ========\n"
	mkdir rRNAFILTER/
	cd rRNAFILTER/
	fq2fa ../../*1P_error_corrected.fastq R1.fa
	fq2fa ../../*2P_error_corrected.fastq R2.fa
	wget http://hulab.ucf.edu/research/projects/rRNAFilter/software/rRNAFilter.zip
	unzip rRNAFilter.zip
	cd rRNAFilter/
	java -jar -Xmx7g rRNAFilter_commandline.jar -i ../R1.fa -r 0
	java -jar -Xmx7g rRNAFilter_commandline.jar -i ../R2.fa -r 0
	cd ..
	rm -r rRNAFilter rRNAFilter.zip
	fasta_to_tab R1.fa_rRNA | cut -f 1 | cut -f1 -d " " > names.txt
	fasta_to_tab R2.fa_rRNA | cut -f 1 | cut -f1 -d " " >> names.txt
	sort -u names.txt > names_sorted.txt
	seqtk subseq R1.fa names_sorted.txt > rRNAFilter_paired_R1.fa
	seqtk subseq R2.fa names_sorted.txt > rRNAFilter_paired_R2.fa
	rm names_sorted.txt names.txt
	cd ..
	echo -e "\n======== rRNAFILTER DONE ========\n"

	echo -e "\n======== RUNNING BARRNAP ========\n"
# code for barrnap
	echo -e "\n======== BARRNAP DONE ========\n"

	echo -e "\n======== MAKING FOLDER UNSORTED/ AND COPYING UNSORTED READS IN THERE TO KEEP THE FOLDER STRUCTURE CONSTANT ========\n"
	mkdir UNSORTED/
	cp ../*1P_error_corrected.fastq ../*2P_error_corrected.fastq UNSORTED/

	echo -e "\n++++++++ FINISHED STEP 2: rRNA SORTING OF TRIMMED READS IN FOLDER $trimming_results ++++++++\n"

	######################### Step 3: Assembly ################################

	rrna_filter_results_list="rRNAFILTER SORTMERNA UNSORTED" # <-- MISSES BARRNAP YET, NEEDS TO BE ADDED
  for rrna_filter_results in $rrna_filter_results_list; do
		if [[ $rrna_filter_results == 'rRNAFILTER' ]] ; then
			R1_sorted='rRNAFILTER/rRNAFilter_paired_R1.fa'
			R2_sorted='rRNAFILTER/rRNAFilter_paired_R2.fa'
  	elif [[ $rrna_filter_results == 'SORTMERNA' ]] ; then
			R1_sorted='SORTMERNA/out/aligned_R1.fq'
			R2_sorted='SORTMERNA/out/aligned_R2.fq'
# 	elif [[ $rrna_filter_results == 'BARRNAP' ]] ; then
#			R1_sorted=
#			R2_sorted=
	  else
			R1_sorted='UNSORTED/*1P_error_corrected.fastq'
			R2_sorted='UNSORTED/*2P_error_corrected.fastq'
		fi

		echo -e "++++++++ START STEP 3: ASSEMBLY OF TRIMMED READS IN FOLDER $trimming_results AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ ++++++++\n"
		mkdir $rrna_filter_results/step_3_assembly/
		cd $rrna_filter_results/step_3_assembly/

		echo -e "\n======== RUNNING SPADES ========\n"
		mkdir SPADES/
		spades.py -1 ../../$R1_sorted -2 ../../$R2_sorted --only-assembler -o SPADES/
		echo -e "\n======== SPADES DONE ========\n"

		echo -e "\n======== RUNNING METASPADES ========\n"
		mkdir METASPADES/
		metaspades.py -1 ../../$R1_sorted -2 ../../$R2_sorted --only-assembler -o METASPADES/
		echo -e "\n======== METASPADES DONE ========\n"

		echo -e "\n======== RUNNING MEGAHIT ========\n"
		megahit --presets meta-large -t 16 -1 ../../$R1_sorted -2 ../../$R2_sorted -o MEGAHIT/
		echo -e "\n======== MEGAHIT DONE ========\n"

		echo -e "\n======== RUNNING IDBA_UD ========\n"
		fq2fa --merge --filter ../../$R1_sorted ../../$R2_sorted idba_ud_input.fa
		idba_ud --num_threads 16 --pre_correction -r idba_ud_input.fa \
    -o IDBA_UD/
		mv idba_ud_input.fa IDBA_UD/
		echo -e "\n======== IDBA_UD DONE ========\n"

		echo -e "\n======== RUNNING RNASPADES ========\n"
		mkdir RNASPADES/
		rnaspades.py -1 ../../$R1_sorted -2 ../../$R2_sorted --only-assembler -o RNASPADES/
		echo -e "\n======== RNASPADES DONE ========\n"

		echo -e "\n======== RUNNING IDBA_TRAN ========\n"
		fq2fa --merge --filter ../../$R1_sorted ../../$R2_sorted idba_tran_input.fa
		idba_tran --num_threads 16 --pre_correction -r idba_tran_input.fa \
    -o IDBA_TRAN/
		mv idba_tran_input.fa IDBA_TRAN/
		echo -e "\n======== IDBA_TRAN DONE ========\n"

		echo -e "\n======== RUNNING TRINITY ========\n"
    if [[ $rrna_filter_results == "rRNAFILTER" ]]; then
  		Trinity --seqType fa --max_memory 64G --left ../../$R1_sorted --right \
      ../../$R2_sorted --CPU 16 --output TRINITY/
    else
      Trinity --seqType fq --max_memory 64G --left ../../$R1_sorted --right \
      ../../$R2_sorted --CPU 16 --output TRINITY/
    fi
    cat TRINITY/Trinity.fasta \
    | sed 's/ len/_len/g' > TRINITY/Trinity_with_length.fasta
		echo -e "\n======== TRINITY DONE ========\n"

		echo -e "\n======== RUNNING TRANSABYSS ========\n"
    transabyss --pe ../../$R1_sorted ../../$R2_sorted --threads 16 \
    --outdir TRANSABYSS/
    sed 's/ /_/g' TRANSABYSS/transabyss-final.fa > TRANSABYSS/transabyss-final_edited.fa
		echo -e "\n======== TRANSABYSS DONE ========\n"

		echo -e "\n++++++++ FINISHED STEP 3: ASSEMBLY OF TRIMMED READS IN FOLDER $trimming_results AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ ++++++++\n"

		assembly_results_list="SPADES METASPADES MEGAHIT IDBA_UD RNASPADES IDBA_TRAN TRINITY TRANSABYSS"
		for assembly_results in $assembly_results_list; do
			if [[ $assembly_results == 'SPADES' ]] ; then
				scaffolds='SPADES/scaffolds.fasta'
			elif [[ $assembly_results == 'METASPADES' ]] ; then
				scaffolds='METASPADES/scaffolds.fasta'
			elif [[ $assembly_results == 'MEGAHIT' ]] ; then
				scaffolds='MEGAHIT/final.contigs.fa'
			elif [[ $assembly_results == 'IDBA_UD' ]] ; then
				scaffolds='IDBA_UD/scaffold.fa'
			elif [[ $assembly_results == 'RNASPADES' ]] ; then
				scaffolds='RNASPADES/transcripts.fasta'
			elif [[ $assembly_results == 'IDBA_TRAN' ]] ; then
				scaffolds='IDBA_TRAN/contig.fa'
			elif [[ $assembly_results == 'TRINITY' ]] ; then
				scaffolds='TRINITY/Trinity_with_length.fasta'
			else
				scaffolds='TRANSABYSS/transabyss-final_edited.fa'
			fi


			echo -e "++++++++ START STEP 4: MAPPING OF TRIMMED READS IN FOLDER $trimming_results AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ AND ASSEMBLY IN FOLDER $assembly_results/ ++++++++\n"
			mkdir $assembly_results/step_4_mapping/
			cd $assembly_results/step_4_mapping/

      for mapper in BWA BOWTIE2; do

        mkdir $mapper
        cd $mapper

        if [[ $mapper == 'BWA' ]] ; then
          echo -e "\n======== starting bwa index ========\n"

          bwa index -p bwa_index ../../../$scaffolds

          echo -e "\n======== bwa index complete. Starting bwa ========\n"

          bwa mem -t 10 bwa_index ../../../../../../*1P_error_corrected.fastq ../../../../../../*2P_error_corrected.fastq > ${mapper}_output.sam

          rm bwa_index*

  			else

          echo -e "\n======== bwa complete. Starting bowtie2 index ========\n"

    			bowtie2-build -f ../../../$scaffolds bowtie_index

    			echo -e "\n======== bowtie2 index complete. Starting bowtie2 ========\n"

    			bowtie2 -q -x bowtie_index -1 ../../../../../../*1P_error_corrected.fastq -2 ../../../../../../*2P_error_corrected.fastq -S ${mapper}_output.sam

    			rm bowtie_index*

          echo -e "\n======== bowtie2 complete ========\n"
        fi

  			# Output file (.sam) - edit
  			samtools view -F 4 ${mapper}_output.sam > mapped_reads_${mapper}.sam
  			samtools view -f 4 ${mapper}_output.sam > unmapped_reads_${mapper}.sam
  			cat mapped_reads_${mapper}.sam > mapped_reads_${mapper}.txt
  			cat unmapped_reads_${mapper}.sam > unmapped_reads_${mapper}.txt
  			cut -f3 mapped_reads_${mapper}.txt > mapped_column3_reads_${mapper}.txt
  			cut -f3 unmapped_reads_${mapper}.txt > unmapped_column3_reads_${mapper}.txt
  			sort mapped_column3_reads_${mapper}.txt | uniq -c > sorted_mapped_column3_reads_${mapper}.txt
  			sort unmapped_column3_reads_${mapper}.txt | uniq -c > sorted_unmapped_column3_reads_${mapper}.txt
  			column -t sorted_mapped_column3_reads_${mapper}.txt > aligned_mapped_${mapper}.txt
  			column -t sorted_unmapped_column3_reads_${mapper}.txt > aligned_unmapped_${mapper}.txt
  			sed 's/  */\t/g' aligned_mapped_${mapper}.txt > out_mapped_${mapper}.txt
  			sed 's/  */\t/g' aligned_unmapped_${mapper}.txt > out_unmappped_${mapper}.txt
  			echo -e "counts\tcontig_number" > merge_input_mapped_${mapper}.txt && cat out_mapped_${mapper}.txt >> merge_input_mapped_${mapper}.txt
  			echo -e "counts\tcontig_number" > merge_input_unmapped_${mapper}.txt && cat out_unmappped_${mapper}.txt >> merge_input_unmapped_${mapper}.txt

  			rm *_index* mapped_reads_${mapper}.sam unmapped_reads_${mapper}.sam mapped_reads_${mapper}.txt unmapped_reads_${mapper}.txt mapped_column3_reads_${mapper}.txt unmapped_column3_reads_${mapper}.txt sorted_mapped_column3_reads_${mapper}.txt sorted_unmapped_column3_reads_${mapper}.txt aligned_mapped_${mapper}.txt aligned_unmapped_${mapper}.txt out_mapped_${mapper}.txt out_unmappped_${mapper}.txt
        cd ..
      done

      cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)

			echo -e "++++++++ FINISHED STEP 4: MAPPING FOR TRIMMED READS IN FOLDER $trimming_results AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ AND ASSEMBLY IN FOLDER $assembly_results ========\n"

			echo -e "++++++++ START STEP 5 AND 6: CLASSIFICATION OF ASSEMBLED SCAFFOLDS FROM TRIMMED READS IN FOLDER $trimming_results AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/ AND ASSEMBLY IN FOLDER $assembly_results ========\n"
			mkdir $assembly_results/step_5_reference_DB/
			cd $assembly_results/step_5_reference_DB/

			ref_DB_list="SILVA NCBI_NT"
			for DB in $ref_DB_list; do
				if [[ $DB == "SILVA" ]] ; then
					krakenDB="/hdd1/databases/kraken2_SILVA_DB"
					#centrifugeDB=XX
					blastDB="/hdd1/databases/SILVA_database_mar_2020/SILVA_138_SSURef_NR99_tax_silva.fasta"
				else
					krakenDB="/hdd1/databases/kraken2_nt_DB"
					#centrifugeDB=XX
					blastDB="/hdd1/databases/nt_database_feb_2020_indexed/nt"
				fi

        mkdir $DB
				mkdir $DB/step_6_classification
				cd $DB/step_6_classification

				echo -e "\n======== RUNNING JUSTBLAST WITH DATABASE $DB ========\n"
				justblast ../../../../$scaffolds $blastDB --cpus $threads --evalue 1e-05 \
				--outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend \
				sstart send evalue bitscore staxids" --out_filename blast_output.txt
				echo -e "\n======== JUSTBLAST WITH DATABASE $DB DONE ========\n"

        echo -e "\n======== RUNNING BLAST FIRST HIT ========\n"
				blast_filtering.bash -i blast_output.txt -f blast -t soft -T 16
				cp blast_output.txt blast_filtering_results/
				mv blast_filtering_results/ BLAST_FIRST_HIT/
        echo -e "\n======== BLAST FIRST HIT DONE ========\n"

        echo -e "\n======== RUNNING BLAST FILTERED ========\n"
				blast_filtering.bash -i blast_output.txt -f blast -t strict -T 16
				mv blast_output.txt blast_filtering_results/
				mv blast_filtering_results/ BLAST_FILTERED/
        echo -e "\n======== BLAST FILTERED DONE========\n"

        echo -e "\n======== RUNNING KRAKEN2 WITH DATABASE $DB ========\n"
        mkdir KRAKEN2/
        cd KRAKEN2/
        # Run kraken2
        kraken2 --db $krakenDB --threads 16 \
        ../../../../../$scaffolds > kraken_output.txt

        if [[ $DB == "SILVA" ]] ; then
          # Now we're gonna edit the output so that is has the same format as CREST output,
          # since we already have a script to deal with SILVA CREST output
          # Extract the taxids column of the standard kraken output
          cut -f3 kraken_output.txt > kraken_taxids.txt
          # Access the SILVA taxonomy file and generate a file containing one column for
          # each SILVA taxid and one column for the respective SILVA taxonomy path:
          tail -n +2 /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/taxmap_slv_ssu_ref_nr_138.txt \
          | cut -f 4,6 | sort -u > SILVA_paths_and_taxids.txt
          # Kraken2 spits out the taxid 0 when no hit is found, but 0 doesn't exist in
          # the SILVA taxonomy, so manually add taxid 0 with path “No hits” to the SILVA
          # path file:
          echo -e "No hits;\t0" > tmp && cat SILVA_paths_and_taxids.txt >> tmp \
          && mv tmp SILVA_paths_and_taxids.txt
          # Merge your kraken taxids with the SILVA path file to assign a SILVA taxonomy
          # path to every kraken hit
          mergeFilesOnColumn.pl SILVA_paths_and_taxids.txt kraken_taxids.txt 2 1 > merged.txt
          cut -f -2 merged.txt | sed 's/;\t/\t/g' > merged_edit.txt # Edit the output
          # Extract the sequence names from the kraken output and generate a final file
          # with sequence name, taxid, and SILVA path
          cut -f 3 kraken_output.txt > names.txt
          paste names.txt merged_edit.txt | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2}' \
          > kraken_SILVA_formatted.txt
          # This file has now the same format as the output of CREST and can be translated
          # into NCBI taxonomy the same way as CREST output

          assign_NCBI_staxids_to_CREST_v4.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt \
          /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt \
          kraken_SILVA_formatted.txt kraken_SILVA_formatted_with_NCBI_taxids.txt
          sed -i '1d' kraken_SILVA_formatted_with_NCBI_taxids.txt # Remove header
          mergeFilesOnColumn.pl kraken_SILVA_formatted_with_NCBI_taxids.txt \
          kraken_SILVA_formatted.txt 1 1 > merged_final.txt # Merge SILVA output with taxids
          cut -f3 merged_final.txt > NCBItaxids.txt # Extract taxids
          assign_taxonomy_NCBI_staxids.sh -b NCBItaxids.txt -c 1 -e ~/.etetoolkit/taxa.sqlite
          sed -i '1d' NCBItaxids_with_taxonomy.txt # Remove header
          cut -f2 kraken_output.txt > contig_names.txt # Get contig names from original kraken2 output
          paste contig_names.txt NCBItaxids_with_taxonomy.txt \
          > contigs_with_NCBItaxids_and_taxonomy.txt # Add contig names to taxonomy file
          echo -e "sequence\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
          > kraken_final.txt && cat contigs_with_NCBItaxids_and_taxonomy.txt \
          >> kraken_final.txt # Add header

          # Sort files
          mkdir intermediate_files
          mv kraken_output.txt kraken_taxids.txt SILVA_paths_and_taxids.txt merged* \
          names.txt kraken_SILVA_formatted* NCBItaxids* contig* intermediate_files/

        else
          cut -f 2-3 kraken_output.txt > kraken_output_contig_taxid.txt # Isolate contig names and taxids
          assign_taxonomy_NCBI_staxids.sh -b kraken_output_contig_taxid.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
          sed -i '1d' kraken_output_contig_taxid_with_taxonomy.txt # Remove header
          echo -e "sequence\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" \
          > kraken_final.txt && cat kraken_output_contig_taxid_with_taxonomy.txt \
          >> kraken_final.txt # Add header

          # Sort files
          mkdir intermediate_files
          mv kraken_output* intermediate_files/
        fi

        cd ..
        echo -e "\n======== KRAKEN2 WITH DATABASE $DB DONE========\n"

				# centrifuge code

				classification_tool="BLAST_FILTERED BLAST_FIRST_HIT KRAKEN2"
					cd $classification_tool
					mkdir FINAL_FILES

					# Add assembly sequence to BWA and BOWTIE file
					if [ $assembly_results == 'SPADES' ] || [ $assembly_results == 'METASPADES' ] || [ $assembly_results == 'IDBA_UD' ] || [ $assembly_results == 'RNASPADES' ] || [ $assembly_results == 'TRANSABYSS' ] ; then

							fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} > ./FINAL_FILES/${assembly_results}_tab_to_merge.txt
							mergeFilesOnColumn.pl ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt ./FINAL_FILES/${assembly_results}_tab_to_merge.txt 2 1 > ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt
							cut -f1,2,4 ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt > ./FINAL_FILES/important_column_3_${mapper}_${assembly_results}.txt
							awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./FINAL_FILES/important_column_3_${mapper}_${assembly_results}.txt > ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt
							echo -e "sequence\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_final_${mapper}_merge_ready.txt && cat ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt >> ./FINAL_FILES/${assembly_results}_final_${mapper}_merge_ready.txt
							rm ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt ./FINAL_FILES/important_column_3_${mapper}_${assembly_results}.txt ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt

					elif [[ $assembly_results == 'MEGAHIT' ]] ; then

							fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} > ./FINAL_FILES/fasta_to_tabbed.txt
							sed 's/ /\t/g' ./FINAL_FILES/fasta_to_tabbed.txt| cut -f1,4,5 | sed 's/len=//g' > ./FINAL_FILES/${assembly_results}_tab_to_merge.txt
							rm ./FINAL_FILES/fasta_to_tabbed.txt

							mergeFilesOnColumn.pl ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt ./FINAL_FILES/${assembly_results}_tab_to_merge.txt 2 1 > ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt
							cut -f1,2,4,5 ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt > ./FINAL_FILES/important_column_4_${mapper}_${assembly_results}.txt
							awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $1, $4}' ./FINAL_FILES/important_column_4_${mapper}_${assembly_results}.txt > ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt
							echo -e "sequence\tcontig_length\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_final_${mapper}_merge_ready.txt && cat ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt >> ./FINAL_FILES/${assembly_results}_final_${mapper}_merge_ready.txt
							rm ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt ./FINAL_FILES/important_column_4_${mapper}_${assembly_results}.txt ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt

					elif [[ $assembly_results == 'IDBA_TRAN' ]] ; then

							fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} > ./FINAL_FILES/fasta_to_tabbed.txt
							sed 's/_/\t/2' ./FINAL_FILES/fasta_to_tabbed.txt | sed 's/_/\t/3' | sed 's/ /\t/g' | cut -f1,3,5,6 > ./FINAL_FILES/${assembly_results}_tab_to_merge.txt
							rm ./FINAL_FILES/fasta_to_tabbed.txt


							mergeFilesOnColumn.pl ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt ./FINAL_FILES/${assembly_results}_tab_to_merge.txt 2 1 > ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt
							cut -f1,2,4,5,6 ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt > ./FINAL_FILES/important_columns_${mapper}_${assembly_results}.txt
							awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $3, $4, $1, $5}' ./FINAL_FILES/important_columns_${mapper}_${assembly_results}.txt > ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt
							echo -e "sequence\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_final_${mapper}_merge_ready.txt && cat ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt >> ./FINAL_FILES/${assembly_results}_final_${mapper}_merge_ready.txt
							rm ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt ./FINAL_FILES/important_columns_${mapper}_${assembly_results}.txt ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt

					else

							fasta_to_tab ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${scaffolds} > ./FINAL_FILES/fasta_to_tabbed.txt
							sed 's/ /\t/1' ./FINAL_FILES/fasta_to_tabbed.txt | cut -f1,3 > ./FINAL_FILES/${assembly_results}_tab_to_merge.txt
							rm ./FINAL_FILES/fasta_to_tabbed.txt

							mergeFilesOnColumn.pl ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_4_mapping/${mapper}/merge_input_mapped_${mapper}.txt ./FINAL_FILES/${assembly_results}_tab_to_merge.txt 2 1 > ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt
							cut -f1,2,4 ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt > ./FINAL_FILES/important_column_3_${mapper}_${assembly_results}.txt
							awk 'BEGIN {FS="\t";OFS="\t"} {print $2, $1, $3}' ./FINAL_FILES/important_column_3_${mapper}_${assembly_results}.txt > ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt
							echo -e "sequence\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_final_${mapper}_merge_ready.txt && cat ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt >> ./FINAL_FILES/${assembly_results}_final_${mapper}_merge_ready.txt
							rm ./FINAL_FILES/merged_original_${mapper}_${assembly_results}.txt ./FINAL_FILES/important_column_3_${mapper}_${assembly_results}.txt ./FINAL_FILES/final_order_${mapper}_${assembly_results}.txt
					fi

					# Add classification data to BWA and BOWTIE2
				for tool in $classification_tool; do
					if [[ $classification_tool == 'BLAST_FILTERED' ]] ; then
						tool='BLAST_FILTERED/blast_output_with_taxonomy_and_bitscore_threshold_and_bitscore_filter_and_pident_cutoff_and_LCA.txt'

						merge_mapped_reads_and_contigs.py ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_5_reference_DB/${ref_DB_list}/step_6_classification/${tool} ./FINAL_FILES/${assembly_results}_final_${mapper}_merge_ready.txt ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt

							if [ $assembly_results == 'SPADES' ] || [ $assembly_results == 'METASPADES' ] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' | sed 's/length_//g' | sed 's/cov_//g' > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'IDBA_UD' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | cut -f2-16 > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'MEGAHIT' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | cut -f2-17 > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcontig_length\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'RNASPADES' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' | sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' | cut -f1,2,3,5-21 > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'IDBA_TRAN' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt |	sed 's/_/\t/1' | cut -f2-18 > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'TRINITY' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							else #TRANSABYSS

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | sed 's/_/\t/1' > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tlowest_hit\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							fi

					elif [[ $classification_tool == 'BLAST_FIRST_HIT' ]] ; then
 	 					tool='BLAST_FIRST_HIT/blast_output_with_taxonomy_and_best_hit.txt'

							if [ $assembly_results == 'SPADES' ] || [ $assembly_results == 'METASPADES' ] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' | sed 's/length_//g' | sed 's/cov_//g' > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'IDBA_UD' ]] ; then

									sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | cut -f2-15 > ./FINAL_FILES/new.txt
									echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'MEGAHIT' ]] ; then

									sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | cut -f2-16 > ./FINAL_FILES/new.txt
									echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'RNASPADES' ]] ; then

									sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' | sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' | cut -f1,2,3,5-20 > ./FINAL_FILES/new.txt
									echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'IDBA_TRAN' ]] ; then

									sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | cut -f2-17 > ./FINAL_FILES/new.txt
									echo -e "sequence_number\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'TRINITY' ]] ; then

									sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > ./FINAL_FILES/new.txt
									echo -e "sequence_number\tcontig_length\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							else #TRANSABYSS

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | sed 's/_/\t/1' > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tcontig_coverage\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							fi

					else
						tool='KRAKEN2/kraken_final.txt'

							if [ $assembly_results == 'SPADES' ] || [ $assembly_results == 'METASPADES' ] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/2' | sed 's/_/\t/3' | sed 's/NODE_//g' | sed 's/length_//g' | sed 's/cov_//g' > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'IDBA_UD' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | cut -f2-20 > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'MEGAHIT' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' ./FINAL_FILES/new.txt | cut -f2-19 > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'RNASPADES' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/2' |sed 's/_/\t/3' | sed 's/NODE_//g' | sed 's/length_//g' | sed 's/cov_//g' | sed 's/_/\t/1' | cut -f1,2,3,5-20 > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'IDBA_TRAN' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | cut -f2-20 > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							elif [[ $assembly_results == 'TRINITY' ]] ; then

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/5' | sed 's/len=//g' | sed 's/TRINITY_//g' > ./FINAL_FILES/new.txt
								echo -e "contig_number\tcontig_length\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							else #TRANSABYSS

								sed '1d' ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_merged.txt | sed 's/_/\t/1' | sed 's/_/\t/1' > ./FINAL_FILES/new.txt
								echo -e "sequence_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt && cat ./FINAL_FILES/new.txt >> ./FINAL_FILES/${assembly_results}_${ref_DB_list}_${classification_tool}_${mapper}_final.txt

							fi
					fi

					cd ..
					echo -e "\n========DONE========\n"

					cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_5_reference_DB/${ref_DB_list}/step_6_classification/)
				done
				cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/${assembly_results}/step_5_reference_DB/)
			done
			cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)
		done
		cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/)
	done
	cd $(realpath --relative-to=$(pwd) $base_directory)
done

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"
