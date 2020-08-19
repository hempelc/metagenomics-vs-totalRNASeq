#################################SPADES

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/BLAST_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_nt_with_taxonomy_merge.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/BLAST_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/BLAST_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./SPADES/BLAST_output_nt_with_taxonomy_merge.txt ./SPADES/final_bwa_merge_ready.txt ./SPADES/BLAST_nt_BWA.txt

		# Move all files - easy to find!!!
		mkdir ./SPADES/MERGE_FILES
		mv ./SPADES/final_bwa_merge_ready.txt ./SPADES/final_bowtie2_merge_ready.txt ./SPADES/CREST_merge.txt ./SPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./SPADES/BLAST_output_nt_with_taxonomy_merge.txt ./SPADES/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./SPADES/CREST_BWA.txt > ./SPADES/new.txt &&  sed 's/_/\t/3' ./SPADES/new.txt > ./SPADES/new2.txt
		sed 's/NODE_//g' ./SPADES/new2.txt > ./SPADES/new3.txt && sed 's/length_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/cov_//g' ./SPADES/new4.txt > ./SPADES/new5.txt
		sed '1d' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/CREST_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/CREST_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/CREST_BWA.txt

		sed 's/_/\t/2' ./SPADES/CREST_BOWTIE2.txt > ./SPADES/new.txt &&  sed 's/_/\t/3' ./SPADES/new.txt > ./SPADES/new2.txt
		sed 's/NODE_//g' ./SPADES/new2.txt > ./SPADES/new3.txt && sed 's/length_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/cov_//g' ./SPADES/new4.txt > ./SPADES/new5.txt
		sed '1d' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/CREST_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/CREST_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./SPADES/BLAST_nt_BWA.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_nt_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_nt_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_nt_BWA.txt

		sed '1d' ./SPADES/BLAST_nt_BOWTIE2.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_nt_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_nt_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./SPADES/BLAST_SILVA_BWA.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_SILVA_BWA_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_SILVA_BWA_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_SILVA_BWA.txt

		sed '1d' ./SPADES/BLAST_SILVA_BOWTIE2.txt > ./SPADES/new.txt
		sed 's/_/\t/2' ./SPADES/new.txt > ./SPADES/new2.txt &&  sed 's/_/\t/3' ./SPADES/new2.txt > ./SPADES/new3.txt
		sed 's/NODE_//g' ./SPADES/new3.txt > ./SPADES/new4.txt && sed 's/length_//g' ./SPADES/new4.txt > ./SPADES/new5.txt && sed 's/cov_//g' ./SPADES/new5.txt > ./SPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./SPADES/BLAST_SILVA_BOWTIE2_final.txt && cat ./SPADES/new6.txt >> ./SPADES/BLAST_SILVA_BOWTIE2_final.txt
		rm ./SPADES/new*.txt ./SPADES/BLAST_SILVA_BOWTIE2.txt

		mkdir ./SPADES/FINAL_FILES_SPADES
		mv ./SPADES/CREST_BWA_final.txt ./SPADES/CREST_BOWTIE2_final.txt ./SPADES/BLAST_SILVA_BOWTIE2_final.txt ./SPADES/BLAST_nt_BOWTIE2_final.txt ./SPADES/BLAST_SILVA_BWA_final.txt ./SPADES/BLAST_nt_BWA_final.txt ./SPADES/FINAL_FILES_SPADES/

#################################METASPADES

		# Prepare CREST file - remove 1st/3rd column, add a header
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./METASPADES/CLASSIFICATION/CREST/otus.csv ./METASPADES/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./METASPADES/CLASSIFICATION/CREST/CREST_output.txt > ./METASPADES/CREST_seperated.txt
		cut -f2,4 ./METASPADES/CREST_seperated.txt > ./METASPADES/CREST_header.txt
		sed '1d' ./METASPADES/CREST_header.txt > ./METASPADES/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./METASPADES/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./METASPADES/
		sed '1d' ./METASPADES/CREST_tax_ready_with_taxonomy.txt > ./METASPADES/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./METASPADES/CREST_merge.txt && cat ./METASPADES/CREST_tax_ready_with_taxonomy_noheader.txt >> ./METASPADES/CREST_merge.txt

		rm ./METASPADES/CREST_seperated.txt ./METASPADES/CREST_header.txt ./METASPADES/CREST_tax_ready.txt ./METASPADES/CREST_tax_ready_with_taxonomy.txt ./METASPADES/CREST_tax_ready_with_taxonomy_noheader.txt

		echo "CREST file is ready."

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./METASPADES/CREST_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/CREST_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/BLAST_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_nt_with_taxonomy_merge.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/BLAST_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/BLAST_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./METASPADES/BLAST_output_nt_with_taxonomy_merge.txt ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/BLAST_nt_BWA.txt

		# Move all files - easy to find!!!
		mkdir ./METASPADES/MERGE_FILES
		mv ./METASPADES/final_bwa_merge_ready.txt ./METASPADES/final_bowtie2_merge_ready.txt ./METASPADES/CREST_merge.txt ./METASPADES/BLAST_output_SILVA_with_taxonomy_merge.txt ./METASPADES/BLAST_output_nt_with_taxonomy_merge.txt ./METASPADES/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./METASPADES/CREST_BWA.txt > ./METASPADES/new.txt &&  sed 's/_/\t/3' ./METASPADES/new.txt > ./METASPADES/new2.txt
		sed 's/NODE_//g' ./METASPADES/new2.txt > ./METASPADES/new3.txt && sed 's/length_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/cov_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt
		sed '1d' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/CREST_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/CREST_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/CREST_BWA.txt

		sed 's/_/\t/2' ./METASPADES/CREST_BOWTIE2.txt > ./METASPADES/new.txt &&  sed 's/_/\t/3' ./METASPADES/new.txt > ./METASPADES/new2.txt
		sed 's/NODE_//g' ./METASPADES/new2.txt > ./METASPADES/new3.txt && sed 's/length_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/cov_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt
		sed '1d' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/CREST_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/CREST_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./METASPADES/BLAST_nt_BWA.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_nt_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_nt_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_nt_BWA.txt

		sed '1d' ./METASPADES/BLAST_nt_BOWTIE2.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_nt_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_nt_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./METASPADES/BLAST_SILVA_BWA.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_SILVA_BWA_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_SILVA_BWA_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_SILVA_BWA.txt

		sed '1d' ./METASPADES/BLAST_SILVA_BOWTIE2.txt > ./METASPADES/new.txt
		sed 's/_/\t/2' ./METASPADES/new.txt > ./METASPADES/new2.txt &&  sed 's/_/\t/3' ./METASPADES/new2.txt > ./METASPADES/new3.txt
		sed 's/NODE_//g' ./METASPADES/new3.txt > ./METASPADES/new4.txt && sed 's/length_//g' ./METASPADES/new4.txt > ./METASPADES/new5.txt && sed 's/cov_//g' ./METASPADES/new5.txt > ./METASPADES/new6.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./METASPADES/BLAST_SILVA_BOWTIE2_final.txt && cat ./METASPADES/new6.txt >> ./METASPADES/BLAST_SILVA_BOWTIE2_final.txt
		rm ./METASPADES/new*.txt ./METASPADES/BLAST_SILVA_BOWTIE2.txt

		mkdir ./METASPADES/FINAL_FILES_METASPADES
		mv ./METASPADES/CREST_BWA_final.txt ./METASPADES/CREST_BOWTIE2_final.txt ./METASPADES/BLAST_SILVA_BOWTIE2_final.txt ./METASPADES/BLAST_nt_BOWTIE2_final.txt ./METASPADES/BLAST_SILVA_BWA_final.txt ./METASPADES/BLAST_nt_BWA_final.txt ./METASPADES/FINAL_FILES_METASPADES/

		echo "final METASPADES output generated."

	#################################MEGAHIT

			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./MEGAHIT/CLASSIFICATION/CREST/otus.csv ./MEGAHIT/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./MEGAHIT/CLASSIFICATION/CREST/CREST_output.txt > ./MEGAHIT/CREST_seperated.txt
			cut -f2,4 ./MEGAHIT/CREST_seperated.txt > ./MEGAHIT/CREST_header.txt
			sed '1d' ./MEGAHIT/CREST_header.txt > ./MEGAHIT/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./MEGAHIT/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/
			sed '1d' ./MEGAHIT/CREST_tax_ready_with_taxonomy.txt > ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./MEGAHIT/CREST_merge.txt && cat ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt >> ./MEGAHIT/CREST_merge.txt

			rm ./MEGAHIT/CREST_seperated.txt ./MEGAHIT/CREST_header.txt ./MEGAHIT/CREST_tax_ready.txt ./MEGAHIT/CREST_tax_ready_with_taxonomy.txt ./MEGAHIT/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./MEGAHIT/CREST_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/CREST_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./MEGAHIT/MERGE_FILES
			mv ./MEGAHIT/final_bwa_merge_ready.txt ./MEGAHIT/final_bowtie2_merge_ready.txt ./MEGAHIT/CREST_merge.txt ./MEGAHIT/BLAST_output_SILVA_with_taxonomy_merge.txt ./MEGAHIT/BLAST_output_nt_with_taxonomy_merge.txt ./MEGAHIT/MERGE_FILES/

			# Edit k*_ - in all files
			# CREST files
			sed '1d' ./MEGAHIT/CREST_BWA.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-19 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/CREST_BWA_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/CREST_BWA_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/CREST_BWA.txt

			sed '1d' ./MEGAHIT/CREST_BOWTIE2.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-19 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/CREST_BOWTIE2_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/CREST_BOWTIE2_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./MEGAHIT/BLAST_nt_BWA.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-30 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_BWA_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_nt_BWA_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_nt_BWA.txt

			sed '1d' ./MEGAHIT/BLAST_nt_BOWTIE2.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-30 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_nt_BOWTIE2_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_nt_BOWTIE2_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_nt_BOWTIE2.txt


			# BLAST SILVA files
			sed '1d' ./MEGAHIT/BLAST_SILVA_BWA.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-30 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_BWA_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_SILVA_BWA_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_SILVA_BWA.txt

			sed '1d' ./MEGAHIT/BLAST_SILVA_BOWTIE2.txt > ./MEGAHIT/new.txt
			sed 's/_/\t/1' ./MEGAHIT/new.txt > ./MEGAHIT/new2.txt
			cut -f2-30 ./MEGAHIT/new2.txt > ./MEGAHIT/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcounts\tassembly_sequence" > ./MEGAHIT/BLAST_SILVA_BOWTIE2_final.txt && cat ./MEGAHIT/new3.txt >> ./MEGAHIT/BLAST_SILVA_BOWTIE2_final.txt
			rm ./MEGAHIT/new*.txt ./MEGAHIT/BLAST_SILVA_BOWTIE2.txt

			mkdir ./MEGAHIT/FINAL_FILES_MEGAHIT
			mv ./MEGAHIT/CREST_BWA_final.txt ./MEGAHIT/CREST_BOWTIE2_final.txt ./MEGAHIT/BLAST_SILVA_BOWTIE2_final.txt ./MEGAHIT/BLAST_nt_BOWTIE2_final.txt ./MEGAHIT/BLAST_SILVA_BWA_final.txt ./MEGAHIT/BLAST_nt_BWA_final.txt ./MEGAHIT/FINAL_FILES_MEGAHIT/




#################################IDBA_UD



			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_UD/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./IDBA_UD/BLAST_output_nt_with_taxonomy.txt > ./IDBA_UD/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt && cat ./IDBA_UD/BLAST_output_nt_with_taxonomy_noheader.txt >> ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./IDBA_UD/BLAST_output_SILVA_with_taxonomy.txt > ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./IDBA_UD/BLAST_output_nt_with_taxonomy.txt ./IDBA_UD/BLAST_output_nt_with_taxonomy_noheader.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."

			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA_UD/CLASSIFICATION/CREST/otus.csv ./IDBA_UD/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./IDBA_UD/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA_UD/CREST_seperated.txt
			cut -f2,4 ./IDBA_UD/CREST_seperated.txt > ./IDBA_UD/CREST_header.txt
			sed '1d' ./IDBA_UD/CREST_header.txt > ./IDBA_UD/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_UD/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./IDBA_UD/
			sed '1d' ./IDBA_UD/CREST_tax_ready_with_taxonomy.txt > ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_UD/CREST_merge.txt && cat ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA_UD/CREST_merge.txt

			rm ./IDBA_UD/CREST_seperated.txt ./IDBA_UD/CREST_header.txt ./IDBA_UD/CREST_tax_ready.txt ./IDBA_UD/CREST_tax_ready_with_taxonomy.txt ./IDBA_UD/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./IDBA_UD/CREST_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_UD/CREST_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./IDBA_UD/MERGE_FILES
			mv ./IDBA_UD/final_bwa_merge_ready.txt ./IDBA_UD/final_bowtie2_merge_ready.txt ./IDBA_UD/CREST_merge.txt ./IDBA_UD/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_UD/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_UD/MERGE_FILES/

			# Edit contig-* in all the final files
			# CREST files
			sed '1d' ./IDBA_UD/CREST_BWA.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-20 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/CREST_BWA_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/CREST_BWA_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/CREST_BWA.txt

			sed '1d' ./IDBA_UD/CREST_BOWTIE2.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-20 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/CREST_BOWTIE2_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/CREST_BOWTIE2_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./IDBA_UD/BLAST_nt_BWA.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-31 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_BWA_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_nt_BWA_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_nt_BWA.txt

			sed '1d' ./IDBA_UD/BLAST_nt_BOWTIE2.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-31 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_nt_BOWTIE2_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_nt_BOWTIE2_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_nt_BOWTIE2.txt

			# BLAST SILVA files
			sed '1d' ./IDBA_UD/BLAST_SILVA_BWA.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-31 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_BWA_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_SILVA_BWA_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_SILVA_BWA.txt

			sed '1d' ./IDBA_UD/BLAST_SILVA_BOWTIE2.txt > ./IDBA_UD/new.txt
			sed 's/_/\t/1' ./IDBA_UD/new.txt > ./IDBA_UD/new2.txt
			cut -f2-31 ./IDBA_UD/new2.txt > ./IDBA_UD/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcoverage\tcounts\tassembly_sequence" > ./IDBA_UD/BLAST_SILVA_BOWTIE2_final.txt && cat ./IDBA_UD/new3.txt >> ./IDBA_UD/BLAST_SILVA_BOWTIE2_final.txt
			rm ./IDBA_UD/new*.txt ./IDBA_UD/BLAST_SILVA_BOWTIE2.txt

			mkdir ./IDBA_UD/FINAL_FILES_IDBA_UD
			mv ./IDBA_UD/CREST_BWA_final.txt ./IDBA_UD/CREST_BOWTIE2_final.txt ./IDBA_UD/BLAST_SILVA_BOWTIE2_final.txt ./IDBA_UD/BLAST_nt_BOWTIE2_final.txt ./IDBA_UD/BLAST_SILVA_BWA_final.txt ./IDBA_UD/BLAST_nt_BWA_final.txt ./IDBA_UD/FINAL_FILES_IDBA_UD/

			echo "final IDBA_UD output generated."



######################### Beginning of RNA pipelines ########################


#################################RNASPADES


		# Prepare CREST file - remove 1st/3rd column, add a header
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./rnaSPAdes/CLASSIFICATION/CREST/otus.csv ./rnaSPAdes/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./rnaSPAdes/CLASSIFICATION/CREST/CREST_output.txt > ./rnaSPAdes/CREST_seperated.txt
		cut -f2,4 ./rnaSPAdes/CREST_seperated.txt > ./rnaSPAdes/CREST_header.txt
		sed '1d' ./rnaSPAdes/CREST_header.txt > ./rnaSPAdes/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./rnaSPAdes/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./rnaSPAdes/
		sed '1d' ./rnaSPAdes/CREST_tax_ready_with_taxonomy.txt > ./rnaSPAdes/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./rnaSPAdes/CREST_merge.txt && cat ./rnaSPAdes/CREST_tax_ready_with_taxonomy_noheader.txt >> ./rnaSPAdes/CREST_merge.txt

		rm ./rnaSPAdes/CREST_seperated.txt ./rnaSPAdes/CREST_header.txt ./rnaSPAdes/CREST_tax_ready.txt ./rnaSPAdes/CREST_tax_ready_with_taxonomy.txt ./rnaSPAdes/CREST_tax_ready_with_taxonomy_noheader.txt

		echo "CREST file is ready."


		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./rnaSPAdes/CREST_merge.txt ./rnaSPAdes/final_bwa_merge_ready.txt ./rnaSPAdes/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./rnaSPAdes/CREST_merge.txt ./rnaSPAdes/final_bowtie2_merge_ready.txt ./rnaSPAdes/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_merge.txt ./rnaSPAdes/final_bowtie2_merge_ready.txt ./rnaSPAdes/BLAST_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./rnaSPAdes/BLAST_output_nt_with_taxonomy_merge.txt ./rnaSPAdes/final_bowtie2_merge_ready.txt ./rnaSPAdes/BLAST_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_merge.txt ./rnaSPAdes/final_bwa_merge_ready.txt ./rnaSPAdes/BLAST_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./rnaSPAdes/BLAST_output_nt_with_taxonomy_merge.txt ./rnaSPAdes/final_bwa_merge_ready.txt ./rnaSPAdes/BLAST_nt_BWA.txt


		# Move all files - easy to find!!!
		mkdir ./rnaSPAdes/MERGE_FILES
		mv ./rnaSPAdes/final_bwa_merge_ready.txt ./rnaSPAdes/final_bowtie2_merge_ready.txt ./rnaSPAdes/CREST_merge.txt ./rnaSPAdes/BLAST_output_SILVA_with_taxonomy_merge.txt ./rnaSPAdes/BLAST_output_nt_with_taxonomy_merge.txt ./rnaSPAdes/MERGE_FILES/


		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed 's/_/\t/2' ./rnaSPAdes/CREST_BWA.txt > ./rnaSPAdes/new.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt
		sed 's/NODE_//g' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt && sed 's/length_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/cov_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt
		sed '1d' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-20 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/CREST_BWA_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/CREST_BWA_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/CREST_BWA.txt

		sed 's/_/\t/2' ./rnaSPAdes/CREST_BOWTIE2.txt > ./rnaSPAdes/new.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt
		sed 's/NODE_//g' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt && sed 's/length_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/cov_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt
		sed '1d' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-20 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/CREST_BOWTIE2_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/CREST_BOWTIE2_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/CREST_BOWTIE2.txt


		# BLAST nt files
		sed '1d' ./rnaSPAdes/BLAST_nt_BWA.txt > ./rnaSPAdes/new.txt
		sed 's/_/\t/2' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt
		sed 's/NODE_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/length_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt && sed 's/cov_//g' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-31 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/BLAST_nt_BWA_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/BLAST_nt_BWA_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/BLAST_nt_BWA.txt

		sed '1d' ./rnaSPAdes/BLAST_nt_BOWTIE2.txt > ./rnaSPAdes/new.txt
		sed 's/_/\t/2' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt
		sed 's/NODE_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/length_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt && sed 's/cov_//g' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-31 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/BLAST_nt_BOWTIE2_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/BLAST_nt_BOWTIE2_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/BLAST_nt_BOWTIE2.txt


		# BLAST SILVA files
		sed '1d' ./rnaSPAdes/BLAST_SILVA_BWA.txt > ./rnaSPAdes/new.txt
		sed 's/_/\t/2' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt
		sed 's/NODE_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/length_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt && sed 's/cov_//g' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-31 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/BLAST_SILVA_BWA_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/BLAST_SILVA_BWA_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/BLAST_SILVA_BWA.txt

		sed '1d' ./rnaSPAdes/BLAST_SILVA_BOWTIE2.txt > ./rnaSPAdes/new.txt
		sed 's/_/\t/2' ./rnaSPAdes/new.txt > ./rnaSPAdes/new2.txt &&  sed 's/_/\t/3' ./rnaSPAdes/new2.txt > ./rnaSPAdes/new3.txt
		sed 's/NODE_//g' ./rnaSPAdes/new3.txt > ./rnaSPAdes/new4.txt && sed 's/length_//g' ./rnaSPAdes/new4.txt > ./rnaSPAdes/new5.txt && sed 's/cov_//g' ./rnaSPAdes/new5.txt > ./rnaSPAdes/new6.txt
		sed 's/_/\t/1' ./rnaSPAdes/new6.txt > ./rnaSPAdes/new7.txt
		cut -f1,2,3,5-31 ./rnaSPAdes/new7.txt > ./rnaSPAdes/new8.txt
		echo -e "contig_number\tcontig_length\tcontig_coverage\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./rnaSPAdes/BLAST_SILVA_BOWTIE2_final.txt && cat ./rnaSPAdes/new8.txt >> ./rnaSPAdes/BLAST_SILVA_BOWTIE2_final.txt
		rm ./rnaSPAdes/new*.txt ./rnaSPAdes/BLAST_SILVA_BOWTIE2.txt

		mkdir ./rnaSPAdes/FINAL_FILES_rnaSPAdes
		mv ./rnaSPAdes/CREST_BWA_final.txt ./rnaSPAdes/CREST_BOWTIE2_final.txt ./rnaSPAdes/BLAST_SILVA_BOWTIE2_final.txt ./rnaSPAdes/BLAST_nt_BOWTIE2_final.txt ./rnaSPAdes/BLAST_SILVA_BWA_final.txt ./rnaSPAdes/BLAST_nt_BWA_final.txt ./rnaSPAdes/FINAL_FILES_rnaSPAdes/

		echo "final rnaSPAdes output generated."

#################################IDBA_TRAN

			# Assign BLAST taxonomy - nt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_nt -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_TRAN/BLAST_output_nt_with_taxonomy.txt
			sed '1d' ./IDBA_TRAN/BLAST_output_nt_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt && cat ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_noheader.txt >> ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt

			# Assign BLAST taxonomy - SILVA
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CLASSIFICATION/BLAST/BLAST_output_SILVA -c 13 -e ~/.etetoolkit/taxa.sqlite
			mv _with_taxonomy.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy.txt
			sed '1d' ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy.txt > ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_noheader.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt && cat ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_noheader.txt >> ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt

			rm ./IDBA_TRAN/BLAST_output_nt_with_taxonomy.txt ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_noheader.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_noheader.txt

			echo "taxonomy has been assigned to BLAST files."
			echo "BLAST files are ready."

			# Prepare CREST file - remove 1st/3rd column, add a header
			assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./IDBA_TRAN/CLASSIFICATION/CREST/otus.csv ./IDBA_TRAN/CLASSIFICATION/CREST/CREST_output.txt
			sed 's|:|\t|g' ./IDBA_TRAN/CLASSIFICATION/CREST/CREST_output.txt > ./IDBA_TRAN/CREST_seperated.txt
			cut -f2,4 ./IDBA_TRAN/CREST_seperated.txt > ./IDBA_TRAN/CREST_header.txt
			sed '1d' ./IDBA_TRAN/CREST_header.txt > ./IDBA_TRAN/CREST_tax_ready.txt
			assign_taxonomy_NCBI_staxids.sh -b ./IDBA_TRAN/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
			mv CREST_tax_ready_with_taxonomy.txt ./IDBA_TRAN/
			sed '1d' ./IDBA_TRAN/CREST_tax_ready_with_taxonomy.txt > ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt

			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./IDBA_TRAN/CREST_merge.txt && cat ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt >> ./IDBA_TRAN/CREST_merge.txt

			rm ./IDBA_TRAN/CREST_seperated.txt ./IDBA_TRAN/CREST_header.txt ./IDBA_TRAN/CREST_tax_ready.txt ./IDBA_TRAN/CREST_tax_ready_with_taxonomy.txt ./IDBA_TRAN/CREST_tax_ready_with_taxonomy_noheader.txt

			echo "CREST file is ready."

			# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/CREST_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/CREST_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_SILVA_BOWTIE2.txt
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/BLAST_nt_BOWTIE2.txt

			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_SILVA_BWA.txt
			merge_mapped_reads_and_contigs.py ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/BLAST_nt_BWA.txt

			# Move all files - easy to find!!!
			mkdir ./IDBA_TRAN/MERGE_FILES
			mv ./IDBA_TRAN/final_bwa_merge_ready.txt ./IDBA_TRAN/final_bowtie2_merge_ready.txt ./IDBA_TRAN/CREST_merge.txt ./IDBA_TRAN/BLAST_output_SILVA_with_taxonomy_merge.txt ./IDBA_TRAN/BLAST_output_nt_with_taxonomy_merge.txt ./IDBA_TRAN/MERGE_FILES/

			# Edit contig-* in all the final files
			# CREST files
			sed '1d' ./IDBA_TRAN/CREST_BWA.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-20 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/CREST_BWA_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/CREST_BWA_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/CREST_BWA.txt

			sed '1d' ./IDBA_TRAN/CREST_BOWTIE2.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-20 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence"  > ./IDBA_TRAN/CREST_BOWTIE2_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/CREST_BOWTIE2_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/CREST_BOWTIE2.txt

			# BLAST nt files
			sed '1d' ./IDBA_TRAN/BLAST_nt_BWA.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-31 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_BWA_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_nt_BWA_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_nt_BWA.txt

			sed '1d' ./IDBA_TRAN/BLAST_nt_BOWTIE2.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-31 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_nt_BOWTIE2_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_nt_BOWTIE2_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_nt_BOWTIE2.txt

			# BLAST SILVA files
			sed '1d' ./IDBA_TRAN/BLAST_SILVA_BWA.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-31 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_BWA_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_SILVA_BWA_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_SILVA_BWA.txt

			sed '1d' ./IDBA_TRAN/BLAST_SILVA_BOWTIE2.txt > ./IDBA_TRAN/new.txt
			sed 's/_/\t/1' ./IDBA_TRAN/new.txt > ./IDBA_TRAN/new2.txt
			cut -f2-31 ./IDBA_TRAN/new2.txt > ./IDBA_TRAN/new3.txt
			echo -e "contig_number\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcontig_length\tcontig_read_count\tcounts\tassembly_sequence" > ./IDBA_TRAN/BLAST_SILVA_BOWTIE2_final.txt && cat ./IDBA_TRAN/new3.txt >> ./IDBA_TRAN/BLAST_SILVA_BOWTIE2_final.txt
			rm ./IDBA_TRAN/new*.txt ./IDBA_TRAN/BLAST_SILVA_BOWTIE2.txt

			mkdir ./IDBA_TRAN/FINAL_FILES_IDBA_TRAN
			mv ./IDBA_TRAN/CREST_BWA_final.txt ./IDBA_TRAN/CREST_BOWTIE2_final.txt ./IDBA_TRAN/BLAST_SILVA_BOWTIE2_final.txt ./IDBA_TRAN/BLAST_nt_BOWTIE2_final.txt ./IDBA_TRAN/BLAST_SILVA_BWA_final.txt ./IDBA_TRAN/BLAST_nt_BWA_final.txt ./IDBA_TRAN/FINAL_FILES_IDBA_TRAN/

			echo "final IDBA_TRAN output generated."

##############################TRINITY



		# Prepare CREST file - remove 1st/3rd column, add a header
		assign_NCBI_staxids_to_CREST_v3.py /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_scientific.txt /hdd1/databases/SILVA_database_mar_2020/taxonomy/files_to_make_NCBI_staxids/NCBI_staxids_non_scientific.txt ./TRINITY/CLASSIFICATION/CREST/otus.csv ./TRINITY/CLASSIFICATION/CREST/CREST_output.txt
		sed 's|:|\t|g' ./TRINITY/CLASSIFICATION/CREST/CREST_output.txt > ./TRINITY/CREST_seperated.txt
		cut -f2,4 ./TRINITY/CREST_seperated.txt > ./TRINITY/CREST_header.txt
		sed '1d' ./TRINITY/CREST_header.txt > ./TRINITY/CREST_tax_ready.txt
		assign_taxonomy_NCBI_staxids.sh -b ./TRINITY/CREST_tax_ready.txt -c 2 -e ~/.etetoolkit/taxa.sqlite
		mv CREST_tax_ready_with_taxonomy.txt ./TRINITY/
		sed '1d' ./TRINITY/CREST_tax_ready_with_taxonomy.txt > ./TRINITY/CREST_tax_ready_with_taxonomy_noheader.txt
		echo -e "contig_number\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus" > ./TRINITY/CREST_merge.txt && cat ./TRINITY/CREST_tax_ready_with_taxonomy_noheader.txt >> ./TRINITY/CREST_merge.txt

		rm ./TRINITY/CREST_seperated.txt ./TRINITY/CREST_header.txt ./TRINITY/CREST_tax_ready.txt ./TRINITY/CREST_tax_ready_with_taxonomy.txt ./TRINITY/CREST_tax_ready_with_taxonomy_noheader.txt

		echo "CREST file is ready."

		# Merge files together!!! (BWA with BLAST - nt and SILVA and CREST / BOWTIE2 with BLAST - nt and SILVA and CREST ) SHOULD have 6 output files in total
		merge_mapped_reads_and_contigs.py ./TRINITY/CREST_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/CREST_BWA.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/CREST_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/CREST_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_SILVA_with_taxonomy_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/BLAST_SILVA_BOWTIE2.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_nt_with_taxonomy_merge.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/BLAST_nt_BOWTIE2.txt

		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_SILVA_with_taxonomy_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/BLAST_SILVA_BWA.txt
		merge_mapped_reads_and_contigs.py ./TRINITY/BLAST_output_nt_with_taxonomy_merge.txt ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/BLAST_nt_BWA.txt

		# Move all files - easy to find!!!
		mkdir ./TRINITY/MERGE_FILES
		mv ./TRINITY/final_bwa_merge_ready.txt ./TRINITY/final_bowtie2_merge_ready.txt ./TRINITY/CREST_merge.txt ./TRINITY/BLAST_output_SILVA_with_taxonomy_merge.txt ./TRINITY/BLAST_output_nt_with_taxonomy_merge.txt ./TRINITY/MERGE_FILES/

		# Edit NODE_#_Length_#_Coverage_# in all the final files
		# CREST files
		sed '1d' ./TRINITY/CREST_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/CREST_BWA_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/CREST_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/CREST_BWA.txt

		sed '1d' ./TRINITY/CREST_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/CREST_BOWTIE2_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/CREST_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/CREST_BOWTIE2.txt

		# BLAST nt files
		sed '1d' ./TRINITY/BLAST_nt_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_nt_BWA_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/BLAST_nt_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_nt_BWA.txt

		sed '1d' ./TRINITY/BLAST_nt_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_nt_BOWTIE2_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/BLAST_nt_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_nt_BOWTIE2.txt

		# BLAST SILVA files
		sed '1d' ./TRINITY/BLAST_SILVA_BWA.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_SILVA_BWA_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/BLAST_SILVA_BWA_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_SILVA_BWA.txt

		sed '1d' ./TRINITY/BLAST_SILVA_BOWTIE2.txt > ./TRINITY/new.txt
		sed 's/_/\t/5' ./TRINITY/new.txt > ./TRINITY/new2.txt && sed 's/len=//g' ./TRINITY/new2.txt > ./TRINITY/new3.txt
		sed 's/TRINITY_//g' ./TRINITY/new3.txt > ./TRINITY/new4.txt
		echo -e "contig_number\tcontig_length\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tlowest_rank\tlowest_hit\tsuperkingdom\tkingdom\tphylum\tsubphylum\tclass\tsubclass\torder\tsuborder\tinfraorder\tfamily\tgenus\tcounts\tassembly_sequence" > ./TRINITY/BLAST_SILVA_BOWTIE2_final.txt && cat ./TRINITY/new4.txt >> ./TRINITY/BLAST_SILVA_BOWTIE2_final.txt
		rm ./TRINITY/new*.txt ./TRINITY/BLAST_SILVA_BOWTIE2.txt

		mkdir ./TRINITY/FINAL_FILES_TRINITY
		mv ./TRINITY/CREST_BWA_final.txt ./TRINITY/CREST_BOWTIE2_final.txt ./TRINITY/BLAST_SILVA_BOWTIE2_final.txt ./TRINITY/BLAST_nt_BOWTIE2_final.txt ./TRINITY/BLAST_SILVA_BWA_final.txt ./TRINITY/BLAST_nt_BWA_final.txt ./TRINITY/FINAL_FILES_TRINITY/

		echo "final TRINITY output generated."

	else

		echo "no FINAL merge generated."

	fi

else
	echo "no TRINITY output generated."
fi
