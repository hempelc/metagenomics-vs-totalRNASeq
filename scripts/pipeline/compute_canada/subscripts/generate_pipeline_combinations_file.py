#!/usr/bin/env python

import itertools

def split(a, n):
    """Splits list a into n equally sized chunks"""
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

# Manually enter tools for each step:
trimming=["5", "10", "15", "20"]
sorting=["sortmerna", "barrnap", "rrnafilter", "unsorted"]
assembly=["spades", "metaspades", "idba_ud", "megahit", "rnaspades", "idba_tran", "trinity", "transabyss"]
mapping=["bwa", "bowtie2"]
db=["ncbi_nt", "silva"]
classification=["blast_first_hit", "blast_filtered", "kraken2"]
all=[trimming, sorting, assembly, mapping, db, classification] # list of steps

# Genearte list of tuples of all combinations:
combinations_tuples_list=list(itertools.product(*all))
# Turn that into list of strings:
combinations_str_list = ['-'.join(tool) for tool in combinations_tuples_list]
# Split up list in 3 chunks because we run pipelien on 3 servers:
combinations_chunks=list(split(combinations_str_list, 3))
combinations_graham=combinations_chunks[0]
combinations_beluga=combinations_chunks[1]
combinations_cedar=combinations_chunks[2]

# Write each string in list as line in file for the three servers:
with open("combinations_metagenomics_metatranscriptomics_pipeline_graham.txt", "w") as file:
    for x in combinations_graham:
        file.write(x + "\n")
with open("combinations_metagenomics_metatranscriptomics_pipeline_beluga.txt", "w") as file:
    for x in combinations_beluga:
        file.write(x + "\n")
with open("combinations_metagenomics_metatranscriptomics_pipeline_cedar.txt", "w") as file:
    for x in combinations_cedar:
        file.write(x + "\n")
