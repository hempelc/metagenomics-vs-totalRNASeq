#!/usr/bin/env python

import itertools
import random

def split(lst, n):
    """Splits list into n equally sized chunks"""
    k, m = divmod(len(lst), n)
    return (lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

def split_list(lst, n):
    """Splits list into chunks of size n"""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# Manually enter tools for each step:
trimming=["5", "10", "15", "20"]
sorting=["sortmerna", "barrnap", "rrnafilter", "unsorted"]
assembly=["spades", "metaspades", "idba_ud", "megahit", "rnaspades", "idba_tran", "trinity", "transabyss"]
mapping=["bwa", "bowtie2"]
db=["ncbi_nt", "silva"]
classification=["blast_first_hit", "blast_filtered", "kraken2"]
all=[trimming, sorting, assembly, mapping, db, classification] # list of steps

# Generate list of tuples of all combinations:
combinations_tuples_list=list(itertools.product(*all))
# Turn that into list of strings:
combinations_str_list = ['-'.join(tool) for tool in combinations_tuples_list]
random.shuffle(combinations_str_list)
chunks=list(split_list(combinations_str_list, 16))

for i in range(len(chunks)):
    out="file_chunk" + str(i+1) + ".txt"
    # Write each element of list as line in file:
    with open(out, "w") as file:
        for x in chunks[i]:
            file.write(x + "\n")

# # Write each string in list as line in file:
# with open("combinations_metagenomics_metatranscriptomics_pipeline_all.txt", "w") as file:
#     for x in combinations_str_list:
#         file.write(x + "\n")

# # For three servers:
# # Split up list in 3 chunks because we run pipeline on 3 servers:
# combinations_chunks=list(split(combinations_str_list, 3))
# combinations_graham=combinations_chunks[0]
# combinations_beluga=combinations_chunks[1]
# combinations_cedar=combinations_chunks[2]

# Write each string in list as line in file for the three servers:
# with open("combinations_metagenomics_metatranscriptomics_pipeline_graham.txt", "w") as file:
#     for x in combinations_graham:
#         file.write(x + "\n")
# with open("combinations_metagenomics_metatranscriptomics_pipeline_beluga.txt", "w") as file:
#     for x in combinations_beluga:
#         file.write(x + "\n")
# with open("combinations_metagenomics_metatranscriptomics_pipeline_cedar.txt", "w") as file:
#     for x in combinations_cedar:
#         file.write(x + "\n")
