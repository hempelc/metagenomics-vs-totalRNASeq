import itertools

# Manually enter tools for each step:
trimming=["phred<=5", "phred<=10","phred<=15","phred<=20",]
sorting=["sortmerna", "barrnap", "rrnafilter", "unsorted"]
assembly=["spades", "metaspades", "idba-ud", "megahit", "rnaspades", "idba-trans", "trinitiy", "transabyss"]
mapping=["bwa", "bowtie2"]
db=["ncbi_nt", "silva"]
classification=["blast_first_hit", "blast_filtered", "kraken2"]
all=[trimming, sorting, assembly, mapping, db, classification] # list of steps

# Genearte list of tuples of all combinations:
combinations_tuples_list=list(itertools.product(*all))
# Turn that into list of strings:
combinations_str_list = [','.join(tool) for tool in combinations_tuples_list]
# Write each string n list as line in file:
with open("combinations_metagenomics_metatranscriptomics_pipeline.txt", "w") as file:
    for x in combinations_str_list:
        file.write(x + "\n")
