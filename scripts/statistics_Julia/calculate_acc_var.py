#!/usr/bin/env python3

import pandas as pd
import numpy as np # Not needed for Chris' code portion
import glob
import os

# NOTE: When running the code, ignore the warning - stems from a function that I had to implement because some of the files were wrong ("NA" in counts column), but I fixed that in the pipeline,
# so after I rerun the pipelines, we can take out that portion of code and the warning will disappear

# @Julia, if you want to have a look at one master table, run the following line once you ran the rest of the code (commented out for now):
#master_dfs["M4_RNA"]
# The only thing you have to adjust is the full path to your folder containing the 3 dirs with files (full path from root)
# Note that the first column in the master dfs is always the expected and is type float while all other columns are type int (shouldn't make a difference for the chi square and variance
# calculation I guess)

# FUTURE ADAPTATION: "lowest_hit" needs to be adapted to "species" once we rerun pipelines (can simply literatelly be find-replace'd)

# TO DO (LONG-TERM): making "mega master table" that aggregates on all ranks for each sample, with one column indicating on which rank the row has been aggregated


# Parameters
workdir = "/Users/christopherhempel/Desktop/mock_community_RNA/" # Full path to directory that contains samples, HAS TO END WITH "/"
samples = ["M4_RNA", "M5_RNA", "M6_RNA"] # Samples to include into the analysis (must equal names of directories in workdir that contain each sample's pipeline results)
groupby_rank = "lowest_hit" # Basis for rank to group rows on. Either based on genus (option "genus") or on species (option "lowest_hit" (NOTE later "species"))
rel_abun_basis = "cell" # Basis for relative abundance calculation of expected mock community taxa abundance. Either based on genomic DNA (option "gen") or on cell number (option "cell")

# TO DO: parameter check to make sure that parameters are one of the allowed options and that variable samples is not empty


# Hardcoded variables
num_reads_abs = {"M4_DNA": 310806, "M4_RNA": 211414, "M5_DNA": 146174, "M5_RNA": 322631, "M6_DNA": 114459, "M6_RNA": 269408} # dict with read numbers per sample for easy accessibility
rel_abun_gen = [89.1, 8.9, 0.89, 0.89, 0.089, 0.089, 0.0089, 0.00089, 0.00089, 0.000089] # Relative abundances of mock community taxa in percent based on genomic DNA
rel_abun_cell = [94.9, 4.2, 0.7, 0.12, 0.058, 0.059, 0.015, 0.001, 0.00007, 0.0001] # Relative abundances of mock community taxa in percent based on cell number
master_dfs = {} # Empty dic that will eventually contain all sample's master tables

# Run the code on every sample separately:
for sample in samples:

    # 1. Make expected mock community df

    ## To calculate the absolute expected read count for the taxa in the mock community, we need the hardcoded taxonomic information of each taxon. We save it in a dictionary:
    expected_dic = {}
    expected_dic["L_monocytogenes"] = ["Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Listeriaceae", "Listeria", "Listeria monocytogenes"]
    expected_dic["P_aeruginosa"] = ["Bacteria", "Proteobacteria", "Gammaproteobacteria", "Pseudomonadales", "Pseudomonadaceae", "Pseudomonas", "Pseudomonas aeruginosa"]
    expected_dic["B_subtilis"] = ["Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Bacilliaceae", "Bacillus", "Bacillus subtilis"]
    expected_dic["S_cerevisiae"] = ["Eukaryota", "Ascomycota", "Saccharomycetes" ,"Saccharomycetales", "Saccharomycetaceae", "Saccharomyces", "Saccharomyces cerevisiae"]
    expected_dic["E_coli"] = ["Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia", "Escherichia coli"]
    expected_dic["S_enterica"] = ["Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Salmonella", "Salmonella enterica"]
    expected_dic["L_fermentum"] = ["Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Limosilactobacillus", "Lactobacillus fermentum"]
    expected_dic["E_faecalis"] = ["Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Enterococcaceae", "Enterococcus", "Enterococcus faecalis"]
    expected_dic["C_neoformans"] = ["Eukaryota", "Basidiomycota", "Tremellomycetes", "Tremellales", "Tremellaceae", "Cryptococcus", "Cryptococcus neoformans"]
    expected_dic["S_aureus"] = ["Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus", "Staphylococcus aureus"]


    ## We also need absolute abundances per taxon. Therefore, we multiply the # of reads in sample with the relative abundance of each taxon in the mock community
    ## We calculate the relative read number for each sample depending on the absoulte read number:
    if rel_abun_basis == "gen":
        num_reads_taxa = [rel_abun / 100 * num_reads_abs[sample] for rel_abun in rel_abun_gen]
    elif rel_abun_basis == "cell":
        num_reads_taxa = [rel_abun / 100 * num_reads_abs[sample] for rel_abun in rel_abun_cell]

    ## We add each taxon's read count to its taxonomy list
    rep=0
    for taxon in expected_dic:
        expected_dic[taxon].append(num_reads_taxa[rep])
        rep+=1

    ## We finally read in the expected mock community df as pandas df
    expected_df=pd.DataFrame.from_dict(expected_dic, orient="index", columns=["superkingdom", "phylum", "class", "order", "family", "genus", "lowest_hit", "counts"])


    # 2. Read in all files and make dataframes for each

    ## Make list for all file names
    sample_files = glob.glob("{workdir}{sample}/*.txt".format(workdir=workdir, sample=sample))

    ## Make dict and set first entry to expected community
    sample_dfs = {"expected": expected_df}
    all_taxa = []
    ## For each file
    for file in sample_files:
        ### Read in file as pandas df and fill NaN with "NA"
        df = pd.read_table(file).fillna("NA")
        ### NA check, if NA or NANA in counts, don't do anything, NOTE: note needed later, just temporarily because some pipeline results are wrong      ### TO BE DELETED
        if df.loc[df['counts'] == "NA"].empty and df.loc[df['counts'] == "NANA"].empty:                                                                  ### TO BE DELETED
            ### Cut dfs down to relevant columns
            if groupby_rank == "lowest_hit":
                df_small = df[["superkingdom", "phylum", "class", "order", "family", "genus", "lowest_hit", "counts"]]
            elif groupby_rank == "genus":
                df_small = df[["superkingdom", "phylum", "class", "order", "family", "genus", "counts"]]
            df_small["counts"] =  df_small["counts"].astype(int)
            ### Group similar taxonomy hits and sum counts
            df_agg = df_small.groupby(list(df_small.columns)[:-1]).sum().reset_index()
            ### Add all species to list
            all_taxa.extend(df_agg[groupby_rank].tolist())
            ### Make dict key with keys = filename (without full path and extension) and value = edited df
            sample_dfs["{file}".format(file=file).split("/")[-1].split(".")[-2]] = df_agg
        else:                                                                                                                                            ### TO BE DELETED
            continue                                                                                                                                     ### TO BE DELETED


    # 3. Generate master table

    ## Drop duplicates of species list
    unique_taxa = list(set(all_taxa))

    ## Make master df with taxa from taxa list as row names
    master_dfs[sample] = pd.DataFrame(index=pd.Index(unique_taxa))

    ## For each df in dict with key=pipeline
    for pipeline, df in sample_dfs.items():
        counts=[]
        ### For each taxon in unique taxa list
        for taxon in unique_taxa:
            #### If taxon is in df groupby_rank column
            if (df[groupby_rank] == taxon).any():
                ##### Sum up all counts of that taxon and add it to list "counts"
                counts.append(df.loc[df[groupby_rank] == taxon, 'counts'].sum())
            #### If taxon not in df groupby_rank column, add 0 to list "counts"
            else:
                counts.append(0)
        ### Make a new column in  master table named after the pipeline and add taxon counts for that pipeline
        master_dfs[sample][pipeline]=counts

# TO DO (HIGHEST PRIORITY): CHISQUARE TEST

# TO DO (HIGHEST PRIORITY): VARIATION (note to delete expected column out of master table)

# TO DO: PLOT ACCURACY (FOR NOW AVERAGE BETWEEN REPLICATES) VS. VARIATION FOR ALL PIPELINES


################### JULIA'S CODE
#
#
# taxa = pd.DataFrame(data = [])
#
# for file in sample_files:
#     df = pd.read_table(i, header = 0)
#     taxa = pd.concat([taxa, df['lowest_hit']])
#
# unique = taxa.drop_duplicates()
#
# test = sample_files[0]
# df = pd.read_table(test, header = 0)
#
# unique_list = unique[0].tolist()
# unique_list = [x for x in unique_list if isinstance(x, str)]
# unique_list = [str(i) for i in unique_list]
#
# for i in range(len(unique_list)):
#     unique_list[i]=' '.join(unique_list[i].split()[:2])
#     if len(unique_list[i].split())==2:
#         if unique_list[i][0].isupper():
#             continue
#         else:
#             unique_list[i]="NA"
#     else:
#         unique_list[i]="NA"
#
# unique_list = list(dict.fromkeys(unique_list))
#
#
#
# #print(unique)
# #for i in range(0,10):
# #    print (unique.iloc[i])
# #    print (type(unique.iloc[i]))
#     #(i)
#
# yes = df['lowest_hit'].isin(unique.iloc[0])
# #print(yes)
#
#
#
#
# # unique_genus = file1['genus'].drop_duplicates()
# # print(unique_genus.iloc[0])
# # print(unique_genus.iloc[1])
#
# # for i in unique_genus:
# #     file1
#
#
# # #Coulmns will be expected
#
# # master = pd.DataFrame(data = 0, index = , columns =)
#
# #Populate expected column
#
# #Populate other colunms
#
# #Perform chisquare and save
#
# #Calculate variance
#
# #Plot two values
