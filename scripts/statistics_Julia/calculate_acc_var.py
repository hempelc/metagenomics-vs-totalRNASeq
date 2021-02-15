#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import csv
from scipy.stats import chisquare
from scipy.stats import f_oneway
# NOTE: When running the code, ignore the warning - stems from a function that I had to implement because some of the files were wrong ("NA" in counts column), but I fixed that in the pipeline,
# so after I rerun the pipelines, we can take out that portion of code and the warning will disappear

# @Julia, if you want to have a look at one master df, run the following line once you ran the rest of the code (commented out for now):
#master_dfs["M4_RNA"]
# The only thing you have to adjust is the full path to your folder containing the 3 dirs with files (full path from root)

# FUTURE ADAPTATION: "lowest_hit" needs to be adapted to "species" once we rerun pipelines (can simply literatelly be find-replace'd)

# TO DO (LONG-TERM): making "mega master df" that aggregates on all ranks for each sample, with one column indicating on which rank the row has been aggregated


# Parameters
workdir = "/Users/christopherhempel/Desktop/mock_community_RNA/" # Full path to directory that contains samples, HAS TO END WITH "/"
#workdir = "/Users/julia/Documents/Projects/MicrobeCommunities/Abstract/mock_community_RNA/"
savedir = "/Users/christopherhempel/Desktop/" # Full path to directory where plots should be saved, HAS TO END WITH "/"
samples = ["M4_RNA", "M5_RNA", "M6_RNA"] # 3 replicate samples to include into the analysis (must equal names of directories in workdir that contain each sample's pipeline results)
groupby_rank = "lowest_hit" # Basis for taxa rank to group rows on. Either based on genus (option "genus") or on species (option "lowest_hit" (NOTE later "species"))
rel_abun_basis = "cell" # Basis for relative abundance calculation of expected mock community taxa abundance. Either based on genomic DNA (option "gen") or on cell number (option "cell")

# MAYBE TO DO: implement that script can be run form commadn line and insert parameter check to make sure that
# parameters are one of the allowed options and that variable "samples" is not empty


# Hardcoded variables
num_reads_abs = {"M4_DNA": 310806, "M4_RNA": 211414, "M5_DNA": 146174, "M5_RNA": 322631, "M6_DNA": 114459, "M6_RNA": 269408} # dict with read numbers per sample for easy accessibility
rel_abun_gen = [89.1, 8.9, 0.89, 0.89, 0.089, 0.089, 0.0089, 0.00089, 0.00089, 0.000089] # Relative abundances of mock community taxa in percent based on genomic DNA
rel_abun_cell = [94.9, 4.2, 0.7, 0.12, 0.058, 0.059, 0.015, 0.001, 0.00007, 0.0001] # Relative abundances of mock community taxa in percent based on cell number
master_dfs = {} # Empty dic that will eventually contain all samples' master dfs
all_taxa = [] # Empty list that will eventually contain all taxa that appear in all samples
chi2_var = {} # Empty dic that will eventually contain all chi-squares statistics and ANOVA variance F-statistics for all pipelines

# 1 Read in all pipeline results for every sample and the expected mock communty as df;
#   add all taxa from all dfs to "all_taxa" list, needed to generate master dfs that contain counts of all taxa from all samples and mock community:
for sample in samples:

    ## 1.1 Make expected mock community df

    ### To calculate the absolute expected read count for the taxa in the mock community, we need the hardcoded taxonomic information of each taxon. We save it in a dictionary:
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

    ### We also need absolute abundances per taxon. Therefore, we multiply the absoulte rread number in the sample with the relative abundance of each taxon in the mock community.
    ### We calculate the absolute read number for each taxon based on the absoulte read number of the sample:
    if rel_abun_basis == "gen":
        num_reads_taxa = [rel_abun / 100 * num_reads_abs[sample] for rel_abun in rel_abun_gen]
    elif rel_abun_basis == "cell":
        num_reads_taxa = [rel_abun / 100 * num_reads_abs[sample] for rel_abun in rel_abun_cell]

    ### We add each taxon's read count to its taxonomy list:
    rep=0
    for taxon in expected_dic:
        expected_dic[taxon].append(num_reads_taxa[rep])
        rep+=1

    ### And finally read in the expected mock community dic as pandas df
    expected_df=pd.DataFrame.from_dict(expected_dic, orient="index", columns=["superkingdom", "phylum", "class", "order", "family", "genus", "lowest_hit", "counts"])


    ## 1.2 Read in all files from sample and make dataframes for each pipeline

    ### Make list for all file names:
    sample_files = glob.glob("{workdir}{sample}/*.txt".format(workdir=workdir, sample=sample))

    ### Make dic that will eventually contain all dfs and set first entry to expected community:
    sample_dfs = {"expected": expected_df}
    ### For each file
    for file in sample_files:
        #### Read in file as pandas df and fill NaN with "NA"
        df = pd.read_table(file).fillna("NA")
        #### NA check, if NA or NANA in counts, don't do anything, NOTE: note needed later, just temporarily because some pipeline results are wrong      ### TO BE DELETED
        if df.loc[df['counts'] == "NA"].empty and df.loc[df['counts'] == "NANA"].empty:                                                                  ### TO BE DELETED
            #### Cut dfs down to relevant columns
            if groupby_rank == "lowest_hit":
                df_small = df[["superkingdom", "phylum", "class", "order", "family", "genus", "lowest_hit", "counts"]]
            elif groupby_rank == "genus":
                df_small = df[["superkingdom", "phylum", "class", "order", "family", "genus", "counts"]]
            #### Transform all entries in counts to type int
            #df_small["counts"] =  df_small["counts"].astype(int)
            #### Group similar taxonomy hits and sum counts
            df_agg = df_small.groupby(list(df_small.columns)[:-1]).sum().reset_index()
            #### Add all taxa to list "all taxa"
            all_taxa.extend(df_agg[groupby_rank].tolist())
            #### Make dict key with keys = filename (without full path and extension) and value = edited df
            sample_dfs["{file}".format(file=file).split("/")[-1].split(".")[-2].split("trimmed_at_phred_")[1].split("_pipeline_final")[0]] = df_agg
        else:                                                                                                                                            ### TO BE DELETED
            continue                                                                                                                                     ### TO BE DELETED

# 2 Make a "unique_taxa" list that contains all taxa that appear in all samples;
#   generate a "master_df" df from each sample that contains counts for all taxa in all samples and mock community:

## 2.1 Drop duplicates of "all_taxa" list
unique_taxa = list(set(all_taxa))

## 2.2 Generate master df for every sample and save in dic "master_dfs":
for sample in samples:
    ### Make master df with taxa from "unique_taxa" list as row names
    master_dfs[sample] = pd.DataFrame(index=pd.Index(unique_taxa))

    ### For each df in dict with key=pipeline
    for pipeline, df in sample_dfs.items():
        counts=[]
        #### For each taxon in unique taxa list
        for taxon in unique_taxa:
            ##### If taxon is in df groupby_rank column
            if (df[groupby_rank] == taxon).any():
                ###### Sum up all counts of that taxon and add it to list "counts"
                counts.append(df.loc[df[groupby_rank] == taxon, 'counts'].sum())
            ##### If taxon not in df groupby_rank column, add 0 to list "counts"
            else:
                counts.append(0)
        #### Make a new column in  master df named after the pipeline and add taxon counts for that pipeline
        master_dfs[sample][pipeline]=counts


# 3 Peform tests to determine accuracy and precision

## 3.1 Chi-Squared test (accuracy)
### 3.1.1 To perform a Chi-Squared test on replicates, we summarize columns across the replicates:
master_df_summarized = pd.DataFrame()
### For all columns (note: all samples contain the same column names, so we manually pick one (the first in "samples" list))
for col in master_dfs[samples[0]].columns.tolist():
    #### Make list per column in every sample and save in dir:
    col_dic_chi2 = {}
    for sample in samples:
        col_dic_chi2[sample] = master_dfs[sample][col].tolist()
    #### Summarize the three columns in the three replicates:
    master_df_summarized[col] = [a + b + c for a, b, c in zip(col_dic_chi2[samples[0]], col_dic_chi2[samples[1]], col_dic_chi2[samples[2]])]

### 3.1.2 Now we perform a Chi-Squared test on the summarized master df, for each pipeline against the expected composition
###     (note, this is the DIRTY version where we add the same amount (1) to each value to get rid of zeros in the expected columns, because otherwise the Chi-Squared test is not possible):
for pipeline in master_df_summarized.columns.tolist()[1:]:
    chi2_var[pipeline] = [(chisquare([abun + 1 for abun in master_df_summarized[pipeline].tolist()],[abun + 1 for abun in master_df_summarized['expected'].tolist()])[0])]

## 3.2 One-way ANOVA (precision)
## Finally, we calculate a one-way ANOVA for all pipelines across the three replicates and add their F-statistics to the pipelines in the "chi2_var" dic
for col in master_dfs[samples[0]].columns.tolist()[1:]: # all samples contain the name column names, so we manually pick one (the first in "samples" list)
    col_dic_anova = {}
    for sample in samples:
        col_dic_anova[sample] = master_dfs[sample][col].tolist()
    chi2_var[col].append(f_oneway(col_dic_anova[samples[0]], col_dic_anova[samples[1]], col_dic_anova[samples[2]])[0])

## 3.3 Assign which tools have been used each step in the pipeline, with one column per step
for pipeline in chi2_var:
    pipeline_replace = pipeline.replace("IDBA_", "IDBA-").replace("NCBI_NT", "NCBI-NT").replace("BLAST_FIRST_HIT", "BLAST-FIRST-HIT").replace("BLAST_FILTERED", "BLAST-FILTERED")
    chi2_var[pipeline].extend(pipeline_replace.split("_"))

## 3.4 Save "chi2_var" dic as csv so that it can be plotted using R:
df_save = pd.DataFrame.from_dict(chi2_var, orient="index", columns=["chi-square statistics", "ANOVA F-statistics", "trimmed PHRED", "rRNA filter", "assembler", "mapper", "DB", "classification"])
df_save.to_csv("{savedir}chi2_var.csv".format(savedir=savedir), index_label="pipeline")


# # 4 Plot accuracy vs. precision (NOTE: switched to R for that part)
#
# ## 4.1 Make plot data out of "chi2_var" dic that contains coordinates and pipeline names
# plot_data = {"x":[], "y":[], "pipeline":[]}
# for pipeline, coord in chi2_var.items():
#     plot_data["x"].append(coord[0])
#     plot_data["y"].append(coord[1])
#     plot_data["pipeline"].append(pipeline)
#
# ## 4.2 Display the plot
# plt.figure(figsize=(10,8))
# plt.title('Scatter Plot', fontsize=20)
# plt.xlabel('Chi-Square', fontsize=15)
# plt.ylabel('ANOVA', fontsize=15)
# #plt.xlim(0, 0.1e12)
# #plt.ylim(min(range), max(range))
# #plt.ylim(-0.2e-30, 0.2e-30)
# plt.axis([0, 4e10, -1e-31, 1e-31])
# plt.scatter(plot_data["x"], plot_data["y"], marker = 'o')
# ### Add labels
# for pipeline, x, y in zip(plot_data["pipeline"], plot_data["x"], plot_data["y"]):
#     plt.annotate(pipeline, xy = (x, y), size="xx-small", va="bottom", ha="center", stretch="condensed")
# ### Save plot as png
# plt.savefig("{savedir}chi2_var.png".format(savedir=savedir), dpi=600)
#
#
#
# ##### test area for fisher test
#
# import numpy as np
# import rpy2.robjects.numpy2ri
# from rpy2.robjects.packages import importr
# rpy2.robjects.numpy2ri.activate()
# test_pip = master_dfs[sample]['trimmed_at_phred_5_UNSORTED_SPADES_BWA_SILVA_KRAKEN2_pipeline_final'].tolist()
# test_exp = master_dfs[sample]['expected'].tolist()
# np_array = np.column_stack([test_exp, test_pip])
#
# stats = importr('stats')
# res = stats.fisher_test(np_array, simulate_p_value="TRUE")
# print(res)
#
#
# # Fisher test loop
# for i in master_dfs[sample].columns.tolist()[1:]:
#     np_array = np.column_stack([master_dfs[sample]['expected'].tolist(), master_dfs[sample][i].tolist()])
#     res = stats.fisher_test(np_array, simulate_p_value="TRUE")
#     print(res[0][0])
