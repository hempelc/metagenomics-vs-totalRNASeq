#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import os
import csv
import math
import statistics
import copy
from scipy.stats import chisquare
from scipy.stats import combine_pvalues
# NOTE: When running the code, ignore the warning - stems from a function that I had to implement because some of the files were wrong ("NA" in counts column), but I fixed that in the pipeline,
# so after I rerun the pipelines, we can take out that portion of code and the warning will disappear

# @Julia, if you want to have a look at one master df, run the following line once you ran the rest of the code (commented out for now):
#master_dfs_uniq_abs["M4_RNA"]
# The only thing you have to adjust is the full path to your folder containing the 3 dirs with files (full path from root)

# FUTURE ADAPTATION: "lowest_hit" needs to be adapted to "species" once we rerun pipelines (can simply literatelly be find-replace'd)

# TO DO (LONG-TERM): making "mega master df" that aggregates on all ranks for each sample, with one column indicating on which rank the row has been aggregated

# Functions
## To calculate the distance between two points:
def point_dist (p1, p2):
    dist = math.sqrt( ((p1[0]-p2[0])**2)+((p1[1]-p2[1])**2) )
    return dist

## To normalize value list[x] to range 0-1 based on min/max of list:
def normalize (x,min_x,max_x):
    normalized = (x-min_x)/(max_x-min_x)
    return normalized

# Parameters
workdir = "/Users/christopherhempel/Desktop/mock_community_RNA/" # Full path to directory that contains samples, HAS TO END WITH "/"
#workdir = "/Users/julia/Documents/Projects/MicrobeCommunities/Abstract/mock_community_RNA/"
savedir = "/Users/christopherhempel/Desktop/jupyter_notebook_plots/" # Full path to directory where plots should be saved, HAS TO END WITH "/"
samples = ["M4_RNA", "M5_RNA", "M6_RNA"] # 3 replicate samples to include into the analysis (must equal names of directories in workdir that contain each sample's pipeline results)
groupby_rank = "lowest_hit" # Basis for taxa rank to group rows on. Either based on genus (option "genus") or on species (option "lowest_hit" (NOTE later "species"))
rel_abun_basis = "cell" # Basis for relative abundance calculation of expected mock community taxa abundance. Either based on genomic DNA (option "gen") or on cell number (option "cell")

# MAYBE TO DO: implement that script can be run form command line and insert parameter check to make sure that
# parameters are one of the allowed options and that variable "samples" is not empty


# Hardcoded variables
num_reads_abs = {"M4_DNA": 310806, "M4_RNA": 211414, "M5_DNA": 146174, "M5_RNA": 322631, "M6_DNA": 114459, "M6_RNA": 269408} # dict with read numbers per sample for easy accessibility
rel_abun_gen = [89.1, 8.9, 0.89, 0.89, 0.089, 0.089, 0.0089, 0.00089, 0.00089, 0.000089] # Relative abundances of mock community taxa in percent based on genomic DNA
rel_abun_cell = [94.9, 4.2, 0.7, 0.12, 0.058, 0.059, 0.015, 0.001, 0.00007, 0.0001] # Relative abundances of mock community taxa in percent based on cell number
master_dfs_raw = {} # Empty dic that will eventually contain all samples' raw pipelines output
master_dfs_uniq_abs = {} # Empty dic that will eventually contain all samples' master dfs with absolute counts (dfs with rows = all unique taxa)
master_dfs_uniq_rel = {} # Empty dic that will eventually contain all samples' master dfs with relative counts (dfs with rows = all unique taxa)
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
    ### Save sample_df in dic "master_dfs_raw":
    master_dfs_raw[sample] = sample_dfs


# 2 Make a "unique_taxa" list that contains all taxa that appear in all samples;
#   generate a "master_df" df from each sample that contains counts for all taxa in all samples and mock community:

## 2.1 Drop duplicates of "all_taxa" list
unique_taxa = list(set(all_taxa))


## 2.2 Generate master df with absolute read counts for every sample and save in dic "master_dfs_uniq_abs":
for sample, pipeline_dfs in master_dfs_raw.items():
    ### Make master df with taxa from "unique_taxa" list as row names
    master_dfs_uniq_abs[sample] = pd.DataFrame(index=pd.Index(unique_taxa))

    ### For each df in dict with key=pipeline
    for pipeline, data in pipeline_dfs.items():
        counts=[]
        #### For each taxon in unique taxa list
        for taxon in unique_taxa:
            ##### If taxon is in df groupby_rank column
            if (data[groupby_rank] == taxon).any():
                ###### Sum up all counts of that taxon and add it to list "counts"
                counts.append(data.loc[data[groupby_rank] == taxon, 'counts'].sum())
            ##### If taxon not in df groupby_rank column, add 0 to list "counts"
            else:
                counts.append(0)
        #### Make a new column in  master df named after the pipeline and add taxon counts for that pipeline
        master_dfs_uniq_abs[sample][pipeline]=counts

## 2.3 Generate master df with relative read counts for every sample and save in dic "master_dfs_uniq_rel":
### Based on absolute counts, divide every column entry by the column sum for all columns
for key, value in master_dfs_uniq_abs.items():
    master_dfs_uniq_rel[key] = master_dfs_uniq_abs[key]/master_dfs_uniq_abs[key][master_dfs_uniq_abs[key].columns].sum()

## 2.3 Generate master df with presence/absence for every sample and save in dic "master_dfs_uniq_pa":
### Deepcopy one of the master dfs
master_dfs_uniq_pa = copy.deepcopy(master_dfs_uniq_abs)
### Replace all non-zero values with 1
for key, value in master_dfs_uniq_pa.items():
    value[value != 0] = 1


# 3 Calculate confusion matrix

## 3.1 Cut down master_dfs to only expected taxa
### Make function to cut down master_dfs to only expected taxa:
def cutdown (master_df, type):
    cutdown_dic = {}
    for sample, pipelines in master_df.items():
        if type == "above_zero":
            cutdown_dic[sample] = pipelines.loc[pipelines['expected'] != 0]
        elif type == "equal_zero":
            cutdown_dic[sample] = pipelines.loc[pipelines['expected'] == 0]
    return cutdown_dic

### Apply cutdown function on absolute, relative, and pa master_dfs
abs_cutdown_expected = cutdown(master_dfs_uniq_abs, "above_zero")
rel_cutdown_expected = cutdown(master_dfs_uniq_rel, "above_zero")
pa_cutdown_expected = cutdown(master_dfs_uniq_pa, "above_zero")

abs_cutdown_FP = cutdown(master_dfs_uniq_abs, "equal_zero")
rel_cutdown_FP = cutdown(master_dfs_uniq_rel, "equal_zero")
pa_cutdown_FP = cutdown(master_dfs_uniq_pa, "equal_zero")

## 3.2 Make confusion matrices for all samples for all 3 master_dfs
def confusion_calc (cutdown_dic_expected, cutdown_dic_FP, type):
    confusion_master={}
    for sample, pipelines in cutdown_dic_expected.items():
        expected_list=pipelines['expected'].tolist()
        confusion_dic = {}
        for pipeline, abundances in pipelines.iloc[:, 1:].iteritems():
            FN=0
            FP=0
            TP=0
            for i in range(len(abundances.tolist())):
                if abundances.tolist()[i] - expected_list[i] < 0:
                    FN += (abundances.tolist()[i] - expected_list[i])*-1
                if abundances.tolist()[i] - expected_list[i] > 0:
                    FP += (abundances.tolist()[i] - expected_list[i])
                if abundances.tolist()[i] >= expected_list[i]:
                    TP += expected_list[i]
                if abundances.tolist()[i] < expected_list[i]:
                    TP += abundances.tolist()[i]
            if type == "abs" or type == "rel":
                confusion_values = {"subseed_reads": FN, "exceed_reads": FP, "true_reads":TP}
                confusion_values["false_reads"]=cutdown_dic_FP[sample][pipeline].sum()
            elif type == "pa":
                confusion_values = {"FN": FN, "TP":TP}
                confusion_values["FP"]=cutdown_dic_FP[sample][pipeline].sum()
            confusion_dic[pipeline] = confusion_values
        confusion_master[sample] = confusion_dic
    return confusion_master

### Apply confusion_calc function on absolute, relative, and pa cutdown master_dfs
abs_params_reps = confusion_calc(abs_cutdown_expected, abs_cutdown_FP, "abs")
rel_params_reps = confusion_calc(rel_cutdown_expected, rel_cutdown_FP, "rel")
pa_params_reps = confusion_calc(pa_cutdown_expected, pa_cutdown_FP, "pa")


## 3.3 Calculate the average between replicates:
### Make dics containing each pipeline as key with one list per parameter for each pipeline:
abs_params={}
rel_params={}
pa_params={}

for params_dic in [abs_params, rel_params, pa_params]:
    for sample in ["M4_RNA", "M5_RNA", "M6_RNA"]:                                                                                                      ### TO BE DELETED since all samples will contain the same pipelines eventually
        for pipeline in abs_params_reps[sample].keys():
            if params_dic==abs_params or params_dic==rel_params:
                params_dic[pipeline]={"subseed_reads": [], "exceed_reads": [], "true_reads": [], "false_reads": []}
            elif params_dic==pa_params:
                params_dic[pipeline]={"FN": [], "FP": [], "TP": []}

### Fill in dics created above so that every pipeline contains a dic with
### parameters as keys and each parameter lists its three values across the
### three replicates:
for param_dic_rep in [abs_params_reps, rel_params_reps, pa_params_reps]:
    for sample, pipelines in param_dic_rep.items():
        for pipeline, params in pipelines.items():
            for param, value in params.items():
                if param_dic_rep==abs_params_reps:
                    abs_params[pipeline][param].append(value)
                elif param_dic_rep==rel_params_reps:
                    rel_params[pipeline][param].append(value)
                elif param_dic_rep==pa_params_reps:
                    pa_params[pipeline][param].append(value)

### Calculate the average across the three replicates:
for params_dic in [abs_params, rel_params, pa_params]:
    for pipeline, params in params_dic.items():
        for param, list in params.items():
            params_dic[pipeline][param]=sum(list)/len(list)


# 4 Calculate variance for absolute and realtive data and mismatches between replicates for p/a data


## 4.1  We calculate the variances for each taxon for all pipelines across the three replicates, sum up the variance of all taxa across the
##      replicates for every pipeline, and add the returned summed variance to each pipeline in the respective params dics.
##      NOTE: We only do that for absolute and relative data as it is unapplicable for p/a data.

## NOTE: For now, some pipelines didn't work, so we have to make a list of shared pipelines                                                            ### TO BE DELETED
shared_pipelines_dupl = []                                                                                                                                   ### TO BE DELETED
for sample in samples:                                                                                                                                  ### TO BE DELETED
    shared_pipelines_dupl.extend(master_dfs_uniq_abs[sample].columns.tolist()[1:])                                                                               ### TO BE DELETED
shared_pipelines = []                                                                                                                                  ### TO BE DELETED
for i in shared_pipelines_dupl:                                                                                                                                  ### TO BE DELETED
    if shared_pipelines_dupl.count(i) > 2:                                                                                                                               ### TO BE DELETED
        if i not in shared_pipelines:                                                                                                                                  ### TO BE DELETED
            shared_pipelines.append(i)                                                                                                                                  ### TO BE DELETED

 for params_dic in [abs_params, rel_params]:
    for pipeline in shared_pipelines:                                                                                                                            ### TO BE DELETED
        if params_dic==abs_params:
            master_dfs_uniq=master_dfs_uniq_abs
        elif params_dic==rel_params:
            master_dfs_uniq=master_dfs_uniq_rel
        ### Each pipeline gets a variance sum variable
        var_sum = 0
        pipeline_col = master_dfs_uniq[samples[0]].columns.tolist()[1:].index(pipeline)
        ### For each taxon in the pipelines across replicates, calculate the variance and sum up all variances across all taxa:
        for taxon in range(0,master_dfs_uniq[samples[0]].shape[0]):
            var_sum += statistics.variance([int(master_dfs_uniq[samples[0]].iloc[taxon,pipeline_col]), int(master_dfs_uniq[samples[1]].iloc[taxon,pipeline_col]), int(master_dfs_uniq[samples[2]].iloc[taxon,pipeline_col])])
        ### Add var_sum to chi2_var dic for respective pipeline
        params_dic[pipeline]["variance"]=var_sum

## 4.2 We calculate the mismatches between replicates for each taxon for all pipelines, sum up the mismatches
## for every pipeline, and add the mismatches to each pipeline in the respective params dics.
for pipeline in shared_pipelines:                                                                                                                            ### TO BE DELETED
    ### Each pipeline gets a mismatches variable to count mismatches between replicates
    mismatches = 0
    pipeline_col = master_dfs_uniq_pa[samples[0]].columns.tolist()[1:].index(pipeline)
    ### For each taxon in the pipelines across replicates, count how often 1 (=present) was detected
    ### - if it was detected only 1 one 2 times, pipelines mismatch in their outcome, which we collect in variable mismatches:
    for taxon in range(0,master_dfs_uniq_pa[samples[0]].shape[0]):
        count=[master_dfs_uniq_pa[samples[0]].iloc[taxon,pipeline_col], master_dfs_uniq_pa[samples[1]].iloc[taxon,pipeline_col], master_dfs_uniq_pa[samples[2]].iloc[taxon,pipeline_col]].count(1)
        if count==1 or count==2:
            mismatches+=1
    pa_params[pipeline]["mismatches"]=mismatches

# 5 Make master params df
master_param_df={}
for pipeline in shared_pipelines:                                                                                                                ### TO BE DELETED
    master_param_df[pipeline]={}
    for params_dic in [abs_params, rel_params, pa_params]:
        for param in params_dic[pipeline].keys():
            master_param_df[pipeline][param]=params_dic[pipeline][param]

# TO DO: to df and then PCA

# # This is for all 3 data types, which is inappropriate:
# for params_dic in [abs_params, rel_params, pa_params]:
#     for pipeline in shared_pipelines:                                                                                                                       ### TO BE DELETED
#         if params_dic==abs_params:
#             master_dfs_uniq=master_dfs_uniq_abs
#         elif params_dic==rel_params:
#             master_dfs_uniq=master_dfs_uniq_rel
#         elif params_dic==pa_params:
#             master_dfs_uniq=master_dfs_uniq_pa
#         col = master_dfs_uniq[samples[0]].columns.tolist()[1:].index(pipeline)
#         ### Each pipeline gets a variance sum variable
#         var_sum = 0
#         ### For each row in the pipelines across replicates, calculate the variance and sum them up:
#         for row in range(0,master_dfs_uniq[samples[0]].shape[0]):
#             var_sum += statistics.variance([int(master_dfs_uniq[samples[0]].iloc[row,col]), int(master_dfs_uniq[samples[1]].iloc[row,col]), int(master_dfs_uniq[samples[2]].iloc[row,col])])
#         ### Add var_sum to chi2_var dic for respective pipeline
#         params_dic[pipeline]["variance"]=var_sum


## 3.3 Assign which tools have been used each step in the pipeline, with one column per step,
##     as well as column for coordinates (Chi-Square statistics and summed variance) normalized to range 0-1,
##     and add a column to find the distance between (1|1) (reverted normalized origin = best case) and each normalized pipeline coordinate

### 3.3.1 Calculate min and max for Chi-Square statistics and summed variance;
chi2_min = chi2_var[min(chi2_var.keys(), key=(lambda k: chi2_var[k][0]))][0]
chi2_max = chi2_var[max(chi2_var.keys(), key=(lambda k: chi2_var[k][0]))][0]
var_min = chi2_var[min(chi2_var.keys(), key=(lambda k: chi2_var[k][1]))][1]
var_max = chi2_var[max(chi2_var.keys(), key=(lambda k: chi2_var[k][1]))][1]

### 3.3.2 Add columns for tools, normalized coordinates, and distance
for pipeline in chi2_var:
    #### Add normalized Chi-Square statistics and summed variance coordinates (reversed for accuracy so that 1 is best and 0 is worst)
    chi2_var[pipeline].extend([1-normalize(chi2_var[pipeline][0], chi2_min, chi2_max), 1-normalize(chi2_var[pipeline][1], var_min, var_max)])
    #### Add distance from normalized coordinates to normalized origin
    chi2_var[pipeline].append(point_dist([chi2_var[pipeline][2], chi2_var[pipeline][3]], [1,1]))
    #### Add tool names
    pipeline_replace = pipeline.replace("IDBA_", "IDBA-").replace("NCBI_NT", "NCBI-NT").replace("BLAST_FIRST_HIT", "BLAST-FIRST-HIT").replace("BLAST_FILTERED", "BLAST-FILTERED")
    chi2_var[pipeline].extend(pipeline_replace.split("_"))





## 3.4 Save "chi2_var" dic as csv so that it can be plotted using R:
# df_save = pd.DataFrame.from_dict(chi2_var, orient="index", columns=["chi-square statistics", "summed variance", "normalized chi-square statistics", "normalized summed variance", "distance from origin", "trimmed PHRED", "rRNA filter", "assembler", "mapper", "DB", "classification"])
# df_save.to_csv("{savedir}chi2_var.csv".format(savedir=savedir), index_label="pipeline")
#
#
# # 4. Get abundance for mock community members from "best" pipeline (NOTE: only for talk since we change the approach later)
# dic_abun_obs = {"abun": [abun / 3 for abun in master_df_summarized["10_BARRNAP_IDBA_TRAN_BOWTIE2_NCBI_NT_KRAKEN2"]], "taxa": unique_taxa}
# df_abun_obs = pd.DataFrame.from_dict(dic_abun_obs)
# df_abun_obs.to_csv("{savedir}df_abun_obs.csv".format(savedir=savedir))
#
# dic_abun_exp = {"abun": [abun / 3 for abun in master_df_summarized["expected"]], "taxa": unique_taxa}
# df_abun_exp = pd.DataFrame.from_dict(dic_abun_exp)
# df_abun_exp.to_csv("{savedir}df_abun_exp.csv".format(savedir=savedir))




# OLD CODE:

# ## 3.1 Chi-Squared test (accuracy)
#
# ### 3.1.1 To perform a Chi-Squared test on replicates, we summarize columns across the replicates:
# master_df_summarized = pd.DataFrame()
# ### For all columns (note: all samples contain the same column names, so we manually pick one (the first in "samples" list))
#
# ### NOTE: For now, some pipelines didn't work, so we have to make a list of shared pipelines                                                            ### TO BE DELETED
# shared_pipelines = []                                                                                                                                   ### TO BE DELETED
# for sample in samples:                                                                                                                                  ### TO BE DELETED
#     shared_pipelines.extend(master_dfs_uniq_abs[sample].columns.tolist()[1:])                                                                               ### TO BE DELETED
# shared_pipelines = list(set(i for i in shared_pipelines if shared_pipelines.count(i) > 2))                                                              ### TO BE DELETED
#
# #for col in master_dfs_uniq_abs[samples[0]].columns.tolist():
# shared_pipelines_exp = shared_pipelines[:]                                                                                                              ### TO BE DELETED
# shared_pipelines_exp.append("expected")                                                                                                                 ### TO BE DELETED
# for pipeline in shared_pipelines_exp:                                                                                                                   ### TO BE DELETED
#     #### Make list per column in every sample and save in dir:
#     pipeline_dic_chi2 = {}
#     for sample in samples:
#         pipeline_dic_chi2[sample] = master_dfs_uniq_abs[sample][pipeline].tolist()
#     #### Summarize the three columns in the three replicates:
#     master_df_summarized[pipeline] = [a + b + c for a, b, c in zip(pipeline_dic_chi2[samples[0]], pipeline_dic_chi2[samples[1]], pipeline_dic_chi2[samples[2]])]
#
# ### 3.1.2 Now we perform a Chi-Squared test on the summarized master df, for each pipeline against the expected composition
# ###       (note, this is the DIRTY version where we add the same amount (1) to each value to get rid of zeros in the expected columns, because otherwise the Chi-Squared test is not possible):
#
# #for pipeline in master_df_summarized.columns.tolist()[1:]:
# for pipeline in shared_pipelines:                                                                                                                       ### TO BE DELETED
#     chi2_var[pipeline] = [(chisquare([abun + 1 for abun in master_df_summarized[pipeline].tolist()],[abun + 1 for abun in master_df_summarized['expected'].tolist()])[0])]
#
#
# ## 3.2 Sum of variance for each taxon (precision)
# ##     Finally, we calculate the variance for all pipelines across the three replicates by summing up the variance of each taxa across the
# ##     replicates for every pipeline and add the returned summed variance to the pipelines in the "chi2_var" dic
#
# #for pipeline in master_dfs_uniq_abs[samples[0]].columns.tolist()[1:]:
# for pipeline in shared_pipelines:                                                                                                                       ### TO BE DELETED
#     ### Turn pipelines into indexes
#     col = master_dfs_uniq_abs[samples[0]].columns.tolist()[1:].index(pipeline)
#     ### Each pipeline gets a variance sum variable
#     var_sum = 0
#     ### For each row in the pipelines across replicates, calculate the variance and sum them up:
#     for row in range(0,master_dfs_uniq_abs[samples[0]].shape[0]):
#         var_sum += statistics.variance([int(master_dfs_uniq_abs[samples[0]].iloc[row,col]), int(master_dfs_uniq_abs[samples[1]].iloc[row,col]), int(master_dfs_uniq_abs[samples[2]].iloc[row,col])])
#     ### Add var_sum to chi2_var dic for respective pipeline
#     chi2_var[pipeline].append(var_sum)
#
#
