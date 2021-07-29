#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 16 Jul 2021

# Preprocesses pipeline data and export metrics table

import pandas as pd
import glob
import os
import copy
import pickle

# FUTURE ADAPTATION: "lowest_hit" needs to be adapted to "species" once we rerun pipelines (can simply literatelly be find-replace'd)

# Parameters
workdir = "/Users/christopherhempel/Desktop/mock_community_RNA/" # Full path to directory that contains samples, HAS TO END WITH "/"
sample_names = ["M4", "M5", "M6", "M_Neg", "M_Ext"] # 3 replicate samples to include into the analysis (must equal names of directories in workdir that contain each sample's pipeline results)
groupby_rank = "lowest_hit" # Basis for taxa rank to group rows on. Either based on genus (option "genus") or on species (option "lowest_hit" (NOTE later "species"))
rel_abun_basis = "cell" # Basis for relative abundance calculation of expected mock community taxa abundance. Either based on genomic DNA (option "gen") or on cell number (option "cell")


# Hardcoded variables
rel_abun_gen = [0.891, 0.089, 0.0089, 0.0089, 0.00089, 0.00089, 0.000089, 0.0000089, 0.0000089, 0.00000089] # Relative abundances of mock community taxa based on genomic DNA
rel_abun_cell = [0.949, 0.042, 0.007, 0.0012, 0.00058, 0.00059, 0.00015, 0.00001, 0.0000007, 0.000001] # Relative abundances of mock community taxa based on cell number


# Functions
## Function to calculate the average absolute deviation (AAD) from a central point:
def aad (reps, central):
    aad=(sum([abs(x-central) for x in reps]))/len(reps)
    return aad

## Function to cut down master_dfs to only expected taxa:
def cutdown (master_df, type):
    cutdown_dic = {}
    for sample, pipelines in master_df.items():
        if type == "above_zero":
            cutdown_dic[sample] = pipelines.loc[pipelines['expected'] != 0]
        elif type == "equal_zero":
            cutdown_dic[sample] = pipelines.loc[pipelines['expected'] == 0]
    return cutdown_dic

## Function to calculate metrics for all samples for all 3 master_dfs:
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
                if abundances.tolist()[i] > expected_list[i]:
                    FP += (abundances.tolist()[i] - expected_list[i])
                if abundances.tolist()[i] != 0 and abundances.tolist()[i] >= expected_list[i]:
                    TP += expected_list[i]
                if abundances.tolist()[i] < expected_list[i]:
                    TP += abundances.tolist()[i]
            if type == "rel":
                for i in range(len(abundances.tolist())):
                    if abundances.tolist()[i] < expected_list[i]:
                        FN += (abundances.tolist()[i] - expected_list[i])*-1
                confusion_values = {"subceed_reads": FN, "exceed_reads": FP, "true_reads":TP}
                confusion_values["false_reads"]=cutdown_dic_FP[sample][pipeline].sum()
            elif type == "pa":
                for i in range(len(abundances.tolist())):
                    if abundances.tolist()[i] - expected_list[i] < 0:
                        FN += (abundances.tolist()[i] - expected_list[i])*-1
                confusion_values = {"FN": FN, "TP":TP}
                confusion_values["FP"]=cutdown_dic_FP[sample][pipeline].sum()
                confusion_values["TN"]=len(cutdown_dic_FP[sample][pipeline])-cutdown_dic_FP[sample][pipeline].tolist().count(1)
            confusion_dic[pipeline] = confusion_values
        confusion_master[sample] = confusion_dic
    return confusion_master

# 1 Make expected mock community
expected_dic = {} # Empty dic that will eventually contain expected species adn their abundances

### To calculate the absolute expected read count for the taxa in the mock community, we need the hardcoded taxonomic information of each taxon. We save it in a dictionary:
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


### We also need relative abundances per taxon. Therefore, we use the relative abundance of each taxon in the mock community, either based on SSU genes or cells.
if rel_abun_basis == "gen":
    rel_abun = rel_abun_gen
elif rel_abun_basis == "cell":
    rel_abun = rel_abun_cell

### We add each taxon's read count to its taxonomy list:
rep=0
for taxon in expected_dic:
    expected_dic[taxon].append(rel_abun[rep])
    rep+=1

### And finally read in the expected mock community dic as pandas df
expected_df=pd.DataFrame.from_dict(expected_dic, orient="index", columns=["superkingdom", "phylum", "class", "order", "family", "genus", "lowest_hit", "rel_abun"])


# 2 Read in all pipeline results for every sample as df;
#   add all taxa from all dfs to "all_taxa" list, needed to generate master dfs that
#   contain counts of all taxa from all samples and mock community:
master_dfs_raw = {} # Empty dic that will eventually contain all samples' raw pipelines output
all_taxa = [] # Empty list that will eventually contain all taxa that appear in all samples

for sample_type in ["DNA", "RNA"]:
    samples=[x + "_" + sample_type for x in sample_names]
    for sample in samples:
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
                #### Turn counts into relative abundances and rename
                df_agg["counts"]=df_agg["counts"]/df_agg["counts"].sum()
                df_agg.rename(columns = {'counts':'rel_abun'}, inplace = True)
                #### Add all taxa to list "all taxa"
                all_taxa.extend(df_agg[groupby_rank].tolist())
                #### Edit file name so that pipeline name is displayed properly
                pipeline_name=sample_type + "_" + file.split("/")[-1].split(".")[-2].split("trimmed_at_phred_")[1].split("_pipeline_final")[0].replace("IDBA_", "IDBA-").replace("NCBI_NT", "NCBI-NT").replace("BLAST_FIRST_HIT", "BLAST-FIRST-HIT").replace("BLAST_FILTERED", "BLAST-FILTERED")
                #### Make dict key with keys = piepline_name and value = edited df
                sample_dfs[pipeline_name] = df_agg
            else:                                                                                                                                            ### TO BE DELETED
                continue                                                                                                                                     ### TO BE DELETED
        ### Save sample_df in dic "master_dfs_raw":
        master_dfs_raw[sample] = sample_dfs



# 3 Make a "unique_taxa" list that contains all taxa that appear in all samples;
#   generate a "master_df" df from each sample that contains counts for all taxa in all samples and mock community:

## 3.1 Add expected taxa to list in case they're not picked up by any pipeline and drop duplicates:
all_taxa.extend(expected_df[groupby_rank].tolist())
unique_taxa = list(set(all_taxa))


## 3.2 Generate master df with relative read counts for every sample and save in dic "master_dfs_uniq":
master_dfs_uniq = {} # Empty dic that will eventually contain all samples' master dfs (dfs with rows = all unique taxa)
for sample, pipeline_dfs in master_dfs_raw.items():
    ### Make master df with taxa from "unique_taxa" list as row names
    master_dfs_uniq[sample] = pd.DataFrame(index=pd.Index(unique_taxa))
    ### For each df in dict with key=pipeline
    for pipeline, data in pipeline_dfs.items():
        abun=[]
        #### For each taxon in unique taxa list
        for taxon in unique_taxa:
            ##### If taxon is in df groupby_rank column
            if (data[groupby_rank] == taxon).any():
                ###### Sum up all counts of that taxon and add it to list "counts"
                abun.append(data.loc[data[groupby_rank] == taxon, 'rel_abun'].sum())
            ##### If taxon not in df groupby_rank column, add 0 to list "counts"
            else:
                abun.append(0)
        #### Make a new column in  master df named after the pipeline and add taxon counts for that pipeline
        master_dfs_uniq[sample][pipeline]=abun

### NOTE: For now, some pipelines didn't work, so we have to make a list of shared pipelines                                                            ### TO BE DELETED
shared_pipelines_dupl = []                                                                                                                                   ### TO BE DELETED
for sample in samples:                                                                                                                                  ### TO BE DELETED
    shared_pipelines_dupl.extend(master_dfs_uniq[sample].columns.tolist())                                                                             ### TO BE DELETED
shared_pipelines = []                                                                                                                                  ### TO BE DELETED
for i in shared_pipelines_dupl:                                                                                                                                  ### TO BE DELETED
    if shared_pipelines_dupl.count(i) > 2:                                                                                                                               ### TO BE DELETED
        if i not in shared_pipelines:                                                                                                                                  ### TO BE DELETED
            shared_pipelines.append(i)                                                                                                                                  ### TO BE DELETED
                                                                                                                           ### TO BE DELETED
## And we cut the master df down to these pipelines:                                                                                                                           ### TO BE DELETED
master_dfs_uniq_shared={}                                                                                                                           ### TO BE DELETED
for sample in samples:                                                                                                                           ### TO BE DELETED
    master_dfs_uniq_shared[sample]={}                                                                                                                           ### TO BE DELETED
                                                                                                                               ### TO BE DELETED
for i in shared_pipelines:                                                                                                                           ### TO BE DELETED
    for sample in samples:                                                                                                                           ### TO BE DELETED
        master_dfs_uniq_shared[sample][i]=master_dfs_uniq[sample][i]                                                                                                                           ### TO BE DELETED

master_dfs_uniq_del={}                                                                                                                           ### TO BE DELETED
for sample in samples:                                                                                                                           ### TO BE DELETED
    master_dfs_uniq_del[sample]=pd.DataFrame(master_dfs_uniq_shared[sample])                                                                                                                           ### TO BE DELETED

master_dfs_uniq = copy.deepcopy(master_dfs_uniq_del)                                                                                                                     ### TO BE DELETED

# We'll save that df as "rel"
master_dfs_uniq_rel = master_dfs_uniq

## 3.3 Substract controls from samples
master_dfs_uniq_rel_sub={}
for sample_type in ["DNA", "RNA"]:
    samples=[x + "_" + sample_type for x in ["M4", "M5", "M6"]]
    for sample in samples:
        # We substract twice the reads occuring in the extractio and filtration control
        master_dfs_uniq_rel_sub[sample]=master_dfs_uniq_rel[sample]-(master_dfs_uniq_rel["M_Neg_" + sample_type] + master_dfs_uniq_rel["M_Ext_" + sample_type])*2
        # And since this sustarction overrides the expected, we replace it with the old expected
        master_dfs_uniq_rel_sub[sample]["expected"]=master_dfs_uniq_rel[sample]["expected"]

## 3.4 Generate master df with presence/absence for every sample and save in dic "master_dfs_uniq_pa":
### Deepcopy one of the master dfs
master_dfs_uniq_pa = copy.deepcopy(master_dfs_uniq_rel_sub)
### Replace all non-zero values with 1
for key, value in master_dfs_uniq_pa.items():
    value[value != 0] = 1



# 4 Calculate metrics

## 4.1 Cut down master_dfs to only expected taxa
### Apply cutdown function on absolute, relative, and pa master_dfs
rel_cutdown_expected = cutdown(master_dfs_uniq_rel_sub, "above_zero")
pa_cutdown_expected = cutdown(master_dfs_uniq_pa, "above_zero")

rel_cutdown_FP = cutdown(master_dfs_uniq_rel_sub, "equal_zero")
pa_cutdown_FP = cutdown(master_dfs_uniq_pa, "equal_zero")


## 4.2 Apply confusion_calc function on absolute, relative, and pa cutdown master_dfs
rel_metrics_reps = confusion_calc(rel_cutdown_expected, rel_cutdown_FP, "rel")
pa_metrics_reps = confusion_calc(pa_cutdown_expected, pa_cutdown_FP, "pa")


## 4.3 Calculate the average between replicates:
### Make dics containing each pipeline as key with one list per metric for each pipeline:
rel_metrics={} # Empty dic that will eventually contain all relative abundance-based metrics for PCA
pa_metrics={} # Empty dic that will eventually contain all p/a-based metrics for PCA
for metrics_dic in [rel_metrics, pa_metrics]:
    for sample in samples:                                                                                                      ### TO BE DELETED since all samples will contain the same pipelines eventually
        for pipeline in rel_metrics_reps[sample].keys():
            if  metrics_dic==rel_metrics:
                metrics_dic[pipeline]={"subceed_reads": [], "exceed_reads": [], "true_reads": [], "false_reads": []}
            elif metrics_dic==pa_metrics:
                metrics_dic[pipeline]={"FN": [], "FP": [], "TP": [], "TN": []}

### Fill in dics created above so that every pipeline contains a dic with
### metrics as keys and each metric lists its three values across the
### three replicates:
for metric_dic_rep in [rel_metrics_reps, pa_metrics_reps]:
    for sample, pipelines in metric_dic_rep.items():
        for pipeline, metrics in pipelines.items():
            for metric, value in metrics.items():
                if metric_dic_rep==rel_metrics_reps:
                    rel_metrics[pipeline][metric].append(value)
                elif metric_dic_rep==pa_metrics_reps:
                    pa_metrics[pipeline][metric].append(value)

### Calculate the average across the three replicates:
for metrics_dic in [rel_metrics, pa_metrics]:
    for pipeline, metrics in metrics_dic.items():
        for metric, lst in metrics.items():
            metrics_dic[pipeline][metric]=sum(lst)/len(lst)


## 4.4 Caluclate the AAD across the three replicates for both types of datasets
### 4.4.1 For p/a data
metric_list_pa=["FN", "FP", "TP", "TN"]
# ideal_vals_pa=[0, 0, 9, len(unique_taxa)-len(expected_df.index)] # for average distance from ecpected
for pipeline in pa_metrics.keys():
    for metric in metric_list_pa:
        metric_rep_values_pa=[]
        for sample in samples:
            metric_rep_values_pa.append(float(pa_metrics_reps[sample][pipeline][metric]))
        pa_metrics[pipeline][metric + "_aad"]=aad(metric_rep_values_pa, sum(metric_rep_values_pa)/len(metric_rep_values_pa))
        #pa_metrics[pipeline][metric + "_aad"]=aad(metric_rep_values_pa, ideal_vals_pa[metric_list_pa.index(metric)]) # for average distance from ecpected


### 4.4.2 For relative data
metric_list_rel=["subceed_reads", "exceed_reads", "true_reads", "false_reads"]
# ideal_vals_rel=[0, 0, 1, 0]  # for average distance from ecpected
for pipeline in rel_metrics.keys():
    for metric in metric_list_rel:
        metric_rep_values_rel=[]
        for sample in samples:
            metric_rep_values_rel.append(float(rel_metrics_reps[sample][pipeline][metric]))
        rel_metrics[pipeline][metric + "_aad"]=aad(metric_rep_values_rel, sum(metric_rep_values_rel)/len(metric_rep_values_rel))
        #pa_metrics[pipeline][metric + "_aad"]=aad(metric_rep_values_rel, ideal_vals_rel[metric_list_rel.index(metric)]) # for average distance from ecpected



# 5 Make master metrics df
## Merge separate metric dics into one
master_metrics_df={} # Empty dic that will eventually contain all metrics for PCA
for pipeline in shared_pipelines[1:len(shared_pipelines)]:                                                                                                                ### TO BE DELETED
    master_metrics_df[pipeline]={}
    for metrics_dic in [rel_metrics, pa_metrics]:
        for metric in metrics_dic[pipeline].keys():
            master_metrics_df[pipeline][metric]=metrics_dic[pipeline][metric]

## Add info on tools used in each step for each pipeline
step_list=["trimming_score", "rRNA_sorting_tool", "assembly_tool", "mapper", "database", "classifier"]
for step in step_list:
    for pipeline in master_metrics_df.keys():
        master_metrics_df[pipeline][step]=pipeline.split("_")[step_list.index(step)]


## Make the df
metrics_df=pd.DataFrame(master_metrics_df)
metrics_df=metrics_df.transpose()

## Save the df
statsdir=os.path.join(workdir, "stats_exports")
if not os.path.exists(statsdir):
    os.makedirs(statsdir)
metrics_df.to_csv(os.path.join(statsdir, "metrics_df.csv"), index_label="pipeline")

## Also save a pickle object with some information needed for next step of code
with open(os.path.join(statsdir, "TN.pkl"), 'wb') as f:
    pickle.dump(len(unique_taxa)-len(expected_df), f)
