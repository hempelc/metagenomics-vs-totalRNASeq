#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 3 Dec 2021

# This script processes pipeline data from multiple replicates of environmental
# samples and compares them to clusters based on mock community samples

import pandas as pd
import glob
import os
import copy
import pickle
import logging

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters
## Full path to directory that contains samples
workdir = "/Users/christopherhempel/Desktop/pipeline_results_env_samples/"
## List of DNA and RNA samples, replicates of 3 plus filtration controls (Neg) and
## extraction controls (Ext); must equal names of directories in workdir that
## contain each sample's pipeline results:
samples = ["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA",
    "M_Neg_DNA", "M_Neg_RNA", "M_Ext_DNA", "M_Ext_RNA"]
## Taxonomic rank to group rows on. Either based on genus (option "genus")
## or on species (option "species"):
groupby_rank = "species"
## Output dir
exportdir=os.path.join(workdir, "results_env_samples")
if not os.path.exists(exportdir):
    os.makedirs(exportdir)


# 1 Read in all pipeline results for every sample as df and
#   add all taxa from all dfs to "all_taxa" list (needed to generate master dfs that
#   contain counts of all taxa from all samples):
master_dfs_raw = {} # Empty dic that will eventually contain all samples' raw pipelines output
all_taxa = [] # Empty list that will eventually contain all taxa that appear in all samples
for sample in samples:
    ## Make a list for all file names in sample dic:
    sample_files = glob.glob(os.path.join(workdir, sample, "*.txt"))
    ## Make a dic that will eventually contain all pipeline dfs and set the first entry to expected community:
    sample_dfs = {}
    ## For each file in the sample dic
    for file in sample_files:
        #### Read in file as pandas df, fill NaN with "NA", "Unknown" by "NA", and fix one taxonomic misambiguation
        df = pd.read_table(file).fillna("NA").replace("Lactobacillus",
            r"Limosilactobacillus", regex=True).replace("Unknown",
            "NA").replace("-", r"", regex=True)
        ### Apply a species filter: if a species is not 2 words (contains a space),
        ### replace species value with "NA"
        #### Therefore, first get indices of species not containing a space
        idx=df['species'].str.contains(" ")[df['species'].str.contains(" ") == False].index
        #### And replace them with "NA" in the df
        df.loc[idx,'species'] = "NA"
        ### Cut df down to relevant columns
        if groupby_rank == "species":
            df_small = df[["superkingdom", "phylum", "class", "order", "family",
                "genus", "species", "counts"]]
        elif groupby_rank == "genus":
            df_small = df[["superkingdom", "phylum", "class", "order", "family",
                "genus", "counts"]]
        ### The negative controls often have no sequences = empty dfs, therefore we need to
        ### ignore them in the next step since we get errors if we use groupby on an epmty df:
        if df.empty:
            df_agg = df_small
        else:
            #### Group similar taxonomy hits and sum their counts:
            df_agg = df_small.groupby(list(df_small.columns)[:-1]).sum().reset_index()
            #### Turn counts into relative abundances:
            df_agg["counts"]=df_agg["counts"]/df_agg["counts"].sum()
        ### Rename counts col
        df_agg.rename(columns = {'counts':'rel_abun'}, inplace = True)
        ### Add all taxa to list "all_taxa"
        all_taxa.extend(df_agg[groupby_rank].tolist())
        ### Edit file name so that we can name dfs based on their file name=pipeline
        pipeline_name = file.split("/")[-1].split(".")[-2].split("trimmed_at_phred_")[1].split("_final")[0].replace("idba_",
            "idba-").replace("ncbi_nt", "ncbi-nt").replace("blast_first_hit",
            "blast-first-hit").replace("blast_filtered", "blast-filtered")
        ### Add df_agg to the sample_dfs dic with key=pipeline_name
        sample_dfs[pipeline_name] = df_agg
        ### Save sample_df in dic master_dfs_raw:
        master_dfs_raw[sample] = sample_dfs


# 2 Generate master dfs with relative read counts

## 2.1 Drop duplicates:
unique_taxa = list(set(all_taxa))

## 2.2 Generate master df with relative read counts for every sample and save
##     in dic master_dfs_uniq:

### Empty dic that will eventually contain all samples' master dfs
### (dfs with rows = all unique taxa):
master_dfs_rel = {}
for sample, pipeline_dfs in master_dfs_raw.items():
    ### Make master df with taxa from unique_taxa list as row names:
    master_dfs_rel[sample] = pd.DataFrame(index=pd.Index(unique_taxa))
    ### For each df in dic with key=pipeline:
    for pipeline, data in pipeline_dfs.items():
        #### Make a list for abundances
        abun=[]
        #### For each taxon in unique taxa list:
        for taxon in unique_taxa:
            ##### If taxon is in df groupby_rank column:
            if (data[groupby_rank] == taxon).any():
                ###### Sum up all counts of that taxon and add them to list abun:
                abun.append(data.loc[data[groupby_rank] == taxon, 'rel_abun'].sum())
            ##### If taxon not in df groupby_rank column, add 0 to list abun:
            else:
                abun.append(0)
        #### Make a new column in master df named after the pipeline and add
        #### taxon counts for that pipeline:
        master_dfs_rel[sample][pipeline]=abun

## 2.3 Add DNA or RNA prefixes to pipeline names to distinguish them later
## (the name of expected gets stripped off the prefix after):
master_dfs_prefix={}
for sample_type in ["DNA", "RNA"]:
    for sample in [x + "_" + sample_type for x in ["M4", "M5", "M6", "M_Neg", "M_Ext"]]:
        master_dfs_prefix[sample]=master_dfs_rel[sample].add_prefix(sample_type
            + "_").rename(columns={sample_type + "_expected": "expected"})

## 2.4 Substract controls from samples
master_dfs_rel_sub={}
for sample_type in ["DNA", "RNA"]:
    for sample in [x + "_" + sample_type for x in ["M4", "M5", "M6"]]:
        #### We substract twice the reads occuring in the filtration and extraction
        #### control from the samples, separately for RNA and DNA controls:
        master_dfs_rel_sub[sample]=master_dfs_prefix[sample]-(master_dfs_prefix["M_Neg_"
            + sample_type] + master_dfs_prefix["M_Ext_" + sample_type])*2
### Convert counts below 0 to 0 (happens if negative control contains more reads than original sample):
for key, value in master_dfs_rel_sub.items():
    value[value < 0] = 0

## 2.5 Generate master df with presence/absence data
## (0=not found, 1=found) for every sample and save in dic master_dfs_pa:
### Deepcopy the relative abundance master dfs:
master_dfs_pa = copy.deepcopy(master_dfs_rel_sub)
### Replace all values above 0 with 1:
for key, value in master_dfs_pa.items():
    value[value > 0] = 1



## 3 Calculate the average between replicates for rel and p/a data
for master_dfs in [master_dfs_rel_sub, master_dfs_pa]:
    if master_dfs is master_dfs_rel_sub:
        data_type="rel"
    else:
        data_type="pa"

    ### Concatenate all 3 dfs into one
    concat=pd.DataFrame({})
    for sample in master_dfs.keys():
        concat=pd.concat((concat, master_dfs[sample]), axis=1)

    ### Fill dic with average
    master_dic={}
    for pipeline in list(set(concat.columns)):
        master_dic[pipeline]={}
        for taxon, values in concat[pipeline].iterrows():
            master_dic[pipeline][taxon]=sum(values)/len(values)

    ## Make the df
    master_df=pd.DataFrame(master_dic).transpose()

    ## Save the df
    master_df.to_csv(os.path.join(exportdir, "master_df_env_" + groupby_rank + "_" + data_type + ".csv"), index_label="pipeline")
