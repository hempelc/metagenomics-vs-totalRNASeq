#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 3 Dec 2021

# This script processes pipeline data from multiple replicates of environmental
# samples and compares them to clusters based on mock community samples

import pandas as pd
import glob
import os
import copy
import logging
from sklearn.preprocessing import StandardScaler
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')

# Parameters set manually
## Full path to directory that contains env/mock samples
workdir_env = "/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_env_samples/"
workdir_mock = "/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_mock_samples/"
## List of mock/env DNA and RNA samples, replicates of 3 plus filtration controls (Neg) and
## extraction controls (Ext); must equal names of directories in workdir that
## contain each sample's pipeline results:
samples_mock = ["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA",
    "M_Neg_DNA", "M_Neg_RNA", "M_Ext_DNA", "M_Ext_RNA"]
samples_env = ["F4_DNA", "F4_RNA", "F5_DNA", "F5_RNA", "F6_DNA", "F6_RNA",
    "F_Neg_DNA", "F_Neg_RNA", "F_Ext_DNA", "F_Ext_RNA"]
## Dics containing number of reads per sample
sample_reads_mock={"M4_DNA": 817619, "M4_RNA": 94633,
    "M5_DNA": 644634, "M5_RNA": 78149, "M6_DNA": 669382, "M6_RNA": 120144,
    "M_Ext_DNA": 399, "M_Ext_RNA": 887, "M_Neg_DNA": 640, "M_Neg_RNA": 1551}
sample_reads_env={"F4_DNA": 1355159, "F4_RNA": 1902388, "F5_DNA": 1373310, "F5_RNA": 1099851,
    "F6_DNA": 1571705, "F6_RNA": 773067, "F_Ext_DNA": 99, "F_Ext_RNA": 2685,
    "F_Neg_DNA": 157, "F_Neg_RNA": 1137}
## Indicate if you want to loop over env and mock (True/False)
looping_sample=True
## If you set looping_rank to False, then define what specific sample set you want
## to process (option "mock" or "env"):
sample_set_name="env"
## Indicate if you want to loop over genus and species (True/False)
looping_rank=True
## If you set looping_rank to False, then define what specific rank you want to process:
### Taxonomic rank to group rows on. Either based on genus (option "genus")
### or on species (option "species"):
rank="genus"


# Parameters set automatically
if looping_sample:
    sample_set_name_lst=["mock", "env"]
else:
    sample_set_name_lst=[sample_set_name]
if looping_rank:
    groupby_rank_lst=["genus", "species"]
else:
    groupby_rank_lst=[rank]

for sample_set in sample_set_name_lst:
    if sample_set=="mock":
        workdir = workdir_mock
        samples=samples_mock
    else:
        workdir = workdir_env
        samples=samples_env

    for groupby_rank in groupby_rank_lst:
        ## Set output dir
        exportdir=os.path.join(workdir, "rel_abun_{0}".format(groupby_rank))
        if not os.path.exists(exportdir):
            os.makedirs(exportdir)

        # 1 Read in all pipeline results for every sample as df and
        #   add all taxa from all dfs to "all_taxa" list (needed to generate master dfs that
        #   contain counts of all taxa from all samples):
        master_dfs_raw = {} # Empty dic that will eventually contain all samples' raw pipelines output
        all_taxa = [] # Empty list that will eventually contain all taxa that appear in all samples
        for sample in samples:
            ## Make a list for all file names in sample dic:
            sample_files = glob.glob(os.path.join(workdir, sample, "*.txt*"))
            ## Make a dic that will eventually contain all pipeline dfs:
            sample_dfs = {}
            ## For each file in the sample dic
            for file in sample_files:
                #### Read in file as pandas df, fill NaN with "NA", "Unknown" by "NA", and fix one taxonomic misambiguation
                df = pd.read_table(file).replace("Lactobacillus", r"Limosilactobacillus", regex=True)\
                    .replace("Unknown", "NA").replace("-", r"", regex=True)
                df=df.rename(columns={df.columns[0]: 'sequence_name'})\
                    .dropna(subset = ['sequence_name']).fillna("NA")
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
                    ### Apply a species filter: if a species is not 2 words (contains a space),
                    ### replace species value with "NA"
                    #### Therefore, first get indices of species not containing a space
                    idx=df['species'].str.contains(" ")[df['species'].str.contains(" ") == False].index
                    #### And replace them with "NA" in the df
                    df.loc[idx,'species'] = "NA"
                    ### Cut df down to relevant columns
                    #### Group similar taxonomy hits and sum their counts:
                    df_agg = df_small.groupby(list(df_small.columns)[:-1]).sum().reset_index()
                    #### Turn counts into relative abundances:
                    df_agg["counts"]=df_agg["counts"]/df_agg["counts"].sum()
                ### Rename counts col
                df_agg.rename(columns = {'counts':'rel_abun'}, inplace = True)
                ### Add all taxa to list "all_taxa"
                all_taxa.extend(df_agg[groupby_rank].tolist())
                ### Edit file name so that we can name dfs based on their file name=pipeline
                pipeline_name = file.split("/")[-1].split(".")[0].split("trimmed_at_phred_")[1].split("_final")[0].replace("idba_",
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
            if sample_set=="mock":
                names_and_contr=["M4", "M5", "M6", "M_Neg", "M_Ext"]
            else:
                names_and_contr=["F4", "F5", "F6", "F_Neg", "F_Ext"]
            for sample in [x + "_" + sample_type for x in names_and_contr]:
                master_dfs_prefix[sample]=master_dfs_rel[sample].add_prefix(sample_type
                    + "_").rename(columns={sample_type + "_expected": "expected"})

        ## 2.4 Substract controls from samples
        master_dfs_rel_sub={}
        for sample_type in ["DNA", "RNA"]:
            if sample_set=="mock":
                read_dic=sample_reads_mock
                neg="M_Neg_" + sample_type
                ext="M_Ext_" + sample_type
                names=["M4", "M5", "M6"]
            else:
                read_dic=sample_reads_env
                neg="F_Neg_" + sample_type
                ext="F_Ext_" + sample_type
                names=["F4", "F5", "F6"]
            neg_readnum=read_dic[neg]
            ext_readnum=read_dic[ext]
            for sample in [x + "_" + sample_type for x in names]:
                #### We substract the reads occuring in the filtration and extraction
                #### control from the samples, separately for RNA and DNA controls:
                ##### We're converting counts back to absolute and substract absolute numbers of reads of the negative controls
                readnum=read_dic[sample]
                master_dfs_rel_sub[sample]=master_dfs_prefix[sample]*readnum-(master_dfs_prefix[neg]*neg_readnum \
                    + master_dfs_prefix[ext]*ext_readnum)
                ### Convert counts below 0 to 0 (happens if negative control contains more reads than original sample):
                master_dfs_rel_sub[sample][master_dfs_rel_sub[sample] < 0] = 0
                ### Convert counts back to relative
                master_dfs_rel_sub[sample]=master_dfs_rel_sub[sample]/master_dfs_rel_sub[sample].sum()


        ## 2.5 Generate master df with presence/absence data
        ## (0=not found, 1=found) for every sample and save in dic master_dfs_pa:
        ### Deepcopy the relative abundance master dfs:
        master_dfs_pa = copy.deepcopy(master_dfs_rel_sub)
        ### Replace all values above 0 with 1:
        for key, value in master_dfs_pa.items():
            value[value > 0] = 1



        # 3 Calculate the average between replicates for rel and p/a data
        for master_dfs in [master_dfs_rel_sub, master_dfs_pa]:
            if master_dfs is master_dfs_rel_sub:
                data_type="rel"
            else:
                data_type="pa"
            for rep, df in master_dfs.items():
                ## Make the df
                master_df=pd.DataFrame(df).transpose()
                ## Standardize data
                ### If rel, then replace 0s and determine central log ratio
                if master_dfs is master_dfs_rel_sub:
                    master_df=pd.DataFrame(clr(multiplicative_replacement(master_df)), index=master_df.index, columns=master_df.columns)
                ## Save the df
                master_df.to_csv(os.path.join(exportdir, rep + "_rel_abun_" + groupby_rank + "_" + data_type + ".csv"), index_label="pipeline")
