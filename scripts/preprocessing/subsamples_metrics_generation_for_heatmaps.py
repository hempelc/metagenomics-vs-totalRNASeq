#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 2 Mar 2022

# This script processes pipeline data from multiple replicates of mock community
# samples and exports a metrics table. The output can be found in the workdir:
# for each processed combination, a directory is created called "metrics_"
# followed by the combination, which contains .csv files with processed metrics
# for each replicate. It's a modified version of the script
# "metrics_generation.py" that is adapted for subsampled DNA samples.

import pandas as pd #v1.3.5
import glob
import os
import copy
import logging
from skbio.stats.composition import multiplicative_replacement #v0.5.6
from skbio.stats.composition import clr #v0.5.6

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters set manually
## Full path to directory that contains samples
workdir = "/Users/christopherhempel/Desktop/pipeline_results_coverage/dna_subsamples/"
## Lists of DNA mock community samples, replicates of 3 plus filtration control (Neg) and
## extraction control (Ext); must equal names of directories in workdir that
## contain each sample's pipeline results:
samples = ["M4_DNA_subsample", "M5_DNA_subsample", "M6_DNA_subsample"]
neg_samples=["M_Neg_DNA", "M_Ext_DNA"]
## Dic containing number of reads per sample
sample_reads={"M4_DNA_subsample": 94633, "M5_DNA_subsample": 78149,
    "M6_DNA_subsample": 120144, "M_Neg_DNA": 640, "M_Ext_DNA": 399}
## Indicate if you want to loop over genus and species (True/False)
looping=True
## If you set looping to False, then define what specific rank and abundance
## you want to process:
### Taxonomic rank to group rows on. Either based on genus (option "genus")
### or on species (option "species"):
rank="species"


# Parameters set automatically
if looping:
    groupby_rank_lst=["genus", "species"]
else:
    groupby_rank_lst=[rank]

## Relative abundances of mock community taxa based on genome copy numbers as given by manufacturer:
rel_abun_gen_manuf = [0.948, 0.042, 0.007, 0.0023, 0.00058, 0.00059, 0.00015,
    0.00001, 0.0000015, 0.000001]
## These don't add up to 1 for some reason, so we shift the relative abundances
## so that they sum up to 1:
rel_abun_gen = [x/sum(rel_abun_gen_manuf) for x in rel_abun_gen_manuf]


# Functions
## Function to cut down master_dfs to only expected taxa,
## options for type="expected" if expected taxa are selected or "FP"
## if false positive taxa are selected:
def cutdown (master_df, type):
    cutdown_dic = {}
    for sample, pipelines in master_df.items():
        if type == "expected":
            cutdown_dic[sample] = pipelines.loc[pipelines['expected'] != 0]
        elif type == "FP":
            cutdown_dic[sample] = pipelines.loc[pipelines['expected'] == 0]
    return cutdown_dic



# 1 Make expected mock community df
expected_dic = {} # Empty dic that will eventually contain expected species and their abundances

### To calculate the absolute expected read count for the taxa in the mock community,
### we need the hardcoded taxonomic information of each taxon. We save it in a dictionary:
expected_dic["L_monocytogenes"] = ["Bacteria", "Firmicutes", "Bacilli",
    "Bacillales", "Listeriaceae", "Listeria", "Listeria monocytogenes"]
expected_dic["P_aeruginosa"] = ["Bacteria", "Proteobacteria",
    "Gammaproteobacteria", "Pseudomonadales", "Pseudomonadaceae", "Pseudomonas",
        "Pseudomonas aeruginosa"]
expected_dic["B_subtilis"] = ["Bacteria", "Firmicutes", "Bacilli", "Bacillales",
    "Bacilliaceae", "Bacillus", "Bacillus subtilis"]
expected_dic["S_cerevisiae"] = ["Eukaryota", "Ascomycota", "Saccharomycetes"
    ,"Saccharomycetales", "Saccharomycetaceae", "Saccharomyces", "Saccharomyces cerevisiae"]
expected_dic["E_coli"] = ["Bacteria", "Proteobacteria", "Gammaproteobacteria",
    "Enterobacterales", "Enterobacteriaceae", "Escherichia", "Escherichia coli"]
expected_dic["S_enterica"] = ["Bacteria", "Proteobacteria", "Gammaproteobacteria",
    "Enterobacterales", "Enterobacteriaceae", "Salmonella", "Salmonella enterica"]
expected_dic["L_fermentum"] = ["Bacteria", "Firmicutes", "Bacilli", "Lactobacillales",
    "Lactobacillaceae", "Limosilactobacillus", "Limosilactobacillus fermentum"]
expected_dic["E_faecalis"] = ["Bacteria", "Firmicutes", "Bacilli", "Lactobacillales",
    "Enterococcaceae", "Enterococcus", "Enterococcus faecalis"]
expected_dic["C_neoformans"] = ["Eukaryota", "Basidiomycota", "Tremellomycetes",
    "Tremellales", "Tremellaceae", "Cryptococcus", "Cryptococcus neoformans"]
expected_dic["S_aureus"] = ["Bacteria", "Firmicutes", "Bacilli", "Bacillales",
    "Staphylococcaceae", "Staphylococcus", "Staphylococcus aureus"]

### And  read in the expected mock community dic as pandas df:
expected_df=pd.DataFrame.from_dict(expected_dic, orient="index",
    columns=["superkingdom", "phylum", "class", "order", "family", "genus", "species"])

## We  need absolute abundances of mock community taxa. Therefore, we use the relative
## abundance of each taxon in the mock community. While we're at it, we determine the lowest relative
## abundance, lower it by 3 orders of magnitude, and use the result to impute 0s later
expected_df["rel_abun"] = rel_abun_gen
impute=min(rel_abun_gen)*1e-03



# 2 Read in all pipeline results for every sample as df and
#   add all taxa from all dfs to "all_taxa" list (needed to generate master dfs that
#   contain counts of all taxa from all samples and the expected community):
## Loop over groupby_rank
for groupby_rank in groupby_rank_lst:
    for data_type in ["rel", "pa"]:
        master_dfs_raw = {} # Empty dic that will eventually contain all samples' raw pipelines output
        all_taxa = [] # Empty list that will eventually contain all taxa that appear in all samples

        for sample in samples:
            ## Make a list for all file names in sample dic:
            if "M_" in sample:
                sample_files = glob.glob(os.path.join(workdir, sample, "*.txt*"))
            else:
                sample_files = glob.glob(os.path.join(workdir, sample, "*", groupby_rank + "_" + data_type + "_coverage", "*", "*.txt"))
            ## Make a dic that will eventually contain all pipeline dfs and set the first entry to expected community:
            sample_dfs = {"expected": expected_df}
            ## For each file in the sample dic
            for file in sample_files:
                ### Read in file as pandas df, fill NaN with "NA", "Unknown" by "NA", and fix one taxonomic misambiguation
                df = pd.read_table(file).replace("Lactobacillus", r"Limosilactobacillus", regex=True)\
                    .replace("Unknown", "NA").replace("-", r"", regex=True)
                df=df.rename(columns={df.columns[0]: 'sequence_name'})\
                    .dropna(subset = ['sequence_name']).fillna("NA")
                if "assembly_sequence " in df.columns:
                    df=df.rename(columns={"assembly_sequence ": 'assembly_sequence'})
                if not df.empty:
                    ### We need the scaffold length to determine covered bases, so
                    ### if that info is not available, we have to generate it from the given sequence
                    if "sequence_length" not in df.columns:
                        df["sequence_length"] = df["assembly_sequence"].str.len()
                    ### Apply a species filter: if a species is not 2 words (contains a space),
                    ### replace species value with "NA"
                    #### Therefore, first get indices of species not containing a space
                    idx=df['species'].str.contains(" ")[df['species'].str.contains(" ") == False].index
                    #### And replace them with "NA" in the df
                    df.loc[idx,'species'] = "NA"
                else:
                    df["sequence_length"]=[]
                ### Determine covered bases of scaffolds to aggregate information
                ### across scaffolds with similar taxonomic annotation
                df["covered_bases"]=df["sequence_length"]*df["coverage"]
                ### Cut df down to relevant columns
                if groupby_rank == "species":
                    df_small = df[["superkingdom", "phylum", "class", "order", "family",
                        "genus", "species", "covered_bases", "sequence_length"]]
                elif groupby_rank == "genus":
                    df_small = df[["superkingdom", "phylum", "class", "order", "family",
                        "genus", "covered_bases", "sequence_length"]]
                ### The negative controls often have no sequences = empty dfs, therefore we need to
                ### ignore them in the next step since we get errors if we use groupby on an epmty df:
                if df.empty:
                    df_agg = df_small
                else:
                    #### Group similar taxonomy hits and sum their covered bases and sequence length:
                    df_agg = df_small.groupby(list(df_small.columns)[:-2])["covered_bases", "sequence_length"].sum().reset_index()
                    #### Determine average per-base coverage for each taxon
                    df_agg["per_base_coverage"] = df_agg["covered_bases"]/df_agg["sequence_length"]
                    df_agg=df_agg.drop(["sequence_length", "covered_bases"], axis=1)
                    #### Turn coverages into relative abundances:
                    df_agg["per_base_coverage"]=df_agg["per_base_coverage"]/df_agg["per_base_coverage"].sum()
                ### Rename counts col
                df_agg.rename(columns = {'per_base_coverage':'rel_abun'}, inplace = True)
                ### Add all taxa to list "all_taxa"
                all_taxa.extend(df_agg[groupby_rank].tolist())
                ### Edit file name so that we can name dfs based on their file name=pipeline
                sample_name = file.split("/")[7]
                ### Add df_agg to the sample_dfs dic with key=subsample name
                sample_dfs[sample_name] = df_agg
            ### Save sample_df in dic master_dfs_raw:
            master_dfs_raw[sample] = sample_dfs



        # Repeat the code above for negative control samples - these were not subsampled,
        # so we import all files and select the pipeline that needs to be substracted later
        master_dfs_neg_raw = {} # Empty dic that will eventually contain all samples' raw pipelines output
        for sample in neg_samples:
            ## Make a list for all file names in sample dic:
            if "M_" in sample:
                sample_files = glob.glob(os.path.join(workdir, sample, "*.txt*"))
            else:
                sample_files = glob.glob(os.path.join(workdir, sample, "*", groupby_rank + "_" + data_type, "*", "*.txt"))
            ## For each file in the sample dic
            sample_dfs={}
            for file in sample_files:
                #### Read in file as pandas df, fill NaN with "NA", "Unknown" by "NA", and fix one taxonomic misambiguation
                df = pd.read_table(file).replace("Lactobacillus", r"Limosilactobacillus", regex=True)\
                    .replace("Unknown", "NA").replace("-", r"", regex=True)
                df=df.rename(columns={df.columns[0]: 'sequence_name'})\
                    .dropna(subset = ['sequence_name']).fillna("NA")
                if "assembly_sequence " in df.columns:
                    df=df.rename(columns={"assembly_sequence ": 'assembly_sequence'})
                if not df.empty:
                    ### We need the scaffold length to determine covered bases, so
                    ### if that info is not available, we have to generate it from the given sequence
                    if "sequence_length" not in df.columns:
                        df["sequence_length"] = df["assembly_sequence"].str.len()
                    ### Apply a species filter: if a species is not 2 words (contains a space),
                    ### replace species value with "NA"
                    #### Therefore, first get indices of species not containing a space
                    idx=df['species'].str.contains(" ")[df['species'].str.contains(" ") == False].index
                    #### And replace them with "NA" in the df
                    df.loc[idx,'species'] = "NA"
                else:
                    df["sequence_length"]=[]
                ### Determine covered bases of scaffolds to aggregate information
                ### across scaffolds with similar taxonomic annotation
                df["covered_bases"]=df["sequence_length"]*df["coverage"]
                ### Cut df down to relevant columns
                if groupby_rank == "species":
                    df_small = df[["superkingdom", "phylum", "class", "order", "family",
                        "genus", "species", "covered_bases", "sequence_length"]]
                elif groupby_rank == "genus":
                    df_small = df[["superkingdom", "phylum", "class", "order", "family",
                        "genus", "covered_bases", "sequence_length"]]
                ### The negative controls often have no sequences = empty dfs, therefore we need to
                ### ignore them in the next step since we get errors if we use groupby on an epmty df:
                if df.empty:
                    df_agg = df_small
                else:
                    #### Group similar taxonomy hits and sum their covered bases and sequence length:
                    df_agg = df_small.groupby(list(df_small.columns)[:-2])["covered_bases", "sequence_length"].sum().reset_index()
                    #### Determine average per-base coverage for each taxon
                    df_agg["per_base_coverage"] = df_agg["covered_bases"]/df_agg["sequence_length"]
                    df_agg=df_agg.drop(["sequence_length", "covered_bases"], axis=1)
                    #### Turn coverages into relative abundances:
                    df_agg["per_base_coverage"]=df_agg["per_base_coverage"]/df_agg["per_base_coverage"].sum()
                ### Rename counts col
                df_agg.rename(columns = {'per_base_coverage':'rel_abun'}, inplace = True)
                ### Add all taxa to list "all_taxa"
                all_taxa.extend(df_agg[groupby_rank].tolist())
                ### Edit file name so that we can name dfs based on their file name=pipeline
                pipeline_name = file.lower().split("/")[-1].split(".")[0].split("trimmed_at_phred_")[1].split("_final")[0].replace("_pipeline", "").replace("idba_",
                    "idba-").replace("ncbi_nt", "ncbi-nt").replace("blast_first_hit",
                    "blast-first-hit").replace("blast_filtered", "blast-filtered")
                ### Add df_agg to the sample_dfs dic with key=pipeline_name
                sample_dfs[pipeline_name] = df_agg
            ### Save sample_df in dic master_dfs_raw:
            master_dfs_neg_raw[sample] = sample_dfs



        # 3 Generate master dfs with relative read counts

        ## 3.1 Add all expected taxa to list in case they're not picked up by any pipeline
        ##     and drop duplicates:
        all_taxa.extend(expected_df[groupby_rank].tolist())
        unique_taxa = list(set(all_taxa))

        ## 3.2 Generate master df with relative read counts for every sample and save
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

        # Repeat that for negative controls (could be done together but copy-pasting the code was faster)
        master_dfs_neg_rel = {}
        for sample, pipeline_dfs in master_dfs_neg_raw.items():
            ### Make master df with taxa from unique_taxa list as row names:
            master_dfs_neg_rel[sample] = pd.DataFrame(index=pd.Index(unique_taxa))
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
                master_dfs_neg_rel[sample][pipeline]=abun

        ## 3.4 Substract controls from samples
        master_dfs_rel_sub={}
        neg_readnum=sample_reads["M_Neg_DNA"]
        ext_readnum=sample_reads["M_Ext_DNA"]
        # Define which pipeline from negative controls needs to be substracted from samples
        for sample in samples:
            if sample=="M4_DNA_subsample":
                if groupby_rank + "_" + data_type=="genus_rel":
                    pip="10_barrnap_spades_bowtie2_silva_blast-first-hit"
                elif groupby_rank + "_" + data_type=="genus_pa":
                    pip="20_barrnap_metaspades_bwa_silva_kraken2"
                elif groupby_rank + "_" + data_type=="species_rel":
                    pip="20_unsorted_metaspades_bowtie2_silva_kraken2"
                elif groupby_rank + "_" + data_type=="species_pa":
                    pip="20_barrnap_metaspades_bwa_silva_kraken2"
            elif sample=="M5_DNA_subsample":
                if groupby_rank + "_" + data_type=="genus_rel":
                    pip="20_barrnap_spades_bowtie2_silva_blast-first-hit"
                elif groupby_rank + "_" + data_type=="genus_pa":
                    pip="20_sortmerna_transabyss_bwa_silva_kraken2"
                elif groupby_rank + "_" + data_type=="species_rel":
                    pip="10_unsorted_transabyss_bwa_silva_kraken2"
                elif groupby_rank + "_" + data_type=="species_pa":
                    pip="15_sortmerna_idba-ud_bwa_silva_kraken2"
            elif sample=="M6_DNA_subsample":
                if groupby_rank + "_" + data_type=="genus_rel":
                    pip="20_barrnap_rnaspades_bowtie2_silva_kraken2"
                elif groupby_rank + "_" + data_type=="genus_pa":
                    pip="20_barrnap_transabyss_bwa_silva_kraken2"
                elif groupby_rank + "_" + data_type=="species_rel":
                    pip="10_barrnap_transabyss_bwa_silva_kraken2"
                elif groupby_rank + "_" + data_type=="species_pa":
                    pip="10_barrnap_idba-ud_bwa_silva_kraken2"
            #### We substract the reads occuring in the filtration and extraction
            #### control from the samples
            ##### We're converting counts back to absolute and substract absolute numbers of reads of the negative controls
            readnum=sample_reads[sample]
            readnum_df=master_dfs_rel[sample]*readnum
            ##### Substract specific negative control pipelines
            master_dfs_rel_sub[sample]=readnum_df.sub((master_dfs_neg_rel["M_Neg_DNA"]\
                [pip]*neg_readnum + master_dfs_neg_rel["M_Ext_DNA"][pip]*ext_readnum), axis=0)
            ### Convert counts below 0 to 0 (happens if negative control contains more reads than original sample):
            master_dfs_rel_sub[sample][master_dfs_rel_sub[sample] < 0] = 0
            ### Convert counts back to relative
            master_dfs_rel_sub[sample]=master_dfs_rel_sub[sample]/master_dfs_rel_sub[sample].sum()
            #### And since this substraction overrides the expected, we replace it with the old expected:
            master_dfs_rel_sub[sample]["expected"]=master_dfs_rel[sample]["expected"]


        ## 3.5 Generate master df with presence/absence data
        ## (0=not found, 1=found) for every sample and save in dic master_dfs_pa:
        ### Deepcopy the relative abundance master dfs:
        master_dfs_pa = copy.deepcopy(master_dfs_rel_sub)
        ### Replace all values above 0 with 1:
        for key, value in master_dfs_pa.items():
            value[value > 0] = 1



        # 4 Calculate metrics

        ## 4.1 Cut down master_dfs to only expected taxa and to false positive (FP) taxa
        ### Apply cutdown function on relative, and pa master_dfs
        rel_cutdown_expected = cutdown(master_dfs_rel_sub, "expected")
        pa_cutdown_expected = cutdown(master_dfs_pa, "expected")

        rel_cutdown_FP = cutdown(master_dfs_rel_sub, "FP")
        pa_cutdown_FP = cutdown(master_dfs_pa, "FP")


        ## 4.2 Use the rel_cutdown_expected dfs as metrics for rel abundances of each
        ##     expected taxon and add additional metrics
        metrics_reps=copy.deepcopy(rel_cutdown_expected) # Copy df under different name
        for sample in samples:
            FP_rel=rel_cutdown_FP[sample].sum(axis=0)
            TP=pa_cutdown_expected[sample].sum(axis=0)
            FP=pa_cutdown_FP[sample].sum(axis=0)
            metrics_reps[sample].loc['FP_rel'] = FP_rel
            metrics_reps[sample].loc['TP'] = TP
            metrics_reps[sample].loc['FP'] = FP



        # 5 Make final metrics_df and standardize rel abundances
        for rep in metrics_reps.keys():
            metrics_df=metrics_reps[rep].transpose()
            ## Standardize by replacing 0s by "impute" (3 orders of magnitude lower than lowest taxon) and taking the centered log ratio
            metrics_df_rel_std=pd.DataFrame(clr(multiplicative_replacement(metrics_df.iloc[:, 0:11].to_numpy(dtype="float"), delta=impute)),
                index=metrics_df.index, columns=metrics_df.columns[0:11])
            metrics_df_concat=pd.concat([metrics_df_rel_std, metrics_df.iloc[:, 11:]], axis=1)
            ## Save the df
            metrics_df_concat.to_csv(os.path.join(workdir, "{0}_{1}_{2}_metrics_df.csv".format(rep, groupby_rank, data_type)), index_label="subsample")
