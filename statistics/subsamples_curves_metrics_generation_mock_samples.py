#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 2 Mar 2022

# This script processes pipeline data from multiple replicates of mock community
# samples and exports a metrics table. It's a modified version of the script
# "metrics_generation_mock_samples.py" that is adapted for subsampled DNA and RNA
# samples that were subsampled at different depths.
# Requires a specific directory structure to work.

import pandas as pd
import glob
import os
import copy
import logging
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters set manually
## Full path to directory that contains samples
workdir = "/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_mock_samples_subsamples_curves/"
## Lists of DNA and RNA mock community samples, replicates of 3 plus filtration controls (Neg) and
## extraction controls (Ext); must equal names of directories in workdir that
## contain each sample's pipeline results:
samples = ["M4_DNA", "M5_DNA", "M6_DNA", "M4_RNA", "M5_RNA", "M6_RNA"]
neg_samples=["M_Neg_DNA", "M_Ext_DNA", "M_Neg_RNA", "M_Ext_RNA"]
## Dic containing number of reads per negative control
neg_sample_reads={"M_Neg_DNA": 682, "M_Neg_RNA": 5200,
    "M_Ext_DNA":  445, "M_Ext_RNA": 2672}
## List that indicates at what number of reads was subsampled
subsample_readnums=[1000]
## Indicate if you want to loop over all combinations of genus/species and silva/ncbi and rel/pa (True/False)
looping=False
## If you set looping to False, then define what specific rank, datatype, and database
## you want to process:
rank="genus"
database="silva"
dt_type="rel"


# Parameters set automatically
if looping:
    groupby_rank_lst=["genus", "species"]
    #db_lst=["ncbi", "silva"]
    db_lst=["silva"]
    data_types=["rel", "pa"]
else:
    groupby_rank_lst=[rank]
    db_lst=[database]
    data_types=[dt_type]
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



# 1 Import excel file containing which pipeline was used for each subsample
pipelines_df=pd.read_excel(os.path.join(workdir, "rerun_dna_subsamples.xlsx"))


# 2 Make expected mock community df
expected_dic = {} # Empty dic that will eventually contain expected species and their abundances

### To calculate the  expected reads for the taxa in the mock community,
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
## abundance of each taxon in the mock community, either based on SSU
## genes or cells. While we're at it, we determine the lowest relative
## abundance, lower it 3 orders of magnitude, and use the result to impute 0s later
expected_df["rel_abun"] = rel_abun_gen
impute=min(rel_abun_gen)*1e-03


# 3 Read in all pipeline results for every sample as df and
#   add all taxa from all dfs to "all_taxa" list (needed to generate master dfs that
#   contain counts of all taxa from all samples and the expected community):
## Loop over all combinations
for subsample_readnum in subsample_readnums:
    for groupby_rank in groupby_rank_lst:
        for data_type in data_types:
            for db in db_lst:

                master_dfs_raw = {} # Empty dic that will eventually contain all samples' raw pipelines output
                all_taxa = [] # Empty list that will eventually contain all taxa that appear in all samples

                for sample in samples:
                    subsample_dir=os.path.join(str(subsample_readnum), sample)
                    ## Make a list for all file names in sample dic:
                    sample_files = glob.glob(os.path.join(workdir, subsample_dir, "*", db, groupby_rank + "_" + data_type, "*.txt"))
                    ## Make a dic that will eventually contain all pipeline dfs and set the first entry to expected community:
                    sample_dfs = {"expected": expected_df}
                    ## For each file in the sample dic
                    for file in sample_files:
                        #### Read in file as pandas df, fill NaN with "NA", "Unknown" by "NA", and fix one taxonomic misambiguation
                        df = pd.read_table(file).fillna("NA").replace("Lactobacillus",
                            r"Limosilactobacillus", regex=True).replace("Unknown", "NA").replace("-", r"", regex=True)
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
                        ### Edit file name so that we can name dfs based on their file name=subsample
                        sample_name = file.split("/")[8]
                        ### Add df_agg to the sample_dfs dic with key=subsample
                        sample_dfs[sample_name] = df_agg
                    ### Save sample_df in dic master_dfs_raw:
                    master_dfs_raw[sample] = sample_dfs

                # Repeat the code above for negative control samples - these were not subsampled,
                # so we import all files and select the pipelines that need to be substracted later
                master_dfs_neg_raw = {} # Empty dic that will eventually contain all samples' raw pipelines output
                for neg_sample in neg_samples:
                    ## Make a list for all file names in sample dic:
                    neg_sample_files = glob.glob(os.path.join(workdir, str(subsample_readnum), neg_sample, "*.txt"))
                    ## For each file in the sample dic
                    neg_sample_dfs={}
                    for file in neg_sample_files:
                        #### Read in file as pandas df, fill NaN with "NA", "Unknown" by "NA", and fix one taxonomic misambiguation
                        df = pd.read_table(file).fillna("NA").replace("Lactobacillus",
                            r"Limosilactobacillus", regex=True).replace("Unknown", "NA").replace("-", r"", regex=True)
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
                        neg_sample_dfs[pipeline_name] = df_agg
                    ### Save sample_df in dic master_dfs_raw:
                    master_dfs_neg_raw[neg_sample] = neg_sample_dfs


                # 4 Generate master dfs with relative read counts

                ## 4.1 Add all expected taxa to list in case they're not picked up by any pipeline
                ##     and drop duplicates:
                all_taxa.extend(expected_df[groupby_rank].tolist())
                unique_taxa = list(set(all_taxa))

                ## 4.2 Generate master df with relative read counts for every sample

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

                ## 4.3 Substract controls from samples
                master_dfs_rel_sub={}
                neg_readnum=neg_sample_reads["M_Neg_DNA"]
                ext_readnum=neg_sample_reads["M_Ext_DNA"]
                ### Define which pipeline from negative controls needs to be substracted from samples
                for sample in samples:
                    #### Define used pipeline
                    pip = list(pipelines_df[(pipelines_df["db"] == db) \
                        & (pipelines_df["sample"] == sample)]["{0}_{1}"\
                        .format(groupby_rank, data_type)])[0].replace("-", "_")\
                        .replace("idba_", "idba-").replace("ncbi_nt", "ncbi-nt")\
                        .replace("blast_first_hit", "blast-first-hit")\
                        .replace("blast_filtered", "blast-filtered")
                    ### We substract the reads occuring in the filtration and extraction
                    ### control from the samples:
                    #### We're converting counts back to absolute and substract absolute numbers of reads of the negative controls
                    readnum_df=master_dfs_rel[sample]*subsample_readnum
                    master_dfs_rel_sub[sample]=readnum_df.sub((master_dfs_neg_rel["M_Neg_DNA"]\
                        [pip]*neg_readnum + master_dfs_neg_rel["M_Ext_DNA"][pip]*ext_readnum), axis=0)
                    ### Convert counts below 0 to 0 (happens if negative control contains more reads than original sample):
                    master_dfs_rel_sub[sample][master_dfs_rel_sub[sample] < 0] = 0
                    ### Convert counts back to relative
                    #### Loop over columns and don't convert to relative if column
                    #### sum == 0, otherwise you get an error:
                    for col in master_dfs_rel_sub[sample].columns:
                        if master_dfs_rel_sub[sample][col].sum()==0:
                            continue
                        else:
                            master_dfs_rel_sub[sample][col]=master_dfs_rel_sub[sample][col]\
                                /master_dfs_rel_sub[sample][col].sum()
                    ### And since this substraction overrides the expected, we replace it with the old expected:
                    master_dfs_rel_sub[sample]["expected"]=master_dfs_rel[sample]["expected"]


                ## 4.4 Generate master df with presence/absence data
                ## (0=not found, 1=found) for every sample and save in dic master_dfs_pa:
                ### Deepcopy the relative abundance master dfs:
                master_dfs_pa = copy.deepcopy(master_dfs_rel_sub)
                ### Replace all values above 0 with 1:
                for key, value in master_dfs_pa.items():
                    value[value > 0] = 1



                # 5 Calculate metrics

                ## 5.1 Cut down master_dfs to only expected taxa and to false positive (FP) taxa
                ### Apply cutdown function on relative, and pa master_dfs
                rel_cutdown_expected = cutdown(master_dfs_rel_sub, "expected")
                pa_cutdown_expected = cutdown(master_dfs_pa, "expected")

                rel_cutdown_FP = cutdown(master_dfs_rel_sub, "FP")
                pa_cutdown_FP = cutdown(master_dfs_pa, "FP")


                ## 5.2 Use the rel_cutdown_expected dfs as metrics for rel abundances of each
                ##     expected taxon and add additional metrics
                metrics_reps=copy.deepcopy(rel_cutdown_expected) # Copy df under different name
                for sample in samples:
                    FP_rel=rel_cutdown_FP[sample].sum(axis=0)
                    TP=pa_cutdown_expected[sample].sum(axis=0)
                    FP=pa_cutdown_FP[sample].sum(axis=0)
                    metrics_reps[sample].loc['FP_rel'] = FP_rel
                    metrics_reps[sample].loc['TP'] = TP
                    metrics_reps[sample].loc['FP'] = FP



                # 6 Make final metrics_df and Standardize rel abundances
                ## Add info on tools used in each step for each pipeline
                for rep in metrics_reps.keys():
                    ### Make the df
                    metrics_df=metrics_reps[rep].transpose()
                    ### Standardize by replacing 0s by "impute" (3 orders of
                    ### magnitude lower than lowest taxon) and taking the centered
                    ### log ratio. If entire row sum == 0, then
                    ### multiplicative_replacement doesn't work, in which case
                    ### we drop these:
                    metrics_df_rel_std=pd.DataFrame({}, index=metrics_df.index, columns=metrics_df.columns[0:11])
                    for index,row in metrics_df.iterrows():
                        if row.sum()==0:
                            continue
                            # Or: replace 0s manually
                            #metrics_df_rel_std.loc[index]=clr(row[0:11].replace(0, impute))
                        else:
                            metrics_df_rel_std.loc[index]=clr(multiplicative_replacement(row[0:11].to_numpy(dtype="float"), delta=impute))
                    metrics_df_rel_std=metrics_df_rel_std.dropna()
                    metrics_df_TP_FP=metrics_df.loc[metrics_df_rel_std.index].iloc[:, 11:]
                    metrics_df_concat=pd.concat([metrics_df_rel_std, metrics_df_TP_FP], axis=1)
                    ### Save the df
                    metrics_df_concat.to_csv(os.path.join(workdir, str(subsample_readnum), "{0}_{1}_{2}_{3}_metrics_df.csv".format(rep, db, groupby_rank, data_type)), index_label="subsample")
