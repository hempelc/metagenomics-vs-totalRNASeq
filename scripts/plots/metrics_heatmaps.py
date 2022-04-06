#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 9 Mar 2021

# This script processes pipeline data from multiple replicates of mock community
# samples and compares the respective metrics in heatmaps

import pandas as pd #v1.3.5
import numpy as np #v1.21.3
import math
import plotly.express as px #v5.5.0
import os
import copy
import logging
from scipy.stats import ttest_rel #v1.7.3
from scipy.spatial.distance import euclidean #v1.7.3

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')

# Parameters set manually
## Full path to directory that contains samples:
workdir="/Users/christopherhempel/Desktop/pipeline_results_coverage/"
## List of DNA and RNA mock community samples, replicates of 3; must equal names of directories in workdir that
## contain each sample's pipeline results:
samples = ["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA"]
types=["DNA", "RNA"]
reps=["M4", "M5", "M6"]
# Include silva, ncbi-nt, or both?
db='silva'
## Set if you want to loop over all result combinations of parameters in script
## "processing_and_metrics.py" and all metrics (True or False)
looping=True
## If you set looping to False, then define what specific combination
## and metrics you want to process:
### ("gen_genus", "gen_species")
combinations=["gen_genus"]
metrics=["rel"]

# Parameters set automatically
## Set lists for aggregation, combinations, and metrics for looping:
if looping:
    combinations=["gen_genus", "gen_species"]
    metrics=["rel", "pa"]

## Make plot export directory:
exportdir=os.path.join(workdir, "stats_summary_plots", "abundance_heatmaps")
if not os.path.exists(exportdir):
    os.mkdir(exportdir)


# Function for df normalization
def diff_norm(series):
    diff=abs(series-series["expected"])
    return (diff-diff.min())/(diff.max()-diff.min())

best_pips_master={}
for combination in combinations:
    for metr in metrics:
        master_df=pd.DataFrame()
        best_pips={}
        metrics_df_master=pd.DataFrame()
        for sample in samples:
            # Import metrics df
            file=os.path.join(workdir, "metrics_" + combination, sample + "_metrics_df.csv")
            metrics_df=pd.read_csv(file, index_col=0)
            # Subset df to wanted database
            if db=='silva':
                metrics_df=metrics_df[metrics_df['database']!='ncbi-nt']
            elif db=='ncbi':
                metrics_df=metrics_df[metrics_df['database']!='silva']
            metrics_df["sample"]=[sample]*len(metrics_df)
            metrics_df_master=pd.concat([metrics_df_master, metrics_df])
        # Drop unwanted columns
        if metr=="rel":
            metrics_df_master=metrics_df_master.drop(metrics_df_master.loc[:,'TP':'classifier'], axis=1)
        elif metr=="pa":
            metrics_df_master=metrics_df_master.loc[:, ["TP","FP","sample"]]
            sample_col=metrics_df_master["sample"]
        # Separate expected from df
        metrics_df_master_no_exp=metrics_df_master.drop(["expected"], axis=0)
        exp=metrics_df_master.drop(["sample"], axis=1).loc["expected"].drop_duplicates().loc["expected"]
        #  Euclidean distances between pipelines and expected:
        metrics_df_master_no_exp=metrics_df_master_no_exp.reset_index()
        ## Calculate euc dist for each pipeline
        for index,row in metrics_df_master_no_exp.iterrows():
            metrics_df_master_no_exp.loc[index,'euc_dist'] = euclidean(row.drop(["pipeline", "sample"]), exp)
        metrics_df_master_no_exp=metrics_df_master_no_exp.set_index("pipeline")
        for sample in samples:
            # Pick best pipeline metrics based on Euclidean distance and collect
            # name of best pipelines in best_pips_master
            subdf=metrics_df_master_no_exp[metrics_df_master_no_exp['sample']==sample].drop(['sample'], axis=1)
            best=subdf.sort_values('euc_dist').iloc[0]
            best_dist=subdf.sort_values('euc_dist')['euc_dist'][0]
            best_pipelines=subdf[subdf['euc_dist']==best_dist].index
            master_df=master_df.append(best.rename("_".join(sample.split("_")[::-1])))
            best_pips[sample]=list(best_pipelines)
        best_pips_master[combination + "_" + metr]=best_pips
        # Add expected and sort
        master_df=master_df.append(exp).fillna(0).sort_values('expected', axis=1, ascending=False)
        master_df=master_df.reindex(["DNA_M4", "DNA_M5", "DNA_M6", "expected", "RNA_M4", "RNA_M5", "RNA_M6"])
        master_df=master_df[[c for c in master_df if c!='euc_dist'] + ['euc_dist']]

        # Until here, the code was based on the initial pipeline results based on all reads.
        # But we had to downsample DNA reads to make DNa and RNA comparable.
        # Therefore, we import data from subsampled DNA samples and replace DNA metrics above:
        # subsamples = ["M4_DNA_subsample", "M5_DNA_subsample", "M6_DNA_subsample"]
        # subsample_dir="/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_mock_samples_DNA_subsample"
        # for subsample in subsamples:
        #     subsample_file=os.path.join(subsample_dir, "{0}_{1}_{2}_metrics_df.csv".format(subsample, combination.replace("gen_",""), metr))
        #     subsample_metrics_df=pd.read_csv(subsample_file, index_col=0).drop("expected")
        #     if metr=="rel":
        #         subsample_metrics_df=subsample_metrics_df.drop(subsample_metrics_df.loc[:,'TP':'FP'], axis=1)
        #         subsample_metrics=subsample_metrics_df.mean(axis=0)
        #     elif metr=="pa":
        #         subsample_metrics_df=subsample_metrics_df.loc[:, ["TP","FP"]]
        #         subsample_metrics=subsample_metrics_df.mean(axis=0).round(0)
        #     subsample_metrics["euc_dist"]=euclidean(subsample_metrics, exp)
        #     name="_".join(subsample.replace("_subsample","").split("_")[::-1])
        #     master_df.loc[name]=subsample_metrics
        #
        # # Calculate p-values between DNA and RNA
        # dna_eucdist=master_df[master_df.index.str.contains('DNA')]["euc_dist"]
        # rna_eucdist=master_df[master_df.index.str.contains('RNA')]["euc_dist"]
        # pval=ttest_rel(dna_eucdist, rna_eucdist)[1]


        # Normalize abundances from 0 to 1
        master_df_norm_w_exp=master_df.apply(diff_norm)


        # Plot heatmap
        combination_short=combination.replace("gen_", "")
        if metr=="pa":
            non_normalized_col=[px.colors.sequential.Viridis[x] for x in [0,5,9]]
        elif metr=="rel":
            non_normalized_col=px.colors.sequential.Viridis
        ## With expected
        ### Normalized
        fig=px.imshow(master_df_norm_w_exp, color_continuous_scale=["#d80054", "#51a9ff"], text_auto=".2f", title="p=" + str(round(pval, 3)))
        fig.show()
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_norm_w_exp_" + combination_short + "_" + metr + ".png"))
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_norm_w_exp_" + combination_short + "_" + metr + ".svg"))

        ### Non-normalized
        fig=px.imshow(master_df.drop(['euc_dist'], axis=1), color_continuous_scale=non_normalized_col, text_auto=".2f", title="p=" + str(round(pval, 3)))
        fig.show()
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_w_exp_" + combination_short + "_" + metr + ".png"))
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_w_exp_" + combination_short + "_" + metr + ".svg"))
        #### Plot eucdist separately on different colourscale to make it distinguishable
        eucdist_df=copy.deepcopy(master_df)
        eucdist_df["random"]=np.random.randint(master_df["euc_dist"].max(), size=len(master_df))
        fig=px.imshow(eucdist_df.loc[:, "euc_dist":"random"], color_continuous_scale=px.colors.sequential.Reds_r[1:], text_auto=".2f", title="p=" + str(round(pval, 3)))
        fig.show()
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_w_exp_" + combination_short + "_" + metr + "_euc_dist.png"))
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_w_exp_" + combination_short + "_" + metr + "_euc_dist.svg"))
