#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 21 Sep 2021

# This script processes the output from the script "regression_pca_clustering.py"

import pandas as pd
import plotly.express as px
import logging
import os
import numpy as np

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters manual setting
workdir="/Users/christopherhempel/Desktop/pipeline_results_mock_community"

# Parameters auto setting
levels=["genus_cell", "genus_gen", "species_cell", "species_gen"]
metr_lst=["all", "rel", "pa"]
steps=["type", "trimming_score", "rRNA_sorting_tool", "assembly_tool", "mapper", "database", "classifier"]

# Import dfs and concatenate them, creating columsn to separate them:
master_df=pd.DataFrame()
for level in levels:
    #level=levels[1]
    statsdir=os.path.join(workdir, "results_{0}/stats_exports".format(level))
    for metr in metr_lst:
        #metr=metr_lst[0]
        df_eucdist=pd.read_csv(os.path.join(statsdir, "{0}_euc_dist_steps.csv".format(metr)))
        df_pval=pd.read_csv(os.path.join(statsdir, "{0}_pvalues_tools.csv".format(metr)))
        # Merge dfs
        df_merged=pd.merge(df_eucdist, df_pval, on="tool")
        # Set one column with combinations:
        df_merged["combination"]=[level + "_" + metr]*len(df_merged)
        # Set two columns to "no" by default:
        df_merged["mean_euc_dist_lowest"]=["no"]*len(df_merged)
        df_merged["min_euc_dist_lowest"]=["no"]*len(df_merged)
        df_merged["euc_dist_category"]=["no"]*len(df_merged)
        for step in steps:
            # Replace "mean_euc_dist_lowest" value of lowest mean tool with "yes" for each step
            step_lowest_mean=df_merged[df_merged['tool'].str.contains(step)]["mean_euc_dist"].min()
            df_merged.loc[df_merged["mean_euc_dist"] == step_lowest_mean , 'mean_euc_dist_lowest'] = "yes"
            # Replace "min_euc_dist_lowest" value of lowest min tool with "yes" for each step
            step_lowest_min=df_merged[df_merged['tool'].str.contains(step)]["min_euc_dist"].min()
            df_merged.loc[df_merged["min_euc_dist"] == step_lowest_min , 'min_euc_dist_lowest'] = "yes"
        master_df=pd.concat([master_df, df_merged])
        master_df.loc[master_df["combination"]=="genus_gen_all"]

# Make col for significance of p-value
master_df.loc[master_df["p-value"] <= 0.001 , 'significance'] = "p <= 0.001"
master_df.loc[master_df["p-value"] > 0.001 , 'significance'] = "0.05 > p > 0.01"
master_df.loc[master_df["p-value"] > 0.05, 'significance'] = "p > 0.05"

# Make col for category of euc dist
master_df.loc[(master_df["mean_euc_dist_lowest"] == "yes") & (master_df["min_euc_dist_lowest"] == "yes") , 'euc_dist_category'] = "both"
master_df.loc[(master_df["mean_euc_dist_lowest"] == "no") & (master_df["min_euc_dist_lowest"] == "no") , 'euc_dist_category'] = "none"
master_df.loc[(master_df["mean_euc_dist_lowest"] == "yes") & (master_df["min_euc_dist_lowest"] == "no") , 'euc_dist_category'] = "mean"
master_df.loc[(master_df["mean_euc_dist_lowest"] == "no") & (master_df["min_euc_dist_lowest"] == "yes") , 'euc_dist_category'] = "min"


# Graph v1
fig = px.scatter(master_df, x="combination", y="tool",
	         size="min_euc_dist", color="significance",
             hover_name="min_euc_dist", size_max=28, height=1250, width=1050,
             symbol="min_euc_dist_lowest", symbol_sequence=["circle", "circle-dot"],
             category_orders={"tool": master_df["tool"].to_list()},
             color_discrete_sequence=["red", "lightgrey", "salmon"])
fig.update_traces(marker=dict(line=dict(width=4,color='black')), selector=dict(mode='markers'))
fig.write_image("/Users/christopherhempel/Desktop/bubbleplot.svg")
fig.write_image("/Users/christopherhempel/Desktop/bubbleplot.png")


# Graph v2
# Invert pval for graph:
master_df["p-value"]=1-master_df["p-value"]
fig = px.scatter(master_df, x="combination", y="tool",
	         size="p-value", color="euc_dist_category",
             hover_name="min_euc_dist", size_max=28, height=1350, width=1050,
             #symbol="mean_euc_dist_lowest", symbol_sequence=["circle", "circle-dot"],
             category_orders={"tool": master_df["tool"].to_list()},
             color_discrete_sequence=["lightgrey", "#00ab56", "#dac751", "#737cd0"])
fig.update_traces(marker=dict(line=dict(width=0.5,color='black')), selector=dict(mode='markers'))
fig.show()
fig.write_image("/Users/christopherhempel/Desktop/bubbleplot_col2.svg")


["#b99bcf",
"#eac221",
"#8c3cc2"]
