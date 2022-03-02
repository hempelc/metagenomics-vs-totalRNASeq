#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 2 Mar 2022

# This script generates performance curves for various subsample sizes of DNA and RNA samples

import pandas as pd
import os
import numpy as np
from scipy.spatial.distance import euclidean
import plotly.graph_objs as go


# Parameters set manually
samples = ["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA"]
workdir = "/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_mock_samples_subsamples_curves/"
## Subsample read numbers
subsample_readnums=[20000, 40000]
## Indicate if you want to loop over all 4 combinations of genus/species and cell/gen (True/False)
looping=True
## If you set looping to False, then define what specific rank and datatype and database
## you want to process:
rank="genus"
database="silva"
dt_type="pa"


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


master_df={}
# Loop over all combinations
for sample in samples:
    for groupby_rank in groupby_rank_lst:
        for data_type in data_types:
            for db in db_lst:

                # Empty storage df
                df_mean_sd=pd.DataFrame({}, index=["mean", "mean+sd", "mean-sd"])

                # Loop over subsample size
                for subsample_readnum in subsample_readnums:
                    # Import data
                    file=os.path.join(workdir, str(subsample_readnum), "{0}_{1}_{2}_{3}_metrics_df.csv".format(sample, db, groupby_rank, data_type))
                    df=pd.read_csv(file, index_col=0)
                    # Drop unwanted columns
                    if data_type=="rel":
                        df=df.drop(df.loc[:,'TP':'FP'], axis=1)
                    elif data_type=="pa":
                        df=df.loc[:, ["TP","FP"]]
                    # Separate expected from df
                    df_no_exp=df.drop(["expected"], axis=0)
                    exp=df.loc["expected"]
                    # Calculate euc dist for each pipeline
                    for index,row in df_no_exp.iterrows():
                        df_no_exp.loc[index,'euc_dist'] = euclidean(row, exp)
                    # Calculate mean and sd and transform into list
                    mean = np.mean(df_no_exp['euc_dist'])
                    sd = np.std(df_no_exp['euc_dist'])
                    mean_sd_list=[mean, mean+sd, mean-sd]
                    # Add to df and add df to master dic
                    df_mean_sd[str(subsample_readnum)]=mean_sd_list
                    master_df["{0}_{1}_{2}_{3}".format(sample, db, groupby_rank, data_type)]=df_mean_sd

# Plot (taken from https://plotly.com/python/continuous-error-bars/)
fig=go.Figure()
for sample in master_df.keys():
    mean=master_df[sample].loc["mean"]
    mean_plus_sd=master_df[sample].loc["mean+sd"]
    mean_minus_sd=master_df[sample].loc["mean-sd"]
    fig.add_trace(
        go.Scatter(
            name=sample,
            x=subsample_readnums,
            y=mean,
            line=dict(color='rgb(0,100,80)'),
            mode='lines'))
    fig.add_trace(
        go.Scatter(
            x=subsample_readnums+subsample_readnums[::-1], # x, then x reversed
            y=list(mean_plus_sd)+list(mean_minus_sd)[::-1], # upper, then lower reversed
            fill='toself',
            fillcolor='rgba(0,100,80,0.2)',
            line=dict(color='rgba(255,255,255,0)'),
            hoverinfo="skip",
            showlegend=False))
fig.show()
