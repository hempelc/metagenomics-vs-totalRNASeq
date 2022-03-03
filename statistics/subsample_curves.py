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
subsample_readnums=[20000, 40000, 60000, 78149]
## Indicate if you want to keep replicates separate
sep_reps=False
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
master_df_no_reps={}
# Loop over all combinations
for sample in samples:
    for groupby_rank in groupby_rank_lst:
        for data_type in data_types:
            for db in db_lst:

                # Empty storage df
                df_eucdist=pd.DataFrame()

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
                    # Add to df and add df to master dic
                    df_eucdist[str(subsample_readnum)]=df_no_exp['euc_dist']
                    master_df["{0}_{1}_{2}_{3}".format(sample, db, groupby_rank, data_type)]=df_eucdist

for na in ["RNA", "DNA"]:
    for groupby_rank in groupby_rank_lst:
        for data_type in data_types:
            for db in db_lst:
                eval_lvl="{0}_{1}_{2}_{3}".format(na, db, groupby_rank, data_type)
                reps=[x for x in master_df.keys() if eval_lvl in x]
                rep_df=pd.DataFrame()
                for rep in reps:
                    rep_df=pd.concat([rep_df, master_df[rep]])
                master_df_no_reps["{0}_{1}_{2}_{3}".format(na, db, groupby_rank, data_type)]=rep_df

# Plot (taken from https://plotly.com/python/continuous-error-bars/)
for groupby_rank in groupby_rank_lst:
    for data_type in data_types:
        for db in db_lst:
            fig=go.Figure()
            if sep_reps==True:
                dic=master_df
            else:
                dic=master_df_no_reps
            lvls=[x for x in dic.keys() if "{0}_{1}_{2}".format(db, groupby_rank, data_type) in x]
            for lvl in lvls:
                if "RNA" in lvl:
                    col='rgba(255,0,0,1)'
                    col_err='rgba(255,0,0,0.2)'
                else:
                    col='rgba(0,0,255,1)'
                    col_err='rgba(0,0,255,0.2)'
                # # Calculate mean and sd
                mean=dic[lvl].mean()
                sd=dic[lvl].std()
                mean_plus_sd=mean+sd
                mean_minus_sd=mean-sd
                fig.add_trace(
                    go.Scatter(
                        name=lvl,
                        x=subsample_readnums,
                        y=mean,
                        line=dict(color=col),
                        mode='lines'))
                fig.add_trace(
                    go.Scatter(
                        x=subsample_readnums+subsample_readnums[::-1], # x, then x reversed
                        y=list(mean_plus_sd)+list(mean_minus_sd)[::-1], # upper, then lower reversed
                        fill='toself',
                        fillcolor=col_err,
                        line=dict(color='rgba(255,255,255,0)'),
                        hoverinfo="skip",
                        showlegend=False))
            fig.update_layout(
            title="{0}_{1}_{2}".format(db, groupby_rank, data_type),
            xaxis_title="Number of reads",
            yaxis_title="Euclidean distance to reference",
            legend_title="Sample")
            fig.show()
