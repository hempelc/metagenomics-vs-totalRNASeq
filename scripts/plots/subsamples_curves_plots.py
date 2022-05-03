#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 2 Mar 2022

# This script generates accuracy curves for various subsample sizes of DNA and RNA samples

import os
import pandas as pd #v1.3.5
import numpy as np #v1.21.3
import plotly.graph_objs as go #v5.5.0
from scipy.spatial.distance import euclidean #v1.7.3
import statsmodels.api as sm #v0.13.0
from statsmodels.formula.api import ols #v0.13.0
from statsmodels.stats.anova import anova_lm #v0.13.0
from scipy.stats import ttest_rel #v1.7.3


# Parameters set manually
samples = ["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA"]
workdir = "/Users/christopherhempel/Desktop/pipeline_results_coverage/subsample_curves/"
## Subsample read numbers
subsample_readnums=[1000, 2500, 5000, 10000, 20000, 40000, 60000, 78149, 94633,
    120144, 200000, 300000, 400000, 500000, 600000, 644634, 669382, 817619]
## Indicate if you want to keep replicates separate
sep_reps=False
## If sep_reps=False, indicate if you want to show separate replcates as gray lines
show_reps=True
## Indicate if you don't want to show euc_dist for up to 40k reads
fourtyk=False
## Indicate if you want to loop over all 4 combinations of genus/species and cell/gen (True/False)
looping=True
## If you set looping to False, then define what specific rank and datatype and database
## you want to process:
rank="species"
database="silva"
dt_type="rel"

# Adding regression curves for data until 120,144 reads
def add_reg(mean, colour, reg_reads):
    mean_sub=mean.loc[:reg_reads]
    y_reg=mean_sub
    x_reg=mean_sub.index
    fit_results=sm.OLS(y_reg, sm.add_constant(x_reg), missing="drop").fit()
    const = round(fit_results.params["const"], 2)
    coeff = round(fit_results.params["x1"], 8)
    fun="{0}+{1}*x".format(const, coeff)
    y_pred = fit_results.predict()
    fig.add_trace(
        go.Scatter(
            name=fun,
            x=x_reg,
            y=y_pred,
            line=dict(color=colour, dash='dash'),
            mode='lines',
            showlegend=False))


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

                # List with names of subsamples
                subsamples=[sample + "_subsample_" + str(x) for x in range(1,11)]
                # Empty storage df
                df_eucdist=pd.DataFrame({}, index=subsamples)

                # Loop over subsample size
                for subsample_readnum in subsample_readnums:
                    # Import data
                    file=os.path.join(workdir, "subsamples_" + str(subsample_readnum) + "_coverage", "{0}_{1}_{2}_{3}_metrics_df.csv".format(sample, db, groupby_rank, data_type))
                    if not os.path.isfile(file):
                        continue
                    df=pd.read_csv(file, index_col=0)
                    # Drop unwanted columns
                    if data_type=="rel":
                        df=df.drop(df.loc[:,'TP':'FP'], axis=1)
                    elif data_type=="pa":
                        df=df.loc[:, ["TP","FP"]]
                    # Separate expected from df
                    df_no_exp=df.drop(["expected"], axis=0)
                    # If df is empty, then skip
                    if df_no_exp.empty:
                        continue
                    exp=df.loc["expected"]
                    # Calculate euc dist for each pipeline
                    for index,row in df_no_exp.iterrows():
                        df_no_exp.loc[index,'euc_dist'] = euclidean(row, exp)
                    # Add to df and add df to master dic
                    df_eucdist[subsample_readnum]=df_no_exp['euc_dist']
                master_df["{0}_{1}_{2}_{3}".format(sample, db, groupby_rank, data_type)]=df_eucdist

# Previous code was for every replicate separately, now we concatenate all
# replicates together and store concatenated results separately
for na in ["RNA", "DNA"]:
    for groupby_rank in groupby_rank_lst:
        for data_type in data_types:
            for db in db_lst:
                eval_lvl="{0}_{1}_{2}_{3}".format(na, db, groupby_rank, data_type)
                reps=[x for x in master_df.keys() if eval_lvl in x]
                rep_df=pd.DataFrame()
                for rep in reps:
                    rep_df=pd.concat([rep_df, master_df[rep]])
                master_df["{0}_{1}_{2}_{3}".format(na, db, groupby_rank, data_type)]=rep_df


# Make a plot for each evaluation level and add each sampel as layer in a for loop
# Also store data in dfs for following stats
for groupby_rank in groupby_rank_lst:
    for data_type in data_types:
        for db in db_lst:
            diffs=pd.DataFrame(index=['diff'])
            coeffs=pd.DataFrame(index=['coef'])
            fig=go.Figure()
            lvls=[x for x in master_df.keys() if "{0}_{1}_{2}".format(db, groupby_rank, data_type) in x]
            ## Sort by DNA and RNA
            lvls=[x for x in lvls if "DNA" in x if "M" in x]\
                +[x for x in lvls if "RNA" in x if "M" in x]\
                +[x for x in lvls if "DNA" in x if "M" not in x]\
                +[x for x in lvls if "RNA" in x if "M" not in x]
            if sep_reps:
                lvls=[x for x in lvls if "M" in x]
            ## Make empty df for model comparison
            model_comp_df=pd.DataFrame()
            for lvl in lvls:
                ## Regression calc based on # of reads
                if '4' in lvl:
                    reg_reads = 94633
                elif '5' in lvl:
                    reg_reads = 78149
                elif '6' in lvl:
                    reg_reads = 120144
                ## Cut down df to reads >=40k to exclude unuseful data
                if not fourtyk and "_rel" in lvl:
                    cols=[x for x in master_df[lvl].columns if x >= 40000]
                else:
                    cols=master_df[lvl].columns
                ## Set colours
                if "RNA" in lvl:
                    color='rgba(220,50,32,1)'
                    color_err='rgba(220,50,32,0.1)'
                else:
                    color='rgba(0,90,181,1)'
                    color_err='rgba(0,90,181,0.1)'
                sequencing_type=lvl.replace("_{0}_{1}_{2}".format(db, groupby_rank, data_type), '')
                if sequencing_type=="DNA":
                    name="Metagenomics"
                else:
                    name="Total RNA-Seq"
                ## Calculate mean and sd
                mean=master_df[lvl].loc[:, cols].mean().rename("accuracy")
                sd=pd.Series()
                for col in cols:
                    col_no_nan=master_df[lvl][col].dropna()
                    if len(col_no_nan)==1:
                        sd.loc[col]=0
                    else:
                        sd.loc[col]=col_no_nan.std()
                mean_plus_sd=mean+sd
                mean_minus_sd=mean-sd
                ## Generate regression curves and determine coefficients and constants
                mean_sub=mean.loc[:reg_reads]
                y_reg=mean_sub
                x_reg=mean_sub.index
                fit_results=sm.OLS(y_reg, sm.add_constant(x_reg), missing="drop").fit()
                const = fit_results.params["const"]
                coef = fit_results.params["x1"]
                # Save coefficients of replicates
                if "M" in lvl:
                    coeffs[lvl]=coef
                # Add to df for model comparison
                if "M" in lvl:
                    df_comp_rep=pd.DataFrame({"accuracy": mean_sub})
                    if "DNA" in lvl:
                        df_comp_rep["seqtype"]=[0]*len(df_comp_rep)
                    elif "RNA" in lvl:
                        df_comp_rep["seqtype"]=[1]*len(df_comp_rep)
                    model_comp_df=pd.concat([model_comp_df, df_comp_rep])
                # Add lines to plot (following https://plotly.com/python/continuous-error-bars/)
                if not sep_reps:
                    # Add grey lines
                    if "M" in lvl:
                        if show_reps:
                            fig.add_trace(
                                go.Scatter(
                                    name=name,
                                    x=list(cols),
                                    y=mean,
                                    line=dict(color='rgba(211,211,211,1)'),
                                    mode='lines',
                                    showlegend=False))
                    else:
                        # Add actual lines
                        fig.add_trace(
                            go.Scatter(
                                name=name,
                                x=list(cols),
                                y=mean,
                                line=dict(color=color),
                                mode='lines+markers'))
                        fig.add_trace(
                            go.Scatter(
                                x=list(cols)+list(cols)[::-1], # x, then x reversed
                                y=list(mean_plus_sd)+list(mean_minus_sd)[::-1], # upper, then lower reversed
                                fill='toself',
                                fillcolor=color_err,
                                line=dict(color='rgba(255,255,255,0)'),
                                hoverinfo="skip",
                                showlegend=False))
                        # Add regression line
                        add_reg(mean, "black", reg_reads)
                else:
                    fig.add_trace(
                        go.Scatter(
                            name=name,
                            x=list(cols),
                            y=mean,
                            line=dict(color=color),
                            mode='lines+markers'))
                    fig.add_trace(
                        go.Scatter(
                            x=list(cols)+list(cols)[::-1], # x, then x reversed
                            y=list(mean_plus_sd)+list(mean_minus_sd)[::-1], # upper, then lower reversed
                            fill='toself',
                            fillcolor=color_err,
                            line=dict(color='rgba(255,255,255,0)'),
                            hoverinfo="skip",
                            showlegend=False))
                    # Add regression line
                    add_reg(mean, "black", reg_reads)

            # Calculate p-values between DNA and RNA coefficients
            coeffs=coeffs.transpose()
            dna_coeffs=coeffs[coeffs.index.str.contains('DNA')]["coef"]
            rna_coeffs=coeffs[coeffs.index.str.contains('RNA')]["coef"]
            pval_coef=ttest_rel(dna_coeffs, rna_coeffs)[1]
            if pval_coef < 0.001:
                pval_coef_title="< 0.001"
            else:
                pval_coef_title="= {0}".format(round(pval_coef, 3))

            # Determine model performance with and without sequencing type
            # following https://nathancarter.github.io/how2data/site/how-to-compare-two-nested-linear-models/
            model_comp_df = model_comp_df.reset_index().rename(columns={'index':'readnum'})
            no_seq_model = ols('accuracy ~ readnum', data = model_comp_df).fit()
            seq_model = ols('accuracy ~ readnum + seqtype', data = model_comp_df).fit()
            pval_model=anova_lm(no_seq_model, seq_model)["Pr(>F)"][1]
            if pval_model < 0.001:
                pval_model_title="< 0.001"
            else:
                pval_model_title="= {0}".format(round(pval_model, 3))

            # Update and save figure
            fig.update_layout(title="p_seq " + pval_model_title + ", p_coef " + pval_coef_title,
                xaxis_title="Number of reads",
                yaxis_title="Euclidean distance to reference",
                legend_title="Sequencing type", template="simple_white")
            fig.update_yaxes(autorange="reversed")
            fig.show()
            fig.write_image(os.path.join(workdir, "subsample_curves_{0}_{1}_{2}.png".format(db, groupby_rank, data_type)), height=400, width=800)
            fig.write_image(os.path.join(workdir, "subsample_curves_{0}_{1}_{2}.svg".format(db, groupby_rank, data_type)), height=400, width=800)
