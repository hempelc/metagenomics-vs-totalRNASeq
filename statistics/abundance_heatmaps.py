import pandas as pd
import numpy as np
from scipy.special import softmax
import math
import plotly.express as px
from scipy.spatial.distance import euclidean
import os
import copy
import logging

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')

# Parameters set manually
## Full path to directory that contains samples:
workdir="/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_mock_samples/"
## List of DNA and RNA mock community samples, replicates of 3; must equal names of directories in workdir that
## contain each sample's pipeline results:
samples = ["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA"]
types=["DNA", "RNA"]
reps=["M4", "M5", "M6"]

## Set if you want to loop over all result combinations of parameters in script
## "processing_and_metrics.py" and all metrics (True or False)
looping=True
## If you set looping to False, then define what specific combination
## and metrics you want to process:
### ("gen_genus", "gen_species")
combinations="gen_species"
metrics="rel"

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

for combination in combinations:
    for metr in metrics:
        master_df=pd.DataFrame()
        for sample in samples:
            # Import metrics df
            file=os.path.join(workdir, "metrics_" + combination, sample + "_metrics_df.csv")
            metrics_df=pd.read_csv(file, index_col=0)
            # Drop unwanted columns
            if metr=="rel":
                metrics_df=metrics_df.drop(metrics_df.loc[:,'TP':'classifier'], axis=1)
            elif metr=="pa":
                metrics_df=metrics_df.loc[:, "TP":"FP"]
            # Separate expected from df
            metrics_df_no_exp=metrics_df.drop(["expected"], axis=0)
            exp=metrics_df.loc["expected"]

            #  Euclidean distances between pipelines and expected:
            ## Calculate euc dist for each pipeline
            for index,row in metrics_df_no_exp.iterrows():
                metrics_df_no_exp.loc[index,'euc_dist'] = euclidean(row, exp)

            # Pick best pipeline metrics based on Euclidean distance
            best=metrics_df_no_exp.sort_values('euc_dist').iloc[0]\
                .rename("_".join(sample.split("_")[::-1]))
            master_df=master_df.append(best)
        # Add expected and sort
        master_df=master_df.append(exp).fillna(0).sort_values('expected', axis=1, ascending=False)
        master_df=master_df.reindex(["DNA_M4", "DNA_M5", "DNA_M6", "expected", "RNA_M4", "RNA_M5", "RNA_M6"])
        master_df=master_df[[c for c in master_df if c!='euc_dist'] + ['euc_dist']]

        # Normalize abundances from 0 to 1
        master_df_norm_w_exp=master_df.apply(diff_norm)


        # Plot heatmap
        combination_short=combination.replace("gen_", "")
        if metr=="pa":
            text_auto=True
            non_normalized_col=[px.colors.sequential.Viridis[x] for x in [0,5,9]]
        elif metr=="rel":
            text_auto=".2f"
            non_normalized_col=px.colors.sequential.Viridis
        ## With expected
        ### Normalized
        fig=px.imshow(master_df_norm_w_exp, color_continuous_scale=["#d80054", "#51a9ff"], text_auto=".2f")
        fig.show()
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_norm_w_exp_" + combination_short + "_" + metr + ".png"))
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_norm_w_exp_" + combination_short + "_" + metr + ".svg"))

        ### Non-normalized
        fig=px.imshow(master_df.drop(['euc_dist'], axis=1), color_continuous_scale=non_normalized_col, text_auto=text_auto)
        fig.show()
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_w_exp_" + combination_short + "_" + metr + ".png"))
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_w_exp_" + combination_short + "_" + metr + ".svg"))
        #### Plot eucdist separately on different colourscale to make it distinguishable
        eucdist_df=copy.deepcopy(master_df)
        eucdist_df["random"]=np.random.randint(master_df["euc_dist"].max(), size=len(master_df))
        fig=px.imshow(eucdist_df.loc[:, "euc_dist":"random"], color_continuous_scale=px.colors.sequential.Reds_r[1:], text_auto=".2f")
        fig.show()
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_w_exp_" + combination_short + "_" + metr + "_euc_dist.png"))
        fig.write_image(os.path.join(exportdir, "abundance_heatmap_w_exp_" + combination_short + "_" + metr + "_euc_dist.svg"))
