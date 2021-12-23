#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 22 Jul 2021

# This script processes the output from the script "processing_and_metrics.py"

import pandas as pd
import numpy as np
import plotly.express as px
import statsmodels.api as sm
import os
import copy
import logging
import pickle
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import euclidean
from scipy.stats import pointbiserialr
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from collections import Counter
from sklearn.decomposition import PCA
from numpy import array, linspace
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from scipy.signal import argrelextrema
import statsmodels.api as sm


# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters set manually
## Full path to directory that contains samples:
workdir="/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_mock_samples/"
## List of DNA and RNA mock community samples, replicates of 3; must equal names of directories in workdir that
## contain each sample's pipeline results:
samples = ["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA"]
## Set if you want to loop over all result combinations of parameters in script
## "processing_and_metrics.py" and all metrics (True or False)
looping=True
## Set lists for combinations and metrics for looping:
combinations=["cell_genus", "cell_species", "gen_genus", "gen_species"]
metr_list=["rel", "pa"]
## If you set looping to False, then define what specific combination and metrics
## you want to process:
### ("cell_genus", "cell_species", "gen_genus", "gen_species")
combination="cell_genus"
### ("rel", "pa")
metrics="rel"



# Parameters set automatically
if looping:
    dir_lst=[os.path.join(workdir, "metrics_" + x) \
        for x in combinations]
    metr_list=metr_list
else:
    dir_lst=[os.path.join(workdir, "metrics_" + combination)]
    metr_list=[metrics]

master_dfs={}
for dir in dir_lst:
    for metr in metr_list:
        for sample in samples:
            ## Make plot export directory:
            exportdir=os.path.join(dir, "stats_" + metr)
            if not os.path.exists(exportdir):
                os.mkdir(exportdir)

            # 1 Import data
            ## Import metrics df
            metrics_df=pd.read_csv(os.path.join(dir, sample + "_metrics_df.csv"), index_col=0)
            ## Convert trimming score column type into str
            metrics_df['trimming_score'] = metrics_df['trimming_score'].astype(str)
            ### Drop additional ones specified by parameter "metr" and set expected dummy
            ### (=ideal pipeline results):
            if metr not in ["rel", "pa"]:
                logging.critical("Parameter metr not in rel_metr, or pa_metr: metr=" + metr)
            if metr=="rel":
                print("Just relative metrics are used.")
                metrics_df=metrics_df.drop(["FP", "TP"], axis=1)
            elif metr=="pa":
                print("Just p/a metrics are used.")
                no_drops=["FP", "TP", 'type', 'trimming_score',
                    'rRNA_sorting_tool', 'assembly_tool', 'mapper', 'database', 'classifier']
                metrics_df=metrics_df.drop([x for x in metrics_df.columns if x not in no_drops], axis=1)
                ## Standardize (needed since TP and FP are not standardized)
                metrics_df_TP_FP_std=pd.DataFrame(StandardScaler().fit_transform(metrics_df.iloc[:, 0:2]), index=metrics_df.index, columns=["TP", "FP"])
                metrics_df=pd.concat([metrics_df_TP_FP_std, metrics_df.iloc[:, 2:]], axis=1)
            ## Separate expected from dfs
            metrics_df_no_exp=metrics_df.drop(["expected"], axis=0)
            exp=metrics_df.iloc[:, :-7].loc["expected"]
            only_metrics_df_no_exp=metrics_df.iloc[:, :-7].drop(["expected"], axis=0)
            only_tools_df_no_exp=metrics_df.iloc[:, -7:].drop(["expected"], axis=0)



            # 2 Euclidean distances between pipelines and expected (note: use standardized columns):
            ## Calculate euc dist for each pipeline
            for index,row in only_metrics_df_no_exp.iterrows():
                metrics_df_no_exp.loc[index,'euc_dist'] = euclidean(row, exp)
            #master_df=pd.concat([master_df_eucdist, metrics_df_no_exp])


            # 4 "Clustering", or more correctly, segmentation of euclidean distances
            segmentation_df=metrics_df_no_exp["euc_dist"].array.reshape(-1,1)

            ## 4.1 Evaluation of appropriate number of clusters using Kernel Density estimation
            ##     by performing a grid search to identify best value for bandwidth
            grid = GridSearchCV(KernelDensity(kernel = 'gaussian'),{'bandwidth': np.linspace(0.1, 2, 20)})
            grid_fit=grid.fit(segmentation_df)
            bandw=list(grid_fit.best_params_.values())[0]

            ## 4.2 Run the actual Kernel Density Estimation
            kde = KernelDensity(kernel='gaussian', bandwidth=bandw).fit(segmentation_df)
            s = linspace(min(metrics_df_no_exp["euc_dist"]),max(metrics_df_no_exp["euc_dist"]))
            e = kde.score_samples(s.reshape(-1,1))
            fig_kde=px.line(x=s, y=e)
            fig_kde.show()

            ## 4.3 Define minima and split data into clusters based on minima
            mi=argrelextrema(e, np.less)[0]
            print("For combination {0}_{1}_{2}, the best number of clusters is {3}."\
                .format(sample, dir.split("/")[-1], metr, len(mi)+1))
            clusters={}
            for i in range(0, len(mi)+1):
                if i==0:
                    clusters[i]=list(segmentation_df[segmentation_df < s[mi][i]])
                elif i==len(mi):
                    clusters[i]=list(segmentation_df[segmentation_df >= s[mi][i-1]])
                else:
                    clusters[i]=list(segmentation_df[(segmentation_df >= s[mi][i-1]) * (segmentation_df <= s[mi][i])])

            ## 4.4 Determine average of clusters
            cluster_aves=[]
            for cluster in clusters.keys():
                cluster_aves.append(sum(clusters[cluster])/len(clusters[cluster]))

            ## 4.4 Turn clusters into mergeable dataframe and add them to metrics df
            clusters_inv={}
            for cluster, eucdist_list in clusters.items():
                for eucdist in eucdist_list:
                    clusters_inv[eucdist]=[str(cluster)]
            cluster_df=pd.DataFrame(clusters_inv).transpose().reset_index()
            cluster_df.columns = ['euc_dist', 'cluster']
            metrics_df_no_exp=pd.merge(metrics_df_no_exp.reset_index(), cluster_df, on="euc_dist").set_index('pipeline')

            ## 4.5 Determine closest cluster to expected based on centroid
            ### Determine euclidean distance of each cluster centroid to expected:
            cluster_centroid_dist=[euclidean(i, exp) for i in cluster_aves]
            cluster_dist_df = pd.DataFrame({"cluster": [str(x) for x in list(range(0, len(cluster_centroid_dist)))],
                "centroid_euc_dist":cluster_centroid_dist})
            ### Merge pipeline info with cluster info:
            metrics_df_no_exp = pd.merge(metrics_df_no_exp.reset_index(), cluster_dist_df, on="cluster").set_index('pipeline')
            ### Determine closest cluster and cut df down to pipelines from closest cluster:
            best_cluster=cluster_dist_df.sort_values("centroid_euc_dist").reset_index().iloc[[0]]["cluster"].to_list()[0]
            best_cluster_df=metrics_df_no_exp.loc[metrics_df_no_exp["cluster"] == best_cluster]
            ### Count occurence of each tool for each step and save in dict:
            counts_dic={}
            for col in best_cluster_df.columns[-10:-3]:
                counts_dic[col]=Counter(best_cluster_df[col])
            ### Save the dict as pickle object, which is needed for the
            ### next step of code in the script "stats_summary_plot.py":
            with open(os.path.join(exportdir, sample + "_" + metr + "_closest_cluster_tool_counts.pkl"), 'wb') as f:
                pickle.dump(counts_dic, f)


            ## Calculate average and minimum euc dist for each tool in each step
            mean_euc_dist_lst=[]
            min_euc_dist_lst=[]
            tools=[]
            for step in metrics_df_no_exp.loc[:, 'type':'classifier'].columns:
                for tool in metrics_df_no_exp[step].unique():
                    ### Make sure trimming score 5 doesn't pick up 15 as well
                    if tool=="5":
                        tool="_5"
                    #### Cut down df to rows containing tool and take the averge and min euc_dist:
                    mean_euc_dist=metrics_df_no_exp.reset_index()[metrics_df_no_exp.reset_index()['pipeline']
                        .str.contains(tool)]["euc_dist"].mean()
                    min_euc_dist=metrics_df_no_exp.reset_index()[metrics_df_no_exp.reset_index()['pipeline']
                        .str.contains(tool)]["euc_dist"].sort_values().min()
                    ### Reverse changes
                    if tool=="_5":
                        tool="5"
                    tools.append(step + "_" + tool)
                    mean_euc_dist_lst.append(mean_euc_dist)
                    min_euc_dist_lst.append(min_euc_dist)

            ## Make and save df
            euc_dist_final_df=pd.DataFrame({"mean_euc_dist": mean_euc_dist_lst, "min_euc_dist": min_euc_dist_lst}, index=tools)
            euc_dist_final_df.to_csv(os.path.join(exportdir, sample + "_" + metr + "_euc_dist_steps.csv"), index_label="tool")

            ## Add df to collection
            master_dfs["{0}_{1}_{2}".format(sample, dir.split("/")[-1], metr)]=metrics_df_no_exp


            # 3 Correlation for separate replicates and sample type
            ## Turn tools into dummie variables
            dummies_df=pd.get_dummies(only_tools_df_no_exp)
            ## Calculate p-values using pearson's correlation coefficient
            pval_dic={}
            for col in dummies_df.columns[1:]:
                X=metrics_df_no_exp["euc_dist"].to_numpy()
                y=np.array(dummies_df[col])
                # logit_model_result=sm.Logit(y,X).fit()
                # pval=logit_model_result.summary2().tables[1]['P>|z|'][0]
                pval=pointbiseriar(X,y)[1]
                pval_dic[col]=pval
            pval_df=pd.DataFrame(pval_dic, index=["p-value"]).transpose()
            pval_df.to_csv(os.path.join(exportdir, sample + "_" + metr + "_pvalues_tools.csv"), index_label="tool")


# Correlation for separate replicates with combined sample type:
for dir in dir_lst:
    for metr in metr_list:
        for rep in ["M4", "M5", "M6"]:
            ## Make plot export directory:
            exportdir=os.path.join(dir, "stats_" + metr)
            if not os.path.exists(exportdir):
                os.mkdir(exportdir)

            rep_concat=pd.concat([master_dfs["{0}_DNA_{1}_{2}".format(rep, dir.split("/")[-1], metr)], master_dfs["{0}_RNA_{1}_{2}".format(rep, dir.split("/")[-1], metr)]])
            dummies_df=pd.get_dummies(rep_concat.iloc[:,-10:-3])
            ## Calculate p-values using pearson's correlation coefficient
            pval_dic={}
            for col in dummies_df.columns:
                X=rep_concat["euc_dist"].to_numpy()
                y=np.array(dummies_df[col])
                # logit_model_result=sm.Logit(y,X).fit()
                # pval=logit_model_result.summary2().tables[1]['P>|z|'][0]
                pval=pointbiserialr(X,y)[1]
                pval_dic[col]=pval
            pval_df=pd.DataFrame(pval_dic, index=["p-value"]).transpose()
            pval_df.to_csv(os.path.join(exportdir, "{0}_{1}_{2}_pvalues_tools.csv".format(rep, metr, dir.split("/")[-1])), index_label="tool")


# Correlation for combined replicates with separate sample type
for dir in dir_lst:
    for metr in metr_list:
        for type in ["DNA", "RNA"]:
            ## Make plot export directory:
            exportdir=os.path.join(dir, "stats_" + metr)
            if not os.path.exists(exportdir):
                os.mkdir(exportdir)

            rep_concat=pd.concat([master_dfs["M4_{0}_{1}_{2}".format(type, dir.split("/")[-1], metr)], master_dfs["M5_{0}_{1}_{2}".format(type, dir.split("/")[-1], metr)]\
                , master_dfs["M6_{0}_{1}_{2}".format(type, dir.split("/")[-1], metr)]])
            dummies_df=pd.get_dummies(rep_concat.iloc[:,-10:-3])
            ## Calculate p-values using pearson's correlation coefficient
            pval_dic={}
            for col in dummies_df.columns:
                X=rep_concat["euc_dist"].to_numpy()
                y=np.array(dummies_df[col])
                # logit_model_result=sm.Logit(y,X).fit()
                # pval=logit_model_result.summary2().tables[1]['P>|z|'][0]
                pval=pointbiserialr(X,y)[1]
                pval_dic[col]=pval
            pval_df=pd.DataFrame(pval_dic, index=["p-value"]).transpose()
            pval_df.to_csv(os.path.join(exportdir, "{0}_{1}_{2}_pvalues_tools.csv".format(type, metr, dir.split("/")[-1])), index_label="tool")

# Correlation for combined replicates with combined sample type
for dir in dir_lst:
    for metr in metr_list:
            ## Make plot export directory:
            exportdir=os.path.join(dir, "stats_" + metr)
            if not os.path.exists(exportdir):
                os.mkdir(exportdir)

            rep_concat=pd.concat([master_dfs["M4_DNA_{0}_{1}".format(dir.split("/")[-1], metr)]\
                , master_dfs["M5_DNA_{0}_{1}".format(dir.split("/")[-1], metr)]\
                , master_dfs["M6_DNA_{0}_{1}".format(dir.split("/")[-1], metr)]\
                , master_dfs["M4_RNA_{0}_{1}".format(dir.split("/")[-1], metr)]\
                , master_dfs["M5_RNA_{0}_{1}".format(dir.split("/")[-1], metr)]\
                , master_dfs["M6_RNA_{0}_{1}".format(dir.split("/")[-1], metr)]])
            dummies_df=pd.get_dummies(rep_concat.iloc[:,-10:-3])
            ## Calculate p-values using pearson's correlation coefficient
            pval_dic={}
            for col in dummies_df.columns:
                X=rep_concat["euc_dist"].to_numpy()
                y=np.array(dummies_df[col])
                # logit_model_result=sm.Logit(y,X).fit()
                # pval=logit_model_result.summary2().tables[1]['P>|z|'][0]
                pval=pointbiserialr(X,y)[1]
                pval_dic[col]=pval
            pval_df=pd.DataFrame(pval_dic, index=["p-value"]).transpose()
            pval_df.to_csv(os.path.join(exportdir, "{0}_{1}_pvalues_tools.csv".format(metr, dir.split("/")[-1])), index_label="tool")
