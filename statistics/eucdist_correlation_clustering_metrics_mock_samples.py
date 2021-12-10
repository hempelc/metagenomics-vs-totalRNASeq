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
from scipy.stats import pearsonr
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from collections import Counter
from sklearn.decomposition import PCA

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters set manually
## Full path to directory that contains samples:
workdir="/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_mock_samples/"
## Set if you want to loop over all result combinations of parameters in script
## "processing_and_metrics.py" and all metrics (True or False)
looping=True
## If you set looping to False, then define what specific combination and metrics
## you want to process:
### ("cell_genus", "cell_species", "gen_genus", "gen_species")
combination="cell_genus"
### ("all", "rel", "pa")
metrics="all"



# Parameters set automatically
if looping:
    dir_lst=[os.path.join(workdir, "metrics_" + x) \
        for x in ["cell_genus", "cell_species", "gen_genus", "gen_species"]]
    metr_list=["all", "rel", "pa"]
else:
    dir_lst=[os.path.join(workdir, "metrics_" + combination)]
    metr_list=[metrics]



for dir in dir_lst:
    for metr in metr_list:
        ## Make plot export directory:
        exportdir=os.path.join(dir, "stats_" + metr)
        if not os.path.exists(exportdir):
            os.mkdir(exportdir)

        # 1 Import data
        ## Import metrics df
        metrics_df=pd.read_csv(os.path.join(dir, "metrics_df.csv"), index_col=0)
        ## Convert trimming score column type into str
        metrics_df['trimming_score'] = metrics_df['trimming_score'].astype(str)
        ### Drop additional ones specified by parameter "metr" and set expected dummy
        ### (=ideal pipeline results):
        if metr not in ["all", "rel", "pa"]:
            logging.critical("Parameter metr not in all_metr, rel_metr, or pa_metr: metr=" + metr)
        if metr=="all":
            print("All metrics are used.")
        elif metr=="rel":
            print("Just relative metrics are used.")
            metrics_df=metrics_df.drop(["FP", "TP"], axis=1)
        elif metr=="pa":
            print("Just p/a metrics are used.")
            drops=["FP", "TP", 'type', 'trimming_score',
                'rRNA_sorting_tool', 'assembly_tool', 'mapper', 'database', 'classifier']
            metrics_df=metrics_df.drop([x for x in metrics_df.columns if x not in drops], axis=1)
        ## Separate metrics and tools
        only_metrics_df_with_exp=metrics_df.iloc[:, :-7]
        only_tools_df_with_exp=metrics_df.iloc[:, -7:]
        ## Standardize metrics
        only_metrics_df_with_exp_std=pd.DataFrame(StandardScaler().fit_transform(only_metrics_df_with_exp),
            columns=only_metrics_df_with_exp.columns, index=only_metrics_df_with_exp.index)
        ## Separate expected from dfs
        exp_std=only_metrics_df_with_exp_std.loc["expected"]
        only_metrics_df_no_exp_std=only_metrics_df_with_exp_std.drop(["expected"], axis=0)
        only_tools_df_no_exp=only_tools_df_with_exp.drop(["expected"], axis=0)
        ## Stick standardized metrics and tools back together
        metrics_df_std_no_exp=pd.concat([only_metrics_df_no_exp_std, only_tools_df_no_exp], axis=1)



        # 2 Euclidean distances between pipelines and expected (note: use standardized columns):
        ## Calculate euc dist for each pipeline
        for index,row in only_metrics_df_no_exp_std.iterrows():
            metrics_df_std_no_exp.loc[index,'euc_dist'] = euclidean(row, exp_std)

        ## Calculate average and minimum euc dist for each tool in each step
        mean_euc_dist_lst=[]
        min_euc_dist_lst=[]
        tools=[]
        for step in metrics_df_std_no_exp.loc[:, 'type':'classifier'].columns:
            for tool in metrics_df_std_no_exp[step].unique():
                ### Make sure trimming score 5 doesn't pick up 15 as well
                if tool=="5":
                    tool="_5"
                #### Cut down df to rows containing tool and take the averge and min euc_dist:
                mean_euc_dist=metrics_df_std_no_exp.reset_index()[metrics_df_std_no_exp.reset_index()['pipeline']
                    .str.contains(tool)]["euc_dist"].mean()
                min_euc_dist=metrics_df_std_no_exp.reset_index()[metrics_df_std_no_exp.reset_index()['pipeline']
                    .str.contains(tool)]["euc_dist"].sort_values().min()
                ### Reverse changes
                if tool=="_5":
                    tool="5"
                tools.append(step + "_" + tool)
                mean_euc_dist_lst.append(mean_euc_dist)
                min_euc_dist_lst.append(min_euc_dist)

        ## Make and save df
        euc_dist_final_df=pd.DataFrame({"mean_euc_dist": mean_euc_dist_lst, "min_euc_dist": min_euc_dist_lst}, index=tools)
        euc_dist_final_df.to_csv(os.path.join(exportdir, metr + "_euc_dist_steps.csv"), index_label="tool")



        # 3 Correlation
        ## Turn tools into dummie variables
        dummies_df=pd.get_dummies(only_tools_df_no_exp)
        ## Calculate p-values using pearson's correlation coefficient
        pval_dic={}
        for col in dummies_df.columns:
            pval_dic[col]=pearsonr(dummies_df[col], metrics_df_std_no_exp["euc_dist"])[1]
        pval_df=pd.DataFrame(pval_dic, index=["p-value"]).transpose()
        pval_df.to_csv(os.path.join(exportdir, metr + "_pvalues_tools.csv"), index_label="tool")#



        # 4 k-means clustering based on euclidean distances
        #   (following https://realpython.com/k-means-clustering-python/)
        kmean_df_std=metrics_df_std_no_exp["euc_dist"].array.reshape(-1,1)

        ## 4.1 Evaluation of appropriate k
        k_list=[]
        reps=100
        ## We repeat k evaluation 100 times to account for randomness:
        for rep in range(1, reps):
            ## 4.1.1 Estimate the best kluster number k
            ### Number of k to test
            ktest=20
            ## 4.1.2 Elbow method based on SSE
            ### The list sse holds the SSE values for each k
            sse_lst_est = []
            for k in range(3, ktest):
                kmeans = KMeans(n_clusters=k)
                kmeans.fit(kmean_df_std)
                sse_lst_est.append(kmeans.inertia_)
            ## 4.1.3 Automatically calculate the best k using KneeLocator
            kl = KneeLocator(range(3, ktest), sse_lst_est, curve="convex", direction="decreasing")
            k_list.append(str(kl.elbow))

        ## 4.1.4 Select k
        ### Count occurence of each tool for each step and save in dict:
        count_df=pd.DataFrame(Counter(k_list), index=["count"]).transpose().reset_index().rename(columns={'index': 'k'})
        # ### Plot counts
        # fig = px.bar(count_df, x='k', y='count')
        # fig.show()
        ### Select k with highest counts
        counts=count_df.sort_values("count")["count"].iloc[-1]
        k=int(count_df.sort_values("count")["k"].iloc[-1])
        print("The best k based on KneeLocator is {0} in {1}/{2} repetitions.".format(str(kl.elbow), counts, reps))


        ## 4.2 Perform the actual kmeans with the estimated k
        ### We cluster with the same k 100 times and determine the cluster labels
        ### with the lowest silhouette score
        k_act_dic={}
        for i in range(1,100):
            k_act={}
            kmeans = KMeans(n_clusters=k)
            kmeans.fit(kmean_df_std)
            k_act["silscore"]=silhouette_score(kmean_df_std, kmeans.labels_)
            k_act["labels"]=kmeans.labels_
            k_act["centers"]=kmeans.cluster_centers_
            k_act_dic[i]=k_act
        silscore_lst_act=[]
        for x in k_act_dic.values():
            silscore_lst_act.append(x["silscore"])
        best_silscore=sorted(silscore_lst_act)[-1]
        ### And we extract the labels and centers for the lowest silhouette score
        for x in k_act_dic.values():
            if x["silscore"]==best_silscore:
                labels=(x["labels"])
                centers=(x["centers"])


        ## 4.3 Add clusters to df
        metrics_df_std_no_exp["cluster"]=[str(x) for x in labels.tolist()]

        ## 4.4 Determine closest cluster to expected based on centroid
        ### Determine euclidean distance of each cluster centroid to expected:
        cluster_centroid_dist=[euclidean(i, exp_std) for i in centers]
        cluster_dist_df = pd.DataFrame({"cluster": [str(x) for x in list(range(0, len(cluster_centroid_dist)))],
            "centroid_euc_dist":cluster_centroid_dist})
        ### Merge pipeline info with cluster info:
        metrics_df_std_no_exp = pd.merge(metrics_df_std_no_exp, cluster_dist_df, on="cluster")
        ### Determine closest cluster and cut df down to pipelines from closest cluster:
        best_cluster=cluster_dist_df.sort_values("centroid_euc_dist").reset_index().iloc[[0]]["cluster"].to_list()[0]
        best_cluster_df=metrics_df_std_no_exp.loc[metrics_df_std_no_exp["cluster"] == best_cluster]
        ### Count occurence of each tool for each step and save in dict:
        counts_dic={}
        for col in best_cluster_df.columns[-10:-3]:
            counts_dic[col]=Counter(best_cluster_df[col])
        ### Save the dict as pickle object, which is needed for the
        ### next step of code in the script "stats_summary_plot.py":
        with open(os.path.join(exportdir, metr + "_closest_cluster_tool_counts.pkl"), 'wb') as f:
            pickle.dump(counts_dic, f)
