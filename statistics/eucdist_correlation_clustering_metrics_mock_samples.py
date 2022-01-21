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
types=["DNA", "RNA"]
reps=["M4", "M5", "M6"]
## Set if you want to loop over all result combinations of parameters in script
## "processing_and_metrics.py" and all metrics (True or False)
looping=False
## Set lists for aggregation, combinations, and metrics for looping:
aggs=["agg_reps_agg_type", "agg_reps_sep_type", "sep_reps_agg_type", "sep_reps_sep_type"]
combinations=["cell_genus", "cell_species", "gen_genus", "gen_species"]
metr_list=["rel", "pa"]
## If you set looping to False, then define what specific aggregation, combination,
## and metrics you want to process:
### ("agg_reps_agg_type", "agg_reps_sep_type", "sep_reps_ag_type", "sep_reps_sep_type")
agg="agg_reps_sep_type"
### ("cell_genus", "cell_species", "gen_genus", "gen_species")
comb="gen_species"
### ("rel", "pa")
metrics="pa"


# Parameters set automatically
if looping:
    aggs=aggs
    combinations=combinations
    metr_list=metr_list
else:
    aggs=[agg]
    combinations=[comb]
    metr_list=[metrics]


# Functions
def sum_nested_dics(dic1, dic2):
    summed_dic={}
    for dic in [dic1, dic2]:
        for key in dic.keys():
                for x,y in dic[key].items():
                    if not key in summed_dic.keys():
                        summed_dic[key] = {x:y}
                    else:
                        if not x in summed_dic[key].keys():
                            summed_dic[key].update({x:y})
                        else:
                            summed_dic[key].update({x:summed_dic[key][x] + y})
    return summed_dic


master_master_counts_dic={}

for combination in combinations:
    for metr in metr_list:
        ## Make plot export directory:
        exportdir_level1=os.path.join(workdir, "metrics_" + combination, "stats_" + metr)
        if not os.path.exists(exportdir_level1):
            os.mkdir(exportdir_level1)

        ### Do the following chunk for all samples separately
        master_eucdist_dfs={}
        for sample in samples:
            # 1 Import data
            ## Import metrics df
            metrics_df=pd.read_csv(os.path.join(workdir, "metrics_" + combination, sample + "_metrics_df.csv"), index_col=0)
            ## Convert trimming score column type into str
            metrics_df['trimming_score'] = metrics_df['trimming_score'].astype(str)
            ### Drop additional ones specified by parameter "metr" and set expected dummy
            ### (=ideal pipeline results):
            if metr not in ["rel", "pa"]:
                logging.critical("Parameter metr not in rel_metr, or pa_metr: metr=" + metr)
            if metr=="rel":
                metrics_df=metrics_df.drop(["FP", "TP"], axis=1)
            elif metr=="pa":
                no_drops=["FP", "TP", 'type', 'trimming_score',
                    'rRNA_sorting_tool', 'assembly_tool', 'mapper', 'database', 'classifier']
                metrics_df=metrics_df.drop([x for x in metrics_df.columns if x not in no_drops], axis=1)
            ## Separate expected from dfs
            metrics_df_no_exp=metrics_df.drop(["expected"], axis=0)
            exp=metrics_df.iloc[:, :-7].loc["expected"]
            only_metrics_df_no_exp=metrics_df.iloc[:, :-7].drop(["expected"], axis=0)


            # 2 Euclidean distances between pipelines and expected:
            ## Calculate euc dist for each pipeline
            for index,row in only_metrics_df_no_exp.iterrows():
                metrics_df_no_exp.loc[index,'euc_dist'] = euclidean(row, exp)
            master_eucdist_dfs[sample]=metrics_df_no_exp


        ## Aggregate samples
        for aggregation in aggs:
            master_counts_dic={}
            exportdir_level2=os.path.join(exportdir_level1, aggregation)
            if not os.path.exists(exportdir_level2):
                os.mkdir(exportdir_level2)
            aggregation_dic={}
            if aggregation=="agg_reps_agg_type":
                agg_df=pd.DataFrame()
                for sample in samples:
                    agg_df=pd.concat([agg_df, master_eucdist_dfs[sample]])
                aggregation_dic["agg"]=agg_df
            elif aggregation=="agg_reps_sep_type":
                for type in types:
                    agg_df=pd.DataFrame()
                    for rep in reps:
                        agg_df=pd.concat([agg_df, master_eucdist_dfs[rep + "_" + type]])
                    aggregation_dic[type]=agg_df
            elif aggregation=="sep_reps_agg_type":
                for rep in reps:
                    agg_df=pd.DataFrame()
                    for type in types:
                        agg_df=pd.concat([agg_df, master_eucdist_dfs[rep + "_" + type]])
                    aggregation_dic[rep]=agg_df
            elif aggregation=="sep_reps_sep_type":
                aggregation_dic=master_eucdist_dfs


            for combo, df in aggregation_dic.items():
                ## Calculate average and minimum euc dist for each tool in each step
                mean_euc_dist_lst=[]
                min_euc_dist_lst=[]
                tools=[]
                for step in df.loc[:, 'type':'classifier'].columns:
                    for tool in df[step].unique():
                        ### Make sure trimming score 5 doesn't pick up 15 as well
                        if tool=="5":
                            tool="_5"
                        #### Cut down df to rows containing tool and take the averge and min euc_dist:
                        mean_euc_dist=df.reset_index()[df.reset_index()['pipeline']
                            .str.contains(tool)]["euc_dist"].mean()
                        min_euc_dist=df.reset_index()[df.reset_index()['pipeline']
                            .str.contains(tool)]["euc_dist"].sort_values().min()
                        ### Reverse changes
                        if tool=="_5":
                            tool="5"
                        tools.append(step + "_" + tool)
                        mean_euc_dist_lst.append(mean_euc_dist)
                        min_euc_dist_lst.append(min_euc_dist)

                ## Make df
                euc_dist_final_df=pd.DataFrame({"mean_euc_dist": mean_euc_dist_lst, "min_euc_dist": min_euc_dist_lst}, index=tools).reset_index()



                # 4 "Clustering", or more correctly, segmentation of euclidean distances
                segmentation_df=df["euc_dist"].array.reshape(-1,1)

                ## 4.1 Evaluation of appropriate number of clusters using Kernel Density estimation
                ##     by performing a grid search to identify best value for bandwidth
                grid = GridSearchCV(KernelDensity(kernel = 'gaussian'),{'bandwidth': np.linspace(0.1, 2, 20)})
                grid_fit=grid.fit(segmentation_df)
                bandw=list(grid_fit.best_params_.values())[0]

                ## 4.2 Run the actual Kernel Density Estimation
                kde = KernelDensity(kernel='gaussian', bandwidth=bandw).fit(segmentation_df)
                s = linspace(min(df["euc_dist"]),max(df["euc_dist"]))
                e = kde.score_samples(s.reshape(-1,1))
                fig_kde=px.line(x=s, y=e)
                fig_kde.show()

                ## 4.3 Define minima and split data into clusters based on minima
                mi=argrelextrema(e, np.less)[0]
                print("For combination {0}_{1}_{2}, the best number of clusters is {3}."\
                    .format(combo, combination, metr, len(mi)+1))
                clusters={}
                for i in range(0, len(mi)+1):
                    if i==0:
                        clusters[i]=list(segmentation_df[segmentation_df < s[mi][i]])
                    elif i==len(mi):
                        clusters[i]=list(segmentation_df[segmentation_df >= s[mi][i-1]])
                    else:
                        clusters[i]=list(segmentation_df[(segmentation_df >= s[mi][i-1]) * (segmentation_df <= s[mi][i])])

                ## 4.4 Turn clusters into mergeable dataframe and add them to metrics df
                clusters_inv={}
                for cluster, eucdist_list in clusters.items():
                    for eucdist in eucdist_list:
                        clusters_inv[eucdist]=[str(cluster)]
                cluster_df=pd.DataFrame(clusters_inv).transpose().reset_index()
                cluster_df.columns = ['euc_dist', 'cluster']
                df_cluster=pd.merge(df.reset_index(), cluster_df, on="euc_dist").set_index('pipeline')

                ## 4.5 Determine closest cluster to expected based on centroid
                ### Determine euclidean distance of each cluster centroid to expected:
                best_cluster_df=df_cluster.loc[df_cluster["cluster"] == '0']
                ### Count occurence of each tool for each step in closest cluster (=always cluster 0) and save in dict:
                counts_dic={}
                for col in best_cluster_df.columns[-9:-2]:
                    counts_dic[col]=Counter(best_cluster_df[col])
                master_counts_dic["{0}_{1}_{2}".format(combo, combination, metr)]=counts_dic
                master_master_counts_dic[aggregation + "_" + metr + "_" + combination]=master_counts_dic
                master_master_counts_dic["sep_reps_sep_type_rel_gen_genus"].keys()
                ### Save the count dict as pickle object, which is needed for the script "stats_summary_plot.py":
                with open(os.path.join(exportdir_level2, "closest_cluster_tool_counts.pkl"), 'wb') as f:
                    pickle.dump(master_counts_dic, f)


                # 3 Correlation
                ## Turn tools into dummie variables
                dummies_df=pd.get_dummies(df.iloc[:, 2:-1])
                ## Calculate p-values using pearson's correlation coefficient
                pval_dic={}
                for col in dummies_df.columns:
                    X=df["euc_dist"].to_numpy()
                    y=np.array(dummies_df[col])
                    # logit_model_result=sm.Logit(y,X).fit()
                    # pval=logit_model_result.summary2().tables[1]['P>|z|'][0]
                    pval=pointbiserialr(X,y)[1]
                    pval_dic[col]=pval
                pval_df=pd.DataFrame(pval_dic, index=["p-value"]).transpose().reset_index()


                # Merge dfs
                merged=pd.merge(euc_dist_final_df, pval_df, on="index").set_index("index")
                merged.to_csv(os.path.join(exportdir_level2, "eucdist_pvalues_tools_" + combo + ".csv"), index_label="tool")
