#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 22 Jul 2021

# This script processes the output from the script "processing_and_metrics.py"

import pandas as pd
import numpy as np
import plotly.express as px
import os
import copy
import logging
import pickle
import statsmodels.api as sm
from scipy.spatial.distance import euclidean
from scipy.stats import pointbiserialr
from kneed import KneeLocator
from sklearn.metrics import silhouette_score
from collections import Counter
from sklearn.decomposition import PCA
from numpy import array, linspace
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from scipy.signal import argrelextrema
from sklearn.preprocessing import StandardScaler


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
## Set lists for aggregation, combinations, and metrics for looping:
aggs=["agg_reps_agg_type", "agg_reps_sep_type", "sep_reps_agg_type", "sep_reps_sep_type"]
combinations=["cell_genus", "cell_species", "gen_genus", "gen_species"]
metr_list=["rel", "pa"]
## If you set looping to False, then define what specific aggregation, combination,
## and metrics you want to process:
### ("agg_reps_agg_type", "agg_reps_sep_type", "sep_reps_ag_type", "sep_reps_sep_type")
agg="agg_reps_agg_type"
### ("cell_genus", "cell_species", "gen_genus", "gen_species")
comb="gen_genus"
### ("rel", "pa")
metrics="pa"

master_master_counts_dic={}


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
# Sum of nested dics:
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
                ## Standardize (needed since TP and FP are not standardized)
                metrics_df_TP_FP_std=pd.DataFrame(StandardScaler().fit_transform(metrics_df.iloc[:, 0:2]), index=metrics_df.index, columns=["TP", "FP"])
                metrics_df=pd.concat([metrics_df_TP_FP_std, metrics_df.iloc[:, 2:]], axis=1)
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
                ## 3 Calculate average and minimum euc dist for each tool in each step
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

                # #NEW: maybe split based on either minima or maxima
                # ## 4.3 Define minima and split data into clusters based on minima
                # mi=argrelextrema(e, np.less)[0]
                # ma=argrelextrema(e, np.greater)[0]
                # if mi[0]>ma[0]:
                #     cuts=ma
                # else:
                #     cuts=mi
                # print("For combination {0}_{1}_{2}, the best number of clusters is {3}."\
                #     .format(combo, combination, metr, len(cuts)+2))
                # clusters={}
                # for i in range(0, len(cuts)+1):
                #     if i==0:
                #         clusters[i]=list(segmentation_df[segmentation_df < s[cuts][i]])
                #     elif i==len(cuts):
                #         clusters[i]=list(segmentation_df[segmentation_df >= s[cuts][i-1]])
                #     else:
                #         clusters[i]=list(segmentation_df[(segmentation_df >= s[cuts][i-1]) * (segmentation_df <= s[cuts][i])])

                # ORIGINAL (just minima):
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

                ## 4.5 Determine closest cluster to expected
                ### The best cluster is the first:
                best_cluster_df=df_cluster.loc[df_cluster["cluster"] == '0']
                ### Count occurence of each tool for each step in closest cluster (=always cluster 0) and save in dict:
                counts_dic={}
                for col in best_cluster_df.columns[-9:-2]:
                    counts_dic[col]=Counter(best_cluster_df[col])
                counts_dic["mean_euc_dist"]=best_cluster_df["euc_dist"].mean()
                master_counts_dic["{0}_{1}_{2}".format(combo, combination, metr)]=counts_dic
                master_master_counts_dic[aggregation + "_" + metr + "_" + combination]=master_counts_dic
                ### Save the count dict as pickle object, which is needed for the script "stats_summary_plot.py":
                with open(os.path.join(exportdir_level2, "closest_cluster_tool_counts.pkl"), 'wb') as f:
                    pickle.dump(master_counts_dic, f)


                # 5 Correlation
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



# def get_jenks_breaks(data_list, number_class):
#     data_list.sort()
#     mat1 = []
#     for i in range(len(data_list) + 1):
#         temp = []
#         for j in range(number_class + 1):
#             temp.append(0)
#         mat1.append(temp)
#     mat2 = []
#     for i in range(len(data_list) + 1):
#         temp = []
#         for j in range(number_class + 1):
#             temp.append(0)
#         mat2.append(temp)
#     for i in range(1, number_class + 1):
#         mat1[1][i] = 1
#         mat2[1][i] = 0
#         for j in range(2, len(data_list) + 1):
#             mat2[j][i] = float('inf')
#     v = 0.0
#     for l in range(2, len(data_list) + 1):
#         s1 = 0.0
#         s2 = 0.0
#         w = 0.0
#         for m in range(1, l + 1):
#             i3 = l - m + 1
#             val = float(data_list[i3 - 1])
#             s2 += val * val
#             s1 += val
#             w += 1
#             v = s2 - (s1 * s1) / w
#             i4 = i3 - 1
#             if i4 != 0:
#                 for j in range(2, number_class + 1):
#                     if mat2[l][j] >= (v + mat2[i4][j - 1]):
#                         mat1[l][j] = i3
#                         mat2[l][j] = v + mat2[i4][j - 1]
#         mat1[l][1] = 1
#         mat2[l][1] = v
#     k = len(data_list)
#     kclass = []
#     for i in range(number_class + 1):
#         kclass.append(min(data_list))
#     kclass[number_class] = float(data_list[len(data_list) - 1])
#     count_num = number_class
#     while count_num >= 2:  # print "rank = " + str(mat1[k][count_num])
#         idx = int((mat1[k][count_num]) - 2)
#         # print "val = " + str(data_list[idx])
#         kclass[count_num - 1] = data_list[idx]
#         k = int((mat1[k][count_num] - 1))
#         count_num -= 1
#     return kclass
#
# import numpy as np
# import matplotlib.pyplot as plt
#
# x = list(aggregation_dic["DNA"]["euc_dist"])
# breaks = get_jenks_breaks(x, 7)
#
# for line in breaks:
#     plt.plot([line for _ in range(len(x))], 'k--')
#
# plt.plot(x)
# plt.grid(True)
# plt.show()
#
#
# from sklearn.cluster import KMeans
# k_list=[]
#
# ## We repeat k evaluation "k_test_reps" times to account for randomness:
# for rep in range(0, 3):
#     ## 2.1.1 Estimate the best kluster number k
#     ### Number of k to test
#     ktest=20
#     ## 2.1.2 Elbow method based on SSE (kmeans)
#     ### The dic fit holds the SSE values for each k
#     fit={}
#     for k in range(3, ktest):
#         cluster_alg = KMeans(n_clusters=k)
#         cluster_alg.fit(aggregation_dic["agg"]["euc_dist"].to_numpy().reshape(-1, 1))
#         fit[k]=(cluster_alg.inertia_)
#
#     ### Plot SSEs
#     fig = px.line(x=fit.keys(), y=fit.values())
#     fig.show()
#
#     ## 2.1.3 Automatically calculate the best k using KneeLocator
#     kl = KneeLocator(range(3, ktest), list(fit.values()), curve="convex", direction="decreasing")
#     k_list.append(str(kl.elbow))
#
#
# ## 2.1.4 Select k
# ### Count occurence of each tool for each step and save in dict:
# count_df=pd.DataFrame(Counter(k_list), index=["count"]).transpose().reset_index().rename(columns={'index': 'k'})
# ### Plot counts
# fig = px.bar(count_df, x='k', y='count')
# fig.show()
# ### Select k with highest counts
# counts=count_df.sort_values("count")["count"].iloc[-1]
# k=int(count_df.sort_values("count")["k"].iloc[-1])
# print("The best k based on KneeLocator is {0} in {1}/{2} repetitions.".format(k, counts, 3))
#
# ## 2.2 Perform the actual kmeans with the estimated k
# ### We cluster with the same k 100 times and determine the cluster labels
# ### with the lowest silhouette score
# k_actual_dic={}
# for i in range(0,20):
#     k_act={}
#     cluster_alg = KMeans(n_clusters=k)
#     cluster_alg.fit(aggregation_dic["agg"]["euc_dist"].to_numpy().reshape(-1, 1))
#     k_act["silscore"]=silhouette_score(aggregation_dic["agg"]["euc_dist"].to_numpy().reshape(-1, 1), cluster_alg.labels_)
#     k_act["labels"]=cluster_alg.labels_
#     k_actual_dic[i]=k_act
# silscore_lst_act=[]
# for x in k_actual_dic.values():
#     silscore_lst_act.append(x["silscore"])
# best_silscore=sorted(silscore_lst_act)[0]
# labels=[]
# ### And we extract the labels for the lowest silhouette score
# for x in k_actual_dic.values():
#     if x["silscore"]==best_silscore:
#         labels=x["labels"]
#
# aggregation_dic["agg"]["cluster"]=labels
#
# ## 4.5 Determine closest cluster to expected based on centroid
# ### Determine euclidean distance of each cluster centroid to expected:
# for i in set(labels):
#     print(aggregation_dic["agg"][aggregation_dic["agg"]["cluster"]==i]["euc_dist"].mean())
# aggregation_dic["agg"][aggregation_dic["agg"]["cluster"]==4]
# best_cluster_df=aggregation_dic["agg"][aggregation_dic["agg"]["cluster"]==4]
# ### Count occurence of each tool for each step in closest cluster (=always cluster 0) and save in dict:
# counts_dic={}
# for col in best_cluster_df.columns[-9:-2]:
#     counts_dic[col]=Counter(best_cluster_df[col])
