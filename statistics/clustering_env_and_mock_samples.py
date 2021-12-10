#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 3 Dec 2021

# This script compares clustering of environmental and mock community samples

import pandas as pd
import plotly.express as px
import os
import logging
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import rand_score
from kneed import KneeLocator
from collections import Counter
from sklearn.metrics import rand_score
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_score

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters set manually
## Full path to directory that contains samples
workdir = "/Users/christopherhempel/Desktop/pipeline_results"
## Directory for cluster plots
plotdir = os.path.join(workdir, "cluster_comparison_results")
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
## Set if you want to loop over all result combinations of genus/species and
## rel/pa metrics (True or False)
looping=False
## If you set looping to False, then define what specific rank and metrics
## you want to process:
### ("genus", "species")
rank="genus"
### ("rel", "pa")
metrics="rel"


# Parameters set automatically
if looping:
    groupby_rank_lst=["genus", "species"]
    metr_lst=["rel", "pa"]
else:
    groupby_rank_lst=[rank]
    metr_lst=[metrics]
# Dict that will contain all rand indices
rand={}


for groupby_rank in groupby_rank_lst:
    for metr in metr_lst:
        # Df that will contain clusters
        cluster=pd.DataFrame()
        dfs={}
        for sample_set in ["env_samples", "mock_samples"]:

            # 1 Import rel abundances
            file = os.path.join(workdir, "pipeline_results_" + sample_set, "rel_abun_" + groupby_rank,
                "rel_abun_" + groupby_rank + "_" + metr + ".csv")
            df=pd.read_csv(file, index_col="pipeline")
            dfs[sample_set]=df
            ## Standardize
            df_std=pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns)


            # 2 k-means clustering based on euclidean distances
            #   (following https://realpython.com/k-means-clustering-python/)
            ## 2.1 Evaluation of appropriate k
            k_list=[]
            reps=1000
            ## We repeat k evaluation 100 times to account for randomness:
            for rep in range(1, reps):
                ## 2.1.1 Estimate the best kluster number k
                ### Number of k to test
                ktest=20
                ## 2.1.2 Elbow method based on SSE
                ### The list sse holds the SSE values for each k
                sse = {}
                for k in range(3, ktest):
                    kmeans = KMeans(n_clusters=k)
                    kmeans.fit(df_std)
                    sse[k]=(kmeans.inertia_)
                ## 2.1.3 Automatically calculate the best k using KneeLocator
                kl = KneeLocator(range(3, ktest), list(sse.values()), curve="convex", direction="decreasing")
                k_list.append(str(kl.elbow))

            ## 2.1.4 Select k
            ### Count occurence of each tool for each step and save in dict:
            count_df=pd.DataFrame(Counter(k_list[:-1]), index=["count"]).transpose().reset_index().rename(columns={'index': 'k'})
            ### Plot counts
            fig = px.bar(count_df, x='k', y='count')
            fig.show()
            ### Select k with highest counts
            counts=count_df.sort_values("count")["count"].iloc[-1]
            k=int(count_df.sort_values("count")["k"].iloc[-1])
            print("The best k based on KneeLocator is {0} in {1}/{2} repetitions.".format(k, counts, reps))

            ## 2.2 Perform the actual kmeans with the estimated k
            ### We cluster with the same k 100 times and determine the cluster labels
            ### with the lowest silhouette score
            k_act_dic={}
            for i in range(1,100):
                k_act={}
                kmeans = KMeans(n_clusters=k)
                kmeans.fit(df_std)
                k_act["silscore"]=silhouette_score(df_std, kmeans.labels_)
                k_act["labels"]=kmeans.labels_
                k_act_dic[i]=k_act
            silscore_lst_act=[]
            for x in k_act_dic.values():
                silscore_lst_act.append(x["silscore"])
            best_silscore=sorted(silscore_lst_act)[-1]
            ### And we extract the labels and for the lowest silhouette score
            for x in k_act_dic.values():
                if x["silscore"]==best_silscore:
                    labels=(x["labels"])
            ### And save them for the sample set
            cluster[sample_set]=labels

        # 3 Calculate Rand index
        rand[groupby_rank + "_" + metr + "rand"]=rand_score(cluster["env_samples"], cluster["mock_samples"])
        rand[groupby_rank + "_" + metr + "adj_rand"]=adjusted_rand_score(cluster["env_samples"], cluster["mock_samples"])
        print(rand[groupby_rank + "_" + metr + "rand"])
        print(rand[groupby_rank + "_" + metr + "adj_rand"])




        # 4 PCA
        #### On two components for graphical visualization
        # pcs=10
        # pca_graph = PCA(n_components=pcs)
        # pcs_no_exp = pca_graph.fit_transform(df_std)
        # # pc_no_exp_df = pd.DataFrame(data = pcs_no_exp, columns = ['PC1', 'PC2', "PC3", "PC4"],
        # #     index=df_std.index.to_list())
        # for i in range(0, pcs):
        #     print(str(round(pca_graph.explained_variance_ratio_[i] * 100, 2)))
        #         PC1_graph = str(round(pca_graph.explained_variance_ratio_[0] * 100, 2))
        # PC2_graph = str(round(pca_graph.explained_variance_ratio_[1] * 100, 2))
        # PC3_graph = str(round(pca_graph.explained_variance_ratio_[2] * 100, 2))
        # PC4_graph = str(round(pca_graph.explained_variance_ratio_[3] * 100, 2))
        # print("Graph PC1: {0}%, Graph PC2: {1}%".format(PC1_graph, PC2_graph))
        #
        # # Combine with kmeans clusters
        # kmean_df=pd.DataFrame(pcs_no_exp, columns = ['PC1', 'PC2'],
        #     index=std_df.index.to_list()).rename_axis("pipeline").reset_index()
        # kmean_df['cluster']=[str(x) for x in kmeans.labels_.tolist()]
        # fig_kmean = px.scatter(kmean_df, x="PC1", y="PC2", color="cluster", hover_data=["pipeline"])
        # fig_kmean = px.scatter(kmean_df, x="PC1", y="PC2",
        #     labels={"PC1": "PC1 ({0}%)".format(PC1_graph),
        #     "PC2": "PC2 ({0}%)".format(PC2_graph)}, symbol="cluster",
        #         color="cluster", hover_data=["pipeline"])
        # fig_kmean.show()
