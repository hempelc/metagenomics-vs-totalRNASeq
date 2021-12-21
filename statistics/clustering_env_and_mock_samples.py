#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 3 Dec 2021

# This script compares clustering of environmental and mock community samples

import pandas as pd
import plotly.express as px
import os
import logging
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from kneed import KneeLocator
from collections import Counter
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.metrics import silhouette_score
from scipy.spatial import distance_matrix
from kmodes.kmodes import KModes

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
## Set number of clustering tests to estimate appropriate k
k_test_reps = 10
## Set number of repetitions for actual clustering to identify clustering with
## lowest silhouette score
k_actual_reps=10
## Set number of PCs to test in PCAs
pcs=3
## Set if you want to loop over all result combinations of genus/species and
## rel/pa metrics (True or False)
looping=False
## If you set looping to False, then define what specific rank and metrics
## you want to process:
### ("genus", "species")
rank="genus"
### ("rel", "pa")
metrics="pa"


# Parameters set automatically
if looping:
    groupby_rank_lst=["genus", "species"]
    metr_lst=["rel", "pa"]
else:
    groupby_rank_lst=[rank]
    metr_lst=[metrics]

# Dict that will contain all labels
labels_master={}
# Df that will contain all rand indices
indices=pd.DataFrame()

for groupby_rank in groupby_rank_lst:
    for metr in metr_lst:

        samples_set_dfs={}

        for sample_set in ["env_samples", "mock_samples"]:

            # Dic that will contain labels of each sample set
            labels_sample_set={}

            # 1 Import transformed abundances
            file = os.path.join(workdir, "pipeline_results_" + sample_set, "rel_abun_" + groupby_rank,
                "rel_abun_" + groupby_rank + "_" + metr + ".csv")
            df=pd.read_csv(file, index_col="pipeline")
            samples_set_dfs[sample_set]=df


            # 2. PCA - Test what number (if any) of PCs in a PCA improve clustering
            silhouette_scores = {}
            pc_percents = {}
            for n in range(1, pcs):
                n=1
                ## We include clustering without PCA in the loop by skipping PCAs with n=1 and instead using n=1 for no PCA
                if n==1:
                    cluster_alg_input=df
                    n="No PCA"
                else:
                    pca=PCA(n_components=n)
                    cluster_alg_input = pca.fit_transform(df)
                    ## Get summed variance of axes
                    pc_percent=[]
                    for i in range(1, n):
                        pc_percent.append(round(pca.explained_variance_ratio_[i] * 100, 2))
                    pc_percents[n]=sum(pc_percent)



                # 2 Clustering based on euclidean distances
                #   Following https://realpython.com/k-means-clustering-python/ for k-means
                #   and then adapted code for k-modes for p/a data
                ## 2.1 Evaluation of appropriate k
                k_list=[]
                ## We repeat k evaluation 100 times to account for randomness:
                for rep in range(1, k_test_reps):
                    ## 2.1.1 Estimate the best kluster number k
                    ### Number of k to test
                    ktest=20
                    ## 2.1.2 Elbow method based on SSE (kmeans)/cost(kmodes)
                    ### The dic fit holds the SSE/cost values for each k
                    fit={}
                    for k in range(3, ktest):
                        if metr=="pa":
                            cluster_alg = KModes(n_clusters=k)
                            cluster_alg.fit(cluster_alg_input)
                            fit[k]=cluster_alg.cost_
                        else:
                            cluster_alg = KMeans(n_clusters=k)
                            cluster_alg.fit(cluster_alg_input)
                            fit[k]=(cluster_alg.inertia_)

                    fig = px.line(x=fit.keys(), y=fit.values())
                    fig.show()
                    ## 2.1.3 Automatically calculate the best k using KneeLocator
                    kl = KneeLocator(range(3, ktest), list(fit.values()), curve="convex", direction="decreasing")
                    k_list.append(str(kl.elbow))

import matplotlib.pyplot as plt
# Elbow curve to find optimal K
cost = []
K = range(1,5)
for num_clusters in list(K):
    kmode = KModes(n_clusters=num_clusters, init = "random", n_init = 5, verbose=1)
    kmode.fit_predict(cluster_alg_input)
    cost.append(kmode.cost_)

plt.plot(K, cost, 'bx-')
plt.xlabel('No. of clusters')
plt.ylabel('Cost')
plt.title('Elbow Method For Optimal k')
plt.show()

                ## 2.1.4 Select k
                ### Count occurence of each tool for each step and save in dict:
                count_df=pd.DataFrame(Counter(k_list[:-1]), index=["count"]).transpose().reset_index().rename(columns={'index': 'k'})
                ### Plot counts
                # fig = px.bar(count_df, x='k', y='count')
                # fig.show()
                ### Select k with highest counts
                counts=count_df.sort_values("count")["count"].iloc[-1]
                k=int(count_df.sort_values("count")["k"].iloc[-1])
                print("The best k based on KneeLocator is {0} in {1}/{2} repetitions.".format(k, counts, k_test_reps))

                ## 2.2 Perform the actual kmeans with the estimated k
                ### We cluster with the same k 100 times and determine the cluster labels
                ### with the lowest silhouette score
                k_actual_dic={}
                for i in range(1,k_actual_reps):
                    k_act={}
                    cluster_alg = KMeans(n_clusters=k)
                    cluster_alg.fit(cluster_alg_input)
                    k_act["silscore"]=silhouette_score(cluster_alg_input, cluster_alg.labels_)
                    k_act["labels"]=cluster_alg.labels_
                    k_actual_dic[i]=k_act
                silscore_lst_act=[]
                for x in k_actual_dic.values():
                    silscore_lst_act.append(x["silscore"])
                best_silscore=sorted(silscore_lst_act)[-1]
                silhouette_scores[n]=best_silscore
                ### And we extract the labels and for the lowest silhouette score
                for x in k_actual_dic.values():
                    if x["silscore"]==best_silscore:
                        labels=x["labels"]
                ### And save them for the sample set
                labels_sample_set[n]=labels
                labels_master[sample_set]=labels_sample_set

                # 3 Plot PCA with 2 axes
                if n==2:
                    pc1_graph = str(round(pca.explained_variance_ratio_[0] * 100, 2))
                    pc2_graph = str(round(pca.explained_variance_ratio_[1] * 100, 2))
                    pca_graph_df=pd.DataFrame(cluster_alg_input, columns = ['PC1', 'PC2'], index=df.index)
                    pca_graph_df=pca_graph_df.reset_index()
                    pca_graph_df["cluster"]=[str(x) for x in labels]
                    fig_kmean = px.scatter(pca_graph_df, x="PC1", y="PC2",
                        labels={"PC1": "PC1 ({0}%)".format(pc1_graph),
                        "PC2": "PC2 ({0}%)".format(pc2_graph)}, symbol="cluster",
                            color="cluster",
                            title="{0} {1} {2}".format(sample_set, groupby_rank, metr),
                            hover_data=['pipeline'])
                    fig_kmean.show()


        # 3 Calculate Rand index
        ami_lst=[]
        ari_lst=[]
        for i in labels_master["env_samples"].keys():
            ami_lst.append(adjusted_mutual_info_score(labels_master["env_samples"][i], labels_master["mock_samples"][i]))
            ari_lst.append(adjusted_rand_score(labels_master["env_samples"][i], labels_master["mock_samples"][i]))
        indices[groupby_rank + "_" + metr + "_ami"]=ami_lst
        indices[groupby_rank + "_" + metr + "_ari"]=ari_lst
        indices.index=list(labels_master["env_samples"].keys())


        # 4 Calculate concordance of distance matrics
        base = importr('base')
        ape = importr('ape')
        ## 4.1 Bring dfs in similar shape
        mock_sample_df=samples_set_dfs["mock_samples"]
        env_sample_df=samples_set_dfs["env_samples"]
        for taxon in set(mock_sample_df.columns).difference(env_sample_df.columns):
            env_sample_df[taxon]=[0]*len(env_sample_df)
        for taxon in set(env_sample_df.columns).difference(mock_sample_df.columns):
            if taxon not in mock_sample_df.columns:
                mock_sample_df[taxon]=[0]*len(mock_sample_df)
        mock_sample_df=mock_sample_df.reindex(sorted(mock_sample_df.columns), axis=1)
        env_sample_df=env_sample_df.reindex(sorted(env_sample_df.columns), axis=1)

        ## 4.2 Generate distance matrices
        with localconverter(ro.default_converter + pandas2ri.converter):
            dist_matrix_mock_r=base.as_data_frame(pd.DataFrame(distance_matrix(mock_sample_df, mock_sample_df)))
            dist_matrix_env_r=base.as_data_frame(pd.DataFrame(distance_matrix(env_sample_df, env_sample_df)))
            dist_matrix_combined_r=base.rbind(dist_matrix_mock_r, dist_matrix_env_r)

            ## 4.3 Do CADM in R
            cadm=ape.CADM_global(dist_matrix_combined_r, 2, len(dist_matrix_combined_r.columns))
            # Kendall's coefficient of concordance
            W=cadm[0][0]
            # Friedman's chi-square statistic
            Chi2=cadm[0][1]
