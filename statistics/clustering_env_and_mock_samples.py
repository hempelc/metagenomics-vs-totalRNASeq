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
from sklearn.manifold import MDS
from kneed import KneeLocator
from collections import Counter
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.metrics import silhouette_score
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from sklearn_extra.cluster import KMedoids
from skbio.stats.ordination import pcoa

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
k_test_reps = 2
## Set number of repetitions for actual clustering to identify clustering with
## lowest silhouette score
k_actual_reps=2
## Set number of PCs to test in PCAs
pcs=2
## Set if you want to loop over all result combinations of genus/species and
## rel/pa metrics (True or False)
looping=True
## If you set looping to False, then define what specific rank and metrics
## you want to process:
### ("genus", "species")
rank="genus"
### ("rel", "pa")
metrics="pa"
samples = ["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA"]


# Parameters set automatically
if looping:
    groupby_rank_lst=["genus", "species"]
    metr_lst=["rel", "pa"]
else:
    groupby_rank_lst=[rank]
    metr_lst=[metrics]
# Dict that will contain all labels and distance matrices
labels_master={}
pairds_master={}


for metr in metr_lst:

    if metr=="rel":
        distance_metric_kmedoids='euclidean'
        distance_metric_paird="euclidean"
    elif metr=='pa':
        distance_metric_kmedoids='precomputed'
        distance_metric_paird="jaccard"

    for groupby_rank in groupby_rank_lst:
        for sample in samples:


            #for sample_set in ["env_samples", "mock_samples"]:
            for sample_set in ["env_samples"]:

                # Dic that will contain labels of each sample set
                labels_sample_set={}

                # 1 Import transformed abundances
                file = os.path.join(workdir, "pipeline_results_" + sample_set, "rel_abun_" + groupby_rank,
                    sample + "_rel_abun_" + groupby_rank + "_" + metr + ".csv")
                df=pd.read_csv(file, index_col="pipeline")

                # Calculate pairwise distance
                paird=squareform(pdist(df, metric=distance_metric_paird))
                pairds_master["{0}_{1}_{2}".format(metr, groupby_rank, sample_set)]=paird


                # 2 Clustering based on euclidean distances
                #   Following https://realpython.com/k-means-clustering-python/ for k-means
                #   and then adapted code for k-modes for p/a data
                ## 2.1 Evaluation of appropriate k
                k_list=[]


                ## We repeat k evaluation 100 times to account for randomness:
                for rep in range(0, k_test_reps):
                    ## 2.1.1 Estimate the best kluster number k
                    ### Number of k to test
                    ktest=10
                    ## 2.1.2 Elbow method based on SSE (kmeans)/cost(kmodes)
                    ### The dic fit holds the SSE/cost values for each k
                    fit={}
                    for k in range(3, ktest):
                        cluster_alg = KMedoids(n_clusters=k, method="pam", max_iter=1000, metric=distance_metric_kmedoids)
                        cluster_alg.fit(paird)
                        fit[k]=(cluster_alg.inertia_)

                    # fig = px.line(x=fit.keys(), y=fit.values())
                    # fig.show()
                    ## 2.1.3 Automatically calculate the best k using KneeLocator
                    kl = KneeLocator(range(3, ktest), list(fit.values()), curve="convex", direction="decreasing")
                    k_list.append(str(kl.elbow))


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
                for i in range(0,k_actual_reps):
                    k_act={}
                    cluster_alg = KMedoids(n_clusters=k, method="pam", max_iter=1000, metric=distance_metric_kmedoids)
                    cluster_alg.fit(paird)
                    k_act["silscore"]=silhouette_score(paird, cluster_alg.labels_)
                    k_act["labels"]=cluster_alg.labels_
                    k_actual_dic[i]=k_act
                silscore_lst_act=[]
                for x in k_actual_dic.values():
                    silscore_lst_act.append(x["silscore"])
                best_silscore=sorted(silscore_lst_act)[-1]
                ### And we extract the labels for the lowest silhouette score
                for x in k_actual_dic.values():
                    if x["silscore"]==best_silscore:
                        labels=x["labels"]
                ### And save them for the sample set
                labels_master["{0}_{1}_{2}".format(metr, groupby_rank, sample_set)]=labels

                # 3 Do PCA with 2 axes
                # PCoA
                pcoa_model = pcoa(paird, number_of_dimensions=3)
                pca_graph_df=pcoa_model.samples[['PC1', 'PC2', 'PC3']]
                pca_graph_df['pipeline']=df.index
                pc1_graph = str(round(pcoa_model.proportion_explained[0] * 100, 2))
                pc2_graph = str(round(pcoa_model.proportion_explained[1] * 100, 2))
                pc3_graph = str(round(pcoa_model.proportion_explained[2] * 100, 2))
                # # Alternatively PCA
                # pca=PCA(n_components=3)
                # pca_input = pca.fit_transform(df)
                # pc1_graph = str(round(pca.explained_variance_ratio_[0] * 100, 2))
                # pc2_graph = str(round(pca.explained_variance_ratio_[1] * 100, 2))
                # pc3_graph = str(round(pca.explained_variance_ratio_[2] * 100, 2))
                # pca_graph_df=pd.DataFrame(pca_input, columns = ['PC1', 'PC2', 'PC3'], index=df.index)
                # pca_graph_df=pca_graph_df.reset_index()
                pca_graph_df["cluster"]=[str(x) for x in labels]
                fig_kmean = px.scatter_3d(pca_graph_df, x="PC1", y="PC2", z="PC3",
                    labels={"PC1": "PC1 ({0}%)".format(pc1_graph),
                    "PC2": "PC2 ({0}%)".format(pc2_graph),
                    "PC3": "PC3 ({0}%)".format(pc3_graph)}, symbol="cluster",
                        color="cluster",
                        title="{0} {1} {2}".format(sample_set, groupby_rank, metr),
                        hover_data=['pipeline'])
                fig_kmean.show()


#LOOP OVER SAMPLES UNSURE

# 3 Calculate Rand index

## Define function to calculate pairwise indices
def dfun_ami(u, v):
    return adjusted_mutual_info_score(u, v)
def dfun_ari(u, v):
    return adjusted_rand_score(u, v)

pairwise_ami, pairwise_ari = pdist(list(labels_master.values()), dfun_ami), pdist(list(labels_master.values()), dfun_ari)

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



# Hierarchial clustering
                # from scipy.cluster.hierarchy import dendrogram, linkage
                # from sklearn.cluster import AgglomerativeClustering
                # from sklearn.metrics import silhouette_score
                # import matplotlib.pyplot as plt
#
#
#                 links = linkage(df, 'single')
#                 points = range(1, len(df))
#                 plt.figure(figsize=(8, 6))
#                 dendrogram(links, orientation='top', labels=df.index, distance_sort='descending', show_leaf_counts=True)
#                 plt.show()
#
#                 silhouette_coefficients = []
#                 for k in range(3, 20):
#                     cluster_alg = AgglomerativeClustering(n_clusters=k, affinity='euclidean', linkage='ward')
#                     cluster_alg.fit_predict(paird)
#                     score = silhouette_score(paird, cluster_alg.labels_)
#                     silhouette_coefficients.append(score)
#                 ### Plot the scores and manually pick the k with the highest score close to
#                 ### the k estimated with SSE
#                 fig_k_silhouette = px.line(pd.DataFrame(list(zip(range(2, 20),
#                     silhouette_coefficients)), columns =['Number of Clusters', 'Silhouette Coefficients']),
#                     x="Number of Clusters", y="Silhouette Coefficients")
#                 fig_k_silhouette.show()
#
#                 groups = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
#                 groups .fit_predict(iris_data)
