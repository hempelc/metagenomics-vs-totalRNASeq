#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 27 Jan 2022

# This script calculates the congruence among distance matrices based on environmental and mock community samples

import pandas as pd
import numpy as np
import plotly.express as px
import os
import logging
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from itertools import product
base = importr('base')
ape = importr('ape')


# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters set manually
## Full path to directory that contains samples
workdir = "/Users/christopherhempel/Desktop/pipeline_results"
# CADM already generated?
cadm_done=True
## Set if you want to loop over all result combinations of genus/species and
## rel/pa metrics (True or False)
looping=True
## If you set looping to False, then define what specific rank and metrics and sample
## you want to process:
### ("genus", "species")
rank="genus"
### ("rel", "pa")
metrics="pa"
### ("env_samples", "mock_samples)
sets="env_samples"


# Parameters set automatically
## Directory for cluster plots
outdir = os.path.join(workdir, "cluster_comparison_results")
if not os.path.exists(outdir):
    os.mkdir(outdir)
if looping:
    groupby_rank_lst=["genus", "species"]
    metr_lst=["rel", "pa"]
    sample_sets=["env_samples", "mock_samples"]

else:
    groupby_rank_lst=[rank]
    metr_lst=[metrics]
    sample_sets=[sets]

# Dict for distance to use
distance_metric_paird={"rel": "euclidean", "pa": "jaccard"}
# Dict that will contain all distance matrices
pairds_master={}


# Define functions:
## Calculate paired CADM using a dictionary of dissimilarity marices
def paird_cadm(dic):
    w_lst=[]
    chi2_lst=[]
    vals=dic.values()
    idx=list(dic.keys())
    len_vals=len(vals)
    with localconverter(ro.default_converter + pandas2ri.converter):
        for v, u in product(vals, vals):
            u_r=base.as_data_frame(pd.DataFrame(u))
            v_r=base.as_data_frame(pd.DataFrame(v))
            combined_r=base.rbind(u_r, v_r)
            ## 4.3 Do CADM in R
            cadm=ape.CADM_global(combined_r, 2, len(combined_r.columns))
            # Kendall's coefficient of concordance
            w_lst.append(cadm[0][0])
            # Friedman's chi-square statistic
            chi2_lst.append(cadm[0][1])
    w_lst_reshaped=np.reshape(w_lst, [len_vals, len_vals])
    return pd.DataFrame(w_lst_reshaped, index=idx, columns=idx)

if cadm_done==False:
    for metr in metr_lst:
        for groupby_rank in groupby_rank_lst:
            for sample_set in sample_sets:
                if sample_set=="mock_samples":
                    samples=["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA"]
                else:
                    samples=["F4_DNA", "F4_RNA", "F5_DNA", "F5_RNA", "F6_DNA", "F6_RNA"]
                for sample in samples:
                        # Import transformed abundances
                        file = os.path.join(workdir, "pipeline_results_" + sample_set, "rel_abun_" + groupby_rank,
                            sample + "_rel_abun_" + groupby_rank + "_" + metr + ".csv")
                        df=pd.read_csv(file, index_col="pipeline")

                        # Calculate pairwise distance
                        paird=squareform(pdist(df, metric=distance_metric_paird[metr]))
                        pairds_master["{0}_{1}_{2}_{3}".format(sample, metr, groupby_rank, sample_set)]=paird

    #  Calculate congruence among distance matrices
    paird_cadm_df=paird_cadm(pairds_master)
else:
    paird_cadm_df=pd.read_csv(os.path.join(outdir, "cadm_df.csv"), index_col=0)


# Organize df
category_orders_comb_M=['M4_DNA_rel_genus','M4_DNA_rel_species',
    'M5_DNA_rel_genus','M5_DNA_rel_species','M6_DNA_rel_genus',
    'M6_DNA_rel_species', 'M4_DNA_pa_genus','M4_DNA_pa_species',
    'M5_DNA_pa_genus','M5_DNA_pa_species','M6_DNA_pa_genus',
    'M6_DNA_pa_species', 'M4_RNA_rel_genus', 'M4_RNA_rel_species',
    'M5_RNA_rel_genus', 'M5_RNA_rel_species', 'M6_RNA_rel_genus',
    'M6_RNA_rel_species', 'M4_RNA_pa_genus', 'M4_RNA_pa_species',
    'M5_RNA_pa_genus', 'M5_RNA_pa_species', 'M6_RNA_pa_genus',
    'M6_RNA_pa_species']
category_orders_comb_F=[x.replace("M", "F") + "_env_samples" for x in category_orders_comb_M]
category_orders_comb_M=[x + "_mock_samples" for x in category_orders_comb_M]
category_orders_comb_both=category_orders_comb_M+category_orders_comb_F
paird_cadm_df=paird_cadm_df.reindex(index = category_orders_comb_both)
paird_cadm_df=paird_cadm_df.transpose()
paird_cadm_df=paird_cadm_df.reindex(index = category_orders_comb_both)

# Plot
fig_cadm=px.imshow(paird_cadm_df, height=1200, width=1200)
fig_cadm.show()
fig_cadm.write_image(os.path.join(outdir, "cadm_heatmap.png"))

# Save df
paird_cadm_df.to_csv(os.path.join(outdir, "cadm_df.csv"))



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
