#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 27 Jan 2022

# This script calculates the ARI and AMI indices between clusters of environmental and mock community samples

import pandas as pd
import numpy as np
import plotly.express as px
import os
import pickle
import logging
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import adjusted_mutual_info_score
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

# Set workdir
workdir = "/Users/christopherhempel/Desktop/pipeline_results/cluster_comparison_results/"

# Define functions:
## Define pdist functions for pairwise ARI and AMI indices
def dfun_ami(u, v):
    return adjusted_mutual_info_score(u, v)
def dfun_ari(u, v):
    return adjusted_rand_score(u, v)


# Import data
with open(os.path.join(workdir, "labels_master.pkl"), 'rb') as f:
    labels_master = pickle.load(f)

# Calculate AMI and ARI
idx=list(labels_master.keys())
pairwise_ami= pd.DataFrame(squareform(pdist(list(labels_master.values()), dfun_ami)), index=idx, columns=idx).replace(0, 1)
pairwise_ari = pd.DataFrame(squareform(pdist(list(labels_master.values()), dfun_ari)), index=idx, columns=idx).replace(0, 1)

# Organize dfs
category_orders_comb=['M4_DNA_gen_rel_genus','M4_DNA_gen_rel_species',
    'M5_DNA_gen_rel_genus','M5_DNA_gen_rel_species','M6_DNA_gen_rel_genus',
    'M6_DNA_gen_rel_species', 'M4_DNA_gen_pa_genus','M4_DNA_gen_pa_species',
    'M5_DNA_gen_pa_genus','M5_DNA_gen_pa_species','M6_DNA_gen_pa_genus',
    'M6_DNA_gen_pa_species', 'M4_RNA_gen_rel_genus', 'M4_RNA_gen_rel_species',
    'M5_RNA_gen_rel_genus', 'M5_RNA_gen_rel_species', 'M6_RNA_gen_rel_genus',
    'M6_RNA_gen_rel_species', 'M4_RNA_gen_pa_genus', 'M4_RNA_gen_pa_species',
    'M5_RNA_gen_pa_genus', 'M5_RNA_gen_pa_species', 'M6_RNA_gen_pa_genus',
    'M6_RNA_gen_pa_species']
category_orders_comb=[x.replace("gen_", "") for x in category_orders_comb]

pairwise_ami.index=pairwise_ami.index.str.replace("_env_samples", "")
pairwise_ami=pairwise_ami.reindex(index = category_orders_comb)
pairwise_ami=pairwise_ami.transpose()
pairwise_ami.index=pairwise_ami.index.str.replace("_env_samples", "")
pairwise_ami=pairwise_ami.reindex(index = category_orders_comb)

pairwise_ari.index=pairwise_ari.index.str.replace("_env_samples", "")
pairwise_ari=pairwise_ari.reindex(index = category_orders_comb)
pairwise_ari=pairwise_ari.transpose()
pairwise_ari.index=pairwise_ari.index.str.replace("_env_samples", "")
pairwise_ari=pairwise_ari.reindex(index = category_orders_comb)


# Plot
fig_ami=px.imshow(pairwise_ami, height=800, width=800)
fig_ami.show()
fig_ami.write_image(os.path.join(workdir, "ami_heatmap.png"))

fig_ari=px.imshow(pairwise_ari, height=800, width=800)
fig_ari.show()
fig_ari.write_image(os.path.join(workdir, "ari_heatmap.png"))
