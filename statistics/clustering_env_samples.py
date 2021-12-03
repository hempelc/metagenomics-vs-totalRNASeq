#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 3 Dec 2021

# This script compares pipeline data from environmental samples to clusters
# based on mock community samples

import pandas as pd
import os
import logging
from sklearn.preprocessing import StandardScaler

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters
## Full path to directory that contains samples
workdir = "/Users/christopherhempel/Desktop/pipeline_results_mock_community/"
## List of DNA and RNA samples, replicates of 3 plus filtration controls (Neg) and
## extraction controls (Ext); must equal names of directories in workdir that
## contain each sample's pipeline results:
samples = ["M4_DNA", "M4_RNA", "M5_DNA", "M5_RNA", "M6_DNA", "M6_RNA",
    "M_Neg_DNA", "M_Neg_RNA", "M_Ext_DNA", "M_Ext_RNA"]
## Taxonomic rank to group rows on. Either based on genus (option "genus")
## or on species (option "species"):
groupby_rank = "genus"
## Output dir
exportdir=os.path.join(workdir, "results_env_samples_" + groupby_rank)
if not os.path.exists(exportdir):
    os.makedirs(exportdir)





### Standardize columns
std_arr=StandardScaler().fit_transform(master_df)
std_df=pd.DataFrame(std_arr, index=master_df.index, columns=master_df.columns)
