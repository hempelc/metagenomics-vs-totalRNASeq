import pandas as pd
import numpy as np
from scipy.special import softmax
import math
import plotly.express as px
from scipy.spatial.distance import euclidean
import os

# Input directory
input_dir='~/Desktop/pipeline_results/pipeline_results_mock_samples/metrics_gen_genus/'

# Define functions
## Reversing multiplicative replacement
def reverse_multiplicative_replacement(lst, known_sum):
    repl=round(1/len(lst)**2, 8)
    rounded=[round(x, 8) for x in lst]
    non_repl_lst=[x-repl if x==repl else x for x in rounded]
    n_zero=len([x for x in non_repl_lst if x==0])
    div=known_sum/(1-repl*n_zero)
    div_lst=[x*div for x in non_repl_lst]
    return np.array(div_lst)
## Apply function for pd.apply (softmax reverses clr transformation)
def reverse(x):
    return reverse_multiplicative_replacement(softmax(x),1)


# Define best pipelines
dic_best_pip={"M4_DNA": "DNA_5_unsorted_metaspades_bwa_ncbi-nt_kraken2",
    "M5_DNA": "DNA_15_unsorted_transabyss_bwa_ncbi-nt_blast-filtered",
    "M6_DNA": "DNA_20_unsorted_spades_bwa_ncbi-nt_kraken2",
    "M4_RNA": "RNA_20_unsorted_rnaspades_bwa_silva_blast-first-hit",
    "M5_RNA": "RNA_5_rrnafilter_transabyss_bwa_silva_kraken2",
    "M6_RNA": "RNA_5_barrnap_rnaspades_bwa_silva_kraken2"}
# And empty dics
dic_reverse={}
dic_clr_mr={}
dic_best_abun={}
dic_best_abun_clr_mr={}


# Loop over replicates
for rep in ["M4", "M5", "M6"]:
    for NA in ["DNA", "RNA"]:
        # Import data
        sam=rep + "_" + NA
        df=pd.read_csv(os.path.join(input_dir, "{0}_metrics_df.csv".format(sam)), index_col="pipeline")
        df_crop=df.iloc[:, 0:11].transpose()
        # Save transformed df
        dic_clr_mr[sam]=df_crop.transpose()
        # Reverse transformation and save df
        df_rev=df_crop.apply(reverse, axis=0).transpose()
        exp=df_rev.loc["expected"]
        exp_clr_mr=df_crop.transpose().loc["expected"]
        dic_reverse[sam]=df_rev


# Extract best pipeline abundances for both transformed and non-transformed abundances
for rep in ["M4", "M5", "M6"]:
    for NA in ["DNA", "RNA"]:
        sam=rep + "_" + NA
        dic_best_abun[NA+ "_" +rep]=dic_reverse[sam].loc[dic_best_pip[sam]]
for rep in ["M4", "M5", "M6"]:
    for NA in ["DNA", "RNA"]:
        sam=rep + "_" + NA
        dic_best_abun_clr_mr[NA+ "_" +rep]=dic_clr_mr[sam].loc[dic_best_pip[sam]]


# Add expected and normalize abundances from 0 to 1
## For non-transformed abundances
fin=pd.DataFrame(dic_best_abun)
fin = fin.reindex(sorted(fin.columns), axis=1)
fin["expected"]=exp
fin=fin.transpose()

norm_abs=pd.DataFrame()
for i in fin.columns:
    max=abs(fin[i]-fin[i]["expected"]).max()
    norm_abs[i]=abs(fin[i]-fin[i]["expected"])/max

## For transformed abundances
fin_clr_mr=pd.DataFrame(dic_best_abun_clr_mr)
fin_clr_mr["expected"]=exp_clr_mr
fin_clr_mr = fin_clr_mr.reindex(sorted(fin_clr_mr.columns), axis=1).transpose()
fin_clr_mr["euc_dist"]=fin_clr_mr.apply(lambda row : euclidean(row, exp_clr_mr), axis=1)


norm_fin_clr_mr=pd.DataFrame()
for i in fin_clr_mr.columns:
    max=abs(fin_clr_mr[i]-fin_clr_mr[i]["expected"]).max()
    norm_fin_clr_mr[i]=abs(fin_clr_mr[i]-fin_clr_mr[i]["expected"])/max


# Heatmap of transformed and non-transformed abundaces before and after normalization
fig=px.imshow(fin, color_continuous_scale=px.colors.sequential.Blues)
fig.show()
fig=px.imshow(norm_abs, color_continuous_scale=px.colors.sequential.Blues)
fig.show()

fig=px.imshow(fin_clr_mr, color_continuous_scale=px.colors.sequential.Blues)
fig.show()
fig.write_image("/Users/christopherhempel/Desktop/heatmap_abun_cell_genus.png")
fig=px.imshow(norm_fin_clr_mr, color_continuous_scale=px.colors.sequential.Blues)
fig.show()
fig.write_image("/Users/christopherhempel/Desktop/heatmap_norm_abun_cell_genus.png")

# Print euclidean distance to expected for transformed and non-transformed abundaces
fin.apply(lambda row : euclidean(row, exp), axis=1)
fin_clr_mr.apply(lambda row : euclidean(row, exp_clr_mr), axis=1)
