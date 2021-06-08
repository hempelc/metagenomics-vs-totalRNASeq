#!/usr/bin/env python3
import pandas as pd

#Read back in each master dataframe
M4_RNA = pd.read_csv("M4_RNA.csv")
M5_RNA = pd.read_csv("M5_RNA.csv")
M6_RNA = pd.read_csv("M6_RNA.csv")

#Identify which pipelines ran successfully for every replicate
print(M4_RNA.shape)
print(M5_RNA.shape)
print(M6_RNA.shape)

pipelines = list(M4_RNA.columns[2:]) + list(M5_RNA.columns[2:]) + list(M6_RNA.columns[2:])
#Here is all pipeline options, now will only keep ones present in every
unique_pipelines = list(set(pipelines))
print(unique_pipelines[:4])
print(len(unique_pipelines))

common_pipelines = []
for name in unique_pipelines:
    if name in M4_RNA.columns and name in M5_RNA.columns and name in M6_RNA.columns:
        common_pipelines.append(name)

print(len(common_pipelines))

#Using set rearranges th eorder of the column names, order is not relavent except for the first two columns,
#so they were excluded from this intitally filtering and will be manually included into the new dataframes to ensure
#they are present and located where expected

M4_RNA_filtered = M4_RNA.iloc[:,0:2].join(M4_RNA[common_pipelines])
M5_RNA_filtered = M5_RNA.iloc[:,0:2].join(M5_RNA[common_pipelines])
M6_RNA_filtered = M6_RNA.iloc[:,0:2].join(M6_RNA[common_pipelines])

print(M4_RNA_filtered.shape)
print(M5_RNA_filtered.shape)
print(M6_RNA_filtered.shape)




