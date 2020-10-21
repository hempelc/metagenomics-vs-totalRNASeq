#!/usr/bin/python3

# A simple script to merge two files on columns that exist in both files, keeps
# all rows, and adds 0 for read counts if a contig has no mapped reads.
# Needs to have a column 'counts' with readcounts in readcount_file.
# If the first file is empty, the script merges on the colulm sequence_name and
# keeps the column order that way.

# usage: ./merge_files.py file2 file2 output_file_name

import pandas as pd, sys

df1_name=sys.argv[1]
df2_name=sys.argv[2]
output_name=sys.argv[3]

df1 = pd.read_csv(df1_name, sep='\t')
df2 = pd.read_csv(df2_name, sep='\t')
if df1.empty==True:
    df3 = pd.merge(df1, df2, how='outer', right_on="sequence_name", left_index=True)
    if "sequence_name_x" in df3.columns:
        df3 = df3.drop(["sequence_name_x", "sequence_name_y"], axis = 1)
else:
    df3 = pd.merge(df1, df2, how='outer')
df3['counts'].fillna(0, inplace=True)
df3.to_csv(output_name, sep='\t', na_rep='NA', index=False)
