#!/usr/bin/python

# usage: ./merge_files.py taxonomy_file readcount_file output_file_name
# merges two files (mapped reads and contigs with taxonomy) on columns that exist in both files, keeps all rows, adds 0 for read counts if a contig has no mapped reads
# needs to have a column 'counts' with readcounts in readcount_file

import pandas as pd, sys

df1_name=sys.argv[1]
df2_name=sys.argv[2]
output_name=sys.argv[3]

df1 = pd.read_csv(df1_name, sep='\t')
df2 = pd.read_csv(df2_name, sep='\t')
df3 = pd.merge(df1, df2, how='outer')
df3['counts'].fillna(0, inplace=True)
df3.to_csv(output_name, sep='\t', na_rep='NA', index=False)
