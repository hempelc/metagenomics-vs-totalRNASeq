#!/usr/bin/env python3

import pandas as pd

df1_name="kraken2_translate_result_edited.txt"
df2_name="SILVA_paths_and_taxids.txt"
df1 = pd.read_csv(df1_name, sep='\t', names=["tax"])
df2 = pd.read_csv(df2_name, sep='\t', names=["tax", "ID"])
df3 = pd.merge(df1, df2, left_on="tax", right_on="tax", how='left')
df3['ID'].fillna(0, inplace=True)
df3=df3.astype({'ID': 'int32'})
df4=df3[["ID", "tax"]]
df4.to_csv("merged.txt", sep='\t', header=False, index=False)
