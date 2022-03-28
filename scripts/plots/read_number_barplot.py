#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 2 Mar 2022

# This script generates a barplot for read numbers and % of duplicates

import pandas as pd #v1.3.5
import plotly.express as px #v5.5.0

# Specify number of reads
reads = {"DNA 1 R1": 817619, "DNA 1 R2": 817619, "RNA 1 R1": 94633,
    "RNA 1 R2": 94633, "DNA 2 R1": 644634, "DNA 2 R2": 644634, "RNA 2 R1": 78149,
    "RNA 2 R2": 78149, "DNA 3 R1": 669382, "DNA 3 R2": 669382, "RNA 3 R1": 120144,
    "RNA 3 R2": 120144, "DNA ExtCon R1": 399, "DNA ExtCon R2": 399,
    "RNA ExtCon R1": 887, "RNA ExtCon R2": 887, "DNA FilCon R1": 640,
    "DNA FilCon R2": 640, "RNA FilCon R1": 1551, "RNA FilCon R2": 1551}

# Specify percentage of duplictes
percentage = {"DNA 1 R1": 6.0, "DNA 1 R2": 5.5, "RNA 1 R1": 74.1,
    "RNA 1 R2": 74.5, "DNA 2 R1": 4.8, "DNA 2 R2": 4.3, "RNA 2 R1": 70.8,
    "RNA 2 R2": 68.0, "DNA 3 R1": 4.6, "DNA 3 R2": 4.3, "RNA 3 R1": 77.1,
    "RNA 3 R2": 77.0, "DNA ExtCon R1": 1.5, "DNA ExtCon R2": 0.3,
    "RNA ExtCon R1": 64.9, "RNA ExtCon R2": 5.3, "DNA FilCon R1": 12.3,
    "DNA FilCon R2": 7.0, "RNA FilCon R1": 74.7, "RNA FilCon R2": 8.5}

# Make df
reads_df=pd.DataFrame([reads, percentage], index=["reads", "dupl_perc"]).transpose().reset_index()
reads_df["dupl_reads"]=reads_df["dupl_perc"]/100*reads_df["reads"]
reads_df["reads_left"]=reads_df["reads"]-reads_df["dupl_reads"]

# Melt df for stacked barplot
reads_df_melted=pd.melt(reads_df, id_vars="index", value_vars=["dupl_reads", "reads_left"]).round(0)
reads_df_melted=pd.concat([reads_df_melted[reads_df_melted["variable"]=="reads_left"], reads_df_melted[reads_df_melted["variable"]=="dupl_reads"]])
reads_df_melted["reads"]=2*list(reads_df["reads"])

# Plot
fig=px.bar(reads_df_melted, x="index", y="value", text="reads", color="variable",
    labels={"value": "Reads", "index": "Sample"})
fig.update_traces(textposition='outside')
fig.write_image("/Users/christopherhempel/Desktop/readnum.svg", height=800, width=1500)
