#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 2 Mar 2022

# This script generates a barplot for read numbers and % of duplicates.

import os
import pandas as pd #v1.3.5
import plotly.express as px #v5.5.0

# Set output dir:
outdir="/Users/christopherhempel/Desktop/"

# Specify number of reads
reads = {"Mock DNA 1 R1": 817619, "Mock DNA 1 R2": 817619, "Mock RNA 1 R1": 94633,
    "Mock RNA 1 R2": 94633, "Mock DNA 2 R1": 644634, "Mock DNA 2 R2": 644634, "Mock RNA 2 R1": 78149,
    "Mock RNA 2 R2": 78149, "Mock DNA 3 R1": 669382, "Mock DNA 3 R2": 669382, "Mock RNA 3 R1": 120144,
    "Mock RNA 3 R2": 120144, "Mock DNA ExtCon R1": 399, "Mock DNA ExtCon R2": 399,
    "Mock RNA ExtCon R1": 887, "Mock RNA ExtCon R2": 887, "Mock DNA FilCon R1": 640,
    "Mock DNA FilCon R2": 640, "Mock RNA FilCon R1": 1551, "Mock RNA FilCon R2": 1551,
    "Aquarium DNA 1 R1": 1355159, "Aquarium DNA 1 R2": 1355159, "Aquarium RNA 1 R1": 1902388,
    "Aquarium RNA 1 R2": 1902388, "Aquarium DNA 2 R1": 1373310, "Aquarium DNA 2 R2": 1373310, "Aquarium RNA 2 R1": 1099851,
    "Aquarium RNA 2 R2": 1099851, "Aquarium DNA 3 R1": 1571705, "Aquarium DNA 3 R2": 1571705, "Aquarium RNA 3 R1": 773067,
    "Aquarium RNA 3 R2": 773067, "Aquarium DNA ExtCon R1": 99, "Aquarium DNA ExtCon R2": 99,
    "Aquarium RNA ExtCon R1": 2685, "Aquarium RNA ExtCon R2": 2685, "Aquarium DNA FilCon R1": 157,
    "Aquarium DNA FilCon R2": 157, "Aquarium RNA FilCon R1": 1137, "Aquarium RNA FilCon R2": 1137}

# Specify percentage of duplictes
percentage = {"Mock DNA 1 R1": 6.0, "Mock DNA 1 R2": 5.5, "Mock RNA 1 R1": 74.1,
    "Mock RNA 1 R2": 74.5, "Mock DNA 2 R1": 4.8, "Mock DNA 2 R2": 4.3, "Mock RNA 2 R1": 70.8,
    "Mock RNA 2 R2": 68.0, "Mock DNA 3 R1": 4.6, "Mock DNA 3 R2": 4.3, "Mock RNA 3 R1": 77.1,
    "Mock RNA 3 R2": 77.0, "Mock DNA ExtCon R1": 1.5, "Mock DNA ExtCon R2": 0.3,
    "Mock RNA ExtCon R1": 64.9, "Mock RNA ExtCon R2": 5.3, "Mock DNA FilCon R1": 12.3,
    "Mock DNA FilCon R2": 7.0, "Mock RNA FilCon R1": 74.7, "Mock RNA FilCon R2": 8.5,
    "Aquarium DNA 1 R1": 2.9, "Aquarium DNA 1 R2": 2.6, "Aquarium RNA 1 R1": 58.6,
    "Aquarium RNA 1 R2": 55.9, "Aquarium DNA 2 R1": 2.3, "Aquarium DNA 2 R2": 2.1, "Aquarium RNA 2 R1": 55.2,
    "Aquarium RNA 2 R2": 52.2, "Aquarium DNA 3 R1": 2.8, "Aquarium DNA 3 R2": 2.5, "Aquarium RNA 3 R1": 53.6,
    "Aquarium RNA 3 R2": 51.4, "Aquarium DNA ExtCon R1": 1, "Aquarium DNA ExtCon R2": 0,
    "Aquarium RNA ExtCon R1": 47.6, "Aquarium RNA ExtCon R2": 6.3, "Aquarium DNA FilCon R1": 3.8,
    "Aquarium DNA FilCon R2": 0.6, "Aquarium RNA FilCon R1": 74.6, "Aquarium RNA FilCon R2": 7.7}

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
fig.write_image(os.path.join(outdir, "readnum.svg"), height=800, width=1500)
