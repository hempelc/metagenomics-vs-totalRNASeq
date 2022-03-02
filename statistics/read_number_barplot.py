#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 2 Mar 2022

# This script generates a barplot for read numbers and % of duplicates

import pandas as pd
import plotly.express as px

# Specify number of reads
reads={"DNA 1 aquarium": 1355159, "DNA 2 aquarium": 1373310, "DNA 3 aquarium": 1571705, "RNA 1 aquarium": 1902388,
"RNA 2 aquarium": 1099851, "RNA 3 aquarium": 773067, "DNA FilCon aquarium": 99, "RNA FilCon aquarium": 2685,
"DNA ExtCon aquarium": 157, "RNA ExtCon aquarium": 1137, "DNA 1 mock": 817619, "DNA 2 mock": 644634,
"DNA 3 mock": 669382, "RNA 1 mock": 94633, "RNA 2 mock": 78149, "RNA 3 mock": 120144,
"DNA FilCon mock": 640, "RNA FilCon mock": 1551, "DNA ExtCon mock": 399, "RNA ExtCon mock": 887}

# Specify % of duplicated reads
percentage=[2.75,2.2,2.65,75.25,53.7,52.5,2.2,41.15,0.5,26.95,5.75,4.55,4.45,74.3,69.4,77.05,9.65,41.6,0.9,35.1]

# Make df
reads_df=pd.DataFrame(reads, index=["reads"]).transpose().reset_index()
reads_df["dupl_perc"]=percentage
reads_df["dupl_reads"]=reads_df["dupl_perc"]/100*reads_df["reads"]
reads_df["reads_left"]=reads_df["reads"]-reads_df["dupl_reads"]

# Melt df for stacked barplot
reads_df=pd.melt(reads_df, id_vars="index", value_vars=["dupl_reads", "reads_left"]).round(0)
reads_df=pd.concat([reads_df[reads_df["variable"]=="reads_left"], reads_df[reads_df["variable"]=="dupl_reads"]])

# Plot
fig=px.bar(reads_df, x="index", y="value", text="value", color="variable", labels={"value": "Reads", "index": "Sample"})
fig.update_xaxes(tickangle=45)
fig.update_traces(textposition='outside')
fig.update_layout(width=1500, height=800)
