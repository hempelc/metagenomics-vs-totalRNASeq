#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 21 Sep 2021

# This script processes the output from the script "regression_pca_clustering.py"

import pandas as pd
import plotly.express as px
import logging
import os
import pickle
import numpy as np
import glob


# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters manual setting
workdir="/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_mock_samples/"
## Should cell reference be included? (True/False)
cell=False

# Parameters auto setting
## List of DNA and RNA mock community samples, replicates of 3; must equal names of directories in workdir that
## contain each sample's pipeline results:
aggs=["agg_reps_agg_type", "agg_reps_sep_type", "sep_reps_agg_type", "sep_reps_sep_type"]
levels=["cell_genus", "cell_species", "gen_genus", "gen_species"]
metr_lst=["rel", "pa"]
steps=["type", "trimming_score", "rRNA_sorting_tool", "assembly_tool", "mapper", "classifier", "database"]
tools=['type_DNA', 'type_RNA',
    'trimming_score_5', 'trimming_score_10', 'trimming_score_15', 'trimming_score_20',
    'rRNA_sorting_tool_barrnap', 'rRNA_sorting_tool_rrnafilter', 'rRNA_sorting_tool_sortmerna',
    'rRNA_sorting_tool_unsorted', 'assembly_tool_idba-ud', 'assembly_tool_idba-tran',
    'assembly_tool_megahit', 'assembly_tool_metaspades', 'assembly_tool_rnaspades',
    'assembly_tool_spades', 'assembly_tool_transabyss', 'mapper_bowtie2', 'mapper_bwa',
    'database_ncbi-nt', 'database_silva',
    'classifier_blast-filtered', 'classifier_blast-first-hit', 'classifier_kraken2']
plotdir_level1=os.path.join(workdir, "stats_summary_plots")
if not os.path.exists(plotdir_level1):
    os.mkdir(plotdir_level1)

# 1 Import dfs to evaluate tools and concatenate them, creating columns to separate
#   them, and import tools count dics to evaluate clusters:

for agg in aggs:
    tool_eval_master_df=pd.DataFrame()
    cluster_counts_all=pd.DataFrame({'tool':tools})
    cluster_counts_master_df=pd.DataFrame()
    plotdir_level2=os.path.join(workdir, "stats_summary_plots", agg)
    if not os.path.exists(plotdir_level2):
        os.mkdir(plotdir_level2)

    test=pd.DataFrame()


    for level in levels:
        for metr in metr_lst:

            with open(os.path.join(workdir, "metrics_" + level, "stats_" + metr, agg, "closest_cluster_tool_counts.pkl"), 'rb') as f:
                counts_dic = pickle.load(f)

            eucdist_pvalue_files = glob.glob(os.path.join(workdir, "metrics_" + level, "stats_" + metr, agg, "*.csv"))

            for eucdist_pvalue_file in eucdist_pvalue_files:
                df=pd.read_csv(eucdist_pvalue_file)
                combo=eucdist_pvalue_file.split("/")[-1].replace("eucdist_pvalues_tools_", "").replace(".csv", "")
                ## 1.1  Process the dfs
                #### Set one column with combinations:
                df["combination"]=["{0}_{1}_{2}".format(combo, level, metr)]*len(df)
                ### Set two columns to "no" by default:
                df["mean_euc_dist_lowest"]=["no"]*len(df)
                df["min_euc_dist_lowest"]=["no"]*len(df)
                df["euc_dist_category"]=["no"]*len(df)
                for step in steps:
                    #### Replace "mean_euc_dist_lowest" value of lowest mean tool with "yes" for each step
                    step_lowest_mean=df[df['tool'].str.contains(step)]["mean_euc_dist"].min()
                    df.loc[df["mean_euc_dist"] == step_lowest_mean , 'mean_euc_dist_lowest'] = "yes"
                    #### Replace "min_euc_dist_lowest" value of lowest min tool with "yes" for each step
                    step_lowest_min=df[df['tool'].str.contains(step)]["min_euc_dist"].min()
                    df.loc[df["min_euc_dist"] == step_lowest_min , 'min_euc_dist_lowest'] = "yes"
                tool_eval_master_df=pd.concat([tool_eval_master_df, df])
                tool_eval_master_df.loc[tool_eval_master_df["combination"]=="genus_gen_all"]

            ## 1.2 Process the count dics
            ### Resolve structure
            for combo, counts in counts_dic.items():
                counts_dic_res={}
                for step, dic in counts.items():
                    for tool, count in dic.items():
                        counts_dic_res[step + "_" + tool]=count
                counts_df=pd.DataFrame(counts_dic_res, index=["counts_abs"]).transpose()\
                    .reset_index().rename(columns={'index': 'tool'})
                ### Bring in order and add missing tools:
                counts_df=pd.merge(cluster_counts_all, counts_df, how='outer').fillna(0)
                ### Add combinations column:
                counts_df["combination"]=[combo]*len(counts_df)
                ### Generate relative counts for each step
                counts_df_rel=pd.DataFrame()
                for step in steps:
                    sub_df=counts_df[counts_df['tool'].str.contains(step)]
                    sub_df["counts_rel"]=sub_df["counts_abs"]/sub_df["counts_abs"].sum()
                    #### And add steps column for colors of bubbleplot later
                    sub_df["step"]=[step]*len(sub_df)
                    counts_df_rel=pd.concat([counts_df_rel, sub_df])
                ### Add to master df
                #cluster_counts_master_df = pd.concat([cluster_counts_master_df, counts_df_rel])
                test = pd.concat([test, counts_df_rel])


    ## 1.3 Further process tools evaluation master df:
    ### Drop rows with NaN
    tool_eval_master_df=tool_eval_master_df.dropna()
    #### Make col for significance category of p-value
    tool_eval_master_df.loc[tool_eval_master_df["p-value"] <= 0.001 , 'significance_cat'] = "p <= 0.001"
    tool_eval_master_df.loc[tool_eval_master_df["p-value"] > 0.001 , 'significance_cat'] = "0.05 > p > 0.01"
    tool_eval_master_df.loc[tool_eval_master_df["p-value"] > 0.05, 'significance_cat'] = "p > 0.05"

    #### Make col for size level p-value
    tool_eval_master_df.loc[tool_eval_master_df["p-value"] <= 0.001 , 'size_lvl'] = 1
    tool_eval_master_df.loc[tool_eval_master_df["p-value"] > 0.001 , 'size_lvl'] = 0.45
    tool_eval_master_df.loc[tool_eval_master_df["p-value"] > 0.05, 'size_lvl'] = 0.03

    #### Make col for category of lowest euc dist
    tool_eval_master_df.loc[(tool_eval_master_df["mean_euc_dist_lowest"] == "yes") & (tool_eval_master_df["min_euc_dist_lowest"] == "yes") , 'euc_dist_category'] = "both"
    tool_eval_master_df.loc[(tool_eval_master_df["mean_euc_dist_lowest"] == "no") & (tool_eval_master_df["min_euc_dist_lowest"] == "no") , 'euc_dist_category'] = "none"
    tool_eval_master_df.loc[(tool_eval_master_df["mean_euc_dist_lowest"] == "yes") & (tool_eval_master_df["min_euc_dist_lowest"] == "no") , 'euc_dist_category'] = "mean"
    tool_eval_master_df.loc[(tool_eval_master_df["mean_euc_dist_lowest"] == "no") & (tool_eval_master_df["min_euc_dist_lowest"] == "yes") , 'euc_dist_category'] = "min"



    # 2 Graph
    category_orders_tool_eval_tool=[x for x in tools if x in tool_eval_master_df["tool"].unique()]
    category_orders_cluster_counts=[x for x in tools if x in cluster_counts_master_df["tool"].unique()]



    if agg=="agg_reps_agg_type":
        width=1050
        category_orders_comb=tool_eval_master_df["combination"].unique()
    elif agg=="agg_reps_sep_type":
        width=1250
        category_orders_comb=sorted(tool_eval_master_df["combination"].unique())
    elif agg=="sep_reps_agg_type":
        width=1600
        category_orders_comb=[
            'M4_cell_genus_rel','M5_cell_genus_rel','M6_cell_genus_rel','M4_cell_genus_pa','M5_cell_genus_pa','M6_cell_genus_pa','M4_cell_species_rel','M5_cell_species_rel',
            'M6_cell_species_rel', 'M4_cell_species_pa','M5_cell_species_pa', 'M6_cell_species_pa', 'M4_gen_genus_rel','M5_gen_genus_rel',
            'M6_gen_genus_rel', 'M4_gen_genus_pa','M5_gen_genus_pa', 'M6_gen_genus_pa', 'M4_gen_species_rel','M5_gen_species_rel',
            'M6_gen_species_rel', 'M4_gen_species_pa','M5_gen_species_pa', 'M6_gen_species_pa']
    elif agg=="sep_reps_sep_type":
        width=2600
        category_orders_comb=[
            'M4_DNA_cell_genus_rel', 'M5_DNA_cell_genus_rel','M6_DNA_cell_genus_rel', 'M4_DNA_cell_genus_pa', 'M5_DNA_cell_genus_pa','M6_DNA_cell_genus_pa','M4_DNA_cell_species_rel','M5_DNA_cell_species_rel',
            'M6_DNA_cell_species_rel','M4_DNA_cell_species_pa','M5_DNA_cell_species_pa', 'M6_DNA_cell_species_pa', 'M4_DNA_gen_genus_rel','M5_DNA_gen_genus_rel','M6_DNA_gen_genus_rel',
            'M4_DNA_gen_genus_pa','M5_DNA_gen_genus_pa','M6_DNA_gen_genus_pa', 'M4_DNA_gen_species_rel','M5_DNA_gen_species_rel','M6_DNA_gen_species_rel', 'M4_DNA_gen_species_pa','M5_DNA_gen_species_pa','M6_DNA_gen_species_pa',
            'M4_RNA_cell_genus_rel','M5_RNA_cell_genus_rel','M6_RNA_cell_genus_rel','M4_RNA_cell_genus_pa','M5_RNA_cell_genus_pa','M6_RNA_cell_genus_pa','M4_RNA_cell_species_rel','M5_RNA_cell_species_rel','M6_RNA_cell_species_rel',
            'M4_RNA_cell_species_pa','M5_RNA_cell_species_pa','M6_RNA_cell_species_pa','M4_RNA_gen_genus_rel', 'M5_RNA_gen_genus_rel', 'M6_RNA_gen_genus_rel', 'M4_RNA_gen_genus_pa', 'M5_RNA_gen_genus_pa','M6_RNA_gen_genus_pa',
            'M4_RNA_gen_species_rel', 'M5_RNA_gen_species_rel','M6_RNA_gen_species_rel','M4_RNA_gen_species_pa', 'M5_RNA_gen_species_pa','M6_RNA_gen_species_pa']

    if cell==False:
        tool_eval_master_df=tool_eval_master_df[~tool_eval_master_df.combination.str.contains("cell")]
        cluster_counts_master_df=cluster_counts_master_df[~cluster_counts_master_df.combination.str.contains("cell")]
        category_orders_comb=[ x for x in category_orders_comb if "cell" not in x ]
        width=width/1.5

    ## 2.1 Tool evaluation:
    ### Invert pval for graph (otherwise bubbles with good p-value are not visible):
    #tool_eval_master_df["p-value"]=1-tool_eval_master_df["p-value"]
    fig_tool_eval = px.scatter(tool_eval_master_df, x="combination", y="tool",
    	         size="size_lvl", color="euc_dist_category",
                 hover_name="min_euc_dist", size_max=28, height=1350, width=width,
                 #symbol="mean_euc_dist_lowest", symbol_sequence=["circle", "circle-dot"],
                 category_orders={"tool": category_orders_tool_eval_tool, "combination": category_orders_comb},
                 color_discrete_map={'none': "lightgrey", 'min': "#0187b3", 'mean': "#e3c284", "both": "#cc0035"})
    fig_tool_eval.update_traces(marker=dict(line=dict(width=0,color='black')), selector=dict(mode='markers'))
    fig_tool_eval.write_image(os.path.join(plotdir_level2, agg + "_bubbleplot_tool_eval.svg"))
    fig_tool_eval.write_image(os.path.join(plotdir_level2, agg + "_bubbleplot_tool_eval.png"))

    ## 2.2.1 Cluster counts abs & rel:
    fig_cluster_counts_abs_rel = px.scatter(cluster_counts_master_df, x="combination", y="tool",
    	         size="counts_rel", color="counts_abs",
                 hover_name="counts_abs", size_max=30, height=1350, width=width,
                 #symbol="mean_euc_dist_lowest", symbol_sequence=["circle", "circle-dot"],
                 category_orders={"tool": category_orders_cluster_counts, "combination": category_orders_comb}, color_continuous_scale=["#7c8e00", "#5bdf72", "#edbf0a"])
                 #color_continuous_scale=px.colors.sequential.haline_r)
                 #color_discrete_sequence=["lightgrey", "#00ab56", "#dac751", "#737cd0"])
    fig_cluster_counts_abs_rel.update_traces(marker=dict(line=dict(width=0,color='black')), selector=dict(mode='markers'))
    fig_cluster_counts_abs_rel.write_image(os.path.join(plotdir_level2, agg + "_bubbleplot_cluster_counts_abs_rel.svg"))
    fig_cluster_counts_abs_rel.write_image(os.path.join(plotdir_level2, agg + "_bubbleplot_cluster_counts_abs_rel.png"))

    ## 2.2.2 Cluster counts rel:
    fig_cluster_counts_rel = px.scatter(cluster_counts_master_df, x="combination", y="tool",
    	         size="counts_rel", color="step",
                 hover_name="counts_rel", size_max=30, height=1350, width=width,
                 #symbol="mean_euc_dist_lowest", symbol_sequence=["circle", "circle-dot"],
                 category_orders={"tool": category_orders_cluster_counts, "combination": category_orders_comb}, color_discrete_sequence=["#3f47ac",
                    "#ff75b6", "#467600", "#cc0035", "#0187b3", "#01bd7c", "#e3c284"])
    fig_cluster_counts_rel.update_traces(marker=dict(line=dict(width=0,color='black')), selector=dict(mode='markers'))
    fig_cluster_counts_rel.write_image(os.path.join(plotdir_level2, agg + "_bubbleplot_cluster_counts_rel.svg"))
    fig_cluster_counts_rel.write_image(os.path.join(plotdir_level2, agg + "_bubbleplot_cluster_counts_rel.png"))
