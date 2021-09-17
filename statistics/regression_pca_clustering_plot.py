#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 22 Jul 2021

# This script processes the output from the script "processing_and_metrics.py"

import pandas as pd
import numpy as np
import plotly.express as px
import statsmodels.api as sm
import rpy2.robjects as ro
import os
import copy
import logging
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.spatial.distance import euclidean
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# Activate logging for debugging
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s')


# Parameters
## Loop over all dirs and all metrics (True or False)
## Note: section 6 (clustering) requires manual adjustment for k so these have
## to be run manually:
looping=True

## If false, set your working dir and target metrics here:
### Full path to directory containing stats_exports from script "processing_and_metrics.py":
statsdir_lst=["/Users/christopherhempel/Desktop/pipeline_results_mock_community/results_genus_cell/stats_exports"]
### Metrics to use, options: "all", "rel", or "pa"
metr_list=["pa"]

if looping:
    statsdir_lst=["/Users/christopherhempel/Desktop/pipeline_results_mock_community/results_{0}/stats_exports".format(x) for x in ["genus_cell", "genus_gen", "species_cell", "species_gen"]]
    metr_list=["all", "rel", "pa"]

for statsdir in statsdir_lst:
    for metr in metr_list:
        # Imports
        ## Import value for expected TN
        with open(os.path.join(statsdir, "TN.pkl"), 'rb') as f:
            TN = pickle.load(f)

        ## Import metics df
        metrics_df=pd.read_csv(os.path.join(statsdir, "metrics_df.csv"), index_col=0)
        ### Convert trimming score column type into str
        metrics_df['trimming_score'] = metrics_df['trimming_score'].astype(str)
        ### Drop metrics that aren't used
        #### Default: drop TN and FN (collinear with TP and FP):
        metrics_df=metrics_df.drop(["FN", "TN", "FN_aad", "TN_aad"], axis=1)
        ### Drop additional ones specified by parameter "metr" and set expected dummy
        #### (=ideal pipeline results):
        if metr not in ["all", "rel", "pa"]:
            logging.critical("Parameter metr not in all_metr, rel_metr, or pa_metr: metr=" + metr)
        if metr=="all":
            print("All metrics are used.")
        elif metr=="rel":
            print("Just relative metrics are used.")
            metrics_df=metrics_df.drop(["FP", "TP", "FP_aad", "TP_aad"], axis=1)
        elif metr=="pa":
            print("Just p/a metrics are used.")
            drops=["FP", "TP", "FP_aad", "TP_aad", 'type', 'trimming_score',
                'rRNA_sorting_tool', 'assembly_tool', 'mapper', 'database', 'classifier']
            metrics_df=metrics_df.drop([x for x in metrics_df.columns if x not in drops], axis=1)
        ### Add expected
        expected_dummy=metrics_df.loc["expected"]
        metrics_df=metrics_df.drop(["expected"], axis=0)

        # 1 PCA (following https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60)
        ## 1.1 Perform initial PCA
        ### Grab columns we need
        pca_df=metrics_df.iloc[:, :-7]
        ### Save features=metrics in variable
        pca_features=pca_df.columns
        ### Standardize data
        pca_df_std=StandardScaler().fit_transform(pca_df)

        ### Perform PCA
        #### On one component for linear regression (later in script)
        pca_lr = PCA(n_components=1)
        pc_lr = pca_lr.fit_transform(pca_df_std)
        PC1_lr=str(pca_lr.explained_variance_ratio_[0])
        print("Linear regression PC1: " + PC1_lr)

        #### On two components for graphical visualization
        pca_graph = PCA(n_components=2)
        pcs_no_exp = pca_graph.fit_transform(pca_df_std)
        pc_no_exp_df = pd.DataFrame(data = pcs_no_exp, columns = ['PC1', 'PC2'],
            index=pca_df.index.to_list())
        PC1_graph = str(round(pca_graph.explained_variance_ratio_[0] * 100, 2))
        PC2_graph = str(round(pca_graph.explained_variance_ratio_[1] * 100, 2))
        print("Graph PC1: {0}%, Graph PC2: {1}%".format(PC1_graph, PC2_graph))
        ##### Save the loadings=arrowhead coordinates of each feature in variable
        loadings = pca_graph.components_.T * np.sqrt(pca_graph.explained_variance_)

        ## 1.2 Add a dummy variable respresenting an ideal pipeline and predict the
        ##     coordinates of all pipelines and the expected dummy in the PCA
        #### Add a dummy column with expected outcome (the expected outcome for TN is
        #### the number of all taxa - the number of expected taxa, generated in previous script)
        metrics_df_exp=metrics_df.transpose()
        metrics_df_exp["expected_dummy"]=expected_dummy
        metrics_df_exp=metrics_df_exp.transpose()

        ### Add a column to indicate which pipeline is the expected dummy
        col=["pipeline"] * (len(metrics_df_exp)-1)
        col.append("expected")
        metrics_df_exp["exp"]=col

        ### Predict the coordinates in the PCA for all pipelines plus the expected
        ### (that way, the expected doesn't affect the PCA itself)
        pca_df_with_exp=metrics_df_exp.iloc[:, 0:-8]
        pca_df_with_exp_std=StandardScaler().fit_transform(pca_df_with_exp)
        # Extract the expected standardized values for calculation of euc dist:
        exp_std=pca_df_with_exp_std[len(pca_df_with_exp_std)-1]
        pcs_with_exp=pca_graph.transform(pca_df_with_exp_std)
        pc_df = pd.DataFrame(data = pcs_with_exp, columns = ['PC1', 'PC2'],
            index=pca_df_with_exp.index.to_list())



        # # 2 Linear regression (following https://datatofish.com/statsmodels-linear-regression/)
        # ## 2.1 Set up X and y
        # y = pc_lr
        # ### Parameter to include or exclude type variables ("False" or "True")
        # type_included=True
        # if type_included:
        #     X_no_dummies=metrics_df_exp.drop("expected_dummy").drop("exp", axis=1).iloc[:, -7:]
        #     X_no_intercept=pd.get_dummies(X_no_dummies)
        # else:
        #     X_no_dummies=metrics_df_exp.drop("expected_dummy").drop(["exp", "type"], axis=1).iloc[:, -6:]
        #     X_no_intercept=pd.get_dummies(X_no_dummies)
        # ### Add a constant as a column which is needed to calculate the intercept of the model
        #
        # X = sm.add_constant(X_no_intercept)
        #
        # ## 2.2 Fit an ordinary least squares regression model
        # lr_model = sm.OLS(y, X).fit()
        #
        # ## 2.3 Check out the coefficients and p-values and save them
        # ### Set pd so that it only shows the last 3 decimal digits in the editor, which
        # ### is easier for manual inspection
        # pd.set_option('display.float_format', lambda x: '%.3f' % x)
        # ### Use the summary function, which automatically summarizes the output
        # summary = lr_model.summary2().tables[1][['Coef.', 'P>|t|']]
        # summary.rename(columns={"Coef.": "Coefficient", "P>|t|": "p-value"})
        # summary.to_csv(os.path.join(statsdir,  metr + "_lr_coeff_pval.csv"), index_label="tool")



        # 3 R envfit function - One very helpful function in R, envfit, is not
        # reproducible in python. Therefore, we use the rpy2 module to run that R
        # function in python on our PCA and metrics

        ## import R packages
        base = importr('base')
        vegan = importr('vegan')

        ## Transform metrics tools to dummy variables:
        metrics_df_dummies=pd.get_dummies(metrics_df_exp.drop("expected_dummy")\
            .drop("exp", axis=1).iloc[:, -7:])

        ## Open localconverter to be able to use pandas df in R function
        with localconverter(ro.default_converter + pandas2ri.converter):
            ### Run R function envfit on our PCA and metrics
            fit_steps=vegan.envfit(pc_no_exp_df, metrics_df.iloc[:, -7:])
            fit_tools=vegan.envfit(pc_no_exp_df, metrics_df_dummies.astype(str))

        ## Get p-values for tools from envfit output as R objects
        pvalues_steps_r=base.as_data_frame(base.as_list(fit_steps[1][3]))

        ## Get scores=arrowhead coordinates and p-values for tools from envfit output as R objects
        scores_tools_r=base.as_data_frame(vegan.scores(fit_steps, "factors"))
        pvalues_tools_r=base.as_data_frame(base.as_list(fit_tools[1][3]))

        ## Translate R objects into pandas dfs
        with localconverter(ro.default_converter + pandas2ri.converter):
          scores_tools=ro.conversion.rpy2py(scores_tools_r)
          pvalues_steps=ro.conversion.rpy2py(pvalues_steps_r).transpose().rename(columns = {'1':'p-value'})
          pvalues_tools=ro.conversion.rpy2py(pvalues_tools_r).transpose().rename(columns = {'1':'p-value'})

        ## Save p-value dfs
        pvalues_steps.to_csv(os.path.join(statsdir, metr + "_pvalues_steps.csv"), index_label="step")
        pvalues_tools.to_csv(os.path.join(statsdir, metr + "_pvalues_tools.csv"), index_label="tool")



        # 4 Euclidean distances between pipelines and expected (note: use standardized columns):
        euc_dist_df = copy.deepcopy(metrics_df)
        ## Calculate euc dist for each pipeline
        for index,row in pd.DataFrame(pca_df_std, index=pca_df.index).iterrows():
            euc_dist_df.loc[index,'euc_dist'] = euclidean(row, exp_std)
        #euc_dist_df.sort_values(by=['euc_dist'])
        mean_euc_dist_lst=[]
        tools=[]
        ## Calculate average euc dist for each tool in each step
        for step in euc_dist_df.loc[:, 'type':'classifier'].columns:
            for tool in euc_dist_df[step].unique():
                #### Cut down df to rows containing tool and take the averge euc_dist:
                mean_euc_dist=euc_dist_df.reset_index()[euc_dist_df.reset_index()['pipeline']
                    .str.contains(tool)]["euc_dist"].mean()
                tools.append(tool)
                mean_euc_dist_lst.append(mean_euc_dist)
        ## Make and save df
        mean_euc_dist_df=pd.DataFrame({"mean_euc_dist": mean_euc_dist_lst}, index=tools)
        mean_euc_dist_df.to_csv(os.path.join(statsdir, metr + "_mean_euc_dist_steps.csv"), index_label="tools")



        # 5 Plots
        ## Make saving directory:
        plotdir=os.path.join(os.path.abspath(os.path.join(statsdir, os.pardir)), "plots_" + metr)
        if not os.path.exists(plotdir):
            os.mkdir(plotdir)
        ## Make df from PC1+2 coordinates and tools per pipeline to plot:
        plot_df=pd.concat([pc_df, metrics_df_exp.iloc[:, -8:]],
            axis = 1).rename_axis("pipeline").reset_index()

        ## Pick scaling multiplier for arrows in PCA biplot to make them more visible:
        scaling=4
        ## Make plot for each tool column and save plots (biplot part taken from
        ## https://plotly.com/python/pca-visualization; arrow part taken from
        ## https://community.plotly.com/t/arrow-heads-at-the-direction-of-line-arrows/32565/3):
        for col in plot_df.columns[3:]:
            ### We show different arrows based on the column we colour on:
            if col=="exp":
                fig = px.scatter(plot_df, x="PC1", y="PC2",
                    labels={"PC1": "PC1 ({0}%)".format(PC1_graph),
                    "PC2": "PC2 ({0}%)".format(PC2_graph)}, color=col, hover_data=["pipeline"])
                #### Loop over features:
                for i, feature in enumerate(pca_features):
                    ##### Add lines:
                    fig.add_shape(type='line', x0=0, y0=0, x1=loadings[i, 0]*scaling, y1=loadings[i, 1]*scaling)
                    ##### Add arrowheads:
                    fig.add_annotation(x=loadings[i, 0]*scaling, y=loadings[i, 1]*scaling, ax=0, ay=0, xref='x',
                        yref='y', axref='x', ayref='y', text='', showarrow=True, arrowhead=2,
                        arrowsize=1.3, arrowwidth=1.3, arrowcolor='black')
                    ##### Add labels:
                    fig.add_annotation(x=loadings[i, 0]*scaling, y=loadings[i, 1]*scaling, ax=0, ay=0,
                        xanchor="center", yanchor="bottom", text=feature)
                    # #### Highlight one specific point:
                    # idx=plot_df.loc[plot_df['pipeline'] == 'RNA_10_rrnafilter_transabyss_bwa_ncbi-nt_kraken2'].index
                    # fig.data[0].update(selectedpoints=idx,
                    #        selected=dict(marker=dict(color='red')),#color of selected points
                    #        unselected=dict(marker=dict(color='rgb(200,200, 200)',#color of unselected pts
                    #                        opacity=0.9)));
            else:
                #### Extract p value and average euc dist for each tool from p value df
                #### and euc dist dfand replace tool names by tool names + p-value + dist
                #### so that p-values and dist are visible in legend
                plot_df_repl=copy.deepcopy(plot_df)
                for tool in plot_df[col].unique()[:-1]:
                    tool_rep=tool.replace("-", ".")
                    pval=pvalues_tools.reset_index()[pvalues_tools.reset_index()['index']\
                        .str.contains("_" + tool_rep)]["p-value"].to_list()[0]
                    eucdist=round(mean_euc_dist_df.loc[tool].to_list()[0], 2)
                    plot_df_repl=plot_df_repl.replace(tool, "{0} (p={1}, d={2})".format(tool, pval, eucdist))
                #### Plot figure with replaced tool names
                fig = px.scatter(plot_df_repl, x="PC1", y="PC2",
                    labels={"PC1": "PC1 ({0}%)".format(PC1_graph),
                    "PC2": "PC2 ({0}%)".format(PC2_graph), col: "{0} (p={1})".format(col,
                        pvalues_steps.loc[col].to_list()[0])}, color=col, hover_data=["pipeline"])
                #### Cut down scores_tools df to rows per col:
                tools_df=scores_tools.reset_index()[scores_tools.reset_index()['index']
                    .str.contains(col)]
                #### Loop over steps as above:
                for i, tool in enumerate(tools_df["index"]):
                                ##### Add lines:
                                fig.add_shape(type='line', x0=0, y0=0, x1=tools_df.iloc[i, 1]*scaling,
                                    y1=tools_df.iloc[i, 2]*scaling)
                                ##### Add arrowheads:
                                fig.add_annotation(x=tools_df.iloc[i, 1]*scaling, y=tools_df.iloc[i, 2]*scaling, ax=0, ay=0, xref='x',
                                    yref='y', axref='x', ayref='y', text='', showarrow=True, arrowhead=2,
                                    arrowsize=1.3, arrowwidth=1.3, arrowcolor='black')
                                ##### Add labels:
                                fig.add_annotation(x=tools_df.iloc[i, 1]*scaling, y=tools_df.iloc[i, 2]*scaling, ax=0, ay=0,
                                    xanchor="center", yanchor="bottom", text=tool)
            fig.show()
            fig.write_image(os.path.join(plotdir, "PCA_{0}_{1}.png".format(metr, col)))



        # 6 k-means clustering (following https://realpython.com/k-means-clustering-python/)
        kmean_df_std = pca_df_std

        ## 6.1 Estimate the best kluster number k
        ## 6.1.1 Elbow method based on SSE
        #### The list sse holds the SSE values for each k
        sse = []
        for k in range(1, 20):
            kmeans = KMeans(n_clusters=k)
            kmeans.fit(kmean_df_std)
            sse.append(kmeans.inertia_)

        ### Plot SSEs to choose the best k
        fig_k_sse = px.line(pd.DataFrame(list(zip(range(1, 20), sse)),
            columns = ['Number of Clusters', 'SSE']), x="Number of Clusters", y="SSE")
        fig_k_sse.show()
        fig_k_sse.write_image(os.path.join(plotdir, "KMeans_" + metr + "_sse.png"))

        ### 6.1.2 Or automatically calculate the best k using KneeLocator
        kl = KneeLocator(range(1, 20), sse, curve="convex", direction="decreasing")
        print("The best k based on KneeLocator is " + str(kl.elbow))

        ## 6.1.3 Silhouette coefficient
        silhouette_coefficients = []
        ### Start at k=2 when using the silhouette coefficient
        for k in range(2, 20):
            kmeans = KMeans(n_clusters=k)
            kmeans.fit(kmean_df_std)
            score = silhouette_score(kmean_df_std, kmeans.labels_)
            silhouette_coefficients.append(score)

        ### Plot the scores and manually pick the k with the highest score close to
        ### the k estimated with SSE
        fig_k_silhouette = px.line(pd.DataFrame(list(zip(range(2, 20),
            silhouette_coefficients)), columns =['Number of Clusters', 'Silhouette Coefficients']),
            x="Number of Clusters", y="Silhouette Coefficients")
        fig_k_silhouette.show()
        fig_k_silhouette.write_image(os.path.join(plotdir, "KMeans_"  + metr + "_silhouette.png"))

        ## 6.1.2 Based on both tests, we manually define k
        k=10

        ## 6.2 Perform the actual kmeans with the estimated k
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(kmean_df_std)

        # Plot some stats of the model
        print("The lowest SSE value of the KMeans model was " + str(kmeans.inertia_))
        print("The number of iterations required to converge the model was " + str(kmeans.n_iter_))

        # 6.3 Make kmeans df and plot
        kmean_df=pd.DataFrame(pcs_no_exp, columns = ['PC1', 'PC2'],
            index=pca_df.index.to_list()).rename_axis("pipeline").reset_index()
        kmean_df['cluster']=[str(x) for x in kmeans.labels_.tolist()]
        exp_kmeans=pc_df.iloc[-1:].reset_index()
        exp_kmeans.rename(columns = {'index':'pipeline'}, inplace = True)
        exp_kmeans['cluster']=["expected"]
        kmean_df=kmean_df.append(exp_kmeans)
        fig_kmean = px.scatter(kmean_df, x="PC1", y="PC2", color="cluster", hover_data=["pipeline"])
        fig_kmean.show()
        fig_kmean.write_image(os.path.join(plotdir, "KMeans_"  + metr + "_plot.png"))
