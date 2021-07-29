#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 22 Jul 2021

import pandas as pd
import numpy as np
import plotly.express as px
import glob
import os
import csv
import copy
import pickle
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# Set dirs
workdir = "/Users/christopherhempel/Desktop/mock_community_RNA/" # Full path to directory that contains samples, HAS TO END WITH "/"
statsdir=os.path.join(workdir, "stats_exports")

# Import metics df
metrics_df=pd.read_csv(os.path.join(statsdir, "metrics_df.csv"), index_col=0)

# Import value for expected TN
with open(os.path.join(statsdir, "TN.pkl"), 'rb') as f:
    TN = pickle.load(f)

# 1 PCA (from https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60)
## 1.1 Perform initial PCA
###Grab columns we need
pca_df=metrics_df.iloc[:, 0:-6]
### Drop subceed and exceed metrics (decided to not use them)
pca_df=pca_df.drop(["subceed_reads", "exceed_reads", "subceed_reads_aad", "exceed_reads_aad"], axis=1)
### Standardize data
pca_df_std=StandardScaler().fit_transform(pca_df)
### Perform PCA
#### On one component for linear regression
pca_lr = PCA(n_components=1)
pc_lr = pca_lr.fit_transform(pca_df_std)
print("PC1: " + str(pca_lr.explained_variance_ratio_[0]))

#### On two components for graphical visualization
pca_graph = PCA(n_components=2)
pcs_no_exp = pca_graph.fit_transform(pca_df_std)
print("PC1: " + str(pca_graph.explained_variance_ratio_[0]) + ", PC2: " + str(pca_graph.explained_variance_ratio_[1]))


## 1.2 Add dummy and predict coordinates in PCA with all pipelines and expected dummy
#### Add a dummy column with expected outcome (the expected outcome for TN is the number of all taxa - the number of expected taxa)
metrics_df_exp=metrics_df.drop(["subceed_reads", "exceed_reads", "subceed_reads_aad", "exceed_reads_aad"], axis=1).transpose()
metrics_df_exp["expected_dummy"]=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, float(TN), 0.0, 0.0, 0.0, 0.0, "expected", "expected", "expected", "expected", "expected", "expected"]
metrics_df_exp=metrics_df_exp.transpose()

### Add column to indicate which pipeline is the expected dummy
col=["pipeline"] * (len(metrics_df_exp)-1)
col.append("expected")
metrics_df_exp["exp"]=col

### Predict the coordinates in the PCA for all pipelines plus the expected (that way, the expected doesn't affect the PCA itself)
pca_df_with_exp=metrics_df_exp.iloc[:, 0:-7]
pca_df_with_exp_std=StandardScaler().fit_transform(pca_df_with_exp)
pcs_with_exp=pca_graph.transform(pca_df_with_exp_std)
pc_df = pd.DataFrame(data = pcs_with_exp, columns = ['PC1', 'PC2'], index=pca_df_with_exp.index.to_list())


# 2 Plot
## Make saving directory
plotdir=os.path.join(workdir, "plots")
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
## Make df from PC1+2 coordinates and tools per pipeline to plot
plot_df=pd.concat([pc_df, metrics_df_exp.iloc[:, -7:]], axis = 1).rename_axis("pipeline").reset_index()

fig1 = px.scatter(plot_df, x="PC1", y="PC2", color="exp", hover_data=["pipeline"])
fig1.show()
fig1.write_image(os.path.join(plotdir, "raw.png"))

fig2 = px.scatter(plot_df, x="PC1", y="PC2", color="trimming_score", hover_data=["pipeline"])
fig2.show()
fig2.write_image(os.path.join(plotdir, "trimming_score.png"))

fig3 = px.scatter(plot_df, x="PC1", y="PC2", color="rRNA_sorting_tool", hover_data=["pipeline"])
fig3.show()
fig3.write_image(os.path.join(plotdir, "rRNA_sorting_tool.png"))

fig4 = px.scatter(plot_df, x="PC1", y="PC2", color="assembly_tool", hover_data=["pipeline"])
fig4.show()
fig4.write_image(os.path.join(plotdir, "assembly_tool.png"))

fig5 = px.scatter(plot_df, x="PC1", y="PC2", color="mapper", hover_data=["pipeline"])
fig5.show()
fig5.write_image(os.path.join(plotdir, "mapper.png"))

fig6 = px.scatter(plot_df, x="PC1", y="PC2", color="database", hover_data=["pipeline"])
fig6.show()
fig6.write_image(os.path.join(plotdir, "database.png"))

fig7 = px.scatter(plot_df, x="PC1", y="PC2", color="classifier", hover_data=["pipeline"])
fig7.show()
fig7.write_image(os.path.join(plotdir, "classifier.png"))


# 3 Linear regression (following https://datatofish.com/statsmodels-linear-regression/)
## 3.1 Set up X and y
y = pc_lr
X_no_dummies=metrics_df_exp.drop("expected_dummy").drop("exp", axis=1).iloc[:, -6:]
X_no_intercept=pd.get_dummies(data=X_no_dummies)
X = sm.add_constant(X_no_intercept) # This adds a constant as a column which is needed to calculate the intercept of the model


## 3.2 Fit an ordinary least squares regression model
lr_model = sm.OLS(y, X).fit()

# 3.3 Check out the coefficients and p-values
pd.set_option('display.float_format', lambda x: '%.3f' % x) # Only show last 3 decimal digits
summary = lr_model.summary2().tables[1][['Coef.', 'P>|t|']]
summary.rename(columns={"Coef.": "Coefficient", "P>|t|": "p-value"})
summary.to_csv(os.path.join(statsdir, "lr_coeff_pval.csv"), index_label="tool")



# 4 k-means clustering (following https://realpython.com/k-means-clustering-python/)
kmean_table_std = pca_df_std

## 4.1 Estimate best kluster number k
## 4.1.1 Elbow method based on SSE
#### List sse holds the SSE values for each k
sse = []
for k in range(1, 20):
    kmeans = KMeans(n_clusters=k)
    kmeans.fit(kmean_table_std)
    sse.append(kmeans.inertia_)

### Plot SSEs to choose best k
fig_k_sse = px.line(pd.DataFrame(list(zip(range(1, 20), sse)), columns =['Number of Clusters', 'SSE']), x="Number of Clusters", y="SSE")
fig_k_sse.show()
fig_k_sse.write_image(os.path.join(plotdir, "kmean_sse.png"))

### 4.1.2 Or calculate using KneeLocator
kl = KneeLocator(range(1, 20), sse, curve="convex", direction="decreasing")
print(kl.elbow)

## 4.1.3 Silhouette coefficient
silhouette_coefficients = []
### start at 2 clusters for silhouette coefficient
for k in range(2, 20):
    kmeans = KMeans(n_clusters=k)
    kmeans.fit(kmean_table_std)
    score = silhouette_score(kmean_table_std, kmeans.labels_)
    silhouette_coefficients.append(score)

fig_k_silhouette = px.line(pd.DataFrame(list(zip(range(2, 20), silhouette_coefficients)), columns =['Number of Clusters', 'Silhouette Coefficients']), x="Number of Clusters", y="Silhouette Coefficients")
fig_k_silhouette.show()
fig_k_silhouette.write_image(os.path.join(plotdir, "kmean_silhouette.png"))

# Based on both tests, we define k
k=4

## 4.2 Actual kmeans
kmeans = KMeans(n_clusters=k)
kmeans.fit(kmean_table_std)

# Plot soem stats of the model
print(kmeans.inertia_) # The lowest SSE value
print(kmeans.n_iter_) # The number of iterations required to converge

# 4.3 Make kmeans df and plot
kmean_df=pd.DataFrame(pcs_no_exp, columns = ['PC1', 'PC2'], index=pca_df.index.to_list()).rename_axis("pipeline").reset_index()
kmean_df['cluster']=[str(x) for x in kmeans.labels_.tolist()]
exp_kmeans=pc_df.iloc[-1:].reset_index()
exp_kmeans.rename(columns = {'index':'pipeline'}, inplace = True)
exp_kmeans['cluster']=["expected"]
kmean_df=kmean_df.append(exp_kmeans)
fig_kmean = px.scatter(kmean_df, x="PC1", y="PC2", color="cluster", hover_data=["pipeline"])
fig_kmean.show()
fig_kmean.write_image(os.path.join(plotdir, "kmean_plot.png"))


######## INSERT: PCA with rel abun

## Calculate the average across the three replicates for rel abundances:
rel_df = pd.DataFrame()
for pipeline in master_dfs_uniq_rel[sample].columns:
    ave_df_rel=pd.DataFrame()
    for sample in samples:
        ave_df_rel[sample]=master_dfs_uniq_rel[sample][pipeline]
    rel_df[pipeline]=ave_df_rel.mean(axis=1)
rel_df = rel_df.transpose().add_suffix('_rel')

## Calculate the average across the three replicates for rel abundances:
pa_df = pd.DataFrame()
for pipeline in master_dfs_uniq_pa[sample].columns:
    ave_df_pa=pd.DataFrame()
    for sample in samples:
        ave_df_pa[sample]=master_dfs_uniq_pa[sample][pipeline]
    pa_df[pipeline]=ave_df_pa.mean(axis=1)
pa_df = pa_df.transpose().add_suffix('_pa')

combined_df = pd.concat([rel_df, pa_df], axis=1)
## PCA (from https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60)
## 1.1 Perform initial PCA
###Grab columns we need
pca_df=combined_df.iloc[1:,:]
### Standardize data
pca_df_std=StandardScaler().fit_transform(pca_df)
### Perform PCA
#### On one component for linear regression
pca_lr = PCA(n_components=1)
pc_lr = pca_lr.fit_transform(pca_df_std)
print("PC1: " + str(pca_lr.explained_variance_ratio_[0]))

#### On two components for graphical visualization
pca_graph = PCA(n_components=2)
pcs_no_exp = pca_graph.fit_transform(pca_df_std)
print("PC1: " + str(pca_graph.explained_variance_ratio_[0]) + ", PC2: " + str(pca_graph.explained_variance_ratio_[1]))

## 1.2 Add dummy and predict coordinates in PCA with all pipelines and expected dummy
#### Add a dummy column with expected outcome (the expected outcome for TN is the number of all taxa - the number of expected taxa)
pca_df_t=pca_df.transpose()
pca_df_t["expected_dummy"]=combined_df.iloc[0,:]
pca_df_with_exp=pca_df_t.transpose()

# ### Add column to indicate which pipeline is the expected dummy
# col=["pipeline"] * (len(pca_df_with_exp)-1)
# col.append("expected")
# pca_df_with_exp["exp"]=col

### Predict the coordinates in the PCA for all pipelines plus the expected (that way, the expected doesn't affect the PCA itself)
pca_df_with_exp_std=StandardScaler().fit_transform(pca_df_with_exp)
pcs_with_exp=pca_graph.transform(pca_df_with_exp_std)
pc_df = pd.DataFrame(data = pcs_with_exp, columns = ['PC1', 'PC2'], index=pca_df_with_exp.index)

plot_df=pd.concat([pc_df, metrics_df.iloc[:, 16:23]], axis = 1).rename_axis("pipeline").reset_index()

fig1 = px.scatter(plot_df, x="PC1", y="PC2", color="exp", hover_data=["pipeline"])
fig1.show()

fig2 = px.scatter(plot_df, x="PC1", y="PC2", color="trimming_score", hover_data=["pipeline"])
fig2.show()

fig3 = px.scatter(plot_df, x="PC1", y="PC2", color="rRNA_sorting_tool", hover_data=["pipeline"])
fig3.show()

fig4 = px.scatter(plot_df, x="PC1", y="PC2", color="assembly_tool", hover_data=["pipeline"])
fig4.show()

fig5 = px.scatter(plot_df, x="PC1", y="PC2", color="mapper", hover_data=["pipeline"])
fig5.show()

fig6 = px.scatter(plot_df, x="PC1", y="PC2", color="database", hover_data=["pipeline"])
fig6.show()

fig7 = px.scatter(plot_df, x="PC1", y="PC2", color="classifier", hover_data=["pipeline"])
fig7.show()


########## INSERT END




# Exports
summary.to_csv('lr.csv')
# Export tables to csv
metrics_df.iloc[:, 16:22].to_csv('tools.csv')
metrics_df.iloc[:, 0:16].to_csv('metrics.csv')

###################### FOR PCA WITHOUT AAD METRICS:

## Make the df
metrics_df=pd.DataFrame(master_metrics_df)
metrics_df=metrics_df.transpose()
metrics_df=metrics_df.drop(["subceed_reads", "exceed_reads", "true_reads", "false_reads", "FN", "FP", "TP", "TN"], axis=1)
############################################ Loop over DNA and RNA samples would go until here

# 6 PCA (from https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60)
## 2.1 Perform initial PCA
###Grab columns we need
pca_df=metrics_df.iloc[:, 0:8]
### Standardize data
pca_df_std=StandardScaler().fit_transform(pca_df)
### Perform PCA
pca = PCA(n_components=2)
pcs_no_exp = pca.fit_transform(pca_df_std)
print("PC1: " + str(pca.explained_variance_ratio_[0]) + ", PC2: " + str(pca.explained_variance_ratio_[1]))


## 2.2 Add dummy and predict coordinates in PCA with all pipelines and expected dummy
#### Add a dummy column with expected outcome (the expected outcome for TN is the number of all taxa - the number of expected taxa)
metrics_df=metrics_df.transpose()
metrics_df["expected_dummy"]=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, len(unique_taxa)-len(expected_df.index), "expected", "expected", "expected", "expected", "expected", "expected"]
metrics_df=metrics_df.transpose()

### Add column to indicate which pipeline is the expected dummy
col=["pipeline"] * (len(metrics_df)-1)
col.append("expected")
metrics_df["exp"]=col

### Predict the coordinates in the PCA for all pipelines plus the expected (that way, the expected doesn't affekt the PCA itself)
pca_df_with_exp=metrics_df.iloc[:, 0:8]
pca_df_with_exp_std=StandardScaler().fit_transform(pca_df_with_exp)
pcs_with_exp=pca.transform(pca_df_with_exp_std)
pc_df = pd.DataFrame(data = pcs_with_exp, columns = ['PC1', 'PC2'], index=pca_df_with_exp.index)



# 7 Plot
## Make df from PC1+2 coordinates and tools per pipeline to plot
plot_df=pd.concat([pc_df, metrics_df.iloc[:, 8:15]], axis = 1).rename_axis("pipeline").reset_index()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="exp", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="trimming_score", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="rRNA_sorting_tool", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="assembly_tool", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="mapper", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="database", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="classifier", hover_data=["pipeline"])
fig.show()


###################### FOR PCA WITHOUT ACCURACY METRICS:

## Make the df
metrics_df=pd.DataFrame(master_metrics_df)
metrics_df=metrics_df.transpose()
metrics_df=metrics_df.drop(["subceed_reads_aad", "exceed_reads_aad", "true_reads_aad", "false_reads_aad", "FN_aad", "FP_aad", "TP_aad", "TN_aad"], axis=1)
############################################ Loop over DNA and RNA samples would go until here

# 6 PCA (from https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60)
## 2.1 Perform initial PCA
###Grab columns we need
pca_df=metrics_df.iloc[:, 0:8]
### Standardize data
pca_df_std=StandardScaler().fit_transform(pca_df)
### Perform PCA
pca = PCA(n_components=2)
pcs_no_exp = pca.fit_transform(pca_df_std)
print("PC1: " + str(pca.explained_variance_ratio_[0]) + ", PC2: " + str(pca.explained_variance_ratio_[1]))


## 2.2 Add dummy and predict coordinates in PCA with all pipelines and expected dummy
#### Add a dummy column with expected outcome (the expected outcome for TN is the number of all taxa - the number of expected taxa)
metrics_df=metrics_df.transpose()
metrics_df["expected_dummy"]=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, len(unique_taxa)-len(expected_df.index), "expected", "expected", "expected", "expected", "expected", "expected"]
metrics_df=metrics_df.transpose()

### Add column to indicate which pipeline is the expected dummy
col=["pipeline"] * (len(metrics_df)-1)
col.append("expected")
metrics_df["exp"]=col

### Predict the coordinates in the PCA for all pipelines plus the expected (that way, the expected doesn't affekt the PCA itself)
pca_df_with_exp=metrics_df.iloc[:, 0:8]
pca_df_with_exp_std=StandardScaler().fit_transform(pca_df_with_exp)
pcs_with_exp=pca.transform(pca_df_with_exp_std)
pc_df = pd.DataFrame(data = pcs_with_exp, columns = ['PC1', 'PC2'], index=pca_df_with_exp.index)



# 7 Plot
## Make df from PC1+2 coordinates and tools per pipeline to plot
plot_df=pd.concat([pc_df, metrics_df.iloc[:, 8:15]], axis = 1).rename_axis("pipeline").reset_index()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="exp", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="trimming_score", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="rRNA_sorting_tool", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="assembly_tool", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="mapper", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="database", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="classifier", hover_data=["pipeline"])
fig.show()

###################### FOR PCA WITHOUT RELATIVE READS:

## Make the df
metrics_df=pd.DataFrame(master_metrics_df)
metrics_df=metrics_df.transpose()
metrics_df=metrics_df.drop(["subceed_reads_aad", "exceed_reads_aad", "true_reads_aad", "false_reads_aad", "subceed_reads", "exceed_reads", "true_reads", "false_reads"], axis=1)
############################################ Loop over DNA and RNA samples would go until here

# 6 PCA (from https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60)
## 2.1 Perform initial PCA
###Grab columns we need
pca_df=metrics_df.iloc[:, 0:8]
### Standardize data
pca_df_std=StandardScaler().fit_transform(pca_df)
### Perform PCA
pca = PCA(n_components=2)
pcs_no_exp = pca.fit_transform(pca_df_std)
print("PC1: " + str(pca.explained_variance_ratio_[0]) + ", PC2: " + str(pca.explained_variance_ratio_[1]))


## 2.2 Add dummy and predict coordinates in PCA with all pipelines and expected dummy
#### Add a dummy column with expected outcome (the expected outcome for TN is the number of all taxa - the number of expected taxa)
metrics_df=metrics_df.transpose()
metrics_df["expected_dummy"]=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, len(unique_taxa)-len(expected_df.index), "expected", "expected", "expected", "expected", "expected", "expected"]
metrics_df=metrics_df.transpose()

### Add column to indicate which pipeline is the expected dummy
col=["pipeline"] * (len(metrics_df)-1)
col.append("expected")
metrics_df["exp"]=col

### Predict the coordinates in the PCA for all pipelines plus the expected (that way, the expected doesn't affekt the PCA itself)
pca_df_with_exp=metrics_df.iloc[:, 0:8]
pca_df_with_exp_std=StandardScaler().fit_transform(pca_df_with_exp)
pcs_with_exp=pca.transform(pca_df_with_exp_std)
pc_df = pd.DataFrame(data = pcs_with_exp, columns = ['PC1', 'PC2'], index=pca_df_with_exp.index)



# 7 Plot
## Make df from PC1+2 coordinates and tools per pipeline to plot
plot_df=pd.concat([pc_df, metrics_df.iloc[:, 8:15]], axis = 1).rename_axis("pipeline").reset_index()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="exp", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="trimming_score", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="rRNA_sorting_tool", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="assembly_tool", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="mapper", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="database", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="classifier", hover_data=["pipeline"])
fig.show()

###################### FOR PCA WITHOUT P/A DATA:

## Make the df
metrics_df=pd.DataFrame(master_metrics_df)
metrics_df=metrics_df.transpose()
metrics_df=metrics_df.drop(["FP_aad", "TP_aad", "FN_aad", "TN_aad", "FP", "TP", "FN", "TN"], axis=1)
############################################ Loop over DNA and RNA samples would go until here

# 6 PCA (from https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60)
## 2.1 Perform initial PCA
###Grab columns we need
pca_df=metrics_df.iloc[:, 0:8]
### Standardize data
pca_df_std=StandardScaler().fit_transform(pca_df)
### Perform PCA
pca = PCA(n_components=2)
pcs_no_exp = pca.fit_transform(pca_df_std)
print("PC1: " + str(pca.explained_variance_ratio_[0]) + ", PC2: " + str(pca.explained_variance_ratio_[1]))


## 2.2 Add dummy and predict coordinates in PCA with all pipelines and expected dummy
#### Add a dummy column with expected outcome (the expected outcome for TN is the number of all taxa - the number of expected taxa)
metrics_df=metrics_df.transpose()
metrics_df["expected_dummy"]=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, len(unique_taxa)-len(expected_df.index), "expected", "expected", "expected", "expected", "expected", "expected"]
metrics_df=metrics_df.transpose()

### Add column to indicate which pipeline is the expected dummy
col=["pipeline"] * (len(metrics_df)-1)
col.append("expected")
metrics_df["exp"]=col

### Predict the coordinates in the PCA for all pipelines plus the expected (that way, the expected doesn't affekt the PCA itself)
pca_df_with_exp=metrics_df.iloc[:, 0:8]
pca_df_with_exp_std=StandardScaler().fit_transform(pca_df_with_exp)
pcs_with_exp=pca.transform(pca_df_with_exp_std)
pc_df = pd.DataFrame(data = pcs_with_exp, columns = ['PC1', 'PC2'], index=pca_df_with_exp.index)



# 7 Plot
## Make df from PC1+2 coordinates and tools per pipeline to plot
plot_df=pd.concat([pc_df, metrics_df.iloc[:, 8:15]], axis = 1).rename_axis("pipeline").reset_index()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="exp", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="trimming_score", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="rRNA_sorting_tool", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="assembly_tool", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="mapper", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="database", hover_data=["pipeline"])
fig.show()

fig = px.scatter(plot_df, x="PC1", y="PC2", color="classifier", hover_data=["pipeline"])
fig.show()






# OLD CODE:

# # This is for all 3 data types, which is inappropriate:
# for metrics_dic in [abs_metrics, rel_metrics, pa_metrics]:
#     for pipeline in shared_pipelines:                                                                                                                       ### TO BE DELETED
#         if metrics_dic==abs_metrics:
#             master_dfs_uniq=master_dfs_uniq_abs
#         elif metrics_dic==rel_metrics:
#             master_dfs_uniq=master_dfs_uniq_rel
#         elif metrics_dic==pa_metrics:
#             master_dfs_uniq=master_dfs_uniq_pa
#         col = master_dfs_uniq[samples[0]].columns.tolist()[1:].index(pipeline)
#         ### Each pipeline gets a variance sum variable
#         var_sum = 0
#         ### For each row in the pipelines across replicates, calculate the variance and sum them up:
#         for row in range(0,master_dfs_uniq[samples[0]].shape[0]):
#             var_sum += statistics.variance([int(master_dfs_uniq[samples[0]].iloc[row,col]), int(master_dfs_uniq[samples[1]].iloc[row,col]), int(master_dfs_uniq[samples[2]].iloc[row,col])])
#         ### Add var_sum to chi2_var dic for respective pipeline
#         metrics_dic[pipeline]["variance"]=var_sum


## 3.3 Assign which tools have been used each step in the pipeline, with one column per step,
##     as well as column for coordinates (Chi-Square statistics and summed variance) normalized to range 0-1,
##     and add a column to find the distance between (1|1) (reverted normalized origin = best case) and each normalized pipeline coordinate

### 3.3.1 Calculate min and max for Chi-Square statistics and summed variance;
# chi2_min = chi2_var[min(chi2_var.keys(), key=(lambda k: chi2_var[k][0]))][0]
# chi2_max = chi2_var[max(chi2_var.keys(), key=(lambda k: chi2_var[k][0]))][0]
# var_min = chi2_var[min(chi2_var.keys(), key=(lambda k: chi2_var[k][1]))][1]
# var_max = chi2_var[max(chi2_var.keys(), key=(lambda k: chi2_var[k][1]))][1]
#
# ### 3.3.2 Add columns for tools, normalized coordinates, and distance
# for pipeline in chi2_var:
#     #### Add normalized Chi-Square statistics and summed variance coordinates (reversed for accuracy so that 1 is best and 0 is worst)
#     chi2_var[pipeline].extend([1-normalize(chi2_var[pipeline][0], chi2_min, chi2_max), 1-normalize(chi2_var[pipeline][1], var_min, var_max)])
#     #### Add distance from normalized coordinates to normalized origin
#     chi2_var[pipeline].append(point_dist([chi2_var[pipeline][2], chi2_var[pipeline][3]], [1,1]))
#     #### Add tool names
#     pipeline_replace = pipeline.replace("IDBA_", "IDBA-").replace("NCBI_NT", "NCBI-NT").replace("BLAST_FIRST_HIT", "BLAST-FIRST-HIT").replace("BLAST_FILTERED", "BLAST-FILTERED")
#     chi2_var[pipeline].extend(pipeline_replace.split("_"))





## 3.4 Save "chi2_var" dic as csv so that it can be plotted using R:
# df_save = pd.DataFrame.from_dict(chi2_var, orient="index", columns=["chi-square statistics", "summed variance", "normalized chi-square statistics", "normalized summed variance", "distance from origin", "trimmed PHRED", "rRNA filter", "assembler", "mapper", "DB", "classification"])
# df_save.to_csv("{savedir}chi2_var.csv".format(savedir=savedir), index_label="pipeline")
#
#
# # 4. Get abundance for mock community members from "best" pipeline (NOTE: only for talk since we change the approach later)
# dic_abun_obs = {"abun": [abun / 3 for abun in master_df_summarized["10_BARRNAP_IDBA_TRAN_BOWTIE2_NCBI_NT_KRAKEN2"]], "taxa": unique_taxa}
# df_abun_obs = pd.DataFrame.from_dict(dic_abun_obs)
# df_abun_obs.to_csv("{savedir}df_abun_obs.csv".format(savedir=savedir))
#
# dic_abun_exp = {"abun": [abun / 3 for abun in master_df_summarized["expected"]], "taxa": unique_taxa}
# df_abun_exp = pd.DataFrame.from_dict(dic_abun_exp)
# df_abun_exp.to_csv("{savedir}df_abun_exp.csv".format(savedir=savedir))


# 4 Calculate variance for absolute and relative data and mismatches between replicates for p/a data

## 4.1  We calculate the variances for each taxon for all pipelines across the three replicates, sum up the variance of all taxa across the
##      replicates for every pipeline, and add the returned summed variance to each pipeline in the respective metrics dics.
##      NOTE: We only do that for absolute and relative data as it is unapplicable for p/a data.

# for metrics_dic in [abs_metrics, rel_metrics]:
#     for pipeline in shared_pipelines[1:len(shared_pipelines)]:                                                                                                                            ### TO BE DELETED
#         if metrics_dic==abs_metrics:
#             master_dfs_uniq=master_dfs_uniq_abs
#         elif metrics_dic==rel_metrics:
#             master_dfs_uniq=master_dfs_uniq_rel
#         ### Each pipeline gets a variance sum variable
#         var_sum = 0
#         pipeline_col = master_dfs_uniq[samples[0]].columns.tolist()[1:].index(pipeline)
#         ### For each taxon in the pipelines across replicates, calculate the variance and sum up all variances across all taxa:
#         for taxon in range(0,master_dfs_uniq[samples[0]].shape[0]):
#             var_sum += statistics.variance([int(master_dfs_uniq[samples[0]].iloc[taxon,pipeline_col]), int(master_dfs_uniq[samples[1]].iloc[taxon,pipeline_col]), int(master_dfs_uniq[samples[2]].iloc[taxon,pipeline_col])])
#         ### Add var_sum to chi2_var dic for respective pipeline
#         metrics_dic[pipeline]["variance"]=var_sum


## 4.2 We calculate the mismatches between replicates for each taxon for all pipelines, sum up the mismatches
## for every pipeline, and add the mismatches to each pipeline in the respective metrics dics.
# for pipeline in shared_pipelines[1:len(shared_pipelines)]:                                                                                                                            ### TO BE DELETED
#     ### Each pipeline gets a mismatches variable to count mismatches between replicates
#     mismatches = 0
#     pipeline_col = master_dfs_uniq_pa[samples[0]].columns.tolist()[1:].index(pipeline)
#     ### For each taxon in the pipelines across replicates, count how often 1 (=present) was detected
#     ### - if it was detected only 1 one 2 times, pipelines mismatch in their outcome, which we collect in variable mismatches:
#     for taxon in range(0,master_dfs_uniq_pa[samples[0]].shape[0]):
#         count=[master_dfs_uniq_pa[samples[0]].iloc[taxon,pipeline_col], master_dfs_uniq_pa[samples[1]].iloc[taxon,pipeline_col], master_dfs_uniq_pa[samples[2]].iloc[taxon,pipeline_col]].count(1)
#         if count==1 or count==2:
#             mismatches+=1
#     pa_metrics[pipeline]["mismatches"]=mismatches



# ## 3.1 Chi-Squared test (accuracy)
#
# ### 3.1.1 To perform a Chi-Squared test on replicates, we summarize columns across the replicates:
# master_df_summarized = pd.DataFrame()
# ### For all columns (note: all samples contain the same column names, so we manually pick one (the first in "samples" list))
#
# ### NOTE: For now, some pipelines didn't work, so we have to make a list of shared pipelines                                                            ### TO BE DELETED
# shared_pipelines = []                                                                                                                                   ### TO BE DELETED
# for sample in samples:                                                                                                                                  ### TO BE DELETED
#     shared_pipelines.extend(master_dfs_uniq_abs[sample].columns.tolist()[1:])                                                                               ### TO BE DELETED
# shared_pipelines = list(set(i for i in shared_pipelines if shared_pipelines.count(i) > 2))                                                              ### TO BE DELETED
#
# #for col in master_dfs_uniq_abs[samples[0]].columns.tolist():
# shared_pipelines_exp = shared_pipelines[:]                                                                                                              ### TO BE DELETED
# shared_pipelines_exp.append("expected")                                                                                                                 ### TO BE DELETED
# for pipeline in shared_pipelines_exp:                                                                                                                   ### TO BE DELETED
#     #### Make list per column in every sample and save in dir:
#     pipeline_dic_chi2 = {}
#     for sample in samples:
#         pipeline_dic_chi2[sample] = master_dfs_uniq_abs[sample][pipeline].tolist()
#     #### Summarize the three columns in the three replicates:
#     master_df_summarized[pipeline] = [a + b + c for a, b, c in zip(pipeline_dic_chi2[samples[0]], pipeline_dic_chi2[samples[1]], pipeline_dic_chi2[samples[2]])]
#
# ### 3.1.2 Now we perform a Chi-Squared test on the summarized master df, for each pipeline against the expected composition
# ###       (note, this is the DIRTY version where we add the same amount (1) to each value to get rid of zeros in the expected columns, because otherwise the Chi-Squared test is not possible):
#
# #for pipeline in master_df_summarized.columns.tolist()[1:]:
# for pipeline in shared_pipelines:                                                                                                                       ### TO BE DELETED
#     chi2_var[pipeline] = [(chisquare([abun + 1 for abun in master_df_summarized[pipeline].tolist()],[abun + 1 for abun in master_df_summarized['expected'].tolist()])[0])]
#
#
# ## 3.2 Sum of variance for each taxon (precision)
# ##     Finally, we calculate the variance for all pipelines across the three replicates by summing up the variance of each taxa across the
# ##     replicates for every pipeline and add the returned summed variance to the pipelines in the "chi2_var" dic
#
# #for pipeline in master_dfs_uniq_abs[samples[0]].columns.tolist()[1:]:
# for pipeline in shared_pipelines:                                                                                                                       ### TO BE DELETED
#     ### Turn pipelines into indexes
#     col = master_dfs_uniq_abs[samples[0]].columns.tolist()[1:].index(pipeline)
#     ### Each pipeline gets a variance sum variable
#     var_sum = 0
#     ### For each row in the pipelines across replicates, calculate the variance and sum them up:
#     for row in range(0,master_dfs_uniq_abs[samples[0]].shape[0]):
#         var_sum += statistics.variance([int(master_dfs_uniq_abs[samples[0]].iloc[row,col]), int(master_dfs_uniq_abs[samples[1]].iloc[row,col]), int(master_dfs_uniq_abs[samples[2]].iloc[row,col])])
#     ### Add var_sum to chi2_var dic for respective pipeline
#     chi2_var[pipeline].append(var_sum)
#
#
