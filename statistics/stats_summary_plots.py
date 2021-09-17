import pandas as pd
import plotly.express as px
import logging
import os

levels=["genus_cell", "genus_gen", "species_cell", "species_gen"]
metr_lst=["all", "rel", "pa"]
master_df=pd.DataFrame()
for level in levels:
    #level=levels[0]
    statsdir="/Users/christopherhempel/Desktop/pipeline_results_mock_community/results_{0}/stats_exports".format(level)
    for metr in metr_lst:
        #metr=metr_lst[0]
        df=pd.read_csv(os.path.join(statsdir, "{0}_mean_euc_dist_steps.csv".format(metr)))
        df["p-value"]=pd.read_csv(os.path.join(statsdir, "{0}_pvalues_tools.csv".format(metr)))["p-value"]
        df["type"]=[level + "_" + metr]*len(df)
        master_df=pd.concat([master_df, df])


import plotly.express as px
df = px.data.gapminder()


fig = px.scatter(master_df, x="type", y="tools",
	         size="mean_euc_dist", color="p-value",
                 hover_name="mean_euc_dist", size_max=10, height=800)
fig.show()



heat_dic={}
for step in df["step"].unique():
    heat_dic[step]=df[df['step'].str.contains(step)]["p-value"].to_list()
heat_df=pd.DataFrame(heat_dic, index=(df["metrics"] + "_" + df["level+cell"]).unique())


##### !!!!!!!! manual part for now
df = pd.read_csv("/Users/christopherhempel/Desktop/test.csv")
import plotly.graph_objects as go

fig = go.Figure(data=go.Scatter3d(
    x=df['step'],
    y=df['level+cell'],
    z=df['metrics'],
    text=df['p-value'],
    mode='markers',
    marker=dict(
        sizemode='diameter',
        sizeref=0.1,
        size=df['p-value'],
        color=df['p-value'],
        colorscale = 'Viridis',
        colorbar_title = 'p-value',
        line_color='rgb(140, 140, 170)'
    )
))
fig.show()

# Radar plot
df_radar=copy.deepcopy(df)
df_radar["comb"] = df_radar["step"] + "_" + df_radar["metrics"] + "_" + df_radar["level+cell"]
df_radar.sort_values(by=['comb'], inplace=True)
fig = px.line_polar(df_radar, r='p-value', theta='comb', line_close=True)
fig.show()

# Heatmap
heat_dic={}
for step in df["step"].unique():
    heat_dic[step]=df[df['step'].str.contains(step)]["p-value"].to_list()
heat_df=pd.DataFrame(heat_dic, index=(df["metrics"] + "_" + df["level+cell"]).unique())
fig_heat = px.imshow(heat_df.sort_index().transpose())
fig_heat.show()
