#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 2 Mar 2022

# This script tests if expected genera and species are found among aquarium samples


import pandas as pd

workdir='/Users/christopherhempel/Desktop/pipeline_results/pipeline_results_env_samples/'

samples=["F4_DNA", "F4_RNA", "F5_DNA", "F5_RNA", "F6_DNA", "F6_RNA"]
genera=["Caridina", "Ostracoda", "Aplocheilus", "Cynodonichthys", "Danio", "Poecillia", "Anubias",
    "Bolbitis", "Cryptocoryne", "Echinodoras", "Limnophila", "Rotala", "Sagittaria",
    "Vesicularia", "Physella", "Planorbarius", "Melanoides", "Neritidae"]
species=["Caridina multidentata", "Aplocheilus lineatus", "Cynodonichthys hildebrandi",
    "Danio rerio", "Poecillia salvatoris", "Anubias barteri", "Bolbitis heudelotii",
    "Limnophila sessiliflora", "Rotala rotundifolia", "Sagittaria subulata",
    "Physella acuta", "Planorbarius corneus", "Melanoides tuberculate"]

df_genus=pd.DataFrame()
df_spe=pd.DataFrame()

gen_dic={}
for sample in samples:
    df=pd.read_csv(workdir + "rel_abun_genus/" + sample + '_rel_abun_genus_pa.csv')
    gen_dic[sample]=df

spe_dic={}
for sample in samples:
    df=pd.read_csv(workdir + "rel_abun_species/" + sample + '_rel_abun_species_pa.csv')
    spe_dic[sample]=df

for genus in genera:
    lst=[]
    for df in gen_dic.values():
        if genus in df.columns:
            if df[genus].sum()!=0:
                lst.append("yes")
            else:
                lst.append("no")
        else:
            lst.append("no")
    df_genus[genus]=lst
df_genus.index=samples

for spe in species:
    lst=[]
    for df in spe_dic.values():
        if spe in df.columns:
            if df[spe].sum()!=0:
                lst.append("yes")
            else:
                lst.append("no")
        else:
            lst.append("no")
    df_spe[spe]=lst
df_spe.index=samples
