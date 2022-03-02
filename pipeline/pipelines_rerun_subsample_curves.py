#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 2 Mar 2022

# This script re-runs specific pipelines for specific samples in a specific folder structure

import pandas as pd
import os

# Set parameters
lvls=["genus_rel", "genus_pa", "species_rel", "species_pa"]
samples=['M4_DNA', 'M5_DNA', 'M6_DNA', 'M4_RNA', 'M5_RNA', 'M6_RNA']
#dbs=["ncbi", "silva"]
dbs=["silva"]

# Import excel file containing info on which pipeline to use for which sample
# pipelines=pd.read_excel("/Users/christopherhempel/Google Drive/PhD UoG/Pilot project/rerun_dna_subsamples.xlsx")
pipelines=pd.read_excel("/hdd1/chempel/pilot_project/subsamples_curves/rerun_dna_subsamples.xlsx")

# Make dic out of excel file
dic={}
for db in dbs:
    dic[db]={}
    df_db=pipelines[pipelines["db"]==db]
    for sample in samples:
        dic[db][sample]={}
        df_db_sample=df_db[df_db["sample"]==sample]
        for lvl in lvls:
            dic[db][sample][lvl]=list(df_db_sample[lvl])[0]

# Run pipelines through the shell, including changing directories and defining filenames 
workdir=os.getcwd()
for sample in samples:
    for sub in range(1,11):
        os.chdir("{0}/{1}".format(sample, sample + "_subsample_" + str(sub)))
        subdir=os.getcwd()
        R1="{0}/{1}_R1_sub{2}.fastq".format(subdir, sample, str(sub))
        R2="{0}/{1}_R2_sub{2}.fastq".format(subdir, sample, str(sub))
        for db in dbs:
            if not os.path.exists(db):
                os.mkdir(db)
            os.chdir(db)
            for lvl in lvls:
                if not os.path.exists(lvl):
                    os.mkdir(lvl)
                os.chdir(lvl)
                pipeline=dic[db][sample][lvl]
                os.system("METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_steinke_server_single_pip.sh -1 {0} -2 {1} -p 32 -P {2}".format(R1, R2, pipeline))
                os.chdir("..")
            os.chdir("..")
        os.chdir(workdir)
