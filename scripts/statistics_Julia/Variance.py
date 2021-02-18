#!/usr/bin/env python3
import pandas as pd
import numpy as np
import statistics

d = {'taxa':["SpeciesA", "SpeciesB", "SpeciesC", "SpeciesD", "SpeciesE", "SpeciesF"], 
        'Pipeline1':[100,34,5,0,0,77], 
        'Pipeline2':[40,30,0,0,0,22],
        'Pipeline3':[400,350,200,70,11,77]}
Rep1 = pd.DataFrame(data = d)
d = {'taxa':["SpeciesA", "SpeciesB", "SpeciesC", "SpeciesD", "SpeciesE", "SpeciesF"], 
        'Pipeline1':[110,40,6,0,1,75], 
        'Pipeline2':[42,36,20,10,0,0],
        'Pipeline3':[430,355,220,73,11,77]}
Rep2 = pd.DataFrame(data = d)
d = {'taxa':["SpeciesA", "SpeciesB", "SpeciesC", "SpeciesD", "SpeciesE", "SpeciesF"], 
        'Pipeline1':[113,36,0,0,0,79], 
        'Pipeline2':[30,43,10,10,0,23],
        'Pipeline3':[420,365,130,30,11,47]}
Rep3 = pd.DataFrame(data = d)

#This is assuming the dataframes are identical in shape and taxa/pipeline order in each replicate, I believe this to be true
#Going to iterate one column (one piepline) at a time through the dataframe, come up with a variance score for tht column, than save it in an output object
VarianceScore = np.repeat(0,(Rep1.shape[1]-1))
#Start with 1, because we want second column, first column always has taxa names
for col in range(1,Rep1.shape[1]):
    #Each column needs its own dynamic score vector
    score = 0
    #For each row calculate variance and sum the up
    for row in range(0,Rep1.shape[0]):
        score = score + statistics.variance([Rep1.iloc[row,col], Rep2.iloc[row,col], Rep3.iloc[row,col]])
    #Once done for every taxa save the finale variance in the output object
    #Note have to index col - 1 in vectore because out count data starts in the second row of our data frame 
    VarianceScore[(col-1)] = score

print(VarianceScore)
        

#Same as above but uses population variance instead of sample variance 
VarianceScore = np.repeat(0,(Rep1.shape[1]-1))
for col in range(1,Rep1.shape[1]):
    score = 0
    for row in range(0,Rep1.shape[0]):
        score = score + statistics.pvariance([Rep1.iloc[row,col], Rep2.iloc[row,col], Rep3.iloc[row,col]])
    VarianceScore[(col-1)] = score

print(VarianceScore)
        