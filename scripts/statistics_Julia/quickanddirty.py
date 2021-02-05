#!/usr/bin/env python

import pandas as pd
import numpy as np
import glob

#Read in expected
expected = pd.read_csv("ExpectedCommunityAbundance.csv", header = 0)

#Read in file
txt_files = glob.glob("mock_community_RNA/M4_RNA/*.txt")

taxa = pd.DataFrame(data = [])

for i in txt_files:
    df = pd.read_table(i, header = 0)
    taxa = pd.concat([taxa, df['lowest_hit']])

print(taxa.shape)
unique = taxa.drop_duplicates()

test = txt_files[0]
df = pd.read_table(test, header = 0)

unique_list = unique.values.tolist()
uniques = [i for j in unique_list for i in j]
print(uniques)




#print((unique_list[0]))
#unique_str = str(unique_list)
#print((unique_str[0]))

# species = []
# for i in unique_list:
#     if len(i.split()) == 2:
#         species + i
        

# print(species)

#rint(unique)
#for i in range(0,10):
#    print (unique.iloc[i])
#    print (type(unique.iloc[i]))
    #(i)

yes = df['lowest_hit'].isin(unique.iloc[0])
#print(yes)




# unique_genus = file1['genus'].drop_duplicates()
# print(unique_genus.iloc[0])
# print(unique_genus.iloc[1])

# for i in unique_genus:
#     file1


# #Coulmns will be expected

# master = pd.DataFrame(data = 0, index = , columns =)

#Populate expected column

#Populate other colunms

#Perform chisquare and save

#Calculate variance

#Plot two values
