#!/usr/bin/env python3

import pandas as pd
import sys
import argparse

# Make class to format helptext of options properly, taken from
# https://www.google.com/search?q=argsparse+recognize+newline+in+help&oq=argsparse+recognize+newline+in+help&aqs=chrome..69i57j33i22i29i30.12450j0j7&sourceid=chrome&ie=UTF-8
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)

# Define arguments
parser = argparse.ArgumentParser(description='Filter BLAST output.',
    formatter_class=SmartFormatter)
parser.add_argument('inputfile',
    help='Input file in BLAST standard output format plus staxids.')
parser.add_argument('filter', choices=['soft', 'strict'],
    help=('R|Type of filtering.\n\nsoft:\nKeeps the best hit (highest bitscore) '
    'for each sequence. If multiple hits have the same highest bitscore, an LCA '
    'approach is applied (assigns the taxonomy to each sequence based on all '
    'taxonomic ranks that are identical in the remaining hits of each sequence).'
    '\n\nstrict:\nPerforms 3 steps:\n  (1) bitscore filtering - keeps all hits with '
    'a bitscore >= argument -b and within argument -p %% of the best bitscore '
    'per sequence.\n  (2) similarity cutoff - only keeps the taxonomy of hits up '
    'to a certain rank, depending on the hits blast %% identity and cutoff '
    'values given in argument -c.\n  (3) LCA approach - assigns the taxonomy to '
    'each sequence based on all taxonomic ranks that are identical in the '
    'remaining hits of each sequence.'))
parser.add_argument('-p', '--percentage', default=2, type=int,
    help=('Percentage threshold to perform bitscore filtering on when choosing '
    'filter option "strict" (default=2).'))
parser.add_argument('-l', '--length', default=100, type=int,
    help=('Alignment length threshold to perform bitscore filtering on when '
    'choosing filter option "strict" (default=100).'))
parser.add_argument('-b', '--bitscore', default=155, type=int,
    help=('Bitscore threshold to perform bitscore filtering on when choosing '
    'filter option "strict" (default=155).'))
parser.add_argument('-c', '--cutoff', metavar="N", nargs=6,
    default=[99, 97, 95, 90, 85, 80], type=int,
    help=('Similarity cutoff per hit based on BLAST pident values when choosing '
    'filter option "strict". pident cutoffs have to be specified via integers '
    'divided by spaces, in the order for ranks species, genus, family, order, '
    'class, phylum. Taxonomy is only kept for a rank of a hit if the hits pident '
    'is >= the respective cutoff (default=99 97 95 90 85 80).'))
parser.add_argument('-o', '--out', default="blast_filtered.txt",
    type=str, help='Name of output file.')
args = parser.parse_args()

# Set arguments
file=args.inputfile
filter=args.filter
percentage=args.percentage/100
length=args.length
bitscore_threshold=args.bitscore
cutoff=args.cutoff
out=args.out

# Define ranks to loop over later
ranks=["superkingdom", "kingdom", "phylum", "subphylum", "class",
    "subclass", "order", "suborder", "infraorder", "family", "genus", "lowest_hit"]
# Read in df
df=pd.read_csv(file, sep="\t", dtype={"qseqid": str}).fillna('NA')
# Make dicts with taxonomy keys and empty lists
tax_dic={"sequence_name":[]}
for x in ranks:
    tax_dic[x]=[]

# For each unique ID
sequence_names=sorted(df["qseqid"].unique())
for sequence_name in sequence_names:
    # Add the id name to the dic
    tax_dic["sequence_name"].append(sequence_name)
    # Subsample the big df to only rows that contain the sequence id
    df_sub=df[df["qseqid"]==sequence_name]
    # Extract the best bitscore
    best_bitscore=sorted(df_sub["bitscore"], reverse=True)[0]
    print("Processing sequence {0} with bitscore {1}".format(sequence_name,
        best_bitscore))

    # Apply bitscore filter depending on used filter
    if filter=="soft":
        # Bitscore filter = best bitscore
        df_sub_filtered=df_sub[df_sub["bitscore"]==best_bitscore]
    elif filter=="strict":
        # 1. Bitscore filter = minimum alignment length of "length" and minimum
        #    bitscore of "bitscore_threshold" and only hits within bitscore range
        #    of the top "percentage" of the best bitscore
        df_sub_filtered=df_sub[(df_sub["bitscore"]>=best_bitscore-best_bitscore*percentage)
            & (df_sub["length"]>=length) & (df_sub["bitscore"]>=bitscore_threshold)]

        # 2. Apply similarity cutoff
        for idx, row in df_sub_filtered.iterrows():
            if row["pident"] < cutoff[0]:
                    df_sub_filtered.at[idx,"lowest_hit"]="NA"
            if row["pident"] < cutoff[1]:
                    df_sub_filtered.at[idx,"genus"]="NA"
            if row["pident"] < cutoff[2]:
                    df_sub_filtered.at[idx,"family"]="NA"
            if row["pident"] < cutoff[3]:
                for rank in ['order', 'suborder', 'infraorder']:
                    df_sub_filtered.at[idx,rank]="NA"
            if row["pident"] < cutoff[4]:
                for rank in ['class', 'subclass']:
                    df_sub_filtered.at[idx,rank]="NA"
            if row["pident"] < cutoff[5]:
                for rank in ['phylum', 'subphylum']:
                    df_sub_filtered.at[idx,rank]="NA"

    # Apply LCA filter
    for rank in ranks:
        if rank=="lowest_hit":
            taxon=df_sub_filtered[rank].str.split(' ').str[:2].str.join(sep=' ')\
                .unique() # Get first two words and each unique name
        else:
             # First two words not needed, should always be one word:
            taxon=df_sub_filtered[rank].unique()
        if len(taxon)!=1: # If not only one unique name
            tax_dic[rank].append("NA")
        elif rank=="lowest_hit": # If one uniqie name but species rank, make more checks
            # If species is two words and if capitalization structure "Genus species":
            if len(taxon[0].split(' '))!=1 and taxon[0][0].isupper() \
                and taxon[0].split(' ')[1][0].islower():
                tax_dic[rank].append(taxon[0])
            else:
                tax_dic[rank].append("NA")
        else: # If oen unique name and not species rank
            tax_dic[rank].append(taxon[0])

# Turn dic to df, cahnge on column name, and save df
tax_df=pd.DataFrame(tax_dic)
tax_df.rename(columns={'lowest_hit': 'species'}, inplace=True)
tax_df.to_csv(out, sep="\t", index=False)
