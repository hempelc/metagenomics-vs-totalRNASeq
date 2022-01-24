#!/usr/bin/env python3

# Written by Chris Hempel (hempelc@uoguelph.ca) on 20 Jan 2022
# A script to filter blast hits with assigned taxonomy

import datetime
import pandas as pd
import sys
import argparse
import warnings

# Set that warnings are not printed to console
warnings.filterwarnings("ignore")

# Define funtion to print datetime and text
def time_print(text):
    datetime_now=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(datetime_now + "  ---  " + text)

# Define function to filter bitscores
def bitscore_cutoff(x):
    min_bitscore=x.max()-x.max()*0.02
    return x[x>=min_bitscore]

# Define class to format helptext of options properly, taken from
# https://www.google.com/search?q=argsparse+recognize+newline+in+help&oq=argsparse+
# recognize+newline+in+help&aqs=chrome..69i57j33i22i29i30.12450j0j7&sourceid=chrome&ie=UTF-8
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

# Define ranks to use
ranks=["superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass",
    "order", "suborder", "infraorder", "family", "genus", "lowest_hit"]
# Define which columns to load in
req_cols=['qseqid', 'pident', 'length', 'bitscore']+ranks


time_print("Reading in file...")
# Only read in columsn we need
df=pd.read_csv(file, usecols=req_cols, sep="\t", dtype={"qseqid": str}).fillna('NA')


time_print('Checking if species are in "Genus species" format...')
## Cut down to first two words
df["lowest_hit"]=df["lowest_hit"].str.split(' ').str[:2].str.join(sep=' ')

## Check if the first word is capitalized. If not then turn species into "NA"
idx_spe=df["lowest_hit"].str[0].str.isupper()
df["lowest_hit"]=[spe if capitalized else "NA" for spe, capitalized in zip(df["lowest_hit"], idx_spe)]


if filter=="strict":
    time_print("Filtering hits based on bitscore and length...")
    df.loc[(df['length'] < 100) & (df['bitscore'] < 155), ranks] = "NA"

    time_print("Grouping qseqids and filtering hits based on bitscore cutoff for each qseqid...")
    idx = df.groupby(['qseqid'])['bitscore'].transform(bitscore_cutoff) == df['bitscore']
    df=df[idx]

    time_print("Applying similarity cutoff...")
    df.loc[df["pident"] < cutoff[0], 'lowest_hit'] = "NA"
    df.loc[df["pident"] < cutoff[1], 'genus'] = "NA"
    df.loc[df["pident"] < cutoff[2], 'family'] = "NA"
    df.loc[df["pident"] < cutoff[3], ['order', 'suborder', 'infraorder']] = "NA"
    df.loc[df["pident"] < cutoff[4], ['class', 'subclass']] = "NA"
    df.loc[df["pident"] < cutoff[5], ['phylum', 'subphylum']] = "NA"

elif filter=="soft":
    time_print("Grouping qseqids and filtering hits based on the highest bitscore of each qseqid...")
    idx=df.groupby(['qseqid'])['bitscore'].transform(max) == df['bitscore']
    df=df[idx]

# Keep only relevant columns and put lowest_hit to last column
df=df[["qseqid"] + ranks]


time_print("Performing LCA filter...")
## Make a df mask: group df by contigs, check if rank has more than one taxon, and if yes, True, else False
lca_mask=df.groupby(["qseqid"]).transform(lambda x: len(set(x)) != 1)

## Replace ranks in df with "NA" based on mask
df=df.mask(lca_mask, "NA")
df["qseqid"]=df['qseqid']

## Drop duplicate rows == aggregate contig info
df.drop_duplicates(inplace=True)


# Change column name and save df
df.rename(columns={'lowest_hit': 'species', 'qseqid': 'sequence_name'}, inplace=True)
df.to_csv(out, sep="\t", index=False)

time_print("Filtering done.")
