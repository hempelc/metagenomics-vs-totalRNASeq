#library(data.table)
library(dplyr)

# Set path
path <- getwd()

# Hardcoded parameters
file1 <- "trimmed_at_phred_5_unsorted_metaspades_bwa_ncbi_nt_kraken2_final.txt" 
num_reads <- 269408 # Number of reads in sample M6_RNA
rel_abun_basis <- "cell" # Basic for relative abundance of mock community taxa. Either based on genomic DNA (option: "gen") or on cell number (option "cell")

## To calculate the absolute expected read count for the taxa in the mock community, we need the hardcoded taxonomic information of each taxon:
L_monocytogenes_tax <- c("Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Listeriaceae", "Listeria", "Listeria monocytogenes")
P_aeruginosa_tax <- c("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Pseudomonadales", "Pseudomonadaceae", "Pseudomonas", "Pseudomonas aeruginosa")
B_subtilis_tax <- c("Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Bacilliaceae", "Bacillus", "Bacillus subtilis")
S_cerevisiae_tax <- c("Eukaryota", "Ascomycota", "Saccharomycetes" ,"Saccharomycetales", "Saccharomycetaceae", "Saccharomyces", "Saccharomyces cerevisiae")
E_coli_tax <- c("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia", "Escherichia coli")
S_enterica_tax <- c("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Salmonella", "Salmonella enterica")
L_fermentum_tax <- c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Limosilactobacillus", "Lactobacillus fermentum")
E_faecalis_tax <- c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Enterococcaceae", "Enterococcus", "Enterococcus faecalis")
C_neoformans_tax <- c("Eukaryota", "Basidiomycota", "Tremellomycetes", "Tremellales", "Tremellaceae", "Cryptococcus", "Cryptococcus neoformans")
S_aureus_tax <- c("Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus", "Staphylococcus aureus")

## We also need absolute abundances per taxon. Therefore, we multiply the # of reads in sample with the relative abundance of each taxon in the mock community
if (rel_abun_basis=="gen") {
  rel_abun <- num_reads*c(89.1, 8.9, 0.89, 0.89, 0.089, 0.089, 0.0089, 0.00089, 0.00089, 0.000089)/100 # Order of rel abundances in vector equals order of taxa above
} else if (rel_abun_basis=="cell") {
  rel_abun <- num_reads*c(94.9, 4.2, 0.7, 0.12, 0.058, 0.059, 0.015, 0.001, 0.00007, 0.0001)/100 # Order of rel abundances in vector equals order of taxa above
}

## Finally, we generate vectors for each taxon with their taxonomy plus the absolute read count
### Note: ugly code but I don't know how to code this pretty in R, spent 1 hour on it trying to make a for loop and gave up lol
L_monocytogenes <- append(L_monocytogenes_tax, rel_abun[1])
P_aeruginosa <- append(P_aeruginosa_tax, rel_abun[2])
B_subtilis <- append(B_subtilis_tax, rel_abun[3])
S_cerevisiae <- append(S_cerevisiae_tax, rel_abun[4])
E_coli <- append(E_coli_tax, rel_abun[5])
S_enterica <- append(S_enterica_tax, rel_abun[6])
L_fermentum <- append(L_fermentum_tax, rel_abun[7])
E_faecalis <- append(E_faecalis_tax, rel_abun[8])
C_neoformans <- append(C_neoformans_tax, rel_abun[9])
S_aureus <- append(S_aureus_tax, rel_abun[10])

# Make a dataframe with expected community composition
expected <- as.data.frame(t(data.frame(L_monocytogenes, P_aeruginosa, B_subtilis, S_cerevisiae, E_coli, S_enterica, L_fermentum, E_faecalis, C_neoformans, S_aureus)))
names(expected) <- c("superkingdom", "phylum", "class", "order", "family", "genus", "lowest_hit", "counts")

# Read in the data 
Inputs <- read.table(paste(path,file1, sep="/"),header = T, sep ="\t")
## We only include the monophyletic ranks superkingdom, phylum, class, order, family, genus, and species, together with their counts
Inputs <- select(Inputs, "superkingdom", "phylum", "class", "order", "family", "genus", "lowest_hit", "counts")
Data <- list(Inputs, expected)
names(Data[1]) <- file1 # Doesn't work
names(Data)

# Set up data frame 
# Eventually may wish to switch to using datatable to speed up indexing, but will rewrite code then to deal with it

# Because expected will always have to be calculated will use that as opening column and then append columns with files stored in data
Master <- data.frame("Expected" = rep(NA, times = 12), 
                     row.names = c("superkingdom", "phylum", "class", "order", "family", "genus", "lowest_hit", "counts"))
head(Master)

Master <- cbind(Master, matrix(NA, nrow = nrow(Master), ncol= length(Data)))
names(Master)[2:ncol(Master)] <- names(Data)

# Accuracy
# Use chi square function
chisq.test()

# Precision
var()

# Heuristic Selection
plot()

# Compare
# F test
var.test()
