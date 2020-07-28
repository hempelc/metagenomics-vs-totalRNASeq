library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
library(vegan)

##### Set parameters

# Decide what taxonomic groups you want to use for taxonomic assignments, data assigned to superkingdom (superkingdom == T) or phylum (superkingdom == F) 
superkingdom <- T
if (superkingdom==TRUE) {taxlist_group <- "superkingdom"} else {taxlist_group <- "phylum"}

# Indicate assemblers, classifiers, and mappers you used
assemblers <- c("IDBA-tran", "rnaSPAdes", "Trinity")
classifiers <- c("BLAST_nt", "BLAST_SILVA", "CREST")
mappers <- c("Bowtie2", "BWA")

# Set workfolder
workfolder <- "/Users/christopherhempel/Desktop/chrisnatjulia/FINAL_RNA"


##### Start analysis

# Setting WD
setwd(workfolder)

# Reading in data from all files in folder
filenames_list <- list.files(pattern = "\\.txt$")
all_list <- lapply(filenames_list, fread) # Creating a list with data frames of all files
names(all_list) <- gsub("\\_final.txt$", "", filenames_list) # Changing the names of the dataframes in list
tax_full_dataframes <- lapply(all_list, function(x) x%>% select(superkingdom:genus)) # Creating a list with full taxanomy of each data frame
counts_dataframes <- lapply(all_list, function(x) x%>% select(counts)) # Creating a list with counts of each data frame

# Create function to create taxlist vector with phylum per row in taxonomy table (or if phylum = NA, then superkingdom: Unkown)
make_tax_phylum_dataframe <- function(tax)
  {
    taxlastrow <- nrow(tax)
    taxlist <- rep(NA,nrow(tax))
    for (i in 1:nrow(tax))
    { if (is.na(tax[i,1]) == TRUE)
      { taxlist[i] <- "NA"
      } else if (tax[i,1] == "Unknown")
      { taxlist[i] <- "Unknown"
      } else if (is.na(tax[i,3]) == TRUE)
      { taxlist[i] <- paste(tax[i,1], ": Unknown", sep="")
      } else
      { taxlist[i] <- paste(tax[i,1], ": ", tax[i,3], sep="") }
    }
    taxdf <- as.data.frame(taxlist)
    return(taxdf)
  }

# Create function to create taxlist vector with superkingdom per row in taxonomy table
make_tax_superkingdom_dataframe <- function(tax)
{
  taxlastrow <- nrow(tax)
  taxlist <- rep(NA,nrow(tax))
  for (i in 1:nrow(tax))
  { if (is.na(tax[i,1]) == TRUE)
      { taxlist[i] <- "NA"
      } else
      { taxlist[i] <- tax[i,1] }
  }
  taxdf <- t(as.data.frame(taxlist))
  return(taxdf)
}

# Apply taxlist fuction to list of dataframes (phylum and superkingdom assignments in separate dataframe lists)
tax_superkingdom_dataframes <- lapply(tax_full_dataframes, make_tax_superkingdom_dataframe)
tax_phylum_dataframes <- lapply(tax_full_dataframes, make_tax_phylum_dataframe)

if (taxlist_group == "superkingdom") {
  tax_group_dataframes <- tax_superkingdom_dataframes
} else {tax_group_dataframes <- tax_phylum_dataframes}

# Combine read count and taxonomic group information and aggregate based on groups
group_counts_dataframes <- mapply(cbind, tax_group_dataframes, counts_dataframes, SIMPLIFY=FALSE)
aggregate_tax_groups <- function(group_counts_dataframes)
  {
    setNames(aggregate(unlist(group_counts_dataframes[[2]])~unlist(group_counts_dataframes[[1]]),data=group_counts_dataframes,FUN=sum), c("Group", "Reads"))
  }
agg_dataframes <- lapply(group_counts_dataframes, aggregate_tax_groups)

# Merge all dataframes in list
merged <- agg_dataframes %>%
  Reduce(function(df1,df2) full_join(df1,df2,by="Group"), .)  %>% 
  mutate_each(funs(replace(., which(is.na(.)), 0)))
colnames(merged) <- c("Group", gsub("\\_final.txt$", "", filenames_list))
merged_ordered <- merged[order(as.character(merged$Group)),] # Order groups alphabetical
merged_ordered$Group <- factor(merged_ordered$Group, levels=merged_ordered$Group) # Command to keep order for ggplot2

# Set color for plots
viridis_colors <- viridis_pal(option = "D")(nrow(merged_ordered)) # Choosing viridis colors (colorblindfriendly), as many colors from gradient as numbers of groups
set.seed(002)
colorvec1 <- sample(viridis_colors) # Randomize color order, otherwise bars next to each other are hard to distinguish
colorvec2 <- c("#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#11c638", "#8dd593", "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7", "#f3e1eb", "#f6c4e1", "#f79cd4", "#4a6fe3", "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a")
colorvec3 <- c("#72d5de","#eca9b0","#7fe1cf","#e1b0dd","#aad9a7","#74aff3","#c6d494","#b9b5f3","#ebc491","#7bc2f1","#dac6a3","#8bd0eb","#94dcba","#b6bee4","#acd8ba","#86bcb1","#afe6db")

# Generate relative abundances
relative_abundances_only <- merged_ordered
relative_abundances_only[,1] <- NULL # Delete column Group
relative_abundances_only <- decostand(relative_abundances_only, "total", 2)
groups_relative_abundances <- cbind(merged_ordered$Group, decostand(relative_abundances_only, "total", 2))
colnames(groups_relative_abundances) <- c("Group", gsub("\\_final.txt$", "", filenames_list))

# Transform dataset for stacked barplots and plot
groups_melted <- melt(groups_relative_abundances,id.vars = "Group", variable.name="Combination", value.name="Reads")

plot1<-ggplot(data=groups_melted, aes(x=Combination, y=Reads, fill=Group))+
  geom_bar(stat="identity")+
  theme_minimal()+
  scale_fill_manual(values=colorvec1)+
  theme(legend.key.size = unit(0.1,"line"))+
  theme(legend.position="bottom")+
  theme(legend.text=element_text(size=1))+
  theme(legend.title=element_blank())
plot1
#ggsave("Plot1.png", plot=plot1,  device="png",  width=240, units="mm")

## Separate bars based on groups

# Generate Assembler, Classifier, and Mapper columns for different combinations
number_of_assemblers <- length(assemblers)
number_of_classifiers <- length(classifiers)
number_of_mappers <- length(mappers)

# Generate assembler column
assembler_column <- NULL 
for (i in 1:number_of_assemblers) 
{
  loopresult <- rep(assemblers[i], nrow(groups_melted)/number_of_assemblers)
  assembler_column <- c(assembler_column,loopresult)
}

# Generate classifier column
classifier_column <- NULL
rep_classifier <- 1
while (rep_classifier<number_of_assemblers+1)
{
  for (i in 1:number_of_classifiers) 
  {
    loopresults <- rep(classifiers[i], nrow(groups_melted)/number_of_assemblers/number_of_classifiers)
    classifier_column <- c(classifier_column,loopresults)
  }
rep_classifier <- rep_classifier+1
}

# Generate mapper column
mapper_column <- NULL
rep_mapper <- 1
while (rep_mapper<number_of_assemblers*number_of_classifiers+1)
{
  for (i in 1:number_of_mappers)
  {
    loopresults <- rep(mappers[i], nrow(groups_melted)/number_of_assemblers/number_of_classifiers/number_of_mappers)
    mapper_column <- c(mapper_column,loopresults)
  }
  rep_mapper <- rep_mapper+1
}

# Attach new columns to groups_melted
groups_melted$Assembler <- assembler_column
groups_melted$Classifier <- classifier_column
groups_melted$Mapper <- mapper_column

# Plot grouped by assembler and classifier
plot2<-ggplot(groups_melted, aes(x = Mapper, y = Reads, fill = Group))+
  #geom_bar(stat = 'identity', position = 'stack', colour="black")+
  geom_bar(stat = 'identity', position = 'stack')+
# facet_grid(Assembler ~ Classifier, scales = "free", space = "free")+
  facet_grid( ~ Assembler + Classifier, scales = "free", space = "free")+
  theme_minimal()+
  scale_fill_manual(values=colorvec1)+
  theme(legend.key.size = unit(0.8,"line"))+
  theme(legend.position="bottom")+
  theme(legend.text=element_text(size=10))+
  theme(legend.title=element_blank())+
  guides(fill=guide_legend(ncol=9,byrow=FALSE))
plot2
#ggsave("Plot2.png", plot=plot2,  device="png",  width=240, units="mm")

# Plot grouped by assembler
plot3<-ggplot(groups_melted, aes(x = Combination, y = Reads, fill = Group))+
  #geom_bar(stat = 'identity', position = 'stack', colour="black")+
  geom_bar(stat = 'identity', position = 'stack')+
  # facet_grid(Assembler ~ Classifier, scales = "free", space = "free")+
  facet_grid( ~ Assembler, scales = "free", space = "free")+
  theme_minimal()+
  scale_fill_manual(values=colorvec1)+
  theme(legend.key.size = unit(0.1,"line"))+
  theme(legend.position="bottom")+
  theme(legend.text=element_text(size=1))+
  theme(legend.title=element_blank())
plot3
#ggsave("Plot3.png", plot=plot3,  device="png",  width=240, units="mm")

# Subsetting merged_ordered for PHYLUMS (only if taxlist_group <- "phylum") that are of low abundance (rowSums < abundance)
if (taxlist_group == "phylum")
{
minimum_abundance <- 100
merged_ordered_subset <- merged_ordered
rownames(merged_ordered_subset) <- merged_ordered_subset$Group
merged_ordered_subset$Group <- NULL
merged_ordered_subset <- subset(merged_ordered_subset, rowSums(merged_ordered_subset) >= minimum_abundance)
}