## LIBRARIES

library(MetaboAnalystR)
library(tidyverse)

file.sources = list.files("C:/Users/vasileioubill95/Desktop/metaboloanalyst/MetaboAnalystR-master/R",
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

rm(list=ls())

## ENRICHMENT ANALYSIS

setwd("C:/Users/vasileioubill95/Desktop/")

# Create vector consisting of compounds for enrichment analysis
tmp.vec <- c("Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine", 
             "Fumaric acid", "Glycine", "Homocysteine", "L-Cysteine", "L-Isolucine", 
             "L-Phenylalanine", "L-Serine", "L-Threonine", "L-Tyrosine", "L-Valine", 
             "Phenylpyruvic acid", "Propionic acid", "Pyruvic acid", "Sarcosine")

# Create mSetObj
mSet<-InitDataObjects("conc", "msetora", FALSE)
#Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, tmp.vec)
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name")

# Example compound name map
mSet$name.map

# Create the mapping results table
mSet<-CreateMappingResultTable(mSet)
# Input the name of the compound without any matches
mSet<-PerformDetailMatch(mSet, "L-Isolucine");
# Create list of candidates to replace the compound
mSet <- GetCandidateList(mSet);
# Identify the name of the compound to replace
mSet<-SetCandidate(mSet, "L-Isolucine", "L-Isoleucine");
# Set the metabolite filter
mSet<-SetMetabolomeFilter(mSet, F);
# Select metabolite set library
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
# Plot the ORA, bar-graph

mSet<-PlotORA(mSetObj = mSet, imgName = "ora_0_", imgOpt = "net",format =  "pdf", dpi = 72, width = 10)


## PATHWAY ENRICHMENT ANALYSIS

rm(list=ls())

# Create vector consisting of compounds for enrichment analysis
tmp.vec <- c("Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine", "Fumaric acid", "Glycine", "Homocysteine", "L-Cysteine", "L-Isolucine", "L-Phenylalanine", "L-Serine", "L-Threonine", "L-Tyrosine", "L-Valine", "Phenylpyruvic acid", "Propionic acid", "Pyruvic acid", "Sarcosine")
# Create mSetObj for storing objects created during your analysis
mSet<-InitDataObjects("conc", "pathora", FALSE)
# Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, tmp.vec);
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name");
# Creates a mapping result table; shows HMDB, KEGG, PubChem, etc. IDs
# Saved as "name_map.csv" or can be found in mSet$dataSet$map.table
# Compounds with no hits will contain NAs across the columns
mSet<-CreateMappingResultTable(mSet);
# From the mapping result table, L-Isolucine has no matches
# Now, perform potential matching with our database against this compound
mSet<-PerformDetailMatch(mSet, "L-Isolucine");
# Get list of candidates for matching
# Results are found in mSet$name.map$hits.candidate.list
mSet<-GetCandidateList(mSet);
# Replace L-Isolucine with selected compound (L-Isoleucine)
mSet<-SetCandidate(mSet, "L-Isolucine", "L-Isoleucine");
# Select the pathway library, ranging from mammals to prokaryotes
# Note the third parameter, where users need to input the KEGG pathway version.

# Use "current" for the latest KEGG pathway library or "v2018" for the KEGG pathway library version prior to November 2019.
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
# Set the metabolite filter
# Default set to false
mSet<-SetMetabolomeFilter(mSet, F);
# Calculate the over representation analysis score, here we selected to use the hypergeometric test (alternative is Fisher's exact test)
# A results table "pathway_results.csv" will be created and found within your working directory
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
# Plot of the Pathway Analysis Overview
mSet<-PlotPathSummary(mSetObj = mSet, show.grid = F, imgName = "path_view_0_",format = "pdf", dpi = 72, width=NA)
# Plot a specific metabolic pathway, in this case "Glycine, serine and threonine metabolism"
mSet<-PlotKEGGPath(mSetObj = mSet, pathName = "Glycine, serine and threonine metabolism" ,width = 528, height =  480, format = "png", dpi =  NULL)

PreparePDFReport(mSet, "My Name")


