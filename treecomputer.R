library("readr")
library("ape")
library("TreeTools")
library("parallel")
library("stringr")
library("seqinr")
library("TreeDist")
source("utilities.R")


args <- commandArgs(trailingOnly = TRUE)

# Computing an array job

serialno <- args[1]
sequencename <- paste("sequence", serialno, ".nex", sep = "")
filepath <- paste(getwd(), "/Sequences/", sequencename, sep = "")

# Computing a directory job

# filename <- args[1]
# trialname <- str_split(filename, "[.]")[[1]][1]
# filepath <- paste(getwd(), "/Cummings_et_al_1995/",filename, sep = "")


sequence <- ReadCharacters(filepath)
# sequence <- as.character.DNAbin(read.dna(filepath))

# Computing a sequence job

start <- Sys.time()
tree <- infotree(sequence)
end <- Sys.time()

# Computing a codon job

# # Preprocess sequence file into codon 

# codon_sequence <- nuc_to_codon(sequence)
# start <- Sys.time()
# tree <- infotree_codon(codon_sequence)
# end <- Sys.time()

# Write to directory (array job)

print(tree)
treename <- paste("tree", serialno, ".txt", sep = "")
targetpath <- paste(getwd(), "/divisive_trees/", treename, sep = "")
write(tree, file = targetpath)

# Write to directory (directory job)

# print(tree)
# treename <- paste(trialname, ".txt", sep = "")
# targetpath <- paste(getwd(), "/1995_trees/", treename, sep = "")
# write(tree, file = targetpath)
# print(end - start)
