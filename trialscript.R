library("readr")
library("ape")
library("TreeTools")
library("parallel")
source("utlities.R")

setwd("/Users/shashanksule/Documents/info_theoretic_phylo/")
sequence <- ReadCharacters("coiii.nex")

tree <- infotree(sequence)
write(tree, file = "computed.txt")
