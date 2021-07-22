library("readr")
library("ape")
library("TreeTools")
library("parallel")
source("utilities.R")

#sequence <- ReadCharacters("coiii.nex")

sequence <- ReadCharacters("press_data.nex")

start <- Sys.time()
tree <- infotree(sequence)
end <- Sys.time()

write(tree, file = "computed.txt")
print(end - start)