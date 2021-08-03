# Vignettes 

We tested our algorithms on the following datasets: 

1. [Cummings et al 1995](https://github.com/ShashankSule/info_theoretic_phylo/blob/main/testing/Cummings-et-al-1995-data.md) 
2. [100 samples of TN93 with 16 species and 1000 sites](https://github.com/ShashankSule/info_theoretic_phylo/blob/main/testing/TN93_Jul21.md)
3. [105 samples of JC69 with 12 species and 2000 sites](https://github.com/ShashankSule/info_theoretic_phylo/blob/main/testing/JC69_Jul28.md)
4. [125 samples of JC69 with 12 species and 2000 sites](https://github.com/ShashankSule/info_theoretic_phylo/blob/main/testing/JC69-diagnostics.md)

We computed the following measures of the algorithms on each dataset: 

a. Mean Robinson-Foulds distance between computed trees and the corresponding original trees. 
b. Proportion of trees whose computed topologies were the same as the original (so checking the proportion of computed trees at zero RF distance from the original ones)
c. For those computed trees which were topologically identical to the original trees, we plotted the branch lengths of the computed trees against the corresponding original trees.
d. Proportion of ultrametric and additive trees.
e. Proportion of computed trees which were rooted the same as the corresponding original ones. 
