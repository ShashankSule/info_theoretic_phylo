Experiments with information theoretic tree generation
================

We’ll generate 100 blocks of sequence data for 16 species where each
block corresponds to a tree randomly generated in the TN93 model.

``` r
trees <- list()
sequences <- list()

for(i in 1:10){
  who_dat(sites = 1000, 
          model = 6, 
          parameters = "5 5", 
          gamma = "0.5 4", 
          equilibrium = "0.25 0.28 0.34 0.13")
  system("paml4.8/src/evolverRandomTree 5 paml4.8/MCbaseRTree.dat") 
  trees[[i]] <- write.tree(read.nexus("mctrees.nex"))
  sequences[[i]] <- ReadCharacters("mc.nex")
}
```

Let’s compute some diagnostics\!

1.  Robinson-Foulds distances:

<!-- end list -->

``` r
#trials <- c(1:100)
divisive <- lapply(lapply(sequences, infotree), function(x) paste(x,";", sep = ""))
agglomerative <- lapply(lapply(sequences, agg_clustering), function(x) paste(x,";", sep = ""))
nj_trees <- sequences %>%
            lapply(as.DNAbin) %>%
            lapply(dist.dna) %>%
            lapply(nj) %>%
            lapply(write.tree)
```

    ## [1] "Average distance of divisive trees from original: 0"

    ## [1] "Average distance of agglomerative trees from original: 0"

    ## [1] "Average distance of NJ trees from original: 0.2"

2.  Low tolerance
    ultrametricity:

<!-- end list -->

``` r
summary(metricity)
```

    ##  as.character.trees. as.character.divisive. as.character.agglomerative.
    ##  Mode:logical        Mode :logical          Mode :logical              
    ##  TRUE:10             FALSE:10               FALSE:10                   
    ##  as.character.nj_trees.
    ##  Mode :logical         
    ##  FALSE:10

3.  High tolerance
    ultrametricity:

<!-- end list -->

``` r
summary(low_metricity)
```

    ##  as.character.trees. as.character.divisive. as.character.agglomerative.
    ##  Mode:logical        Mode :logical          Mode :logical              
    ##  TRUE:10             FALSE:10               FALSE:10                   
    ##  as.character.nj_trees.
    ##  Mode :logical         
    ##  FALSE:10
