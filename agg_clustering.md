Information theoretic agglomerative clustering
================

# A simple agglomerative algorithm

Let’s write the variation of information as
![VI(X,Y) = 2H(X,Y) - H(X) + H(Y)](https://latex.codecogs.com/png.latex?VI%28X%2CY%29%20%3D%202H%28X%2CY%29%20-%20H%28X%29%20%2B%20H%28Y%29 "VI(X,Y) = 2H(X,Y) - H(X) + H(Y)").
We then approximate it through its “algorithmic” counterpart:
![VI(X,Y) \\approx 2H(X \\sqcup Y) - H(X) - H(Y)](https://latex.codecogs.com/png.latex?VI%28X%2CY%29%20%5Capprox%202H%28X%20%5Csqcup%20Y%29%20-%20H%28X%29%20-%20H%28Y%29 "VI(X,Y) \approx 2H(X \sqcup Y) - H(X) - H(Y)")
Then the agglomerative algorithm based on Press et al. is as follows:

1.  Maintain an active `char` list `forests` denoting tips of the tree
2.  Parse through `forests` and find the two closest clusters
3.  Update `forest` by merging the two clusters via tree joining
    `t1 + t2`.

Moreover, when two forests
![X](https://latex.codecogs.com/png.latex?X "X") and
![Y](https://latex.codecogs.com/png.latex?Y "Y") are joined into
![(X,Y)](https://latex.codecogs.com/png.latex?%28X%2CY%29 "(X,Y)") the
branch length from the node to each tip
![X](https://latex.codecogs.com/png.latex?X "X") and
![Y](https://latex.codecogs.com/png.latex?Y "Y") is taken to be
![1/2 VI(X,Y)](https://latex.codecogs.com/png.latex?1%2F2%20VI%28X%2CY%29 "1/2 VI(X,Y)").

``` r
site_info <- function(seq, name1, name2) {
  #Computes variation of information VI(name1, name2) at a given site 
  
  seqx <- seq[name1]
  seqy <- seq[name2]

  pxy_all <- base.freq(as.DNAbin(c(seqx, seqy)), all = TRUE)
  p_xy <- pxy_all[c("a", "c", "g", "t", "-")]

  px_all <- base.freq(as.DNAbin(seqx), all = TRUE)
  p_x <- px_all[c("a", "c", "g", "t", "-")]

  # Computing p(y)
  py_all <- base.freq(as.DNAbin(seqy), all = TRUE)
  p_y <- py_all[c("a", "c", "g", "t", "-")]

  entr_xy <- 0
  entr_x <- 0
  entr_y <- 0
  
  for (i in c(1:5)) {
    if (p_xy[i] != 0) {
      entr_xy <- entr_xy - p_xy[i] * log2(p_xy[i])
    }

    if (p_x[i] != 0) {
      entr_x <- entr_x - p_x[i] * log2(p_x[i])
    }

    if (p_y[i] != 0) {
      entr_y <- entr_y - p_y[i] * log2(p_y[i])
    }
  }
  I_alg_site <- 2 * entr_xy - entr_x - entr_y
  return(I_alg_site)
}
```

``` r
alg_info <- function(seq_matrix, x_names, y_names) {
  I_alg <- sum(apply(seq_matrix, 2, site_info, name1 = x_names, name2 = y_names))
  return(I_alg)
}
```

``` r
agg_clustering <- function(sequence) {
  #inputs:
  # sequence -- aligned dna sequence in phyDat
  #ouput:
  # tree in newick format
  tips <- rownames(sequence)
  forests <- make_newick(tips)
  num_sites = ncol(sequence)
  
  if (length(forests) == 1) {
    # Just one species
    tree_string <- forests[1]
  } else if (length(forests) == 2) {
    # Just two species
    # tree_string <-
    #   make_newick(paste(forests[1], ",", forests[2], sep = ""))
    branch <- alg_info(sequence, tips[1], tips[2])/num_sites
    tree_string <- paste("(", tips[1], ":", branch/2, ", ", tips[2], ":", branch/2, ")", sep = "")
  } else{
    #More than two species
    
    end <- length(forests) - 3
    dist_matrix <- matrix(0, end+3, end+3)

    x_names <- read.tree(text = paste(forests[1], ";", sep = ""))$tip.label
    y_names <- read.tree(text = paste(forests[2], ";", sep = ""))$tip.label
    #print(x_names)
    #print(y_names)

    max_dist <- alg_info(sequence, x_names, y_names)
    max_pair <- c(1, 2)

    #Subroutine for computing the closest two clusters

    for (k in 1:(length(forests) - 1)) {
      for (j in (k + 1):length(forests)) {
        x_names <- read.tree(text = paste(forests[k], ";", sep = ""))$tip.label
        y_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
        # dist <-
        #   alg_info(matrix(sequence[x_names, ], nrow = length(y_names)))
        dist <- alg_info(sequence, x_names, y_names)
        cat("Current pair: ", x_names, "/", y_names, "; IG =", dist,"\n")
        dist_matrix[k,j] <- dist
        if (dist < max_dist) {
          max_dist <- dist
          max_pair <- c(k, j)
        }
      }
    }

    #Subroutine for joining the two forests
    new_branch <- paste("(", forests[max_pair[1]], ",", forests[max_pair[2]], ")", sep = "")
    forests <- forests[-max_pair]
    forests <- c(forests, new_branch)
    new_tip <- paste("(", tips[max_pair[1]], ":", max_dist/(num_sites*2), ",", tips[max_pair[2]], ":", max_dist/(num_sites*2), ")", sep = "")
    #print(new_tip)
    tips <- tips[-max_pair]
    tips <- c(tips, new_tip)
    dist_matrix <- dist_matrix[-max_pair, -max_pair]
    dist_matrix <- rbind(dist_matrix, integer(end+1))
    dist_matrix <- cbind(dist_matrix, integer(end+2))
    
    for(i in c(1:end)){
      l <- length(forests)
      #first calculate the new distances
      for(j in 1:(l-1)) {
        x_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
        y_names <- read.tree(text = paste(forests[l], ";", sep = ""))$tip.label
        dist <- alg_info(sequence, x_names, y_names)
        cat("Current pair: ", x_names, "/", y_names, "; IG =", dist,"\n")
        dist_matrix[j,l] <- dist
      }
      
      #find the minimum distance
      dist_matrix[row(dist_matrix)>=col(dist_matrix)] <- NA
      max_pair <- arrayInd(which.min(dist_matrix), dim(dist_matrix))
      max_dist <- dist_matrix[max_pair]
      
      new_branch <- paste("(", forests[max_pair[1]], ",", forests[max_pair[2]], ")", sep = "")
      forests <- forests[-max_pair]
      forests <- c(forests, new_branch)
      new_tip <- paste("(", tips[max_pair[1]], ":", max_dist/(num_sites*2), ",", tips[max_pair[2]], ":", max_dist/(num_sites*2), ")", sep = "")
      #print(new_tip)
      tips <- tips[-max_pair]
      tips <- c(tips, new_tip)
      dist_matrix <- dist_matrix[-max_pair, -max_pair]
      dist_matrix <- rbind(dist_matrix, integer(l-2))
      dist_matrix <- cbind(dist_matrix, integer(l-1))
    }

    x_names <- read.tree(text = paste(forests[1], ";", sep = ""))$tip.label
    y_names <- read.tree(text = paste(forests[2], ";", sep = ""))$tip.label
    branch <- alg_info(sequence, x_names, y_names)/num_sites
    tree_string <- paste("(", tips[1], ":", branch/2, ",", tips[2], ":", branch/2, ")", sep = "")
  }
  #print(tree_string)
  return(tree_string)
}
```

## Neighbor-joining style agglomerative algorithm

Neighbour-joining (NJ) is an agglomerative algorithm that returns an
unrooted tree given all pairwise distances between tips. In fact, if the
distance is ultrametric then NJ gives the exact tree. Here we design the
following information-theoretic variant of NJ. In particular, the
initial distances are assigned to be
![VI(x,y)](https://latex.codecogs.com/png.latex?VI%28x%2Cy%29 "VI(x,y)")
where ![x](https://latex.codecogs.com/png.latex?x "x") and
![y](https://latex.codecogs.com/png.latex?y "y") are single-OTU
clusters. Moreover, in the step where we recalculate distances between
the new cluster and the existing ones, we use the variation of
information instead of the old distances.

    if (1 sequence){
      return that sequence
    } else if (2 sequence){
      return 2 sequences as two nodes of the tree (what is the branch length here?)
    } else{

      n <- # of sequences
      dist_matrix <- a new n by n matrix with 0's
      fill the dist_matrix with VI between any two nodes
        
      for(1 through n-2){ #each step joins two clusters
        1. calculate the u vector from the dist_matrix
        2. find the minimum of the expression dij - ui - uj
        3. save the names of the two joined clusters to x_names and y_names
        4. join the two clusters and calculate the branch lengths
        5. update the distance matrix
      }
        
      finally, we have two clusters left, join them together with branch length each half of their distance to     }

We write it up in R as follows:

``` r
nj_agg <- function(sequence) {
  #inputs:
  # sequence -- aligned dna sequence in phyDat
  #ouput:
  # tree in newick format
  tips <- rownames(sequence)
  forests <- make_newick(tips)
  num_sites = ncol(sequence)
  
  if (length(forests) == 1) {
    # Just one species
    tree_string <- forests[1]
  } else if (length(forests) == 2) {
    # Just two species
    # tree_string <-
    #   make_newick(paste(forests[1], ",", forests[2], sep = ""))
    branch <- alg_info(sequence, tips[1], tips[2])/num_sites
    tree_string <- paste("(", tips[1], ":", branch/2, ", ", tips[2], ":", branch/2, ")", sep = "")
  } else{
    #More than two species
    
    end <- length(forests) - 2
    dist_matrix <- matrix(0, end+2, end+2)

    #Subroutine for computing the closest two clusters

    #fill in the distance matrix
    for (k in 1:(length(forests) - 1)) {
      for (j in (k + 1):length(forests)) {
        x_names <- read.tree(text = paste(forests[k], ";", sep = ""))$tip.label
        y_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
        # dist <-
        #   alg_info(matrix(sequence[x_names, ], nrow = length(y_names)))
        dist <- alg_info(sequence, x_names, y_names)
        cat("Current pair: ", x_names, "/", y_names, "; IG =", dist,"\n")
        dist_matrix[k,j] <- dist
      }
    }
    
    #Now repeat the process for n-2 times 
    for(i in c(1:end)){
      l <- length(forests)
      u <- integer(l)
      
      #calculate the u vector
      for(k in 1:l){
        sum <- 0
        for(j in 1:l){
          if(k < j){
            sum <- sum + dist_matrix[k,j]
          }
          if(k > j){
            sum <- sum + dist_matrix[j,k]
          }
        }
        u[k] <- sum/(l - 2)
        print(u[k])
      }
  
      #find the minimum of the expression dij - ui - uj
      max_dist <- dist_matrix[1,2] - u[1] - u[2]
      max_pair <- c(1,2)
      for (k in 1:(l - 1)) {
        for (j in (k + 1):l) {
          val <- dist_matrix[k,j] - u[k] - u[j]
          x_names <- read.tree(text = paste(forests[k], ";", sep = ""))$tip.label
          y_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
          cat("Between group", k,"which is:", x_names,  "and group", j, "which is:", y_names, "current value is:", val, "\n")
          if(val < max_dist){
            max_dist <- val
            max_pair <- c(k,j)
          }
        }
      }
      
      #save the names for the combined clusters
      x_names <- read.tree(text = paste(forests[max_pair[1]], ";", sep = ""))$tip.label
      y_names <- read.tree(text = paste(forests[max_pair[2]], ";", sep = ""))$tip.label
      
      #join the two clusters
      new_branch <- paste("(", forests[max_pair[1]], ",", forests[max_pair[2]], ")", sep = "")
      forests <- forests[-max_pair]
      forests <- c(forests, new_branch)
      branch1 <- 0.5*(dist_matrix[max_pair[1], max_pair[2]] + u[max_pair[1]] - u[max_pair[2]])
      branch2 <- 0.5*(dist_matrix[max_pair[1], max_pair[2]] - u[max_pair[1]] + u[max_pair[2]])
      new_tip <- paste("(", tips[max_pair[1]], ":", branch1, ",", tips[max_pair[2]], ":", branch2, ")", sep = "")
      #print(new_tip)
      tips <- tips[-max_pair]
      tips <- c(tips, new_tip)
      
      # now update the distance matrix 
        
      # delete the combined entries and add a new row and colomn
      dist_matrix <- dist_matrix[-max_pair, -max_pair]
      dist_matrix <- rbind(dist_matrix, integer(l-2))
      dist_matrix <- cbind(dist_matrix, integer(l-1))
        
      # we only need to fill the last column
      for(j in 1:(l-2)) {
        # x_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
        # y_names <- read.tree(text = paste(forests[l], ";", sep = ""))$tip.label
        # dist <- alg_info(sequence, x_names, y_names)
        # cat("Current pair: ", x_names, "/", y_names, "; IG =", dist,"\n")
        # dist_matrix[j,l] <- dist
        z_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
        
        # This is the information theoretic step 
        d_xz <- alg_info(sequence, x_names, z_names)
        d_yz <- alg_info(sequence, y_names, z_names)
        d_xy <- alg_info(sequence, x_names, y_names)
        dist <- 0.5*(d_xz + d_yz - d_xy)
        dist_matrix[j,l-1] <- dist
        cat("Distance between the pair updated: ", x_names, y_names, "/", z_names, "; d =", dist,"\n")
      }

    }
    
    tree_string <- paste("(", tips[1], ":", 0.5*(dist_matrix[1,2]), ",", tips[2], ":", 0.5*(dist_matrix[1,2]), ")", sep = "")
  }
  #print(tree_string)
  return(tree_string)
}
```

Now let’s see how it performs on the Co3 dataset from Cummings et al:

![](agg_clustering_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Press et al:

    ## Current pair:  1 / 2 ; IG = 8 
    ## Current pair:  1 / 3 ; IG = 8 
    ## Current pair:  1 / 4 ; IG = 8 
    ## Current pair:  1 / 5 ; IG = 14 
    ## Current pair:  1 / 6 ; IG = 14 
    ## Current pair:  1 / 7 ; IG = 16 
    ## Current pair:  1 / 8 ; IG = 12 
    ## Current pair:  2 / 3 ; IG = 12 
    ## Current pair:  2 / 4 ; IG = 10 
    ## Current pair:  2 / 5 ; IG = 14 
    ## Current pair:  2 / 6 ; IG = 16 
    ## Current pair:  2 / 7 ; IG = 18 
    ## Current pair:  2 / 8 ; IG = 14 
    ## Current pair:  3 / 4 ; IG = 8 
    ## Current pair:  3 / 5 ; IG = 10 
    ## Current pair:  3 / 6 ; IG = 14 
    ## Current pair:  3 / 7 ; IG = 14 
    ## Current pair:  3 / 8 ; IG = 14 
    ## Current pair:  4 / 5 ; IG = 18 
    ## Current pair:  4 / 6 ; IG = 14 
    ## Current pair:  4 / 7 ; IG = 14 
    ## Current pair:  4 / 8 ; IG = 14 
    ## Current pair:  5 / 6 ; IG = 8 
    ## Current pair:  5 / 7 ; IG = 14 
    ## Current pair:  5 / 8 ; IG = 14 
    ## Current pair:  6 / 7 ; IG = 14 
    ## Current pair:  6 / 8 ; IG = 14 
    ## Current pair:  7 / 8 ; IG = 8 
    ## Current pair:  3 / 1 2 ; IG = 9.686217 
    ## Current pair:  4 / 1 2 ; IG = 9.182958 
    ## Current pair:  5 / 1 2 ; IG = 12.52933 
    ## Current pair:  6 / 1 2 ; IG = 13.86266 
    ## Current pair:  7 / 1 2 ; IG = 16.52933 
    ## Current pair:  8 / 1 2 ; IG = 12.02607 
    ## Current pair:  5 / 3 4 ; IG = 12.52933 
    ## Current pair:  6 / 3 4 ; IG = 12.52933 
    ## Current pair:  7 / 3 4 ; IG = 12.52933 
    ## Current pair:  8 / 3 4 ; IG = 12.52933 
    ## Current pair:  1 2 / 3 4 ; IG = 7.490225 
    ## Current pair:  5 / 1 2 3 4 ; IG = 12.18758 
    ## Current pair:  6 / 1 2 3 4 ; IG = 12.98758 
    ## Current pair:  7 / 1 2 3 4 ; IG = 13.9334 
    ## Current pair:  8 / 1 2 3 4 ; IG = 11.99149 
    ## Current pair:  7 / 5 6 ; IG = 13.3594 
    ## Current pair:  8 / 5 6 ; IG = 12.52933 
    ## Current pair:  1 2 3 4 / 5 6 ; IG = 10.6627 
    ## Current pair:  1 2 3 4 / 7 8 ; IG = 11.30469 
    ## Current pair:  5 6 / 7 8 ; IG = 11.86767

![](agg_clustering_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# Other agglomerative methods

We’ll use the ones shown in APER

1.  UPGMA

![](agg_clustering_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

2.  Neighbour joining

![](agg_clustering_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
