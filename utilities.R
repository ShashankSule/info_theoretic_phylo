#---------General utility functions----------------------------------------------

splitset <- function(n){
  #inputs:
  # n -- integer
  #output: 
  # splitset -- a 2^n by n matrix of incident vectors of all possible  of a set of size n 
  
  if(n==1){
    sp <- matrix(c(0,1), nrow = 2, byrow = TRUE)
  } else {
    sp <- matrix(nrow = 2^n, ncol = n)
    sp[1:2^(n-1),ncol(sp)] <- rep(0,nrow(sp)/2)
    sp[(2^(n-1) + 1): nrow(sp),ncol(sp)] <- rep(1,nrow(sp)/2)
    sp[1:2^(n-1),1:(n-1)] <- splitset(n-1)
    sp[(2^(n-1) + 1):2^n, 1:(n-1)] <- splitset(n-1)
  }
  
  return(sp)
}

make_newick <- function(tips){
  # make an n-chotomous tree on an array of strings 
  return(paste("(", tips, ")", sep = ""))
}

DIM <- function(A) {
  if (class(A)[1] == "matrix") {
    size <- dim(A)[1]
  } else {
    size <- 1
  }
  return(size)
}

is_additive <- function(tree){
  ret <- TRUE
  pair_distance <- cophenetic.phylo(tree)
  n <- ncol(pair_distance)
  if(n <= 3){ return(TRUE)
  } else{
    
      #get every subgroup of four
      for(i in c(1:(n-3))){
        for(j in c((i+1):(n-2))){
          for(k in c((j+1):(n-1))){
            for(l in c((k+1):n)){
              d_ij <- pair_distance[i,j]
              d_ik <- pair_distance[i,k]
              d_il <- pair_distance[i,l]
              d_jk <- pair_distance[j,k]
              d_jl <- pair_distance[j,l]
              d_kl <- pair_distance[k,l]
              
              
              four_pt <- sort(c(d_ij + d_kl, d_ik + d_jl, d_il + d_jk), decreasing = TRUE)
              #print(four_pt)
              if( !(all.equal(four_pt[1], four_pt[2])) ){
                ret <- FALSE
              }
            
            }
          }
        }
      }
    
  return(ret)
  }
}  

nuc_to_codon <- function(sequence){
  padding <- ncol(sequence) %% 3 
  #   
  if(padding !=0 ){
    for(i in 1:(3-padding)){
      sequence <- cbind(sequence, rep("-", nrow(sequence)))
    }
  }
  #   
  codon_sequence <- matrix(splitseq(t(sequence)), nrow = nrow(sequence), ncol = ncol(sequence)/3, byrow = TRUE)
  
  rownames(codon_sequence) <- rownames(sequence)
  return(codon_sequence)
  
}

#----------------------------------Divisive Clustering----------------------------------

info_gain_site <- function(sequence, partition) {
  # inputs:
  # partition -- boolean denoting the partitions
  # sequence -- dataframe of type DNAbin or phyDat with each row a single nucleotide
  # pos -- integer denoting the position in the sequence
  # output:
  # I(partition)
  #computing p(x \oplus y)
  
  pxy_all <- base.freq(as.DNAbin(sequence), all = TRUE)
  p_xy <- pxy_all[c("a", "c", "g", "t", "-")]
  
  A <- sequence[partition]
  B <- sequence[!partition]
  
  # Computing p(x)
  px_all <- base.freq(as.DNAbin(A), all = TRUE)
  p_x <- px_all[c("a", "c", "g", "t", "-")]
  
  # Computing p(y)
  py_all <- base.freq(as.DNAbin(B), all = TRUE)
  p_y <- py_all[c("a", "c", "g", "t", "-")]
  
  w_x <- length(A) / length(sequence)
  w_y <- length(B) / length(sequence)
  
  I <- Entropy(p_xy) - w_x*Entropy(p_x) - w_y*Entropy(p_y)
  
}

vi_site <- function(sequence, partition) {
  
  pxy_all <- base.freq(as.DNAbin(sequence), all = TRUE)
  p_xy <- pxy_all[c("a", "c", "g", "t", "-")]
  
  A <- sequence[partition]
  B <- sequence[!partition]
  
  # Computing p(x)
  px_all <- base.freq(as.DNAbin(A), all = TRUE)
  p_x <- px_all[c("a", "c", "g", "t", "-")]
  
  # Computing p(y)
  py_all <- base.freq(as.DNAbin(B), all = TRUE)
  p_y <- py_all[c("a", "c", "g", "t", "-")]
  
  entr_xy <- 0
  entr_x <- 0
  entr_y <- 0
  
  VI <- 2*Entropy(p_xy) - Entropy(p_x) - Entropy(p_y)
  return(VI)
  
}

info_gain <- function(partition, seq) {
  
  # Compute information gain over all sites (assuming all sites are independent)
  part_line <- as.logical(partition)
  I <- c(0,0)
  site_data <- asplit(seq, 2)
  #I <- sum(as.numeric(mclapply(site_data, info_gain_site, partition = part_line)))
  I <- sum(apply(seq, 2, info_gain_site, partition = part_line))
  
  #print(paste("I =", I))
  return(I)
}

vi <- function(partition, seq) {
  # Compute VI over all sites (assuming all sites are independent)
  part_line <- as.logical(partition)
  I <- c(0,0)
  I <- sum(apply(seq, 2, vi_site, partition = part_line))
  
  #print(paste("IG =", I))
  return(I)
}


infotree <- function(sequence, asym = FALSE) {
  #input:
  # sequence -- matrix of characters
  # output:
  # Newick string representing minimum information gain tree
  # if there are only two sequences return dichotomous tree
  
  l = DIM(sequence)
  names = row.names(sequence)
  num_sites = ncol(sequence)
  
  if (l == 1) {
    tree_string <- names[1]
  } else if (l == 2) {
    part_matrix <- splitset(l)[c(2:(2 ^ (l - 1))), ]
    branch <- vi(part_matrix, sequence)/num_sites 
    # asymmetric branch lengths mode
    tree_string <-
      paste("(", names[1], ":", branch/2, ", ", names[2], ":", branch/2, ")", sep = "")
    cat("Done!\n")
  } else{
    # There are more than two sequences so we must find the optimal partition.
    
    cat("Partitioning...")
    part_matrix <- splitset(l)[c(2:(2 ^ (l - 1))), ]
    parts <- asplit(part_matrix,1)
    res <- mclapply(parts, info_gain, seq = sequence)
    max_val <- max(as.numeric(res))
    max_part <- part_matrix[which.max(as.numeric(res)), ]
    # res <- apply(part_matrix, 1, info_gain, seq = sequence)
    # max_val <- max(res)
    # max_part <- part_matrix[which.max(res), ]
    branch <- vi(max_part, sequence)/num_sites
    cur_partition <- as.logical(max_part)
    
    # Partition sequence optimally into left and right sequences 
    
    #print(paste("The partition is ", cur_partition))
    left_sequence <- sequence[cur_partition, , drop = FALSE]
    right_sequence <- sequence[!cur_partition, , drop = FALSE]
    left_string <- infotree(left_sequence, asym = asym)
    right_string <- infotree(right_sequence, asym = asym)
    
    # asymmetric branch length mode
    
    if(asym){
      left_branch <- (nrow(left_sequence)/nrow(sequence))
      right_branch <- (nrow(right_sequence)/nrow(sequence))
      # left_branch <- left_sequence %>% 
      #                asplit(MARGIN = 2) %>%
      #                lapply(table) %>%
      #                lapply(function(x) x/nrow(left_sequence)) %>%
      #                lapply(Entropy) %>%
      #                as.numeric %>%
      #                mean 
      # right_branch <- right_sequence %>%
      #                 asplit(MARGIN = 2) %>%
      #                 lapply(table) %>%
      #                 lapply(function(x) x/nrow(right_sequence)) %>%
      #                 lapply(Entropy) %>%
      #                 as.numeric %>%
      #                 mean
      tree_string <-
        paste("(", left_string, ":", left_branch*branch, ", ", right_string, ":", right_branch*branch, ")", sep = "")  
    } else{
      tree_string <-
        paste("(", left_string, ":", branch/2, ", ", right_string, ":", branch/2, ")", sep = "")
    }
    
    
    
    
  }
  return(tree_string)
}


#---------------------------Divisive Clustering (Codons)--------------------------------------------------------

info_gain_codon_site <- function(sequence, partition) {
  A <- sequence[partition]
  B <- sequence[!partition]
  w_x <- length(A) / length(sequence)
  w_y <- length(B) / length(sequence)
  p_x <- table(A)/length(A)
  p_y <- table(B)/length(B)
  p_xy <- table(sequence)/length(sequence)
  IG <- Entropy(p_xy) - w_x*Entropy(p_x) - w_y*Entropy(p_y)
  return(IG)
}

vi_codon_site <- function(sequence, partition) {
  A <- sequence[partition]
  B <- sequence[!partition]
  w_x <- length(A) / length(sequence)
  w_y <- length(B) / length(sequence)
  p_x <- table(A)/length(A)
  p_y <- table(B)/length(B)
  p_xy <- table(sequence)/length(sequence)
  VI <- 2*Entropy(p_xy) - Entropy(p_x) - Entropy(p_y)
  return(VI)
}

info_gain_codon <- function(sequence, partition) {
  part_line <- as.logical(partition)
  IG <- c(0,0)
  site_data <- asplit(sequence, 2)
  #I <- sum(as.numeric(lapply(site_data, info_gain, partition = part_line)))
  IG <- sum(apply(sequence, 2, info_gain_codon_site, partition = part_line))
  
  #print(paste("I =", I))
  return(IG)
}

vi_codon <- function(sequence, partition) {
  part_line <- as.logical(partition)
  VI <- c(0,0)
  VI <- sum(apply(sequence, 2, vi_codon_site, partition = part_line))
  
  #print(paste("IG =", I))
  return(VI)
}

infotree_codon <- function(sequence) {
  #input:
  # sequence -- matrix of characters
  # output:
  # Newick string representing minimum information gain tree
  # if there are only two sequences return dichotomous tree
  l = DIM(sequence)
  names = row.names(sequence)
  num_sites = ncol(sequence)
  
  if (l == 1) {
    tree_string <- names[1]
  } else if (l == 2) {
    part_matrix <- splitset(l)[c(2:(2 ^ (l - 1))), ]
    branch <- 2*vi_codon(sequence, part_matrix)/num_sites 
    tree_string <-
      paste("(", names[1], ":", branch/2, ", ", names[2], ":", branch/2, ")", sep = "")
    cat("Done!\n")
  } else{
    # There are more than two sequences so we must find the optimal partition.
    
    cat("Partitioning...")
    part_matrix <- splitset(l)[c(2:(2 ^ (l - 1))), ]
    parts <- asplit(part_matrix,1)
    res <- mclapply(parts, info_gain_codon, sequence = sequence)
    max_val <- max(as.numeric(res))
    max_part <- part_matrix[which.max(as.numeric(res)), ]
    #res <- apply(part_matrix, 1, max_info, seq = sequence)
    #max_val <- max(res)
    #max_part <- part_matrix[which.max(res), ]
    branch <- vi_codon(sequence, max_part)/num_sites
    cur_partition <- as.logical(max_part)
    left_sequence <- sequence[cur_partition, , drop = FALSE]
    right_sequence <- sequence[!cur_partition, , drop = FALSE]
    left_string <- infotree_codon(left_sequence)
    right_string <- infotree_codon(right_sequence)
    
    if(asym){
      left_branch <- (nrow(left_sequence)/nrow(sequence))*branch
      right_branch <- (nrow(right_sequence)/nrow(sequence))*branch
      tree_string <-
        paste("(", left_string, ":", left_branch, ", ", right_string, ":", right_branch, ")", sep = "")  
    } else{
      tree_string <-
        paste("(", left_string, ":", branch/2, ", ", right_string, ":", branch/2, ")", sep = "")
    }
    
  }
  return(tree_string)
}

#----------------------------Agglomerative Clustering--------------------------------------------------------------


site_info <- function(seq, name1, name2) {
  #w_x <- length(seqx) / (length(seqx) + length(seqy))
  #w_y <- length(seqy) / (length(seqx) + length(seqy))
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

alg_info <- function(seq_matrix, x_names, y_names) {
  I_alg <- sum(apply(seq_matrix, 2, site_info, name1 = x_names, name2 = y_names))
  return(I_alg)
}

agg_clustering <- function(sequence, asym = FALSE) {
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
      
      if(asym){
        x_opt <- read.tree(text = paste(forests[max_pair[1]], ";", sep = ""))$tip.label
        y_opt <- read.tree(text = paste(forests[max_pair[2]], ";", sep = ""))$tip.label
        left_branch <- (length(x_opt)/nrow(sequence))*(max_dist/num_sites)
        right_branch <- (length(y_opt)/nrow(sequence))*(max_dist/num_sites)
        new_tip <- paste("(", tips[max_pair[1]], ":", left_branch, ",", tips[max_pair[2]], ":", right_branch, ")", sep = "")
      } else{
        new_tip <- paste("(", tips[max_pair[1]], ":", max_dist/(num_sites*2), ",", tips[max_pair[2]], ":", max_dist/(num_sites*2), ")", sep = "")
      }
      
      #print(new_tip)
      forests <- forests[-max_pair]
      forests <- c(forests, new_branch)
      tips <- tips[-max_pair]
      tips <- c(tips, new_tip)
      dist_matrix <- dist_matrix[-max_pair, -max_pair]
      dist_matrix <- rbind(dist_matrix, integer(l-2))
      dist_matrix <- cbind(dist_matrix, integer(l-1))
    }
    
    x_names <- read.tree(text = paste(forests[1], ";", sep = ""))$tip.label
    y_names <- read.tree(text = paste(forests[2], ";", sep = ""))$tip.label
    branch <- alg_info(sequence, x_names, y_names)/num_sites
    
    if(asym){
      left_branch <- (length(x_names)/nrow(sequence))*(max_dist/num_sites)
      right_branch <- (length(y_names)/nrow(sequence))*(max_dist/num_sites)
      tree_string <- paste("(", tips[1], ":", left_branch, ",", tips[2], ":", right_branch, ")", sep = "")
    } else{
      tree_string <- paste("(", tips[1], ":", branch/2, ",", tips[2], ":", branch/2, ")", sep = "")
     }
    
  }
  #print(tree_string)
  return(tree_string)
}

#----------------------------Agglomerative Clustering (Codons)-------------------------------------------------------------

vi_codon_site_names <- function(sequence, name1, name2) {
  A <- sequence[name1]
  B <- sequence[name2]
  p_x <- table(A)/length(A)
  p_y <- table(B)/length(B)
  p_xy <- table(c(A,B))/(length(A) + length(B))
  VI_site <- 2*Entropy(p_xy) - Entropy(p_x) - Entropy(p_y)
  return(VI_site)
}

vi_codon_names <- function(seq_matrix, x_names, y_names) {
  VI <- sum(apply(seq_matrix, 2, vi_codon_site_names, name1 = x_names, name2 = y_names))
  return(VI)
}

agg_clustering_codons <- function(sequence) {
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
    branch <- vi_codon_names(sequence, tips[1], tips[2])/num_sites
    tree_string <- paste("(", tips[1], ":", branch/2, ", ", tips[2], ":", branch/2, ")", sep = "")
  } else{
    #More than two species
    
    end <- length(forests) - 2
    for (i in c(1:end)) {
      #Do this subroutine n-2 times!
      
      x_names <-
        read.tree(text = paste(forests[1], ";", sep = ""))$tip.label
      y_names <-
        read.tree(text = paste(forests[2], ";", sep = ""))$tip.label
      max_dist <- vi_codon_names(sequence, x_names, y_names)
      max_pair <- c(1, 2)
      
      #Subroutine for computing the closest two clusters
      
      for (k in 1:(length(forests) - 1)) {
        for (j in (k + 1):length(forests)) {
          x_names <- read.tree(text = paste(forests[k], ";", sep = ""))$tip.label
          y_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
          # dist <-
          #   alg_info(matrix(sequence[x_names, ], nrow = length(y_names)))
          dist <- vi_codon_names(sequence, x_names, y_names)
          #cat("Current pair: ", x_names, "/", y_names, "; IG =", dist,"\n")
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
    }
    
    x_names <- read.tree(text = paste(forests[1], ";", sep = ""))$tip.label
    y_names <- read.tree(text = paste(forests[2], ";", sep = ""))$tip.label
    branch <- vi_codon_names(sequence, x_names, y_names)/num_sites
    tree_string <- paste("(", tips[1], ":", branch/2, ",", tips[2], ":", branch/2, ")", sep = "")
  }
  #print(tree_string)
  return(tree_string)
}

#----------------------------Neighbor-Joining Agglomerative Algorithm------------------------------------------------------


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
    branch <- alg_info(sequence, tips[1], tips[2])/num_sites
    tree_string <- paste("(", tips[1], ":", branch/2, ", ", tips[2], ":", branch/2, ")", sep = "")
  } else{
    #More than two species

    end <- length(forests) - 2
    dist_matrix <- matrix(0, end+2, end+2)
    
    #fill in the distance matrix
    for (k in 1:(length(forests) - 1)) {
      for (j in (k + 1):length(forests)) {
        x_names <- read.tree(text = paste(forests[k], ";", sep = ""))$tip.label
        y_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
        dist <- alg_info(sequence, x_names, y_names)
        #cat("Current pair: ", x_names, "/", y_names, "; IG =", dist,"\n")
        dist_matrix[k,j] <- dist
      }
    }
    #View(dist_matrix)
    
    #Now repeat the process of joining clusters for n-2 times 
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
        #print(u[k])
      }
      
      #find the minimum of the expression dij - ui - uj
      max_dist <- dist_matrix[1,2] - u[1] - u[2]
      max_pair <- c(1,2)
      for (k in 1:(l - 1)) {
        for (j in (k + 1):l) {
          val <- dist_matrix[k,j] - u[k] - u[j]
          x_names <- read.tree(text = paste(forests[k], ";", sep = ""))$tip.label
          y_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
          #cat("Between group", k,"which is:", x_names,  "and group", j, "which is:", y_names, "current value is:", val, "\n")
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
        z_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
        d_xz <- alg_info(sequence, x_names, z_names)
        d_yz <- alg_info(sequence, y_names, z_names)
        d_xy <- alg_info(sequence, x_names, y_names)
        dist <- 0.5*(d_xz + d_yz - d_xy)
        dist_matrix[j,l-1] <- dist
        #cat("Distance between the pair updated: ", x_names, y_names, "/", z_names, "; d =", dist,"\n")
      }
      
    }
    
    tree_string <- paste("(", tips[1], ":", 0.5*(dist_matrix[1,2]), ",", tips[2], ":", 0.5*(dist_matrix[1,2]), ")", sep = "")
  }
  return(tree_string)
}



#-----------------------------PAML Handling---------------------------------------------------------------------------

who_dat <- function(format = 2,
                    seqs = 5,
                    sites = 10000,
                    reps = 1,
                    birth = 0.1,
                    death = 0.2,
                    sampling= 0.3,
                    mutation = 1.5,
                    model = 3,
                    parameters = "5",
                    gamma = "0 4",
                    equilibrium = "0.1 0.2 0.3 0.4",
                    spit_seed  = FALSE
){
  seed <- 2*sample(1e6:1e7, 1) + 1
  format_string <- paste(format, 
                         "        * 0: paml format (mc.paml); 1:paup format (mc.nex)",sep = "")
  seed_string <- seed   #* random number seed (odd number)
  
  seqs_sites_reps <- paste(seqs, sites, reps, "* <# seqs>  <# nucleotide sites>  <# replicates>", sep = " ")
  
  #"5 10000 1 * <# seqs>  <# nucleotide sites>  <# replicates>"
  
  rates_string <- paste(birth, death, sampling, mutation, sep = " ")
  
  #"0.1 0.2 0.3 1.5   * birth rate, death rate, sampling fraction, and mutation rate (tree height)"
  
  model_string <- model 
  
  #"3          * model: 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV"
  parameter_string <- parameters
  
  #"5 * kappa or rate parameters in model"
  
  gamma_string <- gamma
  
  #"0  4     * <alpha>  <#categories for discrete gamma>"
  
  equilibrium_string <- equilibrium 
  
  #"0.1 0.2 0.3 0.4    * base frequencies
  #  T   C   A   G"
  
  cat(format_string,
      seed_string,
      seqs_sites_reps,
      rates_string,
      model_string,
      parameter_string,
      gamma_string,
      equilibrium_string, 
      file = "MCbaseRTree.dat",
      sep = "\n\n",
      append= FALSE)
  
  if(spit_seed){
    return(seed)
  }
}

ml_tree <- function(){
  system("src/baseml src/baseml.ctl")
  system("tail mlb > tails.txt")
  putfile <- read_lines("tails.txt")
  return(putfile[[length(putfile)]])
}