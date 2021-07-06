#---------all the utility functions----------------------------------------------

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

#----------------------------------Divisive Clustering----------------------------------

mutual_info <- function(partition, sequence, pos) {
  # inputs:
  # partition -- boolean denoting the partitions
  # sequence -- dataframe of type DNAbin or phyDat with each row an aligned sequence
  # pos -- integer denoting the position in the sequence
  # output:
  # I(partition)
  #computing p(x \oplus y)
  
  pxy_all <- base.freq(as.DNAbin(sequence[, pos]), all = TRUE)
  p_xy <- pxy_all[c("a", "c", "g", "t", "-")]
  
  A <- sequence[partition, pos]
  B <- sequence[!partition, pos]
  
  # Computing p(x)
  px_all <- base.freq(as.DNAbin(A), all = TRUE)
  p_x <- px_all[c("a", "c", "g", "t", "-")]
  
  # Computing p(y)
  py_all <- base.freq(as.DNAbin(B), all = TRUE)
  p_y <- py_all[c("a", "c", "g", "t", "-")]
  
  # Computing weight
  #w_x <- length(A)/length(sequence)
  #w_y <- length(B)/length(sequence)
  w_x <- length(A) / DIM(sequence)
  w_y <- length(B) / DIM(sequence)
  
  I <- 0
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
  I <- entr_xy - w_x * entr_x - w_y * entr_y
  #print(paste("I =",I))
  
  #   for(i in c(1:5)){
  # if(p_x[i] !=0 ){
  #   I <- I - w_y*p_y[i]*log2(p_x[i])
  # }
  # if(p_y[i] != 0){
  #   I <- I - w_x*p_x[i]*log2(p_y[i])
  # }
  
  
  return(I)
}

infotree <- function(sequence) {
  #input:
  # sequence -- matrix of characters
  # output:
  # Newick string representing minimum information gain tree
  # if there are only two sequences return dichotomous tree
  l = DIM(sequence)
  names = row.names(sequence)
  
  if (l == 1) {
    tree_string <- names[1]
  } else if (l == 2) {
    tree_string <- paste("(", names[1], ", ", names[2], ")", sep = "")
  } else{
    # There are more than two sequences so we must find the optimal partition.
    # Initialize the data
    
    partition <- as.logical(splitset(l)[2, ])
    I <- 0
    for (j in 1:dim(sequence)[2]) {
      #print(paste("I =",I))
      I <- I + mutual_info(partition, sequence, j)
      
    }
    max_val <- I
    max_part <- partition
    
    for (i in 2:(2 ^ (l - 1))) {
      # Run through all possible partitions
      I <- 0
      # Compute overall mutual information
      #print(paste("computing the ",i,"th partition"))
      partition <- as.logical(splitset(l)[i, ])
      
      for (j in 1:dim(sequence)[2]) {
        I <- I + mutual_info(partition, sequence, j)
      }
      
      print(paste("I =", I))
      
      if (I > max_val) {
        max_val <- I
        max_part <- partition
      }
    }
    print(paste("The partition is ", max_part))
    left_sequence <- sequence[max_part, , drop = FALSE]
    right_sequence <- sequence[!max_part, , drop = FALSE]
    left_string <- infotree(left_sequence)
    right_string <- infotree(right_sequence)
    
    tree_string <-
      paste("(", left_string, ", ", right_string, ")", sep = "")
    
  }
  return(tree_string)
}

#----------------------------Agglomerative Clustering--------------------------------


alg_info <- function(seqx, seqy) {
  n <- dim(seqx)[2] # stores number of sites
  
  w_x <- nrow(seqx) / (nrow(seqx) + nrow(seqy))
  w_y <- nrow(seqy) / (nrow(seqx) + nrow(seqy))
  
  I_alg <- 0
  
  for (j in 1:n) {
    # compute the algorithmic mutual information between two sequences at one site
    
    pxy_all <- base.freq(as.DNAbin(c(seqx[, j], seqy[, j])), all = TRUE)
    p_xy <- pxy_all[c("a", "c", "g", "t", "-")]
    
    px_all <- base.freq(as.DNAbin(seqx[, j]), all = TRUE)
    p_x <- px_all[c("a", "c", "g", "t", "-")]
    
    # Computing p(y)
    py_all <- base.freq(as.DNAbin(seqy[, j]), all = TRUE)
    p_y <- py_all[c("a", "c", "g", "t", "-")]
    
    entr_xy <- 0
    
    entr_x <- 0
    
    entr_y <- 0
    
    I_alg_site <- 0
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
    I_alg <- I_alg + I_alg_site
  }
  
  
  return(I_alg)
}

agg_clustering <- function(sequence) {
  #inputs:
  # sequence -- aligned dna sequence in phyDat
  #ouput:
  # tree in newick format
  forests <- make_newick(rownames(sequence))
  
  if (length(forests) == 1) {
    # Just one species
    tree_string <- forests[1]
  } else if (length(forests) == 2) {
    # Just two species
    tree_string <-
      make_newick(paste(forests[1], ",", forests[2], sep = ""))
  } else{
    #More than two species
    
    end <- length(forests) - 2
    for (i in c(1:end)) {
      #Do this subroutine n-2 times!
      
      x_names <- read.tree(text = paste(forests[1], ";", sep = ""))$tip.label
      y_names <- read.tree(text = paste(forests[2], ";", sep = ""))$tip.label
      print(x_names)
      print(y_names)
      
      max_dist <-
        alg_info(matrix(sequence[x_names, ], nrow = length(x_names)),
                 matrix(sequence[y_names, ], nrow = length(y_names)))
      max_pair <- c(1, 2)
      
      #Subroutine for computing the closest two clusters
      
      for (k in 1:(length(forests) - 1)) {
        for (j in (k + 1):length(forests)) {
          x_names <- read.tree(text = paste(forests[k], ";", sep = ""))$tip.label
          y_names <- read.tree(text = paste(forests[j], ";", sep = ""))$tip.label
          dist <-
            alg_info(matrix(sequence[x_names, ], nrow = length(x_names)),
                     matrix(sequence[y_names, ], nrow = length(y_names)))
          cat("Current pair: ",
              x_names,
              "/",
              y_names,
              "; affinity =",
              dist,
              "\n")
          if (dist < max_dist) {
            max_dist <- dist
            max_pair <- c(k, j)
          }
          
        }
      }
      
      #Just some diagnostics
      
      
      #Subroutine for joining the two forests
      new_branch <-
        paste("(", forests[max_pair[1]], ",", forests[max_pair[2]], ")", sep =
                "")
      forests <- forests[-max_pair]
      forests <- c(forests, new_branch)
    }
    
    tree_string <-
      make_newick(paste(forests[1], ",", forests[2], sep = ""))
  }
  
}

#-----------------------------Sequence Generation--------------------------------------

hit_dat <- function(format = 2,
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
                    equilibrium = "0.1 0.2 0.3 0.4"
){
  seed <- 2*sample(1e6:1e7, 1) + 1
  format_string <- paste(format, 
                         "        * 0: paml format (mc.paml); 1:paup format (mc.nex)",sep = "")
  seed_string <- seed   #* random number seed (odd number)
  
  seqs_sites_reps <- paste(seqs, sites, reps, "<# seqs>  <# nucleotide sites>  <# replicates>", sep = " ")
  
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
      file = "paml4.8/MCbaseRTree.dat",
      sep = "\n\n",
      append= FALSE)
}
