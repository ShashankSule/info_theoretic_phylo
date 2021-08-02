Information gain model for trees
================
14/06/2021

# Introduction

Let ![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0") be a set of
OTU’s and let
![T(x\_0)](https://latex.codecogs.com/png.latex?T%28x_0%29 "T(x_0)") be
a binary tree associated with
![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0"). If
![\|x\_0\| = n](https://latex.codecogs.com/png.latex?%7Cx_0%7C%20%3D%20n "|x_0| = n")
then the number of bifurcations in
![T(\\mathcal{S})](https://latex.codecogs.com/png.latex?T%28%5Cmathcal%7BS%7D%29 "T(\mathcal{S})")
is ![n-1](https://latex.codecogs.com/png.latex?n-1 "n-1") so the task is
to figure out the bifurcations of
![T(\\mathcal{S})](https://latex.codecogs.com/png.latex?T%28%5Cmathcal%7BS%7D%29 "T(\mathcal{S})")
(or more directly, figure out a set of sensible bifurcations
![B\_i](https://latex.codecogs.com/png.latex?B_i "B_i") to make a tree
with
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D "\mathcal{S}")
as tips or leaves).

The information gain model of computing bifurcations/splits/partitions
is as follows: Let
![\\mathcal{P} = x\_1 \\sqcup x\_2](https://latex.codecogs.com/png.latex?%5Cmathcal%7BP%7D%20%3D%20x_1%20%5Csqcup%20x_2 "\mathcal{P} = x_1 \sqcup x_2")
be a partition of
![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0"). If
![x\_1](https://latex.codecogs.com/png.latex?x_1 "x_1") are realizations
of a random variable
![X\_1](https://latex.codecogs.com/png.latex?X_1 "X_1") and
![x\_2](https://latex.codecogs.com/png.latex?x_2 "x_2") are realizations
of a random variable
![X\_2](https://latex.codecogs.com/png.latex?X_2 "X_2"), then
![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0") is the set of
realizations of a random variable
![Y = X\_\\eta](https://latex.codecogs.com/png.latex?Y%20%3D%20X_%5Ceta "Y = X_\eta")
where
![P(\\eta = 1) = p](https://latex.codecogs.com/png.latex?P%28%5Ceta%20%3D%201%29%20%3D%20p "P(\eta = 1) = p")
and
![P(\\eta = 2) = 1-p](https://latex.codecogs.com/png.latex?P%28%5Ceta%20%3D%202%29%20%3D%201-p "P(\eta = 2) = 1-p").
The information gain of
![\\mathcal{P}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BP%7D "\mathcal{P}")
is then

![\\begin{equation}
I(Y, \\eta) = H(Y) - H(X\_\\eta \\mid \\eta) = H(Y) - pH(X\_1) - (1-p)H(X\_2) \\approx H(x\_0) - \\frac{x\_1}{x\_0}H(X\_1) - \\frac{x\_2}{x\_0}H(X\_2)
\\end{equation}](https://latex.codecogs.com/png.latex?%5Cbegin%7Bequation%7D%0AI%28Y%2C%20%5Ceta%29%20%3D%20H%28Y%29%20-%20H%28X_%5Ceta%20%5Cmid%20%5Ceta%29%20%3D%20H%28Y%29%20-%20pH%28X_1%29%20-%20%281-p%29H%28X_2%29%20%5Capprox%20H%28x_0%29%20-%20%5Cfrac%7Bx_1%7D%7Bx_0%7DH%28X_1%29%20-%20%5Cfrac%7Bx_2%7D%7Bx_0%7DH%28X_2%29%0A%5Cend%7Bequation%7D "\begin{equation}
I(Y, \eta) = H(Y) - H(X_\eta \mid \eta) = H(Y) - pH(X_1) - (1-p)H(X_2) \approx H(x_0) - \frac{x_1}{x_0}H(X_1) - \frac{x_2}{x_0}H(X_2)
\end{equation}")

We pick the partition
![\\mathcal{P}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BP%7D "\mathcal{P}")
that maximizes the information gain on the data
![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0").

This strategy of picking the information gain-maximizing partition
closely follows decision tree learning in binary classification, where
at each node the data is classified according to the labels of the
information gain maximizing feature. In such a case the split at a node
is induced by the labels of the optimal feature
![f^\* = \\underset{\\text{feature}}{argmax} I (data, feature)](https://latex.codecogs.com/png.latex?f%5E%2A%20%3D%20%5Cunderset%7B%5Ctext%7Bfeature%7D%7D%7Bargmax%7D%20I%20%28data%2C%20feature%29 "f^* = \underset{\text{feature}}{argmax} I (data, feature)").
Of course, in our case, when each OTU is represented by a nucleotide (or
codon) sequence, we assume that we don’t have feature-level binary data
so we search over all partitions of
![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0") instead.
Moreover, this strategy of making a tree is divisive or “top-down”: at
each stage we take data
![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0") and split it
into two clusters
![x\_1](https://latex.codecogs.com/png.latex?x_1 "x_1") and
![x\_2](https://latex.codecogs.com/png.latex?x_2 "x_2"). We obtain a
tree by repeating this splitting process on
![x\_1](https://latex.codecogs.com/png.latex?x_1 "x_1") and
![x\_2](https://latex.codecogs.com/png.latex?x_2 "x_2") recursively. In
summary the algorithm is as follows:

    ALGORITHM infotree

    1. Start with data $x_0$
    2. If $x_0$ has one or two species, return the trees (OTU 1) or (OTU 1, OTU 2) respectively. 
    3. If $x_0$ has three more more species do: 
       
       Partition $x_0$ into $x_1, x_2 = \text{argmax}_{(a,b)}\,I(x_0;(a,b))$
       
       Compute: 
       
       a. Left tree = infotree(x_1)
       b. Right tree = infotree(x_2)
       
       Return the tree (Left tree, Right tree)

# Computing IG

To implement the above algorithm we would need to compute the
information gain
![I(x\_0; (x\_1,x\_2))](https://latex.codecogs.com/png.latex?I%28x_0%3B%20%28x_1%2Cx_2%29%29 "I(x_0; (x_1,x_2))")
which in turn involves computing the (empirical) entropies
![H(x\_0)](https://latex.codecogs.com/png.latex?H%28x_0%29 "H(x_0)"),
![H(x\_1)](https://latex.codecogs.com/png.latex?H%28x_1%29 "H(x_1)"),
and
![H(x\_2)](https://latex.codecogs.com/png.latex?H%28x_2%29 "H(x_2)"). In
the case where ![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0")
is sequence data of length 1 (essentially all we have is each species
represented by a single nucleotide), we assume that
![X\_i](https://latex.codecogs.com/png.latex?X_i "X_i") is a random
variable which takes values in a four-element set (i.e
![\\{a,c,g,t\\}](https://latex.codecogs.com/png.latex?%5C%7Ba%2Cc%2Cg%2Ct%5C%7D "\{a,c,g,t\}"))
and ![x\_i](https://latex.codecogs.com/png.latex?x_i "x_i") is the set
of realizations of
![X\_i](https://latex.codecogs.com/png.latex?X_i "X_i"). Then we can
compute the entropies
![H(x\_i)](https://latex.codecogs.com/png.latex?H%28x_i%29 "H(x_i)") as

![H(x\_i) = -\\sum\_{j \\in {a,c,g,t}}log\_2\\Big(\\frac{\|j\|}{\|x\_i\|}\\Big)\\frac{\|j\|}{\|x\_i\|}](https://latex.codecogs.com/png.latex?H%28x_i%29%20%3D%20-%5Csum_%7Bj%20%5Cin%20%7Ba%2Cc%2Cg%2Ct%7D%7Dlog_2%5CBig%28%5Cfrac%7B%7Cj%7C%7D%7B%7Cx_i%7C%7D%5CBig%29%5Cfrac%7B%7Cj%7C%7D%7B%7Cx_i%7C%7D "H(x_i) = -\sum_{j \in {a,c,g,t}}log_2\Big(\frac{|j|}{|x_i|}\Big)\frac{|j|}{|x_i|}")

Here ![\|j\|](https://latex.codecogs.com/png.latex?%7Cj%7C "|j|")
represents the number of observations of the outcome
![j](https://latex.codecogs.com/png.latex?j "j") in the data
![x\_i](https://latex.codecogs.com/png.latex?x_i "x_i"). With this, we
can now set up step 3 of `infotree`.

## Sequence data of length N &gt; 1

Most sequence alignment data, however, has more than one nucleotide site
(surprise). In the case when each OTU is represented by a nucleotide
sequence of length ![N](https://latex.codecogs.com/png.latex?N "N"), the
random variables ![X\_i](https://latex.codecogs.com/png.latex?X_i "X_i")
are of the form
![(Y\_{i}^{1}, \\ldots, Y\_{i}^{N})](https://latex.codecogs.com/png.latex?%28Y_%7Bi%7D%5E%7B1%7D%2C%20%5Cldots%2C%20Y_%7Bi%7D%5E%7BN%7D%29 "(Y_{i}^{1}, \ldots, Y_{i}^{N})")
where
![Y\_{i}^{j}](https://latex.codecogs.com/png.latex?Y_%7Bi%7D%5E%7Bj%7D "Y_{i}^{j}")
is the random variable that takes values in
![\\{a,c,g,t\\}](https://latex.codecogs.com/png.latex?%5C%7Ba%2Cc%2Cg%2Ct%5C%7D "\{a,c,g,t\}")
and describes the ![j](https://latex.codecogs.com/png.latex?j "j")th
site in ![X\_i](https://latex.codecogs.com/png.latex?X_i "X_i") for
![i=0,1,2](https://latex.codecogs.com/png.latex?i%3D0%2C1%2C2 "i=0,1,2").
Consequently, the entropies in equation (1) can be written as joint
entropies:

![H(X\_i) = H(Y\_{i}^{1}, \\ldots, Y\_{i}^{N})](https://latex.codecogs.com/png.latex?H%28X_i%29%20%3D%20H%28Y_%7Bi%7D%5E%7B1%7D%2C%20%5Cldots%2C%20Y_%7Bi%7D%5E%7BN%7D%29 "H(X_i) = H(Y_{i}^{1}, \ldots, Y_{i}^{N})").

We assume that all the sites are independent so the joint entropy can be
written as the sum of entropies at each site:

![H(X\_i) = H(Y\_{i}^{1}, \\ldots, Y\_{i}^{N}) = \\sum\_{j=1}^{N} H(Y\_{i}^{j})](https://latex.codecogs.com/png.latex?H%28X_i%29%20%3D%20H%28Y_%7Bi%7D%5E%7B1%7D%2C%20%5Cldots%2C%20Y_%7Bi%7D%5E%7BN%7D%29%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7BN%7D%20H%28Y_%7Bi%7D%5E%7Bj%7D%29 "H(X_i) = H(Y_{i}^{1}, \ldots, Y_{i}^{N}) = \sum_{j=1}^{N} H(Y_{i}^{j})").

Consequently, for sequence data of length
![N](https://latex.codecogs.com/png.latex?N "N") &gt; 1, the entropy
calculations become

![H(x\_i) = -\\sum\_{j=1}^{N}\\sum\_{y \\in {a,c,g,t}}log\_2\\Big(\\frac{\|y(j)\|}{\|x\_i\|}\\Big)\\frac{\|y(j)\|}{\|x\_i\|}](https://latex.codecogs.com/png.latex?H%28x_i%29%20%3D%20-%5Csum_%7Bj%3D1%7D%5E%7BN%7D%5Csum_%7By%20%5Cin%20%7Ba%2Cc%2Cg%2Ct%7D%7Dlog_2%5CBig%28%5Cfrac%7B%7Cy%28j%29%7C%7D%7B%7Cx_i%7C%7D%5CBig%29%5Cfrac%7B%7Cy%28j%29%7C%7D%7B%7Cx_i%7C%7D "H(x_i) = -\sum_{j=1}^{N}\sum_{y \in {a,c,g,t}}log_2\Big(\frac{|y(j)|}{|x_i|}\Big)\frac{|y(j)|}{|x_i|}")

Here
![\|y(j)\|](https://latex.codecogs.com/png.latex?%7Cy%28j%29%7C "|y(j)|")
is the number of observations of the outcome
![y(j)](https://latex.codecogs.com/png.latex?y%28j%29 "y(j)") in site
![j](https://latex.codecogs.com/png.latex?j "j").

## Branch lengths

There are two important aspects to a phylogenetc tree: its topology and
its branch lengths. Thus far our algorithm only computes the topology,
so we must also build in an information theoretic notion of a branch
length. In particular, we use the *Variation of Information* to describe
the distance between the sets
![x\_1](https://latex.codecogs.com/png.latex?x_1 "x_1") and
![x\_2](https://latex.codecogs.com/png.latex?x_2 "x_2"). Recall that the
variation of information (VI) between two random variables is
![VI(X,Y) = 2H(X,Y) - H(X) + H(Y) = H(X\|Y) + H(Y\|X)](https://latex.codecogs.com/png.latex?VI%28X%2CY%29%20%3D%202H%28X%2CY%29%20-%20H%28X%29%20%2B%20H%28Y%29%20%3D%20H%28X%7CY%29%20%2B%20H%28Y%7CX%29 "VI(X,Y) = 2H(X,Y) - H(X) + H(Y) = H(X|Y) + H(Y|X)").
Note that ![VI](https://latex.codecogs.com/png.latex?VI "VI") satisfies
the triangle inequality so it serves as a distance on the space of
random variables. In our case we approximate it through its
“algorithmic” counterpart:
![VI(X,Y) \\approx 2H(X \\cup Y) - H(X) - H(Y)](https://latex.codecogs.com/png.latex?VI%28X%2CY%29%20%5Capprox%202H%28X%20%5Ccup%20Y%29%20-%20H%28X%29%20-%20H%28Y%29 "VI(X,Y) \approx 2H(X \cup Y) - H(X) - H(Y)")
so the
![VI(x\_1,x\_2) = 2H(x\_0) - H(x\_1) - H(x\_2)](https://latex.codecogs.com/png.latex?VI%28x_1%2Cx_2%29%20%3D%202H%28x_0%29%20-%20H%28x_1%29%20-%20H%28x_2%29 "VI(x_1,x_2) = 2H(x_0) - H(x_1) - H(x_2)").
Since the total length between
![x\_1](https://latex.codecogs.com/png.latex?x_1 "x_1") and
![x\_2](https://latex.codecogs.com/png.latex?x_2 "x_2") is
![VI(x\_1,x\_2)](https://latex.codecogs.com/png.latex?VI%28x_1%2Cx_2%29 "VI(x_1,x_2)"),
we assign the branch length between the tips
![x\_i](https://latex.codecogs.com/png.latex?x_i "x_i") and the root
![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0") as
![1/2 VI(x\_1, x\_2)](https://latex.codecogs.com/png.latex?1%2F2%20VI%28x_1%2C%20x_2%29 "1/2 VI(x_1, x_2)").
Note that in this case, both
![x\_1](https://latex.codecogs.com/png.latex?x_1 "x_1") and
![x\_2](https://latex.codecogs.com/png.latex?x_2 "x_2") are equally far
from the root. There is also a way to assign branch lengths
asymmetrically by weighting the
![VI](https://latex.codecogs.com/png.latex?VI "VI") according to the
weight
![p \\approx \|x\_1\|/\|x\_0\|](https://latex.codecogs.com/png.latex?p%20%5Capprox%20%7Cx_1%7C%2F%7Cx_0%7C "p \approx |x_1|/|x_0|")
from Equation 1. With the branch length modification, `infotree` is the
following algorithm:

    ALGORITHM infotree

    1. Start with data $x_0$
    2. If $x_0$ has one or two species, return the trees (OTU 1) or (OTU 1, OTU 2) respectively. 
    3. If $x_0$ has three more more species do: 
       
       Partition $x_0$ into $x_1, x_2 = \text{argmax}_{(a,b)}\,I(x_0;(a,b))$
       
       Compute: 
       
       a. Left tree = infotree(x_1)
       b. Right tree = infotree(x_2)
       
       if branch lengths are asymmetric do
          left branch = |x_1|/|x_0| * VI(x_1,x_2)
          right branch = |x_2|/|x_0| * VI(x_1,x_2)
       else 
         left branch = 1/2 * VI(x_1,x_2)
         right branch = 1/2 * VI(x_1,x_2)
         
       Return the tree (Left tree: left branch, Right tree: right branch)

# Implementation

We implement this algorithm in `R`, taking advantage of the `Entropy`
function in `TreeTools` and the classes `DNAbin` and `phylo` which
enable quick computation of base frequencies and Newick formats.

``` r
info_gain_site <- function(sequence, partition) {
  # inputs:
  # partition -- boolean denoting the partitions
  # sequence -- dataframe of type DNAbin or phyDat with each row a single nucleotide
  # pos -- integer denoting the position in the sequence
  # output:
  # I(partition)
  #computing 
  
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
    
    left_sequence <- sequence[cur_partition, , drop = FALSE]
    right_sequence <- sequence[!cur_partition, , drop = FALSE]
    left_string <- infotree(left_sequence, asym = asym)
    right_string <- infotree(right_sequence, asym = asym)
    
    # asymmetric branch length mode
    
    if(asym){
      left_branch <- (nrow(left_sequence)/nrow(sequence))
      right_branch <- (nrow(right_sequence)/nrow(sequence))
      tree_string <-
        paste("(", left_string, ":", left_branch*branch, ", ", right_string, ":", right_branch*branch, ")", sep = "")  
    } else{
      tree_string <-
        paste("(", left_string, ":", branch/2, ", ", right_string, ":", branch/2, ")", sep = "")
    }
    
    
    
    
  }
  return(tree_string)
}
```

# Examples

## Length 1 sequence

Let’s try this out on a small example. Using `simseq` we’ll make
sequences of length 1 representing the tree

Using `simseq` we’ll make sequences of length 1 representing the tree

![](info_gain_model_files/figure-gfm/unnamed-chunk-4-1.png)<!-- --> Now
let’s look at the tree produced by the information gain-based divisive
algorithm!

``` r
plot(read.tree(text = paste(infotree(as.character.DNAbin(as.DNAbin(seqs))),";", sep = "")))
```

    ## Partitioning...Partitioning...Partitioning...Done!
    ## Partitioning...Done!
    ## Partitioning...Done!

![](info_gain_model_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> This
tree, while not very informative, is the most sensible given the
original sequence! Traversing from the tips, all the OTU’s with the same
nucleotide are put in the same clade before merging with the others.
Note that `infotree` always assumes that the OTUs have a common ancestor
(the root), which may not be the best assumption in the single
nucleotide case. But anyway, the concept seems to work!

## Data from Press et al

We’ll compute a tree from the first half of the nucleotide dataset of 16
OTU’s found in Chapter 16.4 of the Third Edition of Numerical Recipes.

    ## Partitioning...Partitioning...Done!
    ## Done!
    ## Partitioning...Done!
    ## Done!

![](info_gain_model_files/figure-gfm/unnamed-chunk-6-1.png)<!-- --> Note
that this is the same tree as the original one!

Despite its stunning success in all of two examples, `infotree` isn’t
always perfect. Have a look at some more examples in the tests folder.

# Notes on computation

1.  `infotree` can be slow. This is because we optimize the information
    gain over all possible partitions of
    ![x\_0](https://latex.codecogs.com/png.latex?x_0 "x_0") by actually
    evaluating information gains for each partition and then finding the
    maximum. That amounts to
    ![2^{n-1} - 1](https://latex.codecogs.com/png.latex?2%5E%7Bn-1%7D%20-%201 "2^{n-1} - 1")
    computations where
    ![n = \|x\_0\|](https://latex.codecogs.com/png.latex?n%20%3D%20%7Cx_0%7C "n = |x_0|").
    So even for 16 OTUs we make 32767 computations. That may not seem
    like a lot, but computing the entropies over many sites can take a
    while. For reference, here’s how long the `info_gain` computation
    takes for a sequence block of length 1000 with 16 OTU’s for a single
    partition.

``` r
start<- Sys.time()
info_gain(partition, sequence)
```

    ## [1] 128.8878

``` r
end <- Sys.time()
print(end - start)
```

    ## Time difference of 0.1390619 secs

It takes about 0.2 seconds for each run, so over 32767 runs it takes
about 6553 seconds, which is about a 110 minutes, if we do these
computations sequentially!

2.  To quicken the pace of computation, we can parallelize the
    `infotree` procedure by realizing that:

<!-- -->

1.  The information gain computations over all partitions can be done in
    parallel.
2.  The site information gains can also be computed in parallel (since
    the sites are independent)
3.  The left and right trees can be computed in parallel

First of all, all these operations should be vectorized in order to take
advantage of the quick `apply` functions in R. Second, we can
parallelize these using `parallel::mclapply` but we must be careful: if
we parallelize both (a) and (b) via `mclapply`, `R` will end up
distributing only one core to `b`. Empirically, we find that it is
quicker to only use `mclapply` for task (a).

3.  We also have a codon variant of the above algorithm where instead of
    computing over nucleotide sequence data, we compute over codon data.
    The subroutine for computing the information gain in a codon site on
    a partition is the following:

``` r
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
```

Note that the above code is general for any form of data because of the
use of `table` so you can even compute information gain for nucleotide
data using `info_gain_codon_site`. But we still recommend using
`info_gain_site` for nucleotide data because `table` is slower than
`ape`’s `base.freq` at computing base frequencies.

``` r
start <- Sys.time()
info_gain_codon_site(site_data, as.logical(partition))
```

    ## [1] 0.03757964

``` r
end <- Sys.time()
print(end - start)
```

    ## Time difference of 0.04480696 secs

``` r
start <- Sys.time()
info_gain_site(site_data, as.logical(partition))
end <- Sys.time()
print(end - start)
```

    ## Time difference of 0.004918098 secs
