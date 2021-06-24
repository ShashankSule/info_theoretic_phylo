Information gain model for trees
================
Shashank Sule
14/06/2021

# The model

Let
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}") be a set of OTU’s and let
![T(\\mathcal{S})](https://latex.codecogs.com/png.latex?T%28%5Cmathcal%7BS%7D%29
"T(\\mathcal{S})") be a binary tree associated with
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}"). If ![|\\mathcal{S}| =
n](https://latex.codecogs.com/png.latex?%7C%5Cmathcal%7BS%7D%7C%20%3D%20n
"|\\mathcal{S}| = n") then the number of bifurcations in
![T(\\mathcal{S})](https://latex.codecogs.com/png.latex?T%28%5Cmathcal%7BS%7D%29
"T(\\mathcal{S})") is ![n-1](https://latex.codecogs.com/png.latex?n-1
"n-1") so the task is to figure out the bifurcations of
![T(\\mathcal{S})](https://latex.codecogs.com/png.latex?T%28%5Cmathcal%7BS%7D%29
"T(\\mathcal{S})") (or more directly, figure out a set of sensible
bifurcations ![B\_i](https://latex.codecogs.com/png.latex?B_i "B_i") to
make a tree with
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}") as tips or leaves). In APE lingo, these bifurcations are
called “splits”.

The information gain model of bifurcations/splits/partitions is as
follows: Let
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}") be a set of OTU’s and ![\\mathcal{P} = A \\sqcup
B](https://latex.codecogs.com/png.latex?%5Cmathcal%7BP%7D%20%3D%20A%20%5Csqcup%20B
"\\mathcal{P} = A \\sqcup B") any partition. Supposing that
![A](https://latex.codecogs.com/png.latex?A "A") are realizations of a
random process ![X](https://latex.codecogs.com/png.latex?X "X") and
![B](https://latex.codecogs.com/png.latex?B "B") are realizations of a
random process ![Y](https://latex.codecogs.com/png.latex?Y "Y"), then
the information of
![\\mathcal{P}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BP%7D
"\\mathcal{P}") is ![J(\\mathcal{P}) :=
J(X,Y)](https://latex.codecogs.com/png.latex?J%28%5Cmathcal%7BP%7D%29%20%3A%3D%20J%28X%2CY%29
"J(\\mathcal{P}) := J(X,Y)") where
![J](https://latex.codecogs.com/png.latex?J "J") is some meaningful
information theoretic function of
![X](https://latex.codecogs.com/png.latex?X "X") and
![Y](https://latex.codecogs.com/png.latex?Y
"Y").

# Some comments about the choice of ![J](https://latex.codecogs.com/png.latex?J "J")

Of course, this is a very general model but it’s worth spending some
time about the choices of ![J](https://latex.codecogs.com/png.latex?J
"J"). So far we’ve come up with the following:

1.  Information gain

This model is used to make decision trees from given data
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}") containing data from a set of classes (typically it’s a
binary set). Let ![A \\sqcup B =
\\mathcal{S}](https://latex.codecogs.com/png.latex?A%20%5Csqcup%20B%20%3D%20%5Cmathcal%7BS%7D
"A \\sqcup B = \\mathcal{S}") and ![P\_A =
|A|/|\\mathcal{S}|](https://latex.codecogs.com/png.latex?P_A%20%3D%20%7CA%7C%2F%7C%5Cmathcal%7BS%7D%7C
"P_A = |A|/|\\mathcal{S}|"). Then let ![Z \\sim
X\_\\eta](https://latex.codecogs.com/png.latex?Z%20%5Csim%20X_%5Ceta
"Z \\sim X_\\eta") where ![\\eta
= 1](https://latex.codecogs.com/png.latex?%5Ceta%20%3D%201 "\\eta = 1")
with probability ![P\_A](https://latex.codecogs.com/png.latex?P_A "P_A")
and ![2](https://latex.codecogs.com/png.latex?2 "2") with probability
![1 - P\_A](https://latex.codecogs.com/png.latex?1%20-%20P_A "1 - P_A").
Considering
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}") as a set of samples from
![Z](https://latex.codecogs.com/png.latex?Z "Z") we can actually compute
![I(Z;\\eta) = H(Z) - H(Z \\mid \\eta) = H(\\mathcal{S}) - P\_A H(A) -
(1 - P\_{A})
H(B)](https://latex.codecogs.com/png.latex?I%28Z%3B%5Ceta%29%20%3D%20H%28Z%29%20-%20H%28Z%20%5Cmid%20%5Ceta%29%20%3D%20H%28%5Cmathcal%7BS%7D%29%20-%20P_A%20H%28A%29%20-%20%281%20-%20P_%7BA%7D%29%20H%28B%29
"I(Z;\\eta) = H(Z) - H(Z \\mid \\eta) = H(\\mathcal{S}) - P_A H(A) - (1 - P_{A}) H(B)").
This is the corrected version of the formula at the bottom of page 3 in
`info_theory_ideas`. Using the notation of
![x\_i](https://latex.codecogs.com/png.latex?x_i "x_i") from this
document, we may also write it as ![I(x\_0; (x\_1, x\_2)) = H(x\_0) -
\\alpha H(x\_1) -
(1-\\alpha)H(x\_2)](https://latex.codecogs.com/png.latex?I%28x_0%3B%20%28x_1%2C%20x_2%29%29%20%3D%20H%28x_0%29%20-%20%5Calpha%20H%28x_1%29%20-%20%281-%5Calpha%29H%28x_2%29
"I(x_0; (x_1, x_2)) = H(x_0) - \\alpha H(x_1) - (1-\\alpha)H(x_2)")
where ![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha
"\\alpha") is the proportion of taxa in the partition
![x\_1](https://latex.codecogs.com/png.latex?x_1 "x_1"). Here the
expression “![(x\_1,
x\_2)](https://latex.codecogs.com/png.latex?%28x_1%2C%20x_2%29
"(x_1, x_2)")” is to be understood as the *tree* ![(x\_1,
x\_2)](https://latex.codecogs.com/png.latex?%28x_1%2C%20x_2%29
"(x_1, x_2)") in Newick and not the joint distribution ![(x\_1,
x\_2)](https://latex.codecogs.com/png.latex?%28x_1%2C%20x_2%29
"(x_1, x_2)"). We’ll have to be careful about this notation from now on.
But computationally this strategy gives the most sensible trees for the
simple case of ![n](https://latex.codecogs.com/png.latex?n "n")
sequences of length 1.

2.  Symmetrised cross entropy

Let ![p](https://latex.codecogs.com/png.latex?p "p") be the distribution
of ![x\_1](https://latex.codecogs.com/png.latex?x_1 "x_1") and
![q](https://latex.codecogs.com/png.latex?q "q") the distribution of
![x\_2](https://latex.codecogs.com/png.latex?x_2 "x_2"). Then the cross
entropy is ![H(p,q) := -\\sum\_{x}p(x)\\log\_2
q(x)](https://latex.codecogs.com/png.latex?H%28p%2Cq%29%20%3A%3D%20-%5Csum_%7Bx%7Dp%28x%29%5Clog_2%20q%28x%29
"H(p,q) := -\\sum_{x}p(x)\\log_2 q(x)"). Note that ![H(p,q) \\neq
H(q,p)](https://latex.codecogs.com/png.latex?H%28p%2Cq%29%20%5Cneq%20H%28q%2Cp%29
"H(p,q) \\neq H(q,p)") in general so we let ![J(X,Y) = J(x\_1, x\_2) =
H(p,q) +
H(q,p)](https://latex.codecogs.com/png.latex?J%28X%2CY%29%20%3D%20J%28x_1%2C%20x_2%29%20%3D%20H%28p%2Cq%29%20%2B%20H%28q%2Cp%29
"J(X,Y) = J(x_1, x_2) = H(p,q) + H(q,p)"). The cross entropy is a value
that compares how far ![q](https://latex.codecogs.com/png.latex?q "q")
is from ![p](https://latex.codecogs.com/png.latex?p "p") and in source
coding quantifies the minimum expected code length when the underlying
probability distribution of the source is misunderstood as
![q](https://latex.codecogs.com/png.latex?q "q") instead of
![p](https://latex.codecogs.com/png.latex?p "p").

# Algorithm based on Information Gain

In the case where
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}") is aligned molecular sequence data of length 1
(basically all we have is each species represented through a single
nucleotide), we assume that
![X\_1](https://latex.codecogs.com/png.latex?X_1 "X_1") is a random
variable which takes values in a four-element set (i.e
![\\{a,c,g,t\\}](https://latex.codecogs.com/png.latex?%5C%7Ba%2Cc%2Cg%2Ct%5C%7D
"\\{a,c,g,t\\}")) and ![A](https://latex.codecogs.com/png.latex?A "A")
is the set of realizations of
![X\_1](https://latex.codecogs.com/png.latex?X_1 "X_1") (and similarly
for a random variable ![X\_2](https://latex.codecogs.com/png.latex?X_2
"X_2") whose realisations are represented by
![B](https://latex.codecogs.com/png.latex?B "B")). Let
![\\eta](https://latex.codecogs.com/png.latex?%5Ceta "\\eta") represent
a Bernoulli random variable taking values
![1](https://latex.codecogs.com/png.latex?1 "1") and
![2](https://latex.codecogs.com/png.latex?2 "2") with parameter
![\\alpha =
|A|/|\\mathcal{S}|](https://latex.codecogs.com/png.latex?%5Calpha%20%3D%20%7CA%7C%2F%7C%5Cmathcal%7BS%7D%7C
"\\alpha = |A|/|\\mathcal{S}|"). Then letting ![Z =
X\_\\eta](https://latex.codecogs.com/png.latex?Z%20%3D%20X_%5Ceta
"Z = X_\\eta") the mutual between
![Z](https://latex.codecogs.com/png.latex?Z "Z") and the tree
![\\eta](https://latex.codecogs.com/png.latex?%5Ceta "\\eta") can be
estimated as

  
![I(Z, \\eta) \\approx I(\\mathcal{S}; (A,B)) = H(\\mathcal{S}) -
\\frac{|A|}{|\\mathcal{S}|}H(A) -
\\frac{|B|}{|\\mathcal{S}|}H(B)](https://latex.codecogs.com/png.latex?I%28Z%2C%20%5Ceta%29%20%5Capprox%20I%28%5Cmathcal%7BS%7D%3B%20%28A%2CB%29%29%20%3D%20H%28%5Cmathcal%7BS%7D%29%20-%20%5Cfrac%7B%7CA%7C%7D%7B%7C%5Cmathcal%7BS%7D%7C%7DH%28A%29%20-%20%5Cfrac%7B%7CB%7C%7D%7B%7C%5Cmathcal%7BS%7D%7C%7DH%28B%29
"I(Z, \\eta) \\approx I(\\mathcal{S}; (A,B)) = H(\\mathcal{S}) - \\frac{|A|}{|\\mathcal{S}|}H(A) - \\frac{|B|}{|\\mathcal{S}|}H(B)")  

The optimal partition is a solution to

  
![
\\mathcal{P}^\* = \\text{argmax}\_{(A,B)}\\,I(\\mathcal{S};(A,B))
](https://latex.codecogs.com/png.latex?%0A%5Cmathcal%7BP%7D%5E%2A%20%3D%20%5Ctext%7Bargmax%7D_%7B%28A%2CB%29%7D%5C%2CI%28%5Cmathcal%7BS%7D%3B%28A%2CB%29%29%0A
"
\\mathcal{P}^* = \\text{argmax}_{(A,B)}\\,I(\\mathcal{S};(A,B))
")  

The algorithm `infotree` for making a tree
![T(\\mathcal{S})](https://latex.codecogs.com/png.latex?T%28%5Cmathcal%7BS%7D%29
"T(\\mathcal{S})") with input
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}") is as follows:

1)  Compute the optimal partition of
    ![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
    "\\mathcal{S}") as ![\\mathcal{P}^\* = A \\sqcup
    B](https://latex.codecogs.com/png.latex?%5Cmathcal%7BP%7D%5E%2A%20%3D%20A%20%5Csqcup%20B
    "\\mathcal{P}^* = A \\sqcup B").
2)  Run `infotree` on ![A](https://latex.codecogs.com/png.latex?A "A")
    and ![B](https://latex.codecogs.com/png.latex?B "B") as input.

# Implementation

Let’s try this out on a small example. First we make a tree with say 9
tips that evolves by the Yule model where the probability of a
bifurcation is 1/2.

``` r
t <- rtreeshape(1,9,model = "yule")

plot(as.phylo(t[[1]]))
```

![](info_gain_model_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Using `simseq` we’ll make sequences of length 1 representing the tree

``` r
seqs <- simSeq(as.phylo(t[[1]]),l=1, type = "DNA")
```

Now let’s set up the main body of the algorithm

``` r
mutual_info <- function(partition, sequence, pos){
# inputs:
# partition -- boolean denoting the partitions
# sequence -- dataframe of type DNAbin or phyDat with each row an aligned sequence
# pos -- integer denoting the position in the sequence

# output:
# I(partition)

#computing p(x \oplus y)
  
p_xy <- base.freq(as.DNAbin(sequence))

A <- sequence[partition,pos]
B <- sequence[!partition,pos]

# Computing p(x)

p_x <- base.freq(as.DNAbin(A))
p_a <- length(A)/length(sequence)

# Computing p(y)

p_y <- base.freq(as.DNAbin(B))
p_b <- length(B)/length(sequence)

h_xy <- 0 
h_x <- 0 
h_y <- 0 

# Computing p(x,y)

I <- 0

for(i in c(1:4)){
  if(p_xy[i] !=0){
    h_xy <- h_xy -p_xy[i]*log2(p_xy[i])
  }
  if(p_x[i] != 0){
    h_x <- h_x -p_x[i]*log2(p_x[i])
  }
  if(p_y[i] != 0){
    h_y <- h_y -p_y[i]*log2(p_y[i])
  }
}

I <- h_xy - p_a*h_x - p_b*h_y

# for(i in c(1:4)){
#   if(p_x[i] !=0 ){
#     I <- I - p_a*p_y[i]*log2(p_x[i])
#   }
#   if(p_y[i] != 0){
#     I <- I - p_b*p_x[i]*log2(p_y[i])
#   }
# }
return(I)
}
```

``` r
infopart <- function(sequence){
#input: 
# sequence -- aligned sequence in DNAbin or phyDat 
# output: 
# Newick string representing minimum information gain tree


# if there are only two sequences return dichotomous tree
if(length(sequence) == 1){
 tree_string <- names(sequence)[1] 
} else if(length(sequence) == 2){
  tree_string <- paste("(", names(sequence)[1], ", ", names(sequence)[2], ")", sep = "")
} else{
  
  # There are more than two sequences so we must find the optimal partition. 
  # Initialize the data 
  
  partition <- as.logical(splitset(length(sequence))[2,])
  I <- 0   
  for(j in 1:attr(sequence, "nr")){
          I <- I + mutual_info(partition, sequence, j)
  }
  max_val <- I
  max_part <- partition
  
  for(i in 2:(2^(length(sequence))-1)){ # Run through all possible partitions
    I <- 0 
          # Compute overall mutual information 
    #print(paste("computing the ",i,"th partition")) 
    partition <- as.logical(splitset(length(sequence))[i,])
    
        for(j in 1:attr(sequence, "nr")){
            I <- I + mutual_info(partition, sequence, j)
        }
    
    #print(paste("I =",I))
    
    if(I > max_val){
      max_val <- I 
      max_part <- partition
    }
  }
   #print(paste("The partition is ", max_part))
   left_sequence <- sequence[max_part,]
   right_sequence <- sequence[!max_part,]
   left_string <- infopart(left_sequence)
   right_string <- infopart(right_sequence)
   
   tree_string <- paste("(",left_string,", ",right_string,")", sep = "")

   
}

  return(tree_string)
}
```

# Results

## Sequences of length 1

We check the trees produced by information gain and symmetrized cross
entropy for the 9 tips with sequence length
1.

![](info_gain_model_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->![](info_gain_model_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

It’s a totally sensible tree as it is able to differentiate species with
different nucleotides. Clearly there are three distinct classes:
![\\{1,3,4,6\\},
\\{2\\},](https://latex.codecogs.com/png.latex?%5C%7B1%2C3%2C4%2C6%5C%7D%2C%20%5C%7B2%5C%7D%2C
"\\{1,3,4,6\\}, \\{2\\},") and
![\\{5,8,9\\}](https://latex.codecogs.com/png.latex?%5C%7B5%2C8%2C9%5C%7D
"\\{5,8,9\\}"). Starting from the leaves, each tip is first sorted into
its own class before it is merged with tips from other classes. So this
is a sensible dichotomous tree on our data.

## Sequences of length n

We assume that the distributions on different sites are independent.

``` r
seqs2 <- simSeq(as.phylo(t[[1]]),l=12, type = "DNA")
```

Let’s visualize this sequence first:

![](info_gain_model_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Let’s make the IG tree\!

![](info_gain_model_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

We also have the following toy sequence from Press et al, stored in the
file `press_codes.phy`.

``` r
read.dna("press_codes.phy")
```

    ## 16 DNA sequences in binary format stored in a matrix.
    ## 
    ## All sequences of same length: 12 
    ## 
    ## Labels:
    ## 1
    ## 2
    ## 3
    ## 4
    ## 5
    ## 6
    ## ...
    ## 
    ## Base composition:
    ##     a     c     g     t 
    ## 0.109 0.135 0.547 0.208 
    ## (Total: 192 bases)

We’ll convert it to `phyDat` to start crunching trees:

``` r
trial_DNA <- as.phyDat(read.dna("press_codes.phy"))

trial_DNA
```

    ## 16 sequences with 12 character and 12 different site patterns.
    ## The states are a c g t

Let’s make a tree using `infopart`\! Let’s actually avoid computing the
whole tree for the data since we’d have to make ![2^15
- 1](https://latex.codecogs.com/png.latex?2%5E15%20-%201 "2^15 - 1")
computations. Instead, let’s assume the first (and the most costly)
split: ((0:7),(8:15)). Now we’ll run `infopart` on these partitions and
check the results with Figures 16.6.4-16.6.6 in Press et al.

``` r
layout(matrix(c(1,2), 1, 2))
plot(read.tree(text = paste(infopart(trial_DNA[1:8,]), ";", sep="")))
plot(read.tree(text = paste(infopart(trial_DNA[9:16,]), ";", sep="")))
```

![](info_gain_model_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->
