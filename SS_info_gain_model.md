Information gain model for trees
================
Shashank Sule
14/06/2021

# The model

Let’s write the simplest possible model. Let
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
"\\mathcal{P}") is ![I(\\mathcal{P}) :=
I(X,Y)](https://latex.codecogs.com/png.latex?I%28%5Cmathcal%7BP%7D%29%20%3A%3D%20I%28X%2CY%29
"I(\\mathcal{P}) := I(X,Y)") where
![I](https://latex.codecogs.com/png.latex?I "I") is some meaningful
extension of mutual information to a random
process.

# Algorithm based on ![I(\\mathcal{P})](https://latex.codecogs.com/png.latex?I%28%5Cmathcal%7BP%7D%29 "I(\\mathcal{P})")

In the case where
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}") is aligned molecular sequence data, we assume that ![X =
X\_n](https://latex.codecogs.com/png.latex?X%20%3D%20X_n "X = X_n") is a
random process where ![X\_i](https://latex.codecogs.com/png.latex?X_i
"X_i") takes values in a four-element set (i.e
![\\{A,C,G,T\\}](https://latex.codecogs.com/png.latex?%5C%7BA%2CC%2CG%2CT%5C%7D
"\\{A,C,G,T\\}")) and ![A\_j](https://latex.codecogs.com/png.latex?A_j
"A_j"), the set of elements in the
![j](https://latex.codecogs.com/png.latex?j "j")th sites of
![A](https://latex.codecogs.com/png.latex?A "A") is the set of
realizations of ![X\_n](https://latex.codecogs.com/png.latex?X_n "X_n").
We define the mutual information between
![X](https://latex.codecogs.com/png.latex?X "X") and
![Y](https://latex.codecogs.com/png.latex?Y "Y") to be
![\\frac{1}{n}\\sum\_{i=1}^{n}I(X\_n,
Y\_n)](https://latex.codecogs.com/png.latex?%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Bi%3D1%7D%5E%7Bn%7DI%28X_n%2C%20Y_n%29
"\\frac{1}{n}\\sum_{i=1}^{n}I(X_n, Y_n)"). Here   
![I(X\_n, Y\_n) = \\sum\_{(x,y) \\in
(X,Y)}p(x,y)\\log\_{2}(\\frac{p(x,y)}{p(x)p(y)})](https://latex.codecogs.com/png.latex?I%28X_n%2C%20Y_n%29%20%3D%20%5Csum_%7B%28x%2Cy%29%20%5Cin%20%28X%2CY%29%7Dp%28x%2Cy%29%5Clog_%7B2%7D%28%5Cfrac%7Bp%28x%2Cy%29%7D%7Bp%28x%29p%28y%29%7D%29
"I(X_n, Y_n) = \\sum_{(x,y) \\in (X,Y)}p(x,y)\\log_{2}(\\frac{p(x,y)}{p(x)p(y)})")  
. The strategies for estimating these probabilities are:

  - ![p(x,y) = \\frac{\\text{\# Of pairs (x,y) seen in the data } A\_n
    \\times B\_n}{|A\_n \\times
    B\_n|}](https://latex.codecogs.com/png.latex?p%28x%2Cy%29%20%3D%20%5Cfrac%7B%5Ctext%7B%23%20Of%20pairs%20%28x%2Cy%29%20seen%20in%20the%20data%20%7D%20A_n%20%5Ctimes%20B_n%7D%7B%7CA_n%20%5Ctimes%20B_n%7C%7D
    "p(x,y) = \\frac{\\text{# Of pairs (x,y) seen in the data } A_n \\times B_n}{|A_n \\times B_n|}").
    Note that ![|A\_n \\times B\_n| = (N -
    m)m](https://latex.codecogs.com/png.latex?%7CA_n%20%5Ctimes%20B_n%7C%20%3D%20%28N%20-%20m%29m
    "|A_n \\times B_n| = (N - m)m") where ![A\_n =
    m](https://latex.codecogs.com/png.latex?A_n%20%3D%20m "A_n = m"),
    and each sequence has length
    ![N](https://latex.codecogs.com/png.latex?N "N").

  - ![p(x)](https://latex.codecogs.com/png.latex?p%28x%29 "p(x)") is the
    frequency of ![x](https://latex.codecogs.com/png.latex?x "x") in the
    given site in the data ![A](https://latex.codecogs.com/png.latex?A
    "A") (similarly for ![Y](https://latex.codecogs.com/png.latex?Y
    "Y")).

The optimal partition is a solution to

The algorithm `infopart` for making a tree
![T(\\mathcal{S})](https://latex.codecogs.com/png.latex?T%28%5Cmathcal%7BS%7D%29
"T(\\mathcal{S})") with input
![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
"\\mathcal{S}") is as follows:

1)  Compute the optimal partition of
    ![\\mathcal{S}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BS%7D
    "\\mathcal{S}") as ![\\mathcal{P}^\* = A \\sqcup
    B](https://latex.codecogs.com/png.latex?%5Cmathcal%7BP%7D%5E%2A%20%3D%20A%20%5Csqcup%20B
    "\\mathcal{P}^* = A \\sqcup B").
2)  Run `infopart` on ![A](https://latex.codecogs.com/png.latex?A "A")
    and ![B](https://latex.codecogs.com/png.latex?B "B") as input.
