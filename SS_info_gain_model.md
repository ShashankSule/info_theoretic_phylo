The model
=========

Let’s write the simplest possible model. Let 𝒮 be a set of OTU’s and let
*T*(𝒮) be a binary tree associated with 𝒮. If |𝒮| = *n* then the number
of bifurcations in *T*(𝒮) is *n* − 1 so the task is to figure out the
bifurcations of *T*(𝒮) (or more directly, figure out a set of sensible
bifurcations *B*<sub>*i*</sub> to make a tree with 𝒮 as tips or leaves).
In APE lingo, these bifurcations are called “splits”.

The information gain model of bifurcations/splits/partitions is as
follows: Let 𝒮 be a set of OTU’s and 𝒫 = *A* ⊔ *B* any partition.
Supposing that *A* are realizations of a random process *X* and *B* are
realizations of a random process *Y*, then the information of 𝒫 is
*I*(𝒫) := *I*(*X*, *Y*) where *I* is some meaningful extension of mutual
information to a random process.

Algorithm based on *I*(𝒫)
=========================

In the case where 𝒮 is aligned molecular sequence data, we assume that
*X* = *X*<sub>*n*</sub> is a random process where *X*<sub>*i*</sub>
takes values in a four-element set (i.e {*A*, *C*, *G*, *T*}) and
*A*<sub>*j*</sub>, the set of elements in the *j*th sites of *A* is the
set of realizations of *X*<sub>*n*</sub>. We define the mutual
information between *X* and *Y* to be
$\\frac{1}{n}\\sum\_{i=1}^{n}I(X\_n, Y\_n)$. Here
$$I(X\_n, Y\_n) = \\sum\_{(x,y) \\in (X,Y)}p(x,y)\\log\_{2}(\\frac{p(x,y)}{p(x)p(y)})$$
. The strategies for estimating these probabilities are:

-   $p(x,y) = \\frac{\\text{\# Of pairs (x,y) seen in the data } A\_n \\times B\_n}{|A\_n \\times B\_n|}$.
    Note that |*A*<sub>*n*</sub> × *B*<sub>*n*</sub>| = (*N* − *m*)*m*
    where *A*<sub>*n*</sub> = *m*, and each sequence has length *N*.

-   *p*(*x*) is the frequency of *x* in the given site in the data *A*
    (similarly for *Y*).

The optimal partition is a solution to

The algorithm `infopart` for making a tree *T*(𝒮) with input 𝒮 is as
follows:

1.  Compute the optimal partition of 𝒮 as 𝒫<sup>\*</sup> = *A* ⊔ *B*.
2.  Run `infopart` on *A* and *B* as input.
