Information gain model for trees
================
Shashank Sule
14/06/2021

The model
=========

Let's write the simplest possible model. Let $\mathcal{S}$ be a set of OTU's and let $T(\mathcal{S})$ be a binary tree associated with $\mathcal{S}$. If $|\mathcal{S}| = n$ then the number of bifurcations in $T(\mathcal{S})$ is $n-1$ so the task is to figure out the bifurcations of $T(\mathcal{S})$ (or more directly, figure out a set of sensible bifurcations $B_i$ to make a tree with $\mathcal{S}$ as tips or leaves). In APE lingo, these bifurcations are called "splits".

The information gain model of bifurcations/splits/partitions is as follows: Let $\mathcal{S}$ be a set of OTU's and $\mathcal{P} = A\sqcup B$ any partition. Supposing that $A$ are realizations of a random process $X$ and $B$ are realizations of a random process $Y$, then the information of $\mathcal{P}$ is $I(\mathcal{P}):= I(X,Y)$ where $I$ is some meaningful extension of mutual information to a random process.

\# Algorithm based on $I(\mathcal{P})$ In the case where $\mathcal{S}$ is aligned molecular sequence data, we assume that $X = X_n$ is a random process where $X_i$ takes values in a four-element set (i.e $\{A,C,G,T\}$) and $A_j$, the set of elements in the $j$ th sites of $A$ is the set of realizations of $X_n$. We define the mutual information between $X$ and $Y$ to be $\frac{1}{n}\sum_{i=1}^{n}I(X_n, Y_n)$. Here $$I(X_n, Y_n) =\sum_{(x,y)\in (X,Y)}p(x,y)\log_{2}(\frac{p(x,y)}{p(x)p(y)})$$
</pre>
. The strategies for estimating these probabilities are:

-   <pre>$p(x,y)  =\frac{\text{# Of pairs (x,y) seen in the data } A_n\times B_n}{|A_n\times B_n|}$. Note that $|A_n\times B_n| = (N - m)m$     where $A_n = m$, and each sequence has length $N$.

-   <pre>$p(x)$     is the frequency of $x$     in the given site in the data $A$     (similarly for $Y$).

The optimal partition is a solution to
 $$\mathcal{P}^* =\text{argmax}_{\mathcal{P}}\,I(\mathcal{P})$$
</pre>
The algorithm `infopart` for making a tree $T(\mathcal{S})$ with input $\mathcal{S}$ is as follows:

1.  Compute the optimal partition of $\mathcal{S}$     as $\mathcal{P}^* = A\sqcup B$.
2.  Run `infopart` on $A$     and $B$     as input.
