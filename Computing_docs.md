Computing tools: Sequence generation and parallel computing
================

This is extended documentation for our file
[`diagnostics_and_testing.md`](https://github.com/ShashankSule/info_theoretic_phylo/blob/main/diagnostics_and_testing.md)
where we test and explore the divisive and agglomerative algorithms that
make trees from aligned nucleotide sequences. Obviously to test the
accuracy of either algorithm we must compare its result to the “correct”
tree–the tree that actually corresponds to the input sequence. For cases
such as the `coiii.nex` dataset or the dataset from Press et al, the
correct trees are already known. However, we would ideally like to test
our algorithms against several datasets of varying parameters (such as
sequence lengths and number of taxa) to get more detailed information
about how they perform. We work as follows: we make a random tree and
generate a sequence alignment that corresponds to it using
`evolverRandomTree` from PAML. Then given the sequence alignment we
compute the corresponding divisivie and agglomerative tree using
`infotree` and `agg_clustering` We then find the Robinson-Foulds
distances between these computed trees and the tree we used to generate
the sequence alignment. We also compile data on the ultrametricity and
addivity of these computed trees and on the relationship between the
branch lengths of the original tree and those of the computed ones. We
repeat this procedure a large number of times and compute summary
statistics and plots of the data we get.

There are obviously two “generative” steps in the process: 1) generating
the sequence alignment and 2) generating the tree. We accomplish the
latter on the CBCB cluster; but first we document how to do (1) and the
pitfalls we should be wary of.

# Sequence Generation

## Sampling error concerns

We must be very careful when using simulated alignments for evaluating
correctness because, as mentioned on p.186 of Yang 2006, a tree
reconstruction algorithm can go wrong because of sampling error in the
input sequence alignment. In other words, since the sequence alignment
is generated from a stochastic model (such as JC69) on the original
tree, the alignment only corresponds perfectly to the tree as the number
of sites goes to infinity. However, since we always use a finite number
of sites (modern computers don’t yet work with an infinite amount of
data), there is a chance that the simulated alignment corresponds to
more than just one tree. Consequently the reconstructed tree may not
match the original not because of a flaw in the algorithm but simply
because the sequence alignment has too few sites to actually reflect its
phylogeny. So before we even start testing our algorithms we need to
make sure that the data we give it has low enough (though never zero)
sampling error or noise.

There are also parameters other than site number that can affect
sampling error. For example, in the evolver method a large mutation rate
(given all the other parameters are fixed) could cause too many sites to
get substituted–resulting in a sequence alignment that looks almost
random. Conversely, the mutation rate could be too low so the
bifurcation between different groups of taxa would not be accurately
reflected in the simulated sequence alignment.

## Sequence Generation methods

We explored essentially three methods to simulate sequence from randomly
generated trees:

1.  The `ms()` + `seqgen()` approach outlined on p.5 of [Wei-Chen Chen’s
    `phyclust`
    manual](https://snoweye.github.io/phyclust/document/phyclust-guide.pdf)

2.  Using `simSeq()` in `phangorn` in R.

3.  Using `evolver` in PAML outlined in [Ziheng Yang’s
    manual](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf).

We decided to go with option (c) because of the ease with which random
tree generation and sequence generation are integrated in the
`evolverRandomTree` option. We also made some changes in the original
PAML code to write both the original tree and the generated sequence
files. To get this going from scratch you need the following:

1.  Install and set up PAML from the
    [manual](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf). The
    instructions are in Chapter 2.

2.  Open `evolver.c` in `/src` and make the following changes to also
    store the random tree generated in `evolverRandomTree`:

<!-- end list -->

1.  On line 890 change `fixtree = 1` to `fixtree = 0` in `void
    Simulate()`  
2.  On line 887 add `*fseqtree` between `FILE` and `*fin` in `void
    Simulate()`.
3.  After line 1121 add `fseqtree = gfopen("mctrees.nex", "w");`.
4.  On lines 1198-1200 change

<!-- end list -->

    fprintf(fseq,"\nbegin tree;\n   tree true_tree = [&U] "); 
              OutTreeN(fseq,1,1); fputs(";\n",fseq);
              fprintf(fseq,"end;\n\n");

to

    fprintf(fseqtree,"\nbegin trees;\n   tree true_tree = [&R] "); 
              OutTreeN(fseqtree,1,1); fputs(";\n",fseqtree);
              fprintf(fseqtree,"end;\n\n");

5.  On line 1257 add `fclose(fseqtree)`.

<!-- end list -->

3.  Open Terminal and change the directory to `src`. Then run `gcc -O3
    -o evolverRandomTree evolver.c tools.c` if you’re using `gcc` (or
    the corresponding command in `Readme.txt`).

Now the `evolverRandomTree` is ready. To make a random tree and generate
a sequence from it, run `...src/evolverRandomTree 5 #PATHTODATFILE`. A
sample `dat` file is given in `MCbaseRtree.dat`. In this case
`#PATHTODATFILE` will be `paml4.8/MCbaseRtree.dat`. We have also written
an R wrapper that `who_dat()` which writes a `.dat` file based on the
model parameters you give it. Here is a sample, where we generate 100
blocks of sequence data for 10 species where each block corresponds to a
tree randomly generated in the TN93 model.

``` r
setwd("/Users/shashanksule/Documents/info_theoretic_phylo/")
n <- 100
trees <- list(n)
sequences <- list(n)
seeds <- list(n)
for(i in 1:n){
  
  #Initialize the .dat file and store the seed
  seeds[[i]] <- who_dat(seqs = 16,
                sites = 1000, 
                model = 6, 
                parameters = "5 5", 
                gamma = "0.5 4", 
                mutation = 3.5,
                equilibrium = "0.25 0.28 0.34 0.13",
                spit_seed = TRUE)
  
  #Simulate the trees and the sequences
  system("paml4.8/src/evolverRandomTree 5 paml4.8/MCbaseRTree.dat") 
  trees[[i]] <- write.tree(read.nexus("mctrees.nex"))
  sequences[[i]] <- ReadCharacters("mc.nex")
}
```

# Computing the trees

Unfortunately, the divisive algorithm `infotree` is
![O(2^n)](https://latex.codecogs.com/png.latex?O%282%5En%29 "O(2^n)")
where ![n](https://latex.codecogs.com/png.latex?n "n") is the number of
OTU’s so the problem quickly becomes too much for our little laptops to
handle. But worry not, since we have been conveniently allocated some
space on the CBCB cluster\! If you `ssh` into
`#YOURUSERNAME@cbcbsub00.umiacs.umd.edu` with your password, you’ll have
entered the headnode of the cluster, ready to submit jobs to the most
powerful computing henchmen in the world (for those of you who take this
stuff too literally: we’re just kidding).

There are some excellent resources online on how to use the cluster and
submit jobs to it (including UMD-specific resources). Our favourites
(and the ones sufficient for this project) are:

1.  Junaid Merchant’s [introductory
    workshop](https://github.com/UMD-COMBINE/IntroToHPCs) on HPCs at
    UMD.

2.  [Using SLURM with
    R](https://uscbiostats.github.io/slurmr-workshop/): This one is
    especially helpful as it has a simple end-to-end example of a
    computing procedure that uses parallel computing.

3.  Passing variables between different scripts: [bash to
    R](https://www.r-bloggers.com/2015/02/bashr-howto-pass-parameters-from-bash-script-to-r/)
    and [referencing a bash
    variable](https://stackoverflow.com/questions/56227356/how-to-send-a-loop-on-several-nodes-with-slurm).

4.  [Computing with SLURM
    arrays](https://help.rc.ufl.edu/doc/SLURM_Job_Arrays)

The idea is to compute say a 100 divisive trees using the fastest
possible parallel procedure. There are a few things to keep in mind:

1.  One may think parallelizing using `parallel` in R for every command
    is the way to go. But that isn’t necessarily true; for example we
    may want to compute 10 highly intensive parallelizable tasks
    independently and have 10 cores available. But we shouldn’t
    necessarily submit each task to one core. It may be more efficient
    to do them in sequence but dedicate all 10 cores to each task at a
    time (or do two tasks at once, 5 cores to each task and so on). We
    face this sort of problem in `infotree` because it computes
    `max_info` over all possible partitions. Now we have a choice,
    either to compute over the partitions simultaneously while
    dedicating fewer cores to `max_info` or computing over the
    partitions sequentially while dedicating all our cores to
    `max_info`. We find that for low number of taxa (\<10) it is
    actually beneficial to go for the latter strategy.

2.  We can also choose to split a task over multiple nodes (like
    computing the information gain over batches of partitions) but this
    has to be done explicitly in a bash script or via `slurmR` which
    communicates also with the nodes available in the computer.

The workflow is as follows:

1.  Upload the data on the cluster in a directory `sequence_data`, an R
    script `treecomputer.R` that computes a divisive tree from a
    sequence alignment, and a directory `tree_data` that stores the
    tree. The R script should be able to take the file name o the
    sequence alignment and write its output (a divisive tree) to file.
2.  Write a bash script `process.sh` to execute the R script along with
    instructions in `SLURM` over exactly how to compute over the data.
    Note that this is a quick and dirty way to do this. There are many
    more elegant solutions possible since R scripts don’t actually need
    to be wrapped in bash unlike MATLAB. This script is where we pass
    the file name to `treecomputer.R`.
3.  Execute `sbatch --process.sh` and wait for everything to finish\!

So first we write the data produced via the evolver to the disk:

``` r
getwd()
tree_names <- paste(getwd(), "/TN93_Jul21/Trees/tree", 1:100, ".txt", sep = "")
seed_names <- paste(getwd(), "/TN93_Jul21/Seeds/seed", 1:100, ".txt", sep = "")
seqs_names <- paste(getwd(), "/TN93_Jul21/Sequences/sequence", 1:100, ".nex", sep = "")

mapply(function(x,y) write(x,file = y), trees, tree_names)
mapply(function(x,y) write(x, file = y), seeds, seed_names)

mapply(function(x,y) write.nexus.data(x, file = y), sequences, seqs_names)
```

Navigate to the `info_theoretic_phylo` directory in your terminal and
execute the following:

1.  SSH into your directory on the headnode by executing `ssh
    #username@cbcbsub00.umiacs.umd.edu` and create two directories:
    `sequence_data` and `tree_data`. You can do this by the command
    `mkdir #directoryname`.

2.  Enter `exit` and go back to the `info_theoretic_phylo` directory on
    the local machine. Secure copy the `utilities.R`, `treecomputer.R`
    and `process.sh` to your personal directory on the head node and the
    `/#trialname/Sequences` directory to the `.../sequence_data/`
    directory. You can do the former by executing `scp #filepath
    #username@cbcbsub00.umiacs.umd.edu:/cbcbhomes/#username` and the
    latter by executing `scp -r #directorypath #targetpath`

3.  SSH back into the headnode (step 1). Then execute `sbatch
    process.sh`.

Done\! You can check the status of the computation by executing `squeue`
on the head node.
