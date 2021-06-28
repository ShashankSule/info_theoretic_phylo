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