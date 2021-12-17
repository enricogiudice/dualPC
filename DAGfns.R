
### This function gives edges weights between the bounds
# with both positive and negative signs

wFUN <- function(m, lb, ub){ # function for edge weights
  runif(m, lb, ub)*sample(c(-1, 1), m, replace = TRUE)
}


### This function generates Gaussian data from a DAG
# following the topological order

rmvDAG <- function(trueDAGedges, N, standardise = TRUE) {
  trueDAG <- 1*(trueDAGedges != 0) # the edge presence in the DAG
  n <- ncol(trueDAG) # number of variables
  data <- matrix(0, nrow = N, ncol = n) # to store the simulated data
  top_order <- rev(BiDAG:::DAGtopartition(n, trueDAG)$permy) # go down order
  for (jj in top_order) {
    parents <- which(trueDAG[, jj] == 1) # find parents
    lp <- length(parents) # number of parents
    if (lp == 0) { # no parents
      data[, jj] <- 0
    } else if (lp == 1) { # one parent
      data[, jj] <- data[, parents]*trueDAGedges[parents, jj]
    } else { # more than one parent
      data[, jj] <- colSums(t(data[, parents])*trueDAGedges[parents, jj])
    }
    # add random noise
    data[, jj] <- data[, jj] + rnorm(N)
  }
  if(standardise) { # whether to standardise
    scale(data)
  } else {
    data
  }
}

### This function extracts the skeleton from a graph
Gskel <- function(incidence) {
  1*(incidence|t(incidence))
}

### This function compares an estimated graph to the true one

compareGs <- function (estG, trueG) {
  estSkel <- Gskel(estG) # estimated skeleton
  trueSkel <- Gskel(trueG) # true skeleton
  P <- sum(trueSkel)/2 # number of positives
  diffSkel <- estSkel - trueSkel
  extra_edges <- which(diffSkel > 0) # edges in estimated but not true EG
  FP <- length(extra_edges)/2 # count to FPs
  estG[extra_edges] <- 0 # remove them from further comparisons
  missing_edges <- which(diffSkel < 0) # edges in true but not estimated EG
  FN <- length(missing_edges)/2 # count to FNs
  trueG[missing_edges] <- 0 # remove them from further comparisons
  # modified graphs have the same skeletons, so now just need to count mismatches
  mismatches <- 1*(estG != trueG)
  wrong_order <- sum(Gskel(mismatches))/2 # number of wrongly oriented edges
  FP <- FP + wrong_order/2 # include half in FP
  FN <- FN + wrong_order/2 # and half in FN
  SHD <- FP + FN # shd is the sum of errors
  TP <- P - FN # true positives are without false negatives
  # TPR, FPR_P
  if (P == 0) { # true graph is empty
    if (FP >= 0) {
      TPR <- 0
      FPR_P <- 1
    } else {
      TPR <- 1
      FPR_P <- 0
    }
  } else { # true graph is non-empty
    TPR <- TP/P
    FPR_P <- FP/P
  }
  compGs <- c(TP, FP, SHD, TPR, FPR_P, P)
  names(compGs) <- c("TP","FP", "SHD", "TPR", "FPR_P", "P")
  return(compGs)
}
