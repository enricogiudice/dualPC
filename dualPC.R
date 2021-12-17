# This function inverts a matrix and puts the inverse in correlation form
psolve <- function(M) {
  M1 <- solve(M) # should use pseudoinverse
  scale <- 1/sqrt(diag(M1))
  t(M1*scale)*scale 
}

# this function generates all subsets
combinations <- function (n, r, v = 1:n) {
  v0 <- vector(mode(v), 0)
  if (r == 0) 
    v0
  else if (r == 1) 
    matrix(v, n, 1)
  else if (r == n) 
    matrix(v, 1, n)
  else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), 
             Recall(n - 1, r, v[-1]))
}

# this function just generates the next subset
nextSubS <- function(SubS, max) {
  k <- length(SubS)
  if (SubS[k] == max) { # need to shift
    rmost_vec <- which(SubS != 1:k + max - k)
    if (length(rmost_vec) > 0) {
      change_point <- max(which(SubS != 1:k + max - k))
      SubS[change_point:k] <- SubS[change_point] + 1:(k - change_point + 1)
    } else {
      warning("no more sets!")
    }
  } else {
    SubS[k] <- SubS[k] + 1
  }
  SubS
}

# This function performs the Fisher z-test for conditional independence
Fztest <- function(x, Nmk, T_star) {
  T <- abs(sqrt(Nmk - 3) * 0.5 * pcalg:::log.q1pm(x))
  (T > T_star)
}

# This function finds the nodes connected to x
Find_Nbhd <- function(G, x, y) {
  if (G[x, y] == 1) { # if there is an edge to test
    edge_vec <- G[x, ]
    edge_vec[c(x, y)] <- 0 # don't want these
    which(edge_vec == 1)
  } else { # otherwise
    c()
  }
}

# Compute the conditional correlation coefficient
Compute_rho <- function(M) {
  chol_mat <- chol(M[-c(1:2), -c(1:2), drop = FALSE]) # Cholesky decomposition
  x_mat <- backsolve(chol_mat, M[-c(1:2), 1:2, drop = FALSE], transpose = TRUE)
  local_mat <- M[1:2, 1:2] - t(x_mat) %*% x_mat
  rho <- local_mat[1, 2]/sqrt(local_mat[1, 1]*local_mat[2, 2])
  rho
}

# Compare estimated cpdag, skeleton and pattern graph to truth
compare_results <- function(cpdag, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt) {
  method_vec[4] <- "cpdag"
  result_cpdag <- compareGs(cpdag, trueCPDAG)
  result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, result_cpdag))))
  method_vec[4] <- "skeleton"  # Compare skeleton
  skel <- Gskel(cpdag)
  result_skel <- compareGs(skel, trueskel)
  result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, result_skel))))
  method_vec[4] <- "pattern"  # Compare pattern graph
  patt <- pdag2pattern(cpdag)
  result_patt <- compareGs(patt, truepatt)
  result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, result_patt))))
  return(result_df)
}

# Convert a cpdag to a pattern graph
pdag2pattern <- function(G) {
  pattern_graph <- matrix(0, nrow(G), ncol(G))
  ind <- which(G == 1, arr.ind = TRUE)
  for (i in 1:nrow(ind)) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    allZ <- setdiff(which(G[ ,y] == 1), x)  # nodes directed towards y excluding x
    for (z in allZ) {
      if(G[y,x] == 0 && G[y,z] == 0 && G[x,z] == 0 && G[z,x] == 0) {
        pattern_graph[x,y] <- pattern_graph[z,y] <- 1  # save v-structure x -> y <- z
        G[x,y] <- G[z,y] <- 0  # delete v-structure from old pdag
      }
    }
  }
  pdag <- 1 * (G | t(G))
  return(pattern_graph + pdag)
}

# Orient v structures and apply orientation rules
orient_vstructures <- function(G, sepsets, pres_sepsets, solve.confl = FALSE, 
          orientCollider = TRUE, pattern_graph = FALSE) {

  orientConflictCollider <- function(pdag, x, y, z) {
    if (pdag[x, y] == 1) {
      pdag[y, x] <- 0
    }
    else {
      pdag[x, y] <- pdag[y, x] <- 2
    }
    if (pdag[z, y] == 1) {
      pdag[y, z] <- 0
    }
    else {
      pdag[z, y] <- pdag[y, z] <- 2
    }
    pdag
  }
  rule1 <- function(pdag, solve.confl = FALSE) {
    search.pdag <- pdag
    ind <- which(pdag == 1 & t(pdag) == 0, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- which(search.pdag[b, ] == 1 & search.pdag[, b] == 1 & search.pdag[a, ] == 0 & search.pdag[, a] == 0)
      if (length(isC) > 0) {
        for (ii in seq_along(isC)) {
          c <- isC[ii]
          if (!solve.confl | (pdag[b, c] == 1 & pdag[c, b] == 1)) {
              pdag[b, c] <- 1
              pdag[c, b] <- 0
            }
          else if (pdag[b, c] == 0 & pdag[c, b] == 1) {
              pdag[b, c] <- 2
              pdag[c, b] <- 2
          }
        }
      }
      if (!solve.confl) 
        search.pdag <- pdag
    }
    pdag
  }
  rule2 <- function(pdag, solve.confl = FALSE) {
    search.pdag <- pdag
    ind <- which(search.pdag == 1 & t(search.pdag) == 1, 
                 arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- which(search.pdag[a, ] == 1 & search.pdag[, 
                                                       a] == 0 & search.pdag[, b] == 1 & search.pdag[b, 
                                                       ] == 0)
      for (ii in seq_along(isC)) {
        c <- isC[ii]
        if (!solve.confl | (pdag[a, b] == 1 & pdag[b, 
                                                   a] == 1)) {
          pdag[a, b] <- 1
          pdag[b, a] <- 0
        }
        else if (pdag[a, b] == 0 & pdag[b, a] == 1) {
          pdag[a, b] <- 2
          pdag[b, a] <- 2
        }
      }
      if (!solve.confl) 
        search.pdag <- pdag
    }
    pdag
  }
  rule3 <- function(pdag, solve.confl = FALSE) {
    search.pdag <- pdag
    ind <- which(search.pdag == 1 & t(search.pdag) == 1, 
                 arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      c <- which(search.pdag[a, ] == 1 & search.pdag[, 
                                                     a] == 1 & search.pdag[, b] == 1 & search.pdag[b, 
                                                     ] == 0)
      if (length(c) >= 2) {
        cmb.C <- combn(c, 2)
        cC1 <- cmb.C[1, ]
        cC2 <- cmb.C[2, ]
        for (j in seq_along(cC1)) {
          c1 <- cC1[j]
          c2 <- cC2[j]
          if (search.pdag[c1, c2] == 0 && search.pdag[c2, 
                                                      c1] == 0) {
              if (!solve.confl | (pdag[a, b] == 1 & pdag[b, 
                                                         a] == 1)) {
                pdag[a, b] <- 1
                pdag[b, a] <- 0
                if (!solve.confl) 
                  search.pdag <- pdag
                break
              }
              else if (pdag[a, b] == 0 & pdag[b, a] == 
                       1) {
                pdag[a, b] <- pdag[b, a] <- 2
                break
              }
          }
        }
      }
    }
    pdag
  }
  if (sum(G) == 0) 
    return(G)
  p <- nrow(G)
  pdag <- G
  if (orientCollider) {
    ind <- which(G == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(G[y, ] == 1), x)  # nodes adjacent to y excluding x
      for (z in allZ) {
        if (G[x, z] == 0  && pres_sepsets[x,z] == T && !(y %in% sepsets[[x,z]]) ) { 
            if (!solve.confl) {
              pdag[x, y] <- pdag[z, y] <- 1
              pdag[y, x] <- pdag[y, z] <- 0
            }
            else {
              pdag <- orientConflictCollider(pdag, x, y, z)
            }
        }
      }
    }
  }
  repeat {
    old_pdag <- pdag
    if (pattern_graph == F) {
      pdag <- rule1(pdag, solve.confl = solve.confl)
      pdag <- rule2(pdag, solve.confl = solve.confl)
      pdag <- rule3(pdag, solve.confl = solve.confl)
    }
    if (all(pdag == old_pdag)) 
      break
  }
  pdag
}


dual_pc <- function(cor_mat, N, alpha, ord_ind = TRUE, skeleton = FALSE, pattern_graph = FALSE, max_ord = NULL) {
  n <- ncol(cor_mat) # number of variables
  N <- nrow(data) # number of observations
  if (is.null(max_ord)) { # the maximum subset size to test
    max_ord <- n
  }
  if (length(alpha) == 1) {
    T_c <- abs(qnorm(alpha/2)) # p-value cut-off in z-space for normal tests
    T_p <- T_c # p-value cut-off in z-space for dual tests
  } else { # in case we want different thresholds for normal and dual tests
    T_c <- abs(qnorm(alpha[1]/2))
    T_p <- abs(qnorm(alpha[2]/2))
  }
  c_mat <- cor_mat # local correlation matrix 
  ord <- 0 # size of conditioning sets, 0th level is correlation/precision matrix
  # T statistics on correlation space
  # keep edges which are significant, ie delete those which aren't
  Gc_0 <- Fztest(cor_mat*upper.tri(c_mat), N, T_c)
  # precision matrix
  p_mat <- psolve(c_mat)
  # T statistics on precision space
  Gp_0 <- Fztest(p_mat*upper.tri(p_mat), N-n, T_p)
  # keep track of sepsets of precision matrix
  pres_sepsets <- Gp_0 | t(Gp_0)
  # combined 
  G_0 <- (Gc_0 & Gp_0)*upper.tri(Gc_0)
  # we keep track of the skeleton in a symmetric matrix
  G_cur <- 1*(G_0 | t(G_0)) # current graph
  if (ord_ind) { # to track the edges we need to delete
    del_mat <- matrix(1, n, n)
  }
  if (skeleton == FALSE) { # so we work out directions
    sepsets <- as.list(rep(NA, n*n))
    dim(sepsets) <- c(n, n) # this will record separating sets
  }
  done_flag <- FALSE # track when we have performed all tests needed
  while (ord < max_ord && done_flag == FALSE) { # keep track of order
    done_flag <- TRUE
    ord <- ord + 1
    edges <- which(G_cur == 1, arr.ind = TRUE)
    for (ii in 1:nrow(edges)) {
      x <- edges[ii, 1]
      y <- edges[ii, 2]
      S <- Find_Nbhd(G_cur, x, y)
      nbhd_size <- length(S)
      if (ord <= nbhd_size) { # only need dual tests up to when ord is half nbhd_size
        # we could track which larger tests have already been performed and skip those
        c_mat <- cor_mat[c(x, y, S), c(x, y, S)] # local correlation matrix
        p_mat <- psolve(c_mat) # local precision matrix
        # test current subset S
        test_flag <- Fztest(p_mat[1, 2], N-nbhd_size, T_p)
        if (ord == nbhd_size) {
          n_subsets <- 0 # nothing else to test
        } else {
          # only need subsets up to half the size - the dual part takes care of the others!
          subset_size <- min(ord, nbhd_size - ord)
          n_subsets <- choose(nbhd_size, subset_size)
          if (ord < nbhd_size - 1) { # we need another round
            done_flag <- FALSE
          }
        }
        jj <- 0
        while (jj < n_subsets && test_flag == TRUE) {
          jj <- jj + 1 # loop
          if (jj == 1) {
            subset <- 1:subset_size
          } else {
            subset <- nextSubS(subset, nbhd_size)
          }
          cond_set <- subset + 2 # which rows/columns to condition on
          if (ord <= nbhd_size/2) {
            # normal test
            rho <- Compute_rho(c_mat[c(1, 2, cond_set), c(1, 2, cond_set)])
            test_flag <- Fztest(rho, N-ord, T_c)
          } else {
            # normal test, but use dual space to compute more efficiently
            rho <- Compute_rho(p_mat[c(1, 2, cond_set), c(1, 2, cond_set)])
            test_flag <- Fztest(rho, N-ord, T_c)
          }
          if (skeleton == FALSE && test_flag == FALSE) {  # record separating subset
            sepsets[[x, y]] <- sepsets[[y, x]] <- S[subset]
          }
          if (ord < nbhd_size/2 && test_flag == TRUE) {
            # dual test as well
            rho <- Compute_rho(p_mat[c(1, 2, cond_set), c(1, 2, cond_set)])
            test_flag <- Fztest(rho, N-nbhd_size+ord, T_p)
          }
          if (skeleton == FALSE && test_flag == FALSE) {  # record separating subset
            sepsets[[x, y]] <- sepsets[[y, x]] <- S[-subset]
          }
        }
        if (test_flag == FALSE) { # a test was rejected
          if (ord_ind) { # we only delete at the end of each main loop
            del_mat[x, y] <- 0
            del_mat[y, x] <- 0
          } else {
            G_cur[x, y] <- 0 # delete edges
            G_cur[y, x] <- 0
          }
          if (skeleton == FALSE && jj == 0) { # record separating subsets
            sepsets[[x, y]] <- sepsets[[y, x]] <- S
          }
        }
      }
    }
    if (ord_ind) { # we delete edges now instead
      G_cur <- 1*(G_cur & del_mat)
    }
  }
  if (skeleton) {
    G_cur
  } else {
    orient_vstructures(G_cur, sepsets, pres_sepsets, pattern_graph)
  }
  
}


own_pc <- function(cor_mat, N, alpha, ord_ind = TRUE, skeleton = FALSE, prec = FALSE, max_ord = NULL) {
  n <- ncol(cor_mat) # number of variables
  N <- nrow(data) # number of observations
  if (is.null(max_ord)) { # the maximum subset size to test
    max_ord <- n
  }
  T_c <- abs(qnorm(alpha/2)) # p-value cut-off in z-space for normal tests
  c_mat <- cor_mat # local correlation matrix 
  ord <- 0 # first level is correlation/precision matrix
  # T statistics on correlation space
  # keep edges which are significant, ie delete those which aren't
  Gc_0 <- Fztest(cor_mat*upper.tri(c_mat), N, T_c)
  G_cur <- 1*(Gc_0 | t(Gc_0)) # current graph
  pres_sepsets <- matrix(T, ncol = n, nrow = n)
  if (prec) { # also look at precision matrix
    # precision matrix
    p_mat <- psolve(c_mat)
    # T statistics on precision space
    Gp_0 <- Fztest(p_mat*upper.tri(p_mat), N-n, T_c)
    # keep track of sepsets of precision matrix
    pres_sepsets <- Gp_0 | t(Gp_0)
    # combined 
    G_0 <- (Gc_0 & Gp_0)*upper.tri(Gc_0)
    G_cur <- 1*(G_0 | t(G_0)) # current graph
  }
  if (ord_ind) { # to track the edges we need to delete
    del_mat <- matrix(1, n, n)
  }
  if (skeleton == FALSE) { # so we work out directions
    sepsets <- as.list(rep(NA, n*n))
    dim(sepsets) <- c(n, n) # this will record separating sets
  }
  done_flag <- FALSE
  while (ord < max_ord && done_flag == FALSE) { # keep track of order
    done_flag <- TRUE
    ord <- ord + 1
    edges <- which(G_cur == 1, arr.ind = TRUE)
    for (ii in 1:nrow(edges)) {
      x <- edges[ii, 1]
      y <- edges[ii, 2]
      S <- Find_Nbhd(G_cur, x, y)
      nbhd_size <- length(S)
      if (ord <= nbhd_size) {
        c_mat <- cor_mat[c(x, y, S), c(x, y, S)]
        test_flag <- FALSE
        n_subsets <- choose(nbhd_size, ord)
        # and we need another round
        done_flag <- FALSE
        jj <- 0
        while (jj < n_subsets && test_flag == FALSE) {
          jj <- jj + 1 # loop
          if (jj == 1) {
            subset <- 1:ord
          } else {
            subset <- nextSubS(subset, nbhd_size)
          }
          cond_set <- subset + 2 # which rows/columns to condition on
          # normal test
          rho <- Compute_rho(c_mat[c(1, 2, cond_set), c(1, 2, cond_set)])
          if (Fztest(rho, N-ord, T_c) == 0) { # delete edge
            test_flag <- TRUE
            if (ord_ind) { # we only delete at the end of each main loop
              del_mat[x, y] <- 0
              del_mat[y, x] <- 0
            } else {
              G_cur[x, y] <- 0 # delete edges
              G_cur[y, x] <- 0
            }
            if (skeleton == FALSE) { # record separating subsets
              sepsets[[x, y]] <- S[subset]
              sepsets[[y, x]] <- S[subset]
            }
          }
        }
      }
    }
    if (ord_ind) { # we delete edges now instead
      G_cur <- 1*(G_cur & del_mat)
    }
  }
  if (skeleton) {
    G_cur
  } else {
    orient_vstructures(G_cur, sepsets, pres_sepsets, pattern_graph = FALSE)
  }
  
}
