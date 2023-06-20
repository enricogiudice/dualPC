## PC order independent (stable)
method_vec[1:2] <- c("PCoi", "alpha")
print("PC order independent")
for (alpha in pc_alphas[1:levelly[1]]) {
  method_vec[3] <- alpha
  # Start the clock!
  ptm <- proc.time()
  ## estimate cpdag directly
  pcpdag <- pc(suffStat = list(C = cor_mat, n = nrow(data)),
            indepTest = gaussCItest, ## indep.test: partial correlations
            alpha = alpha, labels = paste0("V", 1:ncol(data)),
            skel.method = "stable")
  # Stop the clock
  PCtime <- (proc.time() - ptm)[1]
  # Extract cpdag
  pc_cpdag <- 1*as(pcpdag@graph, "matrix")
  # Compare results to generating model
  time_df <- rbind(time_df, data.frame(t(c(setup_vec, method_vec[-4], PCtime))))
  result_df <- compare_results(pc_cpdag, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt)
}

## PC standard (original)
method_vec[1:2] <- c("PCst", "alpha")
print("PC standard")
for (alpha in pc_alphas[1:levelly[1]]) {
  method_vec[3] <- alpha
  # Start the clock!
  ptm <- proc.time()
  ## estimate cpdag directly
  pcpdag <- pc(suffStat = list(C = cor_mat, n = nrow(data)),
                      indepTest = gaussCItest, ## indep.test: partial correlations
                      alpha = alpha, labels = paste0("V", 1:ncol(data)),
                      skel.method = "original")
  # Stop the clock
  PCtime <- (proc.time() - ptm)[1]
  # Extract cpdag
  pc_cpdag <- 1*as(pcpdag@graph, "matrix")
  # Compare results to generating model
  time_df <- rbind(time_df, data.frame(t(c(setup_vec, method_vec[-4], PCtime))))
  result_df <- compare_results(pc_cpdag, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt)
}

# dual PC, order independent
method_vec[1:2] <- c("dualPCoi", "alpha")
print("dual PC order independent")
for (alpha in pc_alphas[1:levelly[2]]) {
  method_vec[3] <- alpha
  # Start the clock!
  ptm <- proc.time()
  dualpc_cpdag <- dual_pc(cor_mat, nrow(data), alpha, exact = exact, min_ESS = min_ESS)
  # Stop the clock
  dualPCtime <- (proc.time() - ptm)[1]
  # Compare results to generating model
  time_df <- rbind(time_df, data.frame(t(c(setup_vec, method_vec[-4], dualPCtime))))
  result_df <- compare_results(dualpc_cpdag, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt)
}

# dual PC, standard
method_vec[1:2] <- c("dualPCst", "alpha")
print("dual PC standard")
for (alpha in pc_alphas[1:levelly[2]]) {
  method_vec[3] <- alpha
  # Start the clock!
  ptm <- proc.time()
  dualpc_cpdag <- dual_pc(cor_mat, nrow(data), alpha, ord_ind = FALSE, exact = exact, min_ESS = min_ESS)
  # Stop the clock
  dualPCtime <- (proc.time() - ptm)[1]
  # Compare results to generating model
  time_df <- rbind(time_df, data.frame(t(c(setup_vec, method_vec[-4], dualPCtime))))
  result_df <- compare_results(dualpc_cpdag, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt)
}

if (add_ges) {
## GES
method_vec[1:2] <- c("GES", "lambda")
print("GES")
for (lambda_scale in ges_lambdas) {
  method_vec[3] <- lambda_scale 
  # Start the clock!
  ptm <- proc.time()
  ## estimate cpdag directly
  ges_score <- new("GaussL0penObsScore", data, lambda = lambda_scale*log(nrow(data)))
  ges_fit <- ges(ges_score)
  # Stop the clock
  GEStime <- (proc.time() - ptm)[1]
  # Extract cpdag
  gesCPDAG <- 1*as(ges_fit$essgraph, "matrix")
  # Compare results to generating model
  time_df <- rbind(time_df, data.frame(t(c(setup_vec, method_vec[-4], GEStime))))
  result_df <- compare_results(gesCPDAG, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt)
}
}

if (run_own) { # run my own implementation of the pc algorithm, which is slower
# own PC, order independent
method_vec[1:2] <- c("ownPCoi", "alpha")
print("own PC order independent")
for (alpha in pc_alphas[1:levelly[1]]) {
  method_vec[3] <- alpha
  # Start the clock!
  ptm <- proc.time()
  ownpc_cpdag <- own_pc(cor_mat, nrow(data), alpha, exact = exact)
  # Stop the clock
  ownPCtime <- (proc.time() - ptm)[1]
  # Compare results to generating model
  time_df <- rbind(time_df, data.frame(t(c(setup_vec, method_vec[-4], ownPCtime))))
  result_df <- compare_results(ownpc_cpdag, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt)
}

# own PC, standard
method_vec[1:2] <- c("ownPCst", "alpha")
print("own PC standard")
for (alpha in pc_alphas[1:levelly[1]]) {
  method_vec[3] <- alpha
  # Start the clock!
  ptm <- proc.time()
  ownpc_cpdag <- own_pc(cor_mat, nrow(data), alpha, ord_ind = FALSE, exact = exact)
  # Stop the clock
  ownPCtime <- (proc.time() - ptm)[1]
  # Compare results to generating model
  time_df <- rbind(time_df, data.frame(t(c(setup_vec, method_vec[-4], ownPCtime))))
  result_df <- compare_results(ownpc_cpdag, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt)
}
}

if (run_prec) { # run PC also using the precision matrix
# PC with pre matrix, order independent
method_vec[1:2] <- c("precPCoi", "alpha")
print("precision PC order independent")
for (alpha in pc_alphas[1:levelly[3]]) {
  method_vec[3] <- alpha
  # Start the clock!
  ptm <- proc.time()
  ownpc_cpdag <- own_pc(cor_mat, nrow(data), alpha, prec = TRUE)
  # Stop the clock
  ownPCtime <- (proc.time() - ptm)[1]
  # Compare results to generating model
  time_df <- rbind(time_df, data.frame(t(c(setup_vec, method_vec[-4], ownPCtime))))
  result_df <- compare_results(ownpc_cpdag, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt)
}

# PC with pre matrix, standard
method_vec[1:2] <- c("precPCst", "alpha")
print("precision PC standard")
for (alpha in pc_alphas[1:levelly[3]]) {
  method_vec[3] <- alpha
  # Start the clock!
  ptm <- proc.time()
  ownpc_cpdag <- own_pc(cor_mat, nrow(data), alpha, ord_ind = FALSE, prec = TRUE)
  # Stop the clock
  ownPCtime <- (proc.time() - ptm)[1]
  # Compare results to generating model
  time_df <- rbind(time_df, data.frame(t(c(setup_vec, method_vec[-4], ownPCtime))))
  result_df <- compare_results(ownpc_cpdag, setup_vec, method_vec, result_df, trueCPDAG, trueskel, truepatt)
}
}