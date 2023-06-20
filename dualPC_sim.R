
local_run <- TRUE
# to run locally we need to set the seed and network size
# this is done externally on the cluster runs
write_data <- FALSE
run_own <- FALSE # whether to run our implementation of pc alg
run_prec <- FALSE # whether to run the precision pc
exact <- TRUE # whether to use the exact test or Fisher z approximation
add_ges <- FALSE # whether to add GES as a benchmark
t_noise <- FALSE # whether to add student t noise
sig_trans <- FALSE# whether to sigmoid transform the data
cop_correct <- FALSE # whether to map the data to a Gaussian
sparse_sim <- TRUE # whether to run a sparse simulation

if (local_run) {
  kk <- 4
  jj <- 2
  seed_number <- 101 # the seed
}

n <- c(50, 100, 150, 200)[kk] # number of nodes
N <- n*c(25, 50, 100)[jj] # number of observations
exp_parents <- 2 # expected number of parents
nu <- NULL # degrees of freedom of t-Student noise, NULL is Gaussian
sig_scale <- NULL # how much to transform the data, 0 is no transform
min_ESS <- NULL # whether to restrict the dual tests

if (t_noise) { 
  nu <- c(10, 5, 2)[jj]
  N <- n*50
}

if (sig_trans) {
  sig_scale <- c(0.1, 0.25, 1)[jj]
  N <- n*50
}

if (sparse_sim) {
  exp_parents <- c(0.05, 0.1, 0.2)[jj]
  n <- n*10
  N <- n*0.5
  min_ESS <- 20
}

# range of thresholds
pc_alphas <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25)
ges_lambdas <- c(5, 7, 10, 14, 20, 28, 40, 57, 80)
# how far to run each algorithm
# first for pcalg versions
# second for dual pc
# third for pc + precision
levelly <- c(9, 9, 9)

if (sparse_sim) {
  pc_alphas <- pc_alphas/1000
}

# load libraries
library(pcalg)
library(graph)
library(BiDAG)
source("DAGfns.R")
source("dualPC.R")

# store values for later dataframe
setup_vec <- c(n, N, exp_parents, seed_number)
names(setup_vec) <- c("n", "N", "parents", "seed")
if(!is.null(nu)) {
  setup_vec <- append(setup_vec, nu, after = 3)
  names(setup_vec)[4] <- "df"
  if(!is.null(cop_correct)) {
    setup_vec <- append(setup_vec, 1*cop_correct, after = 4)
    names(setup_vec)[5] <- "correct"
  }
}
if(!is.null(sig_scale)) {
  setup_vec <- append(setup_vec, sig_scale, after = 3)
  names(setup_vec)[4] <- "sig"
  setup_vec <- append(setup_vec, 1*cop_correct, after = 4)
  names(setup_vec)[5] <- "correct"
}
if(!is.null(min_ESS)) {
  setup_vec <- append(setup_vec, min_ESS, after = 3)
  names(setup_vec)[4] <- "min_ESS"
  setup_vec <- append(setup_vec, 1*exact, after = 4)
  names(setup_vec)[5] <- "exact"
}


# create a name for the directory to store stuff
f_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
subdir_name <- paste(paste(names(setup_vec[-length(setup_vec)]), 
                           setup_vec[-length(setup_vec)], sep = "_"), collapse = "_")
dir_name <- paste("./sims", subdir_name, sep = "/")


### Generate data
set.seed(seed_number) # set seed
# generate random DAG
trueDAGedges <- as(pcalg::randDAG(n = n, d = 2*exp_parents, 
                                    wFUN = list(runif, min=0.4, max=2)), "matrix")
trueDAG <- 1*(trueDAGedges != 0)
trueCPDAG <- BiDAG:::dagadj2cpadj(trueDAG)
trueskel <- 1*(trueDAG | t(trueDAG))
truepatt <- pdag2pattern(trueCPDAG)

set.seed(seed_number) # set seed
# generate simulated data
data <- rmvDAG(trueDAGedges, N, t_df = nu)

if (!is.null(sig_scale)) { # sigmoid transformation
  data <- Sig_data(data, sig_scale) 
}
if (cop_correct) { # transform to marginal Gaussians
  data <- Gau_data(data)
}

# create directory if none exists
if (!dir.exists("./sims")) { 
  dir.create("./sims")
}
if (!dir.exists(dir_name)) { 
  dir.create(dir_name)
}

if (write_data) {
  if (!file.exists(paste0(dir_name, "/", f_name, "_data.csv"))) {
    write.csv(data, paste0(dir_name, "/", f_name, "_data.csv"), 
                         row.names = FALSE)
  }
}

# to store method name and parameter values
method_vec <- rep(NA, 4)
names(method_vec) <- c("method", "parameter", "value", "comparison")

# fill out missing simulations
if (!file.exists(paste0(dir_name, "/", f_name, ".Rdata"))) {
  result_df <- NULL
  time_df <- NULL
  # correlation matrix  
  cor_mat <- cor(data)
  source("./dualPC_algs.R")
  
  save(result_df, time_df, file = paste0(dir_name, "/", f_name, ".Rdata"))
}
