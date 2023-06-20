
t_noise <- FALSE # whether we have Student-t noise
sig_trans <- FALSE # whether to sigmoid transform the data
cop_correct <- FALSE # whether to map the data to a Gaussian
exact <- TRUE # whether to use the exact test or Fisher z approximation
sparse_sim <- TRUE # whether to run a sparse simulation

### range of seeds
seed_numbers <- 100 + 1:100 # the seeds

### Settings
for (n2 in c(50, 100, 150, 200)) { # number of nodes
  n <- n2
  for (k in 1:3){
    N <- n*c(25, 50, 100)[k] # number of observations
    exp_parents <- 2 # expected number of parents
    nu <- NULL # Gaussian case
    sig_scale <- NULL # no transform
    min_ESS <- NULL # whether to restrict the dual tests
    if (t_noise) {
      nu <- c(10, 5, 2)[k] # degrees of freedom of t_Student noise, NULL is Gaussian
      N <- n*50
    }
    if (sig_trans) {
      sig_scale <- c(0, 0.25, 1)[k] # how much to transform the data, 0 is no transform
      N <- n*50
    }
    if (sparse_sim) {
      exp_parents <- c(0.05, 0.1, 0.2)[k]
      n <- n2*10
      N <- n*0.5
      min_ESS <- 20
    }
# store values for later dataframe
setup_vec <- c(n, N, exp_parents)
names(setup_vec) <- c("n", "N", "parents")
if(!is.null(nu)) {
  setup_vec <- c(setup_vec, nu)
  names(setup_vec)[4] <- "df"
  if(!is.null(cop_correct)) {
    setup_vec <- append(setup_vec, 1*cop_correct, after = 4)
    names(setup_vec)[5] <- "correct"
  }
}
if(!is.null(sig_scale)) {
  setup_vec <- c(setup_vec, sig_scale)
  names(setup_vec)[4] <- "sig"
  setup_vec <- c(setup_vec, 1*cop_correct)
  names(setup_vec)[5] <- "correct"
}
if(!is.null(min_ESS)) {
  setup_vec <- append(setup_vec, min_ESS, after = 3)
  names(setup_vec)[4] <- "min_ESS"
  setup_vec <- append(setup_vec, 1*exact, after = 4)
  names(setup_vec)[5] <- "exact"
}

# create a name for the directory to store stuff
subdir_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
dir_name <- paste("./sims", subdir_name, sep = "/")
dir_name2 <- paste("./sims_collated", subdir_name, sep = "/")

# create directory if none exists
if (!dir.exists("./sims_collated")) { 
  dir.create("./sims_collated")
}

#if (!file.exists(paste0(dir_name, ".Rdata"))) {

results_df <- NULL
times_df <- NULL

for (seed_number in seed_numbers) {
  seed_part <- paste("", "seed", seed_number, sep = "_")
    if (file.exists(paste0(dir_name, "/", subdir_name, seed_part, ".Rdata"))) {
      load(paste0(dir_name, "/", subdir_name, seed_part, ".Rdata"))
      results_df <- rbind(results_df, result_df)
      times_df <- rbind(times_df, time_df)
    }
}

save(results_df, times_df, file = paste0(dir_name2, ".Rdata"))

  }
}
