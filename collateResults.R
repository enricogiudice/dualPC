
### range of seeds
seed_numbers <- 100 + 1:100 # the seeds

### Settings
for (n in c(50, 100, 150, 200)) { # number of nodes
  for (k in c(25, 50, 100)){ # different simulation settings
    N <- n*k # number of observations
    exp_parents <- 2 # expected number of parents

# store values for later dataframe
setup_vec <- c(n, N, exp_parents)
names(setup_vec) <- c("n", "N", "parents")
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
