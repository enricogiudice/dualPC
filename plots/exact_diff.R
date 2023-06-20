
all_results_df <- NULL
all_times_df <- NULL

for (jj in 1:2) {

for (exp_parents in c(0.05, 0.1, 0.2)) { # expected number of parents
  
  for (n in c(50, 100, 150, 200)*10) {
    
      exact <- c(0,1)[jj] # number of observations
      min_ESS <- 20 # sparse case
      N <- n*0.5

      # store values for later dataframe
      setup_vec <- c(n, N, exp_parents)
      names(setup_vec) <- c("n", "N", "parents")
      if(!is.null(min_ESS)) {
        setup_vec <- append(setup_vec, min_ESS, after = 3)
        names(setup_vec)[4] <- "min_ESS"
        setup_vec <- append(setup_vec, 1*exact, after = 4)
        names(setup_vec)[5] <- "exact"
      }
      # create a name for the directory to store stuff
      subdir_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
      dir_name <- paste("../sims_collated", subdir_name, sep = "/")
      
      load(paste0(dir_name, ".Rdata"))
      
      all_results_df <- rbind(all_results_df, results_df)
      all_times_df <- rbind(all_times_df, times_df)
    }
  }
  
}

  results_df <- all_results_df
  
  results_df$TPR <- as.numeric(as.character(results_df$TPR))
  results_df$FPRp <- as.numeric(as.character(results_df$FPR_P))
  results_df$SHD <- as.numeric(as.character(results_df$SHD))
  results_df$value <- as.numeric(as.character(results_df$value))
  
  # note in the labels, st means standard and oi means order independent (which we later label stable!)
  
  results_df$method <- factor(results_df$method, 
                              levels = c("dualPCst", "dualPCoi", "ownPCst", "ownPCoi",
                                         "PCst", "PCoi", "precPCst", "precPCoi", "GES"))
  
  times_df <- all_times_df
  
  times_df$time <- as.numeric(as.character(times_df$user.self))
  
  times_df$method <- factor(times_df$method, 
                            levels = c("dualPCst", "dualPCoi", "ownPCst", "ownPCoi", 
                                       "PCst", "PCoi", "precPCst", "precPCoi", "GES"))
  
  library(tidyverse)
  
  # keep one alpha level
  results_df %>% filter(method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi")) %>%
    filter(value == 5e-5) %>%
    filter(comparison == "cpdag") %>%
    filter(method %in% c("dualPCst", "dualPCoi")) %>%
    mutate(method = recode_factor(method, dualPCst = "dual PC", 
                                  dualPCoi = "dual PC\nstable",
                                  PCst = "PC", 
                                  PCoi = "PC\nstable"),
           SHDn = SHD/as.numeric(as.character(n))) -> result_df
  
  result_df %>% group_by(n, parents, seed, method) %>% 
    summarise(SHDndiff = sum(2*(exact==1)*SHDn) - sum(SHDn)) -> result_diff
  
  # keep one alpha level
  times_df %>% filter(method %in% c("dualPCst", "dualPCoi")) %>%
    filter(value == 5e-5) %>%
    mutate(method = recode_factor(method, dualPCst = "dual PC", 
                                  dualPCoi = "dual PC\nstable",
                                  PCst = "PC", 
                                  PCoi = "PC\nstable")) -> time_df
  
  time_df %>% group_by(n, parents, seed, method) %>% 
    summarise(timediff = (sum(2*(exact==1)*time) - sum(time))/sum(time)) -> time_diff 
  
  add_facet_text <- function (x) {
    paste("n =", x)
  }
  add_facet_text2 <- function (x) {
    paste0("d = ", x)
  }
  
  colrysall <- c("#8350ba", "#6a3d9a", "#512f75",
                 "#2b94db", "#1f78b4", "#185b88")
  
  #labels <- paste0(rep(levels(time_df$method)[1:4], each=3), ", ", "N = ", times_df$k %>% unique, "n")
  labels <- paste0(rep(c("dual PC", "dual PC stable"), each=3), "\n", "d = ", result_diff$parents %>% unique)

  ggplot(result_diff, 
         aes(x = method, y = SHDndiff, colour = interaction(parents, method), fill = interaction(parents, method))) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    geom_jitter(shape = 20, size = 0.75, alpha = 0.25, position=position_jitterdodge()) +
    facet_wrap(~ n, labeller = labeller(n = add_facet_text), nrow = 1) + 
    scale_color_manual(values = colrysall, name= " ", labels = labels) + 
    scale_fill_manual(values = colrysall, name= " ", labels = labels) +
    theme(legend.key.height = unit(0.725, "cm")) +
    labs(y=expression(paste(Delta, " SHDn"))) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(paste0("SHD_exact_diff.pdf"), width=8, height=3.75)  
  
  ggplot(time_diff, 
         aes(x = method, y = timediff, colour = interaction(parents, method), fill = interaction(parents, method))) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    geom_jitter(shape = 20, size = 0.5, alpha = 0.25, position=position_jitterdodge()) +
    facet_wrap(~ n, labeller = labeller(n = add_facet_text), nrow = 1) + 
    scale_color_manual(values = colrysall, name= " ", labels = labels) + 
    scale_fill_manual(values = colrysall, name= " ", labels = labels) +
    theme(legend.key.height = unit(0.725, "cm")) + 
    labs(y=expression(paste(Delta, " time / ", Sigma, " time"))) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(paste0("time_exact_diff.pdf"), width=8, height=3.75)    
