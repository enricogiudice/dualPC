
t_noise <- TRUE

for (exp_parents in c(1.5, 2)) { # expected number of parents

  times_combined_df <- NULL
  
  all_results_df <- NULL
  all_times_df <- NULL
  
  for (n in c(50, 100, 150, 200)) {
    
    for (jj in 1:3) { 
      N <- n*c(25, 50, 100)[jj] # number of observations
      nu <- NULL # Gaussian case
      
      if (t_noise) {
        nu <- c(10, 5, 2)[jj] # degrees of freedom of t_Student noise, NULL is Gaussian
        N <- n*50
      }
      
      # store values for later dataframe
      setup_vec <- c(n, N, exp_parents)
      names(setup_vec) <- c("n", "N", "parents")
      if(!is.null(nu)) {
        setup_vec <- c(setup_vec, nu)
        names(setup_vec)[4] <- "df"
      }
      # create a name for the directory to store stuff
      subdir_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
      dir_name <- paste("../sims_collated", subdir_name, sep = "/")
      
      load(paste0(dir_name, ".Rdata"))
      
      all_results_df <- rbind(all_results_df, results_df)
      all_times_df <- rbind(all_times_df, times_df)
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

  library(tidyverse)
  
  results_df <- mutate(results_df, k = as.numeric(as.character(N))/as.numeric(as.character(n))) 
  
  results_df %>% group_by(n, df, method, comparison, value) %>%
    summarise(TPR = mean(TPR), 
              FPRp = mean(FPRp)) -> summarised_results_df
  
  compy <- "cpdag"
    
    ggplot() + 
      geom_point(data = results_df %>% 
                   filter(comparison == compy, method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi", "GES"))
                 , aes(x = FPRp, y = TPR, color=method), 
                 alpha = 0.15, size = 0.5, shape = 20) + 
      geom_path(data = summarised_results_df %>% 
                  filter(comparison == compy, method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi", "GES"))
                , aes(x = FPRp, y = TPR, colour = method),
                size = 1) + facet_wrap(~ n + df, nrow = 2, 
                     ## labeller=labeller(n = add_facet_text, df = add_facet_text2, .multi_line = FALSE)) + 
      labeller = label_bquote(cols = n==.(as.numeric(as.character(n)))*','~nu==.(as.numeric(as.character(df))))) + 
      coord_fixed(ratio = 1, xlim = c(0, 0.55), ylim = c(0.25, 1), expand = FALSE) +
      scale_x_continuous(breaks=seq(0, 0.4, 0.2)) +
      scale_color_discrete(labels = c("\ndual PC\n", "\ndual PC\nstable\n", "\nPC\n", "\nPC\nstable\n", "\nGES\n"),
                           type = c("#6a3d9a","#1f78b4", "#e31a1c", "#ff7f00", "#27a538")) + 
      theme(panel.background = element_rect(fill = "#f2f2f2"))
    ggsave(paste0("t_noise_parents_", exp_parents, "_all.png"), width=8.5, height=4, dpi = 200)

}
  
 