
exp_parents <- 2

n <- 50

  all_results_df <- NULL
  
 for (jj in 1:3) { 
   nu <- c(10, 5, 2)[jj] # degrees of freedom of t_Student noise, NULL is Gaussian
   N <- n*50
   
   for (cop_correct in 0:1) {
  
      # store values for later dataframe
      setup_vec <- c(n, N, exp_parents)
      names(setup_vec) <- c("n", "N", "parents")
      setup_vec <- c(setup_vec, nu)
      names(setup_vec)[4] <- "df"
      setup_vec <- append(setup_vec, 1*cop_correct, after = 4)
      names(setup_vec)[5] <- "correct"
      
      # create a name for the directory to store stuff
      subdir_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
      dir_name <- paste("../sims_collated", subdir_name, sep = "/")
      
      load(paste0(dir_name, ".Rdata"))
      
      all_results_df <- rbind(all_results_df, results_df)
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
  
  results_df$correct <- factor(results_df$correct, levels = 1:0, labels = c("Copula", "Uncorrected"))

  library(tidyverse)
  
  results_df <- mutate(results_df, k = as.numeric(as.character(N))/as.numeric(as.character(n))) 
  
  results_df %>% group_by(n, df, correct, method, comparison, value) %>%
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
                size = 1) + facet_grid(correct ~ df, 
                     ## labeller=labeller(n = add_facet_text, df = add_facet_text2, .multi_line = FALSE)) + 
      labeller = label_bquote(cols = nu==.(as.numeric(as.character(df))))) + 
      coord_fixed(ratio = 0.65, xlim = c(0, 0.55), ylim = c(0.25, 1), expand = FALSE) +
      scale_x_continuous(breaks=seq(0, 0.4, 0.2)) +
      scale_color_discrete(labels = c("\ndual PC\n", "\ndual PC\nstable\n", "\nPC\n", "\nPC\nstable\n", "\nGES\n"),
                           type = c("#6a3d9a","#1f78b4", "#e31a1c", "#ff7f00", "#27a538")) + 
      theme(panel.background = element_rect(fill = "#f2f2f2"))
    ggsave(paste0("t_copula_parents_", exp_parents, ".png"), width=7, height=3.75, dpi = 200)


  
 