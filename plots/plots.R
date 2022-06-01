t_noise <- FALSE # whether we have Student-t noise

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
  
  times_df <- all_times_df
  
  times_df$time <- as.numeric(as.character(times_df$user.self))
  
  times_df$method <- factor(times_df$method, 
                            levels = c("dualPCst", "dualPCoi", "ownPCst", "ownPCoi", 
                                       "PCst", "PCoi", "precPCst", "precPCoi", "GES"))
  
  
  library(tidyverse)
  
  results_df <- mutate(results_df, k = as.numeric(as.character(N))/as.numeric(as.character(n))) 
  
  if(!t_noise) {
  results_df %>% group_by(n, k, method, comparison, value) %>%
    summarise(TPR = mean(TPR), 
              FPRp = mean(FPRp)) -> summarised_results_df
  
  times_df <- mutate(times_df, k = as.numeric(as.character(N))/as.numeric(as.character(n))) 
  
  #times_df %>% group_by(n, k, method, value) %>% 
  #  summarise(time = mean(time)) %>% filter(value == 0.05) %>% view()
  
  # keep one alpha level
  times_df %>% filter(method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi")) %>%
    filter(value == 0.05) %>%
    mutate(method = recode_factor(method, dualPCst = "dual PC", 
                                  dualPCoi = "dual PC\nstable",
                                  PCst = "PC", 
                                  PCoi = "PC\nstable")) -> time_df
  
  # keep one alpha level
  results_df %>% filter(method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi")) %>%
    filter(value == 0.05) %>%
    mutate(method = recode_factor(method, dualPCst = "dual PC", 
                                  dualPCoi = "dual PC\nstable",
                                  PCst = "PC", 
                                  PCoi = "PC\nstable"),
           SHDn = SHD/as.numeric(as.character(n))) -> result_df
  
  add_facet_text <- function (x) {
    paste("n =", x)
  }
  add_facet_text2 <- function (x) {
    paste0("N = ", x, "n")
  }
  
  for (compy in c("cpdag", "pattern", "skeleton")){
    
    ggplot() + 
      geom_point(data = results_df %>% 
                   filter(comparison == compy, method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi", "GES"))
                 , aes(x = FPRp, y = TPR, color=method), 
                 alpha = 0.15, size = 0.5, shape = 20) + 
      geom_path(data = summarised_results_df %>% 
                  filter(comparison == compy, method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi", "GES"))
                , aes(x = FPRp, y = TPR, colour = method),
                size = 1) + facet_grid(k ~ n, labeller=labeller(k = add_facet_text2, n = add_facet_text)) + 
      coord_fixed(ratio = 1, xlim = c(0, 0.55), ylim = c(0.25, 1), expand = FALSE) + 
      scale_color_discrete(labels = c("\ndual PC\n", "\ndual PC\nstable\n", "\nPC\n", "\nPC\nstable\n", "\nGES\n"),
                           type = c("#6a3d9a","#1f78b4", "#e31a1c", "#ff7f00", "#27a538"))
    
    ggsave(paste0(compy, "_ROC_parents_", exp_parents, ".pdf"), width=8, height=7)
  }
  
  
  colrysall <- c("#8350ba", "#6a3d9a", "#512f75",
                 "#2b94db", "#1f78b4", "#185b88",
                 "#ea4648", "#e31a1c", "#b51516",
                 "#ff9933", "#ff7f00", "#cc6600")
  
  #labels <- paste0(rep(levels(time_df$method)[1:4], each=3), ", ", "N = ", times_df$k %>% unique, "n")
  labels <- paste0(rep(c("dual PC", "dual PC stable", "PC", "PC stable"), each=3), "\n", "N = ", times_df$k %>% unique, "n")
  
  ggplot(time_df, 
         aes(x = method, y = time, colour = interaction(k, method), fill = interaction(k, method))) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    geom_jitter(shape = 20, size = 0.5, alpha = 0.25, position=position_jitterdodge()) + scale_y_log10() +
    facet_wrap(~ n, labeller = labeller(n = add_facet_text), nrow = 1) + 
    scale_color_manual(values = colrysall, name= " ", labels = labels) + 
    scale_fill_manual(values = colrysall, name= " ", labels = labels) +
    theme(legend.key.height = unit(0.725, "cm")) + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(paste0("time_parents_", exp_parents,".pdf"), width=8, height=3.75)  
  
  ggplot(result_df, 
         aes(x = method, y = SHDn, colour = interaction(k, method), fill = interaction(k, method))) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    geom_jitter(shape = 20, size = 0.5, alpha = 0.25, position=position_jitterdodge()) +
    facet_wrap(~ n, labeller = labeller(n = add_facet_text), nrow = 1) + 
    scale_color_manual(values = colrysall, name= " ", labels = labels) + 
    scale_fill_manual(values = colrysall, name= " ", labels = labels) +
    theme(legend.key.height = unit(0.725, "cm")) + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(paste0("SHD_parents_", exp_parents,".pdf"), width=8, height=3.75)  
  
  }
  
  else {  # if we have t-Student noise
    results_df %>% group_by(n, k, method, comparison, value, df) %>%
      summarise(TPR = mean(TPR), 
                FPRp = mean(FPRp)) -> summarised_results_df
    
    n_lev <- 50  # Keep one n level
    ggplot() + 
      geom_point(data = results_df %>% 
                    filter(comparison == "cpdag", method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi", "GES"), 
                           n == n_lev), aes(x = FPRp, y = TPR, color=method), 
                  alpha = 0.12, size = 0.5, shape = 20) + 
      geom_path(data = summarised_results_df %>% 
                  filter(comparison == "cpdag", method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi", "GES"), 
                         n == n_lev), aes(x = FPRp, y = TPR, colour = method), size = 1) +
      facet_grid(~ df, labeller = label_bquote(cols = nu==.(as.numeric(as.character(df))))) +  
      xlim(0, 0.55) + ylim(0.25, 1) +
      scale_color_discrete(labels = c("\ndual PC\n", "\ndual PC\nstable\n", "\nPC\n", "\nPC\nstable\n", "\nGES\n"),
                            type = c("#6a3d9a","#1f78b4", "#e31a1c", "#ff7f00", "#27a538")) +
      theme(panel.background = element_rect(fill = "#f2f2f2")) 
    ggsave(paste0("t_noise_parents_", exp_parents, "_n_", n_lev, ".pdf"), width=7, height=2.5)  
    
  }
}
