
for (jj in 1:2) {

times_combined_df <- NULL

all_results_df <- NULL
all_times_df <- NULL

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
  
  results_df %>% group_by(n, parents, method, comparison, value) %>%
    summarise(TPR = mean(TPR), 
              FPRp = mean(FPRp)) -> summarised_results_df
  
  # keep one alpha level
  times_df %>% filter(method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi")) %>%
    filter(value == 5e-5) %>%
    mutate(method = recode_factor(method, dualPCst = "dual PC", 
                                  dualPCoi = "dual PC\nstable",
                                  PCst = "PC", 
                                  PCoi = "PC\nstable")) -> time_df
  
  # keep one alpha level
  results_df %>% filter(method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi")) %>%
    filter(value == 5e-5) %>%
    filter(comparison == "cpdag") %>%
    mutate(method = recode_factor(method, dualPCst = "dual PC", 
                                  dualPCoi = "dual PC\nstable",
                                  PCst = "PC", 
                                  PCoi = "PC\nstable"),
           SHDn = SHD/as.numeric(as.character(n))) -> result_df
  
  add_facet_text <- function (x) {
    paste("n =", x)
  }
  add_facet_text2 <- function (x) {
    paste0("d = ", x)
  }
  
  for (compy in c("cpdag", "pattern", "skeleton")){
    
    ggplot() + 
      geom_point(data = results_df %>% 
                   filter(comparison == compy, method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi"))
                 , aes(x = FPRp, y = TPR, color=method), 
                 alpha = 0.15, size = 0.5, shape = 20) + 
      geom_path(data = summarised_results_df %>% 
                  filter(comparison == compy, method %in% c("dualPCst", "dualPCoi", "PCst", "PCoi"))
                , aes(x = FPRp, y = TPR, colour = method),
                size = 1) + facet_wrap(~ n + parents, nrow = 2, labeller=labeller(parents = add_facet_text2, n = add_facet_text, .multi_line = FALSE)) + 
      coord_fixed(ratio = 1, xlim = c(0, 0.185), ylim = c(0.815, 1), expand = FALSE) + 
      scale_color_discrete(labels = c("\ndual PC\n", "\ndual PC\nstable\n", "\nPC\n", "\nPC\nstable\n"),
                           type = c("#6a3d9a","#1f78b4", "#e31a1c", "#ff7f00"))
    
    ggsave(paste0(compy, "_ROC_exact_", exact, ".png"), width=8.5, height=3.25, dpi = 200)
  }
  
  
  colrysall <- c("#8350ba", "#6a3d9a", "#512f75",
                 "#2b94db", "#1f78b4", "#185b88",
                 "#ea4648", "#e31a1c", "#b51516",
                 "#ff9933", "#ff7f00", "#cc6600")
  
  #labels <- paste0(rep(levels(time_df$method)[1:4], each=3), ", ", "N = ", times_df$k %>% unique, "n")
  labels <- paste0(rep(c("dual PC", "dual PC stable", "PC", "PC stable"), each=3), "\n", "d = ", times_df$parents %>% unique)
  
  ggplot(time_df, 
         aes(x = method, y = time, colour = interaction(parents, method), fill = interaction(parents, method))) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    geom_jitter(shape = 20, size = 0.5, alpha = 0.25, position=position_jitterdodge()) + scale_y_log10() +
    facet_wrap(~ n, labeller = labeller(n = add_facet_text), nrow = 1) + 
    scale_color_manual(values = colrysall, name= " ", labels = labels) + 
    scale_fill_manual(values = colrysall, name= " ", labels = labels) +
    theme(legend.key.height = unit(0.725, "cm")) + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(paste0("time_exact_", exact,".pdf"), width=12, height=3.75)  
  
  ggplot(result_df, 
         aes(x = method, y = SHDn, colour = interaction(parents, method), fill = interaction(parents, method))) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    geom_jitter(shape = 20, size = 0.5, alpha = 0.25, position=position_jitterdodge()) +
    facet_wrap(~ n, labeller = labeller(n = add_facet_text), nrow = 1) + 
    scale_color_manual(values = colrysall, name= " ", labels = labels) + 
    scale_fill_manual(values = colrysall, name= " ", labels = labels) +
    theme(legend.key.height = unit(0.725, "cm")) + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(paste0("SHD_exact_", exact,".pdf"), width=12, height=3.75)  
  
  }
  
