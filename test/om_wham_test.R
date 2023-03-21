# Purpose: To test OM consistenecy using WHAM
# Date: 1/2/23
# Creator: Matthew LH. Cheng

library(here)
library(reshape2)
library(tidyverse)
library(wham)
library(TMB)
library(crayon)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
source(here(fxn_path, "simulate_data.R"))

# Path to general input biological parameters
# spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")
spreadsheet_path <- here("input", "EBS_Pollock_Inputs.xlsx")

# simulate data
# simulate data
simulate_data(fxn_path = fxn_path, 
              check_equil = FALSE,
              spreadsheet_path = spreadsheet_path, 
              rec_type = "BH",
              Start_F = c(0.01), 
              Fish_Start_yr = c(1), 
              Surv_Start_yr = c(1), 
              max_rel_F_M = c(1,1), 
              desc_rel_F_M = c(0.05), 
              F_type = c("Contrast"),
              yr_chng = c(25), 
              yr_chng_end = 30,
              fish_likelihood = "dirichlet_multinomial",
              srv_likelihood = "multinomial",
              DM_Fish_Param = 100,
              DM_Srv_Param = 100, 
              Input_Fish_N_Max = c(200), 
              Input_Srv_N_Max = c(200),
              fish_CV = c(0.1, 0.1),
              srv_CV = c(0.1), 
              catch_CV = c(0.01), 
              Input_N_Fish_Time = "F_Vary", 
              Input_N_Fish_Fixed = c(30),
              Mort_Time = "Constant", 
              q_Mean_Fish = c(0.05), 
              q_Mean_Surv = 0.01, 
              fish_selex = c("logistic"), 
              srv_selex = c("logistic"), 
              # if switching to a single sex, be sure to change the nrow to the number of sexes,
              # and to make sure the selex parameters for the fleets align n_pars * n_sexes
              # e.g., (7, 0.8, 4, 0.3) for a logistic with two sexes, nrow = 2
              fish_pars = list(Fleet_1_L = matrix(data = c(5, 0.85), 
                                                  nrow = 1, byrow = TRUE)), # fish fleet 2
              srv_pars = list(Fleet_3_SL = matrix(data = c(2, 0.85), 
                                                  nrow = 1, byrow = TRUE)), # survey fleet 1
              f_ratio = 1, m_ratio = 0)


plot_OM(path = here("figs", "Base_OM_Figs"), file_name = "OM_Check.pdf")

# WHAM checks ------------------------------------------------------------- 

# Fix M for now
M <- list(model = "constant", initial_means = Mean_M) # need to specify est_ages for M to be estimated!

# Set modelling structure here
selectivity <- list(Fleet_2 = list(model = rep("logistic", 3), re = rep("none", 3)),
                    Fleet_1 = list(model = rep("logistic", 2), re = rep("none", 2)),
                    Fleet_1_1TB = list(model = rep("logistic", 3), re = rep("none", 3)),
                    Fleet_1_6TB = list(model = rep("logistic", 7), re = rep("none", 7)),
                    Fleet_1_ar1_y = list(model = rep("logistic", 2 ), re = c("ar1_y", "none")),
                    Fleet_1_2dar1 = list(model = rep("logistic", 2), re = c("2dar1", "none")))


# Set catchability
catchability <- list(Fleet_2 = list(re = rep("none", 3), q_upper = rep(0.5, 3)),
                     Fleet_1 = list(re = rep("none", 1), q_upper = rep(0.5, 1)),
                     Fleet_1_1TB = list(re = rep("none", 1), q_upper = rep(0.5, 1)),
                     Fleet_1_6TB = list(re = rep("none", 1), q_upper = rep(0.5, 1)),
                     Fleet_1_ar1_y = list(re = rep("none", 1), q_upper = rep(0.5, 1)),
                     Fleet_1_2dar1 = list(re = rep("none", 1), q_upper = rep(0.5, 1))) 

# Specify single fleet + timeblock
single_fleet <- c(F, rep(T, 5))
time_block <- c(rep(F,2), rep(T,2), rep(F,2))

# Block period sel and index
block_period_sel <- list(NA, NA, list(16, NA),
                         list(c(5, 10, 15, 20, 26), NA),
                         NA, NA)

block_period_idx <- list(NA, NA, list(16, NA),
                         list(c(5, 10, 15, 20, 26), NA),
                         NA, NA)

# N fleets, indicies, and cv specifiation
n_fleets <- c(2, rep(1,5))
n_indices <- c(3, 2, rep(2,4))

# Catch cv
Catch_CV_Val <- list(c(0.01, 0.01), 0.01, 0.01, 0.01, 0.01, 0.01)
units_indices <- list(c(1,1,2), c(1,2), c(1,2), c(1,2), c(1,2), c(1,2))
units_index_paa <- list(c(2,2,2), c(2, 2,2), c(2,2), c(2,2), c(2,2), c(2,2))
noFish_Idx <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

# Save results
ssb_results <- data.frame()
f_results <- data.frame()
dir_par <- vector()

for(i in 2:length(selectivity)) {
  
  # Create objects to store stuff in
  conv_vec <- vector(length = n_sims)
  results_list <- list() # to store results
  sdrep_list <- list()
  wham_mod <- list()
  
  for(sim in 1:n_sims) {
    
    # set.seed(123)
    # Force our inputs into a list - so that it reads into wham
    basic <- make_wham_input(n_fleets = n_fleets[i], 
                             n_indices = n_indices[i], 
                             Catch_CV_Val = Catch_CV_Val[[i]], 
                             catch_error = FALSE,
                             n_sims = sim, bias_obs = TRUE, 
                             bias_process = TRUE,
                             units_indices = units_indices[[i]], 
                             units_index_paa = units_index_paa[[i]],
                             single_fleet = single_fleet[i], 
                             time_block = time_block[i], 
                             block_period_sel = block_period_sel[[i]], 
                             block_period_idx = block_period_sel[[i]], 
                             noFish_Idx = noFish_Idx[i]) 
    # When you use noFish Idx, leave the nFleets and nIndices as is
    
    # Make WHAM inputs here 
    test_wham <- wham::prepare_wham_input(basic_info = basic, 
                                          selectivity = selectivity[[i]], 
                                          recruit_model = 2, M = M, 
                                          catchability = catchability[[i]],
                                          NAA_re = NULL, age_comp = list(fleets = "dir-mult",
                                                                         indices = "multinomial"))
    
    # Fit WHAM here
    tryCatch( {
      em_fit <- wham::fit_wham(input = test_wham, do.fit = T, do.osa = F, do.retro = F,
                               save.sdrep = TRUE, do.check = F, MakeADFun.silent = T)
      
      wham_mod[[sim]] <- em_fit
      
      convergence_check <- check_convergence(em_fit, ret = T) 
      
      results_list[[sim]] <- get_results(em_fit)
      sdrep_list[[sim]] <- em_fit$sdrep
      dir_par[sim] <- em_fit$sdrep$par.fixed[names(em_fit$sdrep$par.fixed ) == "catch_paa_pars"] 
      
      # threp <- em_fit$report()
      # sapply(grep("nll",names(threp),value=T), function(x) sum(threp[[x]]))
      
      # Get convergence status
      if(convergence_check$convergence == 0 &  convergence_check$is_sdrep == TRUE & convergence_check$na_sdrep == FALSE) {
        conv_vec[sim] <- "Converged"
      } else {
        conv_vec[sim] <- "Not Converged"
        
      }
    } , error = function(error) {cat("ERROR :",conditionMessage(error), "\n")})
    
    cat(yellow("### Done with Simulation", sim, "###"))
    
  } # end loop of number of simulations we want to run
  
  
  # Get SSB -----------------------------------------------------------------
  
  # Save simulations into dataframe
  SSB_df <- compare_results(results_list = results_list, n_sims = sim, 
                            EM_variable = "SSB", OM_df = SSB, conv_vec = conv_vec)
  
  # Rename results and bind
  SSB_df <- SSB_df[[2]] %>%  mutate(Par = names(selectivity)[[i]],
                                    Converged = sum(conv_vec == "Converged", na.rm = TRUE))
  
  
  # Get Fs ------------------------------------------------------------------
  
  # Do the same for f
  F_df <- compare_results(results_list = results_list, n_sims = sim, 
                          EM_variable = "F", OM_df = fish_mort, conv_vec = conv_vec)
  
  F_df <- F_df[[2]] %>%  mutate(Par = names(selectivity)[[i]],
                                Converged = sum(sum(conv_vec == "Converged", na.rm = TRUE)))
  
  ssb_results <- rbind(ssb_results, SSB_df)
  f_results <- rbind(f_results, F_df)
  
  cat(red("### Done with Scenario", i, "###"))
  cat(red("### Number of converged models = ", sum(conv_vec == "Converged", na.rm = TRUE)))
  
  # Rudimentary Comparisons -------------------------------------------------
  
  ssb_plot <- ssb_results %>% mutate(Type = "SSB")
  f_plot <- f_results %>% mutate(Type = "F")
  all_res <- rbind(ssb_plot, f_plot)
  
  png(here("figs", "sanity_checks", "diff_F_diff_selex_rapid.png"), width = 1500, height = 800)
  
  # Relevel factors
  all_res <- all_res %>% 
    mutate(Par = factor(Par, levels = c("Fleet_2", "Fleet_1", "Fleet_1_1TB", "Fleet_1_6TB",
                                        "Fleet_1_ar1_y", "Fleet_1_2dar1")))
  
  # Plot!
  print(
    all_res %>% 
      ggplot(aes(x = Year, y = median)) +
      # geom_label(mapping = aes(x = 72, y = Inf, label = Converged), vjust = 2, size = 5) +
      # annotate("rect", fill = "blue", alpha = 0.2, xmin = 86, xmax = 100, ymin = -Inf, ymax = Inf) +
      geom_ribbon(aes(ymin = lwr_80, ymax = upr_80, group = Par), alpha = 0.5, fill = "grey4") +
      geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, group = Par), alpha = 0.3, fill = "grey4") +
      geom_point(shape = 21, colour = "black", fill = "white", size = 3.8, stroke = 1, alpha = 1) +
      geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 1, alpha = 0.85) +
      theme_bw() + 
      facet_grid(~Type, scales = "free") +
      # ylim(-0.3,0.3)+
      labs(x = "Year", y  ="Relative Error")+
      theme(strip.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 13, color = "black")) 
  )
  dev.off()
  
} # end i loop


# Other plots -------------------------------------------------------------

plot_OM(path = here("figs", "Base_OM_Figs"), file_name = "OM_Check.pdf")

# Output plots
# Create directory to ouptut plots to
wham_out <- here("figs", "wham_checks")
# dir.create(wham_out)
plot_wham_output(wham_mod[[1]], dir.main = wham_out)

