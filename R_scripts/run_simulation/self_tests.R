# Purpose: To run self-tests of OM and EM, but with different number of years for data i.e., 5, 10, and Terminal year
# following change in fleet structure
# Note: All EMs run here are 2 fleet models
# Creator: Matthew LH. Cheng
# Date 3/23/23

# Set up ------------------------------------------------------------------
rm(list=ls()) # remove objects prior to running

library(here)
library(tidyverse)
library(TMB)
library(doSNOW)
library(parallel)

ncores <- detectCores() 
cl <- makeCluster(ncores - 2)
registerDoSNOW(cl)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

compile_tmb(wd = here("src"), cpp = "EM.cpp")

# Read in OM and EM Scenarios
om_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios.xlsx"), sheet = "OM")
em_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios.xlsx"), sheet = "EM_Self_Tests")

# Read in spreadsheet for life history parameters
lh_path <- here("input", "Sablefish_Inputs.xlsx")

# Get number of OM and EM scenarios
n_OM_scen <- length(om_scenarios$OM_Scenarios)
n_EM_scen <- length(em_scenarios$EM_Scenario)

for(n_om in 1:n_OM_scen) {
  
  # Load in OM data
  om_path <- here("output", "OM_Scenarios", om_scenarios$OM_Scenarios[n_om])
  load(here(om_path, paste(om_scenarios$OM_Scenarios[n_om],".RData",sep = "")))
  list2env(oms,globalenv()) # output into global environment

  for(n_em in 1:n_EM_scen) {
    
      # Pre-processing for EM
      n_fleets <- em_scenarios$n_fleets[n_em] # Get number of fleets
      years <- 1:em_scenarios$n_years[n_em] # Get vector of years
      fish_selex_opt <- unlist(strsplit(em_scenarios$Selex[n_em], ",")) # get fishery selectivity options
      time_selex <- em_scenarios$Time_Selex[n_em] # Get time-varying selectivity options
      time_selex_npars <- em_scenarios$Time_Selex_Npars[n_em] # Number of time-varying selectivity parameters 
      if(time_selex_npars == "NA") time_selex_npars <- NULL # replace with NULL
      om_application <- unlist(strsplit(em_scenarios$OM_Application[n_em], ",")) # EMs to apply to particular OM cases
      
      # Only run simulations if OM corresponds with an EM
    if(sum(om_application %in% c(om_scenarios$OM_Scenarios[n_om])) >= 1) {
      
      sim_models <- foreach(sim = 1:n_sims, 
                            .packages = c("TMB", "here", "tidyverse")) %dopar% {

      compile_tmb(wd = here("src"), cpp = "EM.cpp")
      
      # Prepare inputs for EM
      input <- prepare_EM_input(years = years,
                                n_fleets = n_fleets, 
                                Fish_Start_yr = 1,
                                catch_cv = rep(0.01, n_fleets),
                                
                                # Selectivity Blocks
                                F_Slx_Blocks_Input = matrix(c(rep(0)),
                                                            nrow = length(years),
                                                            ncol = n_fleets), # fishery blocks
                                S_Slx_Blocks_Input = matrix(c(0), # selectivity blocks
                                                            nrow = length(years), 
                                                            ncol = 1),
                                use_fish_index = FALSE,
                                Fish_Comp_Like_Model = c("multinomial", "multinomial"),
                                Srv_Comp_Like_Model = c("multinomial"),
                                rec_model = "BH", 
                                F_Slx_Model_Input = fish_selex_opt,
                                S_Slx_Model_Input = c("logistic"), # logistic
                                time_selex = time_selex, 
                                n_time_selex_pars = time_selex_npars,
                                fix_pars = c( "ln_SigmaRec", "ln_q_fish", "ln_h", "ln_M"),
                                sim = sim)
      
      # Run EM model here and get sdrep
      model <- run_EM(data = input$data, parameters = input$parameters, 
                                      map = input$map, 
                                      n.newton = 3, 
                                      silent = TRUE, getsdrep = TRUE)
      
      # Check model convergence
      convergence_status <- check_model_convergence(mle_optim = model$mle_optim,
                                                    mod_rep = model$model_fxn,
                                                    sd_rep = model$sd_rep,
                                                    min_grad = 0.01)
      
      # Get quantities
      quants_df <- get_quants(sd_rep = model$sd_rep,
                              model_fxn = model$model_fxn, 
                              sim = sim, 
                              conv = convergence_status$Convergence,
                              n_years = length(years), 
                              fish_selex_opt = fish_selex_opt,
                              ntrue_fish_fleets = dim(Fish_Age_Comps)[3], 
                              nmod_fish_fleets = input$data$n_fleets,
                              ages = ages, 
                              F_x = 0.4)
      
      # Combine objects to save
      all_obj_list <- list(model, quants_df$Par_df, quants_df$TS_df)
      
      all_obj_list
    } # end foreach loop

    # After we're done running EMs, output objects to save in folder
    
    # Pre-processing here
    model_list <- list()
    params <- data.frame()
    time_series <- data.frame()
    
    for(s in 1:n_sims) {
      model_list[[s]] <- sim_models[[s]][[1]] # save model objects
      params <- rbind(params, sim_models[[s]][[2]]) # save parameters
      time_series <- rbind(time_series, sim_models[[s]][[3]]) # save time series
    } # end s loop
    
    # Now, save our results - create directory to store results first
    em_path_res <- here(om_path, em_scenarios$EM_Scenario[n_em])
    dir.create(em_path_res)
    save(model_list, file = here(em_path_res, paste(em_scenarios$EM_Scenario[n_em], ".RData", sep = ""))) # save models
    
    # Differentiate OM and EMs
    params <- params %>% dplyr::mutate(OM_Scenario = om_scenarios$OM_Scenarios[n_om],
                                       EM_Scenario = em_scenarios$EM_Scenario[n_em])
    time_series <- time_series %>% dplyr::mutate(OM_Scenario = om_scenarios$OM_Scenarios[n_om],
                                                 EM_Scenario = em_scenarios$EM_Scenario[n_em])
    
    # Now, output these as csvs
    write.csv(params, here(em_path_res, "Param_Results.csv"))
    write.csv(time_series, here(em_path_res, "TimeSeries_Results.csv"))
    cat(crayon::green("OM", n_om, "out of", n_OM_scen))
    cat(crayon::yellow("EM", n_em, "out of", n_EM_scen))

    } # end if statement for corresponding om and em
  } # end n_em loop
} # end n_om loop

stopCluster(cl) 


# 
ts_plot <- get_RE_precentiles(df = time_series %>%
                                filter(conv == "Converged", type == "Spawning Stock Biomass"),
                              est_val_col = 1, true_val_col = 5,
                              par_name = "Total Recruitment", group_vars = "year")

(est_plot <- plot_RE_ts_ggplot(data = ts_plot, x = year, y = median,
                               lwr_1 = lwr_80, upr_1 = upr_80,
                               lwr_2 = lwr_95, upr_2 = upr_95,
                               facet_name = par_name))
params %>%
  mutate(RE = (mle_val - t) / t) %>%
  filter(!is.na(t)) %>%
  ggplot(aes(y = RE, fill = type)) +
  geom_boxplot() +
  ylim(-0.3 , 0.3) +
  geom_hline(yintercept = 0, col = "red") +
  facet_wrap(~type, scales = "free") +
  theme_bw()

params %>%
  mutate(RE = (mle_val - t) / t) %>%
  group_by(type) %>%
  summarize(median = median(RE))
