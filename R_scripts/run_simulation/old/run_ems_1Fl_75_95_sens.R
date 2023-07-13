# Purpose: To run EMs that are 1 Fleet as a fast time-block. This is a sensitivity test
# to figure out why the performance of the low data quality is better in a fast LL shift relative
# to the low data quality, and why there are differences in the 75% and 95% shift.
# Creator: Matthew LH. Cheng
# Date: 7/3/23


# Set up ------------------------------------------------------------------
rm(list=ls()) # remove objects prior to running

library(here)
library(tidyverse)
library(TMB)
library(doSNOW)
library(parallel)

ncores <- detectCores() 
# Register cluster here
cl <- makeCluster(ncores - 2)
registerDoSNOW(cl)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

compile_tmb(wd = here("src"), cpp = "EM.cpp")

# Read in OM and EM Scenarios
om_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios_v2.xlsx"), sheet = "OM") %>% 
  filter(str_detect(OM_Scenarios, "Const"))
em_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios_v2.xlsx"), sheet = "EM_1Fl_TI_Blk") %>% 
  filter(str_detect(EM_Scenario, "LL_Blk"))

# Read in spreadsheet for life history parameters
lh_path <- here("input", "Sablefish_Inputs.xlsx")

# Get number of OM and EM scenarios
n_OM_scen <- length(om_scenarios$OM_Scenarios)
n_EM_scen <- length(em_scenarios$EM_Scenario)

for(n_om in 1:n_OM_scen) {
  
  # Load OM -----------------------------------------------------------------
  
  # Load in OM data
  om_path <- here("output", "OM_Scenarios", om_scenarios$OM_Scenarios[n_om])
  load(here(om_path, paste(om_scenarios$OM_Scenarios[n_om],".RData",sep = "")))
  list2env(oms,globalenv()) # output into global environment
  
  for(n_em in 1:n_EM_scen) {
    
    # Pre-process EM ----------------------------------------------------------
    
    # Pre-processing to specify EM here
    n_fleets <- em_scenarios$n_fleets[n_em] # Number of fleets to model
    
    # Specify year options - differs if it is Fast vs Slow OM Scenario
    years_opt <- as.numeric(unlist(strsplit(em_scenarios$n_years[n_em], ",")) )
    if(str_detect(om_scenarios$OM_Scenarios[n_om], "Fast")) years <- 1:years_opt[1]
    if(str_detect(om_scenarios$OM_Scenarios[n_om], "Slow")) years <- 1:years_opt[2]
    
    time_selex <- em_scenarios$Time_Selex[n_em] # Get time-varying selectivity options
    time_selex_npars <- em_scenarios$Time_Selex_Npars[n_em] # Number of time-varying selectivity parameters 
    if(time_selex_npars == "NA") time_selex_npars <- NULL # replace with NULL
    time_block_boolean <- em_scenarios$Time_Block[n_em] # whether or not time-blocking is present
    fish_selex_opt <- unlist(strsplit(em_scenarios$Selex[n_em], ",")) # get fishery selectivity options
    
    # Set up time-blocking structure here for fishery
    if(time_block_boolean == FALSE) { # No Time Blocking
      F_Slx_Blocks_Input <- matrix(c(rep(0)), nrow = length(years), ncol = n_fleets)
      
      # munge into matrix format here
      fish_selex_opt <- matrix(c(rep(fish_selex_opt, length(years))), 
                               nrow = length(years), ncol = n_fleets)
    } # end if no time-block
    
    # if we have time blocking
    if(time_block_boolean == TRUE) { 
      # set up blocking structure
      F_Slx_Blocks_Input <- matrix(c(rep(0, 25), rep(1, length(years) - 25)), 
                                   nrow = length(years), ncol = n_fleets)
      
      # set up matrix for specifying selex form at a given block
      fish_selex_opt <- matrix(c(rep(fish_selex_opt[1], 25), 
                                 rep(fish_selex_opt[2], length(years) - 25)), 
                               nrow = length(years), ncol = n_fleets)
    } # end if Time-block true
    
    # Run Simulation ----------------------------------------------------------
    
    sim_models <- foreach(sim = 1:n_sims,
                          .packages = c("TMB", "here", "tidyverse")) %dopar% {
                            
                            compile_tmb(wd = here("src"), cpp = "EM.cpp")
                            
                            # Prepare inputs for EM
                            input <- prepare_EM_input(years = years,
                                                      n_fleets = n_fleets, 
                                                      Fish_Start_yr = 1,
                                                      ages = ages,
                                                      catch_cv = rep(0.01, n_fleets),
                                                      
                                                      # Selectivity Blocks
                                                      F_Slx_Blocks_Input = F_Slx_Blocks_Input, # fishery blocks
                                                      S_Slx_Blocks_Input = matrix(c(0), # selectivity blocks
                                                                                  nrow = length(years), 
                                                                                  ncol = 1),
                                                      use_fish_index = FALSE,
                                                      Fish_Comp_Like_Model = c("multinomial"),
                                                      Srv_Comp_Like_Model = c("multinomial"),
                                                      rec_model = "BH", 
                                                      F_Slx_Model_Input = fish_selex_opt,
                                                      S_Slx_Model_Input = c("logistic"), # logistic
                                                      time_selex = time_selex, 
                                                      n_time_selex_pars = time_selex_npars,
                                                      fix_pars = c( "ln_SigmaRec", 
                                                                    "ln_q_fish", 
                                                                    "ln_h", "ln_M"), sim = sim)
                            
                            # Run EM model here and get sdrep
                            tryCatch(expr = model <- run_EM(data = input$data, parameters = input$parameters, 
                                                            map = input$map, n.newton = 3,
                                                            silent = TRUE, getsdrep = TRUE), error = function(e){e}) 
                            
                            # Check model convergence
                            tryCatch(expr = convergence_status <- check_model_convergence(mle_optim = model$mle_optim,
                                                                                          mod_rep = model$model_fxn,
                                                                                          sd_rep = model$sd_rep,
                                                                                          min_grad = 0.001), error = function(e){e}) 
                            
                            # Get quantities
                            tryCatch(expr = quants_df <- get_quants(sd_rep = model$sd_rep,
                                                                    model_fxn = model$model_fxn, 
                                                                    sim = sim,  
                                                                    n_sex = n_sex, 
                                                                    sex_ratio = 0.5,
                                                                    conv = convergence_status$Convergence,
                                                                    n_years = length(years), 
                                                                    fish_selex_opt = fish_selex_opt,
                                                                    ntrue_fish_fleets = dim(Fish_Age_Comps)[3], 
                                                                    nmod_fish_fleets = input$data$n_fleets,
                                                                    ages = ages, 
                                                                    F_x = 0.4), error = function(e){e}) 
                            
                            # Combine objects to save
                            tryCatch(expr = all_obj_list <- list(model, 
                                                                 quants_df$Par_df, 
                                                                 quants_df$TS_df), error = function(e){e}) 
                            
                          }# end foreach loop
    
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
    
    # Also, get AIC results
    AIC_df <- data.frame()
    # Get convergence summary here
    Convergence = params %>% group_by(sim) %>% summarize(conv = unique(conv))
    
    for(s in 1:n_sims) {
      # Get AIC here
      AIC_df_tmp <- data.frame(
        AIC = TMB_AIC(model_list[[s]]$mle_optim),
        OM_Scenario = om_scenarios$OM_Scenarios[n_om],
        EM_Scenario = em_scenarios$EM_Scenario[n_em],
        sim = s, conv = Convergence$conv[s]
      )
      AIC_df <- rbind(AIC_df, AIC_df_tmp)
    } # end s sim AIC loop
    
    # Now, output these as csvs
    write.csv(AIC_df, here(em_path_res, "AIC_Results.csv"))
    write.csv(params, here(em_path_res, "Param_Results.csv"))
    write.csv(time_series, here(em_path_res, "TimeSeries_Results.csv"))
    
    # Progress
    cat(crayon::green("OM", n_om, "out of", n_OM_scen))
    cat(crayon::yellow("EM", n_em, "out of", n_EM_scen))
    
  } # end n_em loop
} # end n_om loop

stopCluster(cl) # stop cluster


# Plot Results (Time Series) ------------------------------------------------------------

all_results <- get_results(om_scenario_path = here("output", "Sensitivity"))
param_df <- all_results$Parameter_Sum # parameter dataframe
ts_df <- all_results$TimeSeries_Sum # time series dataframe
aic_df <- all_results$AIC_Sum
unique_oms <- unique(param_df$OM_Scenario) # unique oms

# Munging on AIC dataframe
aic_df <- aic_df %>% 
  # Changing names
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
  ), 
  EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'),
  time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End",
                                           "Terminal")))

# Write out time series and aic dataframe
# data.table::fwrite(ts_df, here("output", "TimeSeries_Summary.csv"))
data.table::fwrite(aic_df, here("output", "Sensitivity", "AIC_Convergence_Summary.csv"))

# Parameter Results ------------------------------------------------------
# Create relative error and CV metrics for parameters
param_df <- param_df %>% mutate(RE = (mle_val - t) / t, # Relative Error
                                TE = mle_val - t, # Total Error
                                ARE = abs(mle_val - t)/t) %>%  # Absolute Relative Error
  dplyr::select(-X)

# Differentiate time components
param_df <- param_df %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int") ~ "Fleet Intersect", # Fleet Transition Intersects
  ),
  EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'), # remove preceeding letters
  time_comp = factor(time_comp, levels = c("Terminal", "Fleet Trans End",
                                           "Fleet Intersect"))) %>% 
  data.frame()

data.table::fwrite(param_df, here("output", "Sensitivity", "Parameter_Summary.csv"))

# Time Series Summary -----------------------------------------------------

# Get relative error of time series
ts_pars <- unique(ts_df$type)

# empty dataframes to store values in
ts_re_df <- data.frame() 
ts_te_df <- data.frame()
ts_are_df <- data.frame()

# Loop through to extract time series components
for(i in 1:length(ts_pars)) {
  
  # Relative Error
  ts_comp <- get_RE_precentiles(df = ts_df %>% 
                                  filter(type == ts_pars[i],
                                         conv == "Converged") %>% data.frame(),
                                est_val_col = 2, true_val_col = 6,
                                par_name = ts_pars[i], 
                                group_vars = c("year", "OM_Scenario", "EM_Scenario")) %>% 
    data.frame()
  
  # Total Error
  ts_te_comp <- get_TE_precentiles(df = ts_df %>% 
                                     filter(type == ts_pars[i], 
                                            conv == "Converged") %>% data.frame(),
                                   est_val_col = 2, true_val_col = 6,
                                   par_name = ts_pars[i], 
                                   group_vars = c("year", "OM_Scenario", "EM_Scenario")) %>% 
    data.frame()
  
  # Absolute Relative error
  ts_are_comp <- get_ARE_precentiles(df = ts_df %>% 
                                       filter(type == ts_pars[i],  
                                              conv == "Converged") %>% data.frame(),
                                     est_val_col = 2, true_val_col = 6,
                                     par_name = ts_pars[i], 
                                     group_vars = c("year", "OM_Scenario", "EM_Scenario")) %>% 
    data.frame()
  
  ts_re_df <- rbind(ts_re_df, ts_comp)
  ts_te_df <- rbind(ts_te_df, ts_te_comp)
  ts_are_df <- rbind(ts_are_df, ts_are_comp)
  
} # end i loop

# Residual munging for time series summary plots
ts_re_df <-  ts_re_df %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
  ),
  EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'))

ts_te_df <-  ts_te_df %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
  ),
  EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'))

ts_are_df <-  ts_are_df %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
  ),
  EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'))

data.table::fwrite(ts_re_df, here("output", "Sensitivity", "TimeSeries_RE.csv"))
data.table::fwrite(ts_te_df, here("output", "Sensitivity", "TimeSeries_TE.csv"))
data.table::fwrite(ts_are_df, here("output", "Sensitivity", "TimeSeries_ARE.csv"))



# Time Series Plots -------------------------------------------------------

# Read in time series RE from actual model runs and bind to sensitivity runs for comparison
ts_re_df_og <- read.csv(here("output", "TimeSeries_RE.csv")) %>% 
  filter(OM_Scenario %in% c("Fast_LL_75_High", "Fast_LL_95_High"),
         EM_Scenario %in% c("1Fl_LL_Blk")) %>% 
  select(-X)

ts_re_df <- rbind(ts_re_df_og, ts_re_df)

# Now relevel factor for organizing plot
ts_re_df <- ts_re_df %>% 
  mutate(time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")))

    # Now loop through each time series component and print the plot out
    print(
      ggplot(ts_re_df, aes(x = year, y = median))  +
        geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = time_comp, group = time_comp), alpha = 0.3) +
        geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
        geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
        facet_grid(OM_Scenario~par_name) +
        coord_cartesian(ylim = c(-0.5,0.5)) +
        scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
        scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
        labs(x = "Year", y = "RE", 
             fill = "Time", color = "Time") +
        theme_matt() +
        # theme(aspect.ratio=1)
        theme(legend.position = "top",
              title = element_text(size = 20),
              axis.text = element_text(size = 13), 
              strip.text = element_text(size = 13)) 
    )
    
    
    
param_df_og <- read.csv(here("output", "Parameter_Summary.csv")) %>% 
  filter(OM_Scenario %in% c("Fast_LL_75_High", "Fast_LL_95_High"),
         EM_Scenario %in% c("1Fl_LL_Blk", "1Fl_LL_Blk_3")) %>% 
  select(-X)
    
param_df %>% 
  rbind(param_df_og) %>% 
  filter(type %in% c("ABC")) %>% 
  group_by(OM_Scenario, type, time_comp, EM_Scenario) %>% 
  summarize(median = median(RE))
  