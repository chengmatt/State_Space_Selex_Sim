# Purpose: To process and collate all model outputs
# Creator: Matthew LH. Cheng
# Date: 6/13/23

library(here)
library(tidyverse)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))
compile_tmb(wd = here("src"), cpp = "EM.cpp")

# Paths
om_scenario_path <- here("output", "OM_Scenarios") # path to OM folder
dir_out <- here("output", "Summary_Plots") # path to output folder
dir.create(dir_out)


# Collate Parameter Results -----------------------------------------------

# Get results from self tests
all_results <- get_results(om_scenario_path = om_scenario_path)
param_df <- all_results$Parameter_Sum # parameter dataframe
ts_df <- all_results$TimeSeries_Sum # time series dataframe
aic_df <- all_results$AIC_Sum
unique_oms <- unique(param_df$OM_Scenario) # unique oms

# Write out time series and aic dataframe
write.csv(ts_df, here("output", "TimeSeries_Summary.csv"))
write.csv(aic_df, here("output", "AIC_Convergence_Summary.csv"))

# Parameter Results ------------------------------------------------------

### ABC using the last 5 years of selex ------------------------------------------------

# Preallocate memory for lists
om_ref_pt_list <- vector("list", length(unique_oms))

for(i in 1:length(unique_oms)) {
  
    # List files in folder
    files <- list.files(here(om_scenario_path, unique_oms[i]))
    # Load in OM data
    load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
    # Remove.RData and .pdf
    files <- files[str_detect(files, ".RData|.pdf") == FALSE]
    
    # Pre allocate em list
    em_ref_pt_list <- vector("list", length(files))

  for(j in 1:length(files)) {
    
    # Load in .RData from model runs
    load(here(om_scenario_path, unique_oms[i], files[j], paste(files[j], ".RData", sep = "")))
    
    # Pre allocate memory for list
    ref_pt_list <- vector("list", length = length(model_list))
    
    # Get convergence statistics here
    convergence <- param_df %>% 
      filter(OM_Scenario == unique_oms[i], EM_Scenario == files[j], type %in% c("F_0.4", "ABC"))
    
    for(sim in 1:length(model_list)) {
      
      # Define to speed up process
      sim_model <- model_list[[sim]]
      
      # Define years based on em year lengths
      years <- 1:dim(sim_model$model_fxn$rep$FAA)[1]
      # Number of modelled fishery fleets
      nmod_fish_fleets <- dim(sim_model$model_fxn$rep$F_Slx)[3]
      # Get ages
      ages <- 1:dim(sim_model$model_fxn$rep$NAA)[2]
      # Get number of sexes
      n_sex <- dim(sim_model$model_fxn$rep$NAA)[3]
 
      # Pre-Processing first
      # Get Fs here
      if(nmod_fish_fleets > 1) {
        # Vector of fleets (define to multiply)
        fleet_num <- seq(1, nmod_fish_fleets)
        # Extract Fs here
        est_TermF <- exp(sim_model$sd_rep$par.fixed[names(sim_model$sd_rep$par.fixed) == "ln_Fy"])[c(length(years) * fleet_num)]
      } else{ # only 1 fleet
        # Extract Fs here
        est_TermF <- exp(sim_model$sd_rep$par.fixed[names(sim_model$sd_rep$par.fixed) == "ln_Fy"])[c(length(years))]
      } # else statement
      
      # Get true f40
      true_f40 <- get_Fx_refpt(ages = ages,
                               MortAA = oms$Mort_at_age[length(years),,sim], 
                               SelexAA = t(colMeans(oms$Fish_selex_at_age[(length(years) - 4):length(years),,,1,sim])),  # females selex
                               MatAA = oms$mat_at_age[length(years),,1,sim],  # females
                               WAA = oms$wt_at_age[length(years),,1,sim],  # females
                               Terminal_F = oms$fish_mort[length(years),,sim],
                               F_x = 0.4)
      
      # Get estimated selectivity here (last 5 years)
      est_Selex <- colMeans(sim_model$model_fxn$rep$F_Slx[years[(length(years) - 4):length(years)],,,1]) # selectivity only for females
      
      # Get estimated Fx% value
      Fx_val <- get_Fx_refpt(ages = ages,
                             MortAA = oms$Mort_at_age[length(years),,sim], 
                             SelexAA = t(est_Selex), 
                             MatAA = oms$mat_at_age[length(years),,1,sim], # females
                             WAA = oms$wt_at_age[length(years),,1,sim],  # females
                             Terminal_F = est_TermF, 
                             F_x = 0.4)
      # Get true abc
      true_ABC <- get_ABC_refpt(Fx_proj = true_f40, # F projection value
                                terminal_yr = length(years), # Terminal year of the assessment
                                sex_ratio = 0.5, # Assumed sex ratio
                                n_ages = length(ages),  # Number of ages
                                n_sex = n_sex,  # Number of sexes
                                n_fleets = 2,  # Number of fishery fleets
                                mean_rec = mean(oms$rec_total[1:length(years),sim], na.rm = TRUE), # mean recruitment
                                term_NAA = oms$N_at_age[length(years),,,sim], # terminal year NAA
                                term_F_Slx = array(colMeans(oms$Fish_selex_at_age[(length(years)-4):length(years),,,,sim]), 
                                                   dim = c(length(ages), 2, n_sex)), # last 5 years of selex averaged
                                term_F = oms$fish_mort[length(years),,sim], # terminal year F
                                MortAA = oms$Mort_at_age[length(years),,sim], # natural mortality
                                WAA = oms$wt_at_age[length(years),,,sim]  ) # weight at age
      
      # Get estimated abc
      ABC_val <- get_ABC_refpt(Fx_proj = Fx_val, # F projeciton value
                               terminal_yr = length(years), # Terminal year of the assessment
                               sex_ratio = 0.5, # Assumed sex ratio
                               n_ages = length(ages),  # Number of ages
                               n_sex = n_sex,  # Number of sexes
                               n_fleets = nmod_fish_fleets,  # Number of fishery fleets
                               mean_rec = mean(sim_model$sd_rep$value[names(sim_model$sd_rep$value) == "Total_Rec"]),  # mean recruitment
                               term_NAA = sim_model$model_fxn$rep$NAA[length(years),,], # terminal year NAA
                               term_F_Slx = array(colMeans(sim_model$model_fxn$rep$F_Slx[years[(length(years) - 4):length(years)],,,]), # terminal year selex
                                                  dim = c(length(ages), nmod_fish_fleets, n_sex)), # last 5 years of selex averaged
                               term_F =  matrix(exp(sim_model$sd_rep$par.fixed[names(sim_model$sd_rep$par.fixed) == "ln_Fy"]), # terminal year F
                                                ncol = nmod_fish_fleets, nrow = length(years))[length(years), ], 
                               MortAA = oms$Mort_at_age[length(years),,sim], # natural mortality at age
                               WAA = oms$wt_at_age[length(years),,,sim]  ) # weight at age
      
      # Create empty df for Fx% to rbind to par all dataframe
      F_df <- data.frame(matrix(ncol = length(names(param_df))))
      colnames(F_df) <- names(param_df) # replace colnames in empty df
      ABC_df <- data.frame(matrix(ncol = length(names(param_df))))
      colnames(ABC_df) <- names(param_df) # replace colnames in empty df
      
      # Get convergence status
      conv = convergence$conv[sim]
      
      # Now, input into empty dataframe
      F_df <- F_df %>% dplyr::mutate(mle_val = Fx_val, 
                                     type = paste("F", 0.4, "Last5", sep = "_"),
                                     t = true_f40,
                                     OM_Scenario = unique_oms[i],
                                     EM_Scenario = files[j])
      
      F_df$sim <- sim
      F_df$conv <- conv
      
      # Do the same for ABC estimates
      ABC_df <- ABC_df %>% dplyr::mutate(mle_val = ABC_val, 
                                         type = "ABC_Last5",
                                         t = true_ABC,
                                         OM_Scenario = unique_oms[i],
                                         EM_Scenario = files[j])
      
      ABC_df$sim <- sim
      ABC_df$conv <- conv
      
      
      ref_pt_list[[sim]] <- rbind(ABC_df, F_df)
      cat(crayon::green(sim))
      
    } # end simulation loop
    
    ref_pt_df <- data.table::rbindlist(ref_pt_list) # list to df
    em_ref_pt_list[[j]] <- ref_pt_df # put df into list
    cat(crayon::blue(j))
    
  } # end j loop (ems)
    
    # list to df
    om_ref_pt_df <- data.table::rbindlist(em_ref_pt_list)
    om_ref_pt_list[[i]] <- om_ref_pt_df
    cat(crayon::red(i))
    
} # end i loop (oms)


# Now, put this back into a regular old dataframe
all_ref_pt_df <- data.table::rbindlist(om_ref_pt_list)


# Create relative error and CV metrics for parameters
param_df <- rbind(all_ref_pt_df, param_df)
param_df <- param_df %>% mutate(RE = (mle_val - t) / t, # Relative Error
                                TE = mle_val - t, # Total Error
                                ARE = abs(mle_val - t)/t) %>%  # Absolute Relative Error
  dplyr::select(-X)

# Differentiate time components
param_df <- param_df %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
  ),
  time_comp = factor(time_comp, levels = c("Terminal", "Fleet Trans End",
                                           "Fleet Intersect"))) %>% 
  data.frame()

write.csv(param_df, here("output", "Parameter_Summary.csv"))

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

write.csv(ts_re_df, here("output", "TimeSeries_RE.csv"))
write.csv(ts_te_df, here("output", "TimeSeries_TE.csv"))
write.csv(ts_are_df, here("output", "TimeSeries_ARE.csv"))


# Get Selectivity estimates ---------------------------------------------------------
# Pre-allocate om list here
om_slx_list <- list() # om selectivity estimates
all_em_slx_list <- list() # em selectivity estimates

for(i in 1:length(unique_oms)) {
  
  # List files in folder
  files <- list.files(here(om_scenario_path, unique_oms[i]))
  # Load in OM data
  load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
  # Remove.RData and .pdf
  files <- files[str_detect(files, ".RData|.pdf") == FALSE]
  # Pre-allocate list to store all EMs
  em_slx_list <- list()
  
  for(j in 1:length(files)) {
    
    # Load in .RData from model runs
    load(here(om_scenario_path, unique_oms[i], files[j], paste(files[j], ".RData", sep = "")))
    
    # Pre-allocate list object here
    F_Slx_list <- list()
    
    # Get convergence statistics here
    convergence <- param_df %>% 
      filter(OM_Scenario == unique_oms[i], EM_Scenario == files[j], type == "F_0.4")
    
    # Get total biomass estimate from models here
    for(m in 1:length(model_list)) {
      
      # Get fishery selectivity from model
      F_Slx_df <- reshape2::melt(model_list[[m]]$model_fxn$rep$F_Slx)
      names(F_Slx_df) <- c("Year", "Age", "Fleet", "Sex", "Selex")
      
      # Construct our dataframe
      F_Slx_df_bind <- cbind(F_Slx_df, sim = m, 
                             conv = convergence$conv[m],
                             OM_Scenario = unique_oms[i],
                             EM_Scenario = files[j])
      
      
      # Now bind everything together
      F_Slx_list[[m]] <- F_Slx_df_bind
      
    } # end m 
    
    # Turn from list to dataframe
    F_Slx_mod_df <- data.table::rbindlist(F_Slx_list)
    em_slx_list[[j]] <- F_Slx_mod_df # put biomass dataframes from EMs into list
    print(j)
    
  } # end j
  
  # put selex dataframes from ems into all models list
  all_em_slx_list[[i]] <- data.table::rbindlist(em_slx_list) 
  # Get OM selectivity 
  om_slx_df <- reshape2::melt(oms$Fish_selex_at_age[,,,,1])
  names(om_slx_df) <- c("Year", "Age", "Fleet", "Sex", "Selex")
  om_slx_df$OM_Scenario <- unique_oms[i]
  om_slx_list[[i]] <- om_slx_df
  print(i)
  
} # end i

# Rbind to list all ems 
EM_fish_Slx <- data.table::rbindlist(all_em_slx_list)

# Rbind to list all oms
OM_fish_Slx <- data.table::rbindlist(om_slx_list)
# Now parse our the numbers rom this
OM_fish_Slx <- OM_fish_Slx %>% 
  mutate(Year = parse_number(paste(Year)),
         Age = parse_number(paste(Age)),
         Fleet = parse_number(paste(Fleet)),
         Sex = parse_number(paste(Sex)))

write.csv(EM_fish_Slx, here("output", "EM_Fish_Selex.csv"), row.names = FALSE)
write.csv(OM_fish_Slx, here("output", "OM_Fish_Selex.csv"), row.names = FALSE)


# Get composition fits ----------------------------------------------------
all_em_fish_comps_list <- list() # em pred composition list
om_fish_comps_list <- list() # om composition list

for(i in 1:length(unique_oms)) {
  
  # List files in folder
  files <- list.files(here(om_scenario_path, unique_oms[i]))
  # Load in OM data
  load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
  # Remove.RData and .pdf
  files <- files[str_detect(files, ".RData|.pdf") == FALSE]
  # Pre-allocate list to store all EMs
  em_fish_comps_list <- list()
  
  for(j in 1:length(files)) {
    
    # Load in .RData from model runs
    load(here(om_scenario_path, unique_oms[i], files[j], paste(files[j], ".RData", sep = "")))
    
    # Pre-allocate list object here
    Fish_comps_list <- list()
    
    # Get total biomass estimate from models here
    for(m in 1:length(model_list)) {
      
      # Get fishery selectivity from model
      Fish_Comps_df <- reshape2::melt(model_list[[m]]$model_fxn$rep$pred_fish_age_comps)
      names(Fish_Comps_df) <- c("Year", "Age", "Fleet", "Sex", "Prop")
      
      # Check convergence here
      convergence_status <- aic_df %>% 
          filter(sim == m, OM_Scenario == unique_oms[i], EM_Scenario == files[j])
        
      # Construct our dataframe
      F_Comps_df_bind <- cbind(Fish_Comps_df, sim = m, 
                               conv = convergence_status$conv,
                               OM_Scenario = unique_oms[i], EM_Scenario = files[j])
      
      
      # Now bind everything together
      Fish_comps_list[[m]] <- F_Comps_df_bind
      
    } # end m 
    
    # Turn from list to dataframe
    Fish_Comps_mod_df <- data.table::rbindlist(Fish_comps_list)
    em_fish_comps_list[[j]] <- Fish_Comps_mod_df # put biomass dataframes from EMs into list
    print(j)
    
  } # end j
  
  # put selex dataframes from ems into all models list
  all_em_fish_comps_list[[i]] <- data.table::rbindlist(em_fish_comps_list) 
  
  # OM Munging Comps and CAA ------------------------------------------------
  
  # Empty array to store stuff in for each OM 1 fleet combo
  obs_fish_age_comps <- array(data = 0, dim = c(dim(oms$N_at_age)[1] - 1,  # years
                                                dim(oms$N_at_age)[2], # ages
                                                1, # fleets
                                                dim(oms$Fish_Age_Comps)[4], # sexes
                                                dim(oms$N_at_age)[4])) # simulations
  
  # Get fishery age comps for 1 fleet model
  for(sim in 1:dim(oms$N_at_age)[4]) { # number of simulations loop 
    # Observed catches 
    obs_catches <- as.matrix(apply(as.matrix(oms$Catch_agg[1:50,,sim]), 1, FUN = sum), ncol = 1)
    # Catch weighting here
    catch_sum <- rowSums(as.matrix(oms$Catch_agg[1:50,,sim]))
    catch_weight <- as.matrix(oms$Catch_agg[1:50,,sim]) / catch_sum 
    
    for(s in 1:dim(oms$Fish_Age_Comps)[4]) { # sex loop
      for(f in 1:dim(oms$Fish_Age_Comps)[3]) { # fleet loop
        # Filter to save as an object
        fish_age_comps <-  oms$Fish_Age_Comps[1:(dim(oms$N_at_age)[1] - 1),,f,s,sim] * catch_weight[,f]
        # Increment comps - fixing fleet index to 1 here
        obs_fish_age_comps[,,1,s,sim] <- obs_fish_age_comps[,,1,s,sim] + fish_age_comps
      } # end f loop
    } # end s loop
    
    # Now, apply the proportion function over a single fleet
    for(s in 1:dim(oms$Fish_Age_Comps)[4]) {
      obs_fish_age_comps[,,,s, sim] <- t(apply(obs_fish_age_comps[,,,s,sim], MARGIN = 1, 
                                               FUN=function(x) { x/sum(x) }))
    } # s loop
  } # end sim loop
  
  # one fleet comps
  one_fleet_comps <- reshape2::melt(obs_fish_age_comps)
  names(one_fleet_comps) <- c("Year", "Ages", "Fleets", "Sexes", "Sim", "Prop")
  one_fleet_comps$fleet_type <- "One Fleet"
  
  # Empty array to store stuff in for each OM 2fleet combo
  two_obs_fish_age_comps <- array(data = 0, dim = c(dim(oms$N_at_age)[1] - 1,  # years
                                                    dim(oms$N_at_age)[2], # ages
                                                    dim(oms$Fish_Age_Comps)[3], # fleets
                                                    dim(oms$Fish_Age_Comps)[4], # sexes
                                                    dim(oms$N_at_age)[4])) # simulations
  
  # Get two fleet comps
  for(sim in 1:dim(oms$N_at_age)[4]) {
    for(s in 1:dim(oms$Fish_Age_Comps)[4]) { # sex loop
      for(f in 1:dim(oms$Fish_Age_Comps)[3]) { # fleet loop
        two_obs_fish_age_comps[,,f,s, sim] <- t(
          apply(oms$Fish_Age_Comps[1:(dim(oms$N_at_age)[1] - 1),,f,s,sim], MARGIN = 1, 
                FUN=function(x) { x/sum(x) })
        )
      } # end s loop
    } # end f loop
  } # end sim loop
  
  two_fleet_comps <- reshape2::melt(two_obs_fish_age_comps)
  names(two_fleet_comps) <- c("Year", "Ages", "Fleets", "Sexes", "Sim", "Prop")
  two_fleet_comps$fleet_type <- "Two Fleets"
  
  # Now, bind the 2 fleet and one fleet dataframe together
  fleet_comps <- rbind(one_fleet_comps, two_fleet_comps)
  fleet_comps$OM_Scenario <- unique_oms[i] # name OM
  om_fish_comps_list[[i]] <- fleet_comps # input into our om list
  print(i)
  
} # end i

# Rbind to list all ems 
EM_fish_comps <- data.table::rbindlist(all_em_fish_comps_list)
om_fish_comps_df <- data.table::rbindlist(om_fish_comps_list)

write.csv(EM_fish_comps, here("output", "EM_Fish_Comps.csv"), row.names = FALSE)
write.csv(om_fish_comps_df, here("output", "OM_Fish_Comps.csv"), row.names = FALSE)


