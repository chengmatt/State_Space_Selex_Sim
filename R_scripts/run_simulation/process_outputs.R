# Purpose: To pre-process and collate all model outputs
# Creator: Matthew LH. Cheng
# Date: 4/5/23

library(here)
library(tidyverse)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

# Paths
om_scenario_path <- here("output", "OM_Scenarios") # path to OM folder
dir_out <- here("output", "Summary_Plots") # path to output folder
dir.create(dir_out)


# Get Self Test Results ---------------------------------------------------

# Get results from self tests
self_test_res <- get_results(om_scenario_path = om_scenario_path)
param_df <- self_test_res$Parameter_Sum # parameter dataframe
ts_df <- self_test_res$TimeSeries_Sum # time series dataframe
unique_oms <- unique(param_df$OM_Scenario) # unique oms

# Create relative error and CV metrics for parameters
param_df <- param_df %>% mutate(RE = (mle_val - t) / t,
                                CV = (mle_var / mle_val) * 100) %>% 
  dplyr::select(-X)


# Differentiate time components
param_df <- param_df %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal",
    str_detect(EM_Scenario, "10_") ~ "10 Years Post Fl Chg",
    str_detect(EM_Scenario, "5_") ~ "5 Years Post Fl Chg"
  ),
  time_comp = factor(time_comp, levels = c("Terminal", "10 Years Post Fl Chg",
                                           "5 Years Post Fl Chg"))) %>% 
  data.frame()

# Get time series summaries -----------------------------------------------

# Get relative error of time series
ts_pars <- unique(ts_df$type)
ts_re_df <- data.frame() # empty dataframe to store
ts_te_df <- data.frame()

# Loop through to extract time series components
for(i in 1:length(ts_pars)) {
  # Get time series components here - relative error
  ts_comp <- get_RE_precentiles(df = ts_df %>% filter(type == ts_pars[i],
                                                      conv == "Converged") %>% data.frame(),
                                est_val_col = 2, true_val_col = 6,
                                par_name = ts_pars[i], group_vars = c("year", 
                                                                      "OM_Scenario",
                                                                      "EM_Scenario")) %>% 
    data.frame()
  
  # Get time series components here - total error
  ts_te_comp <- get_TE_precentiles(df = ts_df %>% filter(type == ts_pars[i],
                                                      conv == "Converged") %>% data.frame(),
                                est_val_col = 2, true_val_col = 6,
                                par_name = ts_pars[i], group_vars = c("year", 
                                                                      "OM_Scenario",
                                                                      "EM_Scenario")) %>% 
    data.frame()
  
  ts_re_df <- rbind(ts_re_df, ts_comp)
  ts_te_df <- rbind(ts_te_df, ts_te_comp)
  
} # end i loop

write.csv(ts_re_df, here("output", "TimeSeries_RE.csv"))
write.csv(ts_te_df, here("output", "TimeSeries_TE.csv"))


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
    
    # Get total biomass estimate from models here
    for(m in 1:length(model_list)) {
      
      # Get fishery selectivity from model
      F_Slx_df <- reshape2::melt(model_list[[m]]$model_fxn$rep$F_Slx)
      names(F_Slx_df) <- c("Year", "Age", "Fleet", "Sex", "Selex")

      # Check convergence here
      convergence_status <- check_model_convergence(mle_optim = model_list[[m]]$mle_optim,
                                                    mod_rep = model_list[[m]]$model_fxn,
                                                    sd_rep = model_list[[m]]$sd_rep,
                                                    min_grad = 0.01)
      # Construct our dataframe
      F_Slx_df_bind <- cbind(F_Slx_df, sim = m, 
                             conv = convergence_status$Convergence,
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
      convergence_status <- check_model_convergence(mle_optim = model_list[[m]]$mle_optim,
                                                    mod_rep = model_list[[m]]$model_fxn,
                                                    sd_rep = model_list[[m]]$sd_rep,
                                                    min_grad = 0.01)
      # Construct our dataframe
      F_Comps_df_bind <- cbind(Fish_Comps_df, sim = m, 
                             conv = convergence_status$Convergence,
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
        fish_age_comps <-  oms$Fish_Age_Comps[1:50,,f,s,sim] * catch_weight[,f]
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
          apply(oms$Fish_Age_Comps[1:50,,f,s,sim], MARGIN = 1, 
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


# ABC Calculations --------------------------------------------------------

# empty dataframe to store all om abc values
# all_om_abc_df <- data.frame()
# 
# # Get ABC values
# for(i in 1:length(unique(ts_df$OM_Scenario))) {
#   
#   # Quick set up here
#   # List files in folder
#   files <- list.files(here(om_scenario_path, unique_oms[i]))
#   # Load in OM data
#   load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
#   # Remove.RData and .pdf
#   files <- files[str_detect(files, ".RData|.pdf") == FALSE]
#   
#   # OM ABC calculations (Retroactive) --------------------------------------------------------------
#   
#   # empty dataframe to store values
#   om_abc_df <- data.frame(Terminal = NA, `10_Years_Post_Fl_Chg` = NA, 
#                           `5_Years_Post_Fl_Chg` = NA,
#                           sim = NA, OM_Scenario = NA)
#   
#   for(j in 1:oms$n_sims) {
#     
#     term_yrs <- c((dim(oms$N_at_age)[1] -1), 35, 30) # Get terminal years
#     n_sex <- dim(oms$N_at_age)[3]  # Get n_sexs
#     n_ages <- dim(oms$N_at_age)[2] # Get number of ages
#     
#     # Empty matrices to store projections in
#     N_proj <- array(data = NA, dim = c(n_ages, n_sex, length(term_yrs)))
#     F40_proj <- array(data = NA, dim = c(n_ages, n_sex, length(term_yrs)))
#     Z_proj <- array(data = NA, dim = c(n_ages, n_sex, length(term_yrs)))
#     CAA_proj <- array(data = NA, dim = c(n_ages, n_sex, length(term_yrs)))
#     Catch_proj <- array(data = NA, dim = c(n_ages, n_sex, length(term_yrs)))
#     
#     # Get FABC projection values here - these should be the same for some cases (i.e., fratio is constant) 
#     # but can change if the fleet structure change is slow.
#     F40 <- param_df %>% # extract f40 first
#       filter(OM_Scenario == unique(param_df$OM_Scenario)[i],
#              sim == j, type == "F_0.4") %>% 
#       dplyr::select(t, time_comp) %>% unique() %>% 
#       data.frame()# get F40 estimated for the OM previously
#     
#     # Sort to go from terminal -> 10 yrs -> 15 yrs
#     F40 <- F40[match(c("Terminal", "10 Years Post Fl Chg", "5 Years Post Fl Chg"), F40$time_comp, ),]
#     
#     for(t in 1:length(term_yrs)) { # time component of EMs (x years post fl change)
#       
#       # Extract mean recruitment over modelled period - multiply by sex ratio
#       mean_rec <- mean(oms$rec_total[1:term_yrs[t],j], na.rm = TRUE) * 0.5
#       fratio <- oms$fish_mort[term_yrs[t],,j] / sum(oms$fish_mort[term_yrs[t],,j]) # Get fratios
#       N_proj[1,,t] <- mean_rec # Input mean recruitment into projection for age-1
#       
#       for(s in 1:n_sex) {
#         for(a in 1:n_ages) {
#           
#           # Calculate age, fleet, and sex specific mortality (returns vector fo fleet specific mortalities)
#           fleet_mort <- sum(oms$fish_mort[term_yrs[t],,j] * oms$Fish_selex_at_age[term_yrs[t],a,,s,j], na.rm = TRUE)
#           
#           if(a < n_ages) { # not plus group
#             # Project forward with terminal year survival
#             N_proj[a + 1,s, t] <- oms$N_at_age[term_yrs[t],a,s,j] * 
#               exp(-(oms$Mort_at_age[term_yrs[t], a, j] + fleet_mort))
#           } else{
#             N_proj[a,s,t] <- N_proj[a,s,t] + (oms$N_at_age[term_yrs[t],a,s,j] * 
#                                                 exp(-(oms$Mort_at_age[term_yrs[t], a, j] + fleet_mort)))
#           } # else for plus group
#           
#           # F40 projections for each age here
#           F40_proj[a,s,t] <- (fratio[1] * F40$t[t] * oms$Fish_selex_at_age[term_yrs[t],a,1,s,j]) + # FABC from first fleet
#             (fratio[2] * F40$t[t] * oms$Fish_selex_at_age[term_yrs[t],a,2,s,j]) # FABC from second fleet
#           
#           # Mortality and catch projections
#           Z_proj[a,s,t] <- F40_proj[a,s,t] + oms$Mort_at_age[term_yrs[t], a, j] # Total instantaneous mortality
#           # Now, get catch projections
#           CAA_proj[a,s,t] <- N_proj[a,s,t] * (1 - exp(-Z_proj[a,s,t])) * (F40_proj[a,s,t] / Z_proj[a,s,t]) 
#           Catch_proj[a,s,t] <- CAA_proj[a,s,t] * oms$wt_at_age[term_yrs[t], a, s, j] # Turn to biomass units
#           
#         } # end a loop
#       } # end sex loop
#     } # end t loop (terminal years of different assessment combinations)
#     
#     # Now, get ABC estimates for the three different time periods
#     OM_ABC <- apply(Catch_proj, MARGIN = 3, sum)
#     # munge into dataframe for storage purposes
#     om_abc_df[j,1:3] <- OM_ABC
#     om_abc_df$sim[j] <- j
#     om_abc_df$OM_Scenario[j] <- unique(ts_df$OM_Scenario)[i]
#     print(j)
#     
#   } # end j loop
#   
#   # now, bind all values
#   all_om_abc_df <- rbind(all_om_abc_df, om_abc_df)
#   print(i)
#   
# } # end i loop
# 
# # Quick munging for this dataframe
# names(all_om_abc_df) <- c("Terminal", "10 Years Post Fl Chg", "5 Years Post Fl Chg", "sim", "OM_Scenario")
# 
# # pivot longer
# all_om_abc_df <- all_om_abc_df %>% pivot_longer(names_to = "time_comp",
#                                                 values_to = "t", 
#                                                 cols = c("Terminal", 
#                                                          "10 Years Post Fl Chg", 
#                                                          "5 Years Post Fl Chg"))
# 
# # EM ABC Calculations (Retroactive) -----------------------------------------------------
# 
# # empty dataframe to store all em abc values
# all_em_abc_df <- data.frame()
# 
# # Now, we need to grab ABC calculations from each respective OM EM combination...
# for(i in 1:length(unique(ts_df$OM_Scenario))) {
#   
#   # List files in folder
#   files <- list.files(here(om_scenario_path, unique_oms[i]))
#   # Load in OM data
#   load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
#   # Remove.RData and .pdf
#   files <- files[str_detect(files, ".RData|.pdf") == FALSE]
#   
#   for(j in 1:length(files)) {
#     
#     # Load in .RData from model runs
#     load(here(om_scenario_path, unique_oms[i], files[j], 
#               paste(files[j], ".RData", sep = "")))
#     
#     # Get f40 estimated
#     F40 <- param_df %>% # extract f40 first
#       filter(OM_Scenario == unique(ts_df$OM_Scenario)[i],
#              EM_Scenario == files[j], type == "F_0.4") %>% 
#       data.frame()
#     
#     # Empty dataframe for storage
#     em_abc_df <- setNames(data.frame(matrix(ncol = length(names(param_df)) - 1, 
#                                             nrow = 0)), names(param_df)[names(param_df) != "t"])
#     
#     for(sim in 1:length(model_list)) {
#       
#       # Empty matrices to store projections in
#       N_proj <- array(data = NA, dim = c(n_ages, n_sex))
#       F40_proj <- array(data = NA, dim = c(n_ages, n_sex))
#       Z_proj <- array(data = NA, dim = c(n_ages, n_sex))
#       CAA_proj <- array(data = NA, dim = c(n_ages, n_sex))
#       Catch_proj <- array(data = NA, dim = c(n_ages, n_sex))
#       
#       # Get terminal year of model
#       term_yr <- dim(model_list[[sim]]$model_fxn$rep$NAA)[1] - 1
#       n_ages <- dim(model_list[[sim]]$model_fxn$rep$NAA)[2] # number of ages
#       n_sex <- dim(model_list[[sim]]$model_fxn$rep$NAA)[3] # number of sexes
#       n_fleets <- dim(model_list[[sim]]$model_fxn$rep$FAA)[3] # number of fleets
#       
#       # Get mean recruitment (divide by 0.5 for sex-ratio)
#       mean_rec <- mean(rowSums(model_list[[sim]]$model_fxn$rep$NAA[1:term_yr,1,])) * 0.5
#       NAA <- model_list[[sim]]$model_fxn$rep$NAA[term_yr,,] # Get numbers at age in terminal year
#       F_Selex <- model_list[[sim]]$model_fxn$rep$F_Slx[term_yr,,,] # fishery selectivity
#       FAA <- model_list[[sim]]$model_fxn$rep$FAA[term_yr,,,] # fishing mortality at age by fleet
#       ZAA <- model_list[[sim]]$model_fxn$rep$ZAA[term_yr,,] # Instatneous mortality
#       F_t <-  exp(model_list[[sim]]$sd_rep$par.fixed[names(model_list[[sim]]$sd_rep$par.fixed) == "ln_Fy"]) # get Fs
#       F_t <- matrix(F_t, nrow = length(1:term_yr), ncol = n_fleets) # munge into matrix
#       N_proj[1,] <- mean_rec # input into first age
#       
#       for(s in 1:n_sex) {
#         for(a in 1:n_ages) {
#           
#           # Exponential survival model 
#           if(a < n_ages) { # not plus group
#             # Project forward with terminal year survival
#             N_proj[a + 1,s] <- NAA[a,s] * exp(-(ZAA[a,s]))
#           } else{
#             N_proj[a,s] <- N_proj[a,s] + (NAA[a,s] * exp(-(ZAA[a,s])))
#           } # else for plus group
#           
#           # Get F40 projections here 
#           if(n_fleets == 1) F40_proj[a,s] <- F40$mle_val[sim] * F_Selex[a,s] # single fleet
#           if(n_fleets > 1){
#             fratio <- F_t[term_yr, ] / sum(F_t[term_yr, ]) # getting fratio here
#             # get F40 projections using fratio
#             F40_proj[a,s] <- (fratio[1] * F40$mle_val[sim] * F_Selex[a,1,s]) + # FABC fl 1
#               (fratio[2] * F40$mle_val[sim] * F_Selex[a,2,s]) # FABC fl 2
#           } # if more than one fleet
#           
#           # Get mortality and catch projections
#           Z_proj[a,s] <- F40_proj[a,s] + oms$Mort_at_age[term_yr, a, sim]
#           # Now, get catch projections
#           CAA_proj[a,s] <- N_proj[a,s] * (1 - exp(-Z_proj[a,s])) * (F40_proj[a,s] / Z_proj[a,s]) 
#           Catch_proj[a,s] <- CAA_proj[a,s] * oms$wt_at_age[term_yr, a, s, j] # Turn to biomass units
#           ABC <- sum(Catch_proj) # sum to get abc
#           
#         } # end a loop
#       } # end s loop
#       
#       # Fill in dataframe
#       em_abc_df[sim,] <- NA # initialize - prob a better way to do this ugh
#       em_abc_df$mle_val[sim] <- ABC # estimated abc
#       em_abc_df$sim[sim] <- sim  # simulation number
#       em_abc_df$OM_Scenario[sim] <- unique_oms[i] # OMs
#       em_abc_df$EM_Scenario[sim] <- files[j] # EMs
#       em_abc_df$type[sim] <- "ABC"
#       em_abc_df$conv[sim] <- F40$conv[sim] # convergence taken from f40 subset from param_df
#       
#       # Fill in time_comp to differeniate
#       if(str_detect(files[j], "Term_")) em_abc_df$time_comp[sim] <- "Terminal"
#       if(str_detect(files[j], "10_")) em_abc_df$time_comp[sim] <- "10 Years Post Fl Chg"
#       if(str_detect(files[j], "5_")) em_abc_df$time_comp[sim] <- "5 Years Post Fl Chg"
#       
#     } # end sim loop
#     
#     # store in list to speed up
#     all_em_abc_df <- rbind(em_abc_df, all_em_abc_df)
#     print(j)
#   } # end j loop
#   print(i)
# } # end i loop
# 
# # Left join the om df to the estiamted abc df
# all_em_abc_df <- all_em_abc_df %>% left_join(all_om_abc_df,
#                                              by = c("sim", "OM_Scenario", "time_comp")) %>% 
#   mutate(RE = (mle_val - t) / t) %>% 
#   data.frame()
# 
# # Now, bind this to parameter summary
# param_df <- rbind(all_em_abc_df, param_df)
# write.csv(param_df, here("output", "Parameter_Summary.csv"), row.names = FALSE)

# # Get Total Biomass -------------------------------------------------------
# 
# # Get total biomass estimates
# total_biom_df <- list()
# for(i in 1:length(unique(ts_df$OM_Scenario))) {
#   
#   # List files in folder
#   files <- list.files(here(om_scenario_path, unique_oms[i]))
#   # Load in OM data
#   load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
#   # Remove.RData and .pdf
#   files <- files[str_detect(files, ".RData|.pdf") == FALSE]
#   
#   # Pre-allocate list for EMs within OMs
#   em_biom_list <- list()
#   
#   for(j in 1:length(files)) {
#     
#     # Load in .RData from model runs
#     load(here(om_scenario_path, unique_oms[i], files[j], paste(files[j], ".RData", sep = "")))
#     
#     # Pre-allocate list object here
#     biom_list <- list()
#     
#     # Get total biomass estimate from models here
#     for(m in 1:length(model_list)) {
#       
#       # Get model biomass
#       mod_biom <- model_list[[m]]$sd_rep$value[names(model_list[[m]]$sd_rep$value) == "Total_Biom"]
#       # Get true biomass values
#       true_biom <- rowSums(oms$Biom_at_age[1:length(mod_biom),,,m])
#       
#       # Check convergence here
#       convergence_status <- check_model_convergence(mle_optim = model_list[[m]]$mle_optim,
#                                                     mod_rep = model_list[[m]]$model_fxn,
#                                                     sd_rep = model_list[[m]]$sd_rep,
#                                                     min_grad = 0.01)
#       
#       # Now, stick this into a dataframe
#       mod_biom_df_bind <- data.frame(X = NA, mle_val = mod_biom, mle_se = NA,
#                                      lwr_95 = NA, upr_95 = NA, t = true_biom,
#                                      sim = m, conv = convergence_status$Convergence, 
#                                      year = 1:length(mod_biom), type = "Total Biomass",
#                                      OM_Scenario = unique_oms[i], EM_Scenario = files[j]) %>% 
#         data.frame()
#       
#       # Now bind everything together
#       biom_list[[m]] <- mod_biom_df_bind
#       
#     } # end m loop
#     
#     # Turn from list to dataframe
#     mod_biom_df <- data.table::rbindlist(biom_list)
#     em_biom_list[[j]] <- mod_biom_df # put biomass dataframes from EMs into list
#     print(j)
#   } # end j loop
#   
#   # turn em list from list into dataframe
#   em_biom_df <- data.table::rbindlist(em_biom_list)
#   total_biom_df[[i]] <- em_biom_df # put biomass dataframes from ems into om list
#   print(i)
# } # end i loop
# 
# # Bind together and write out csv
# ts_df <- rbind(data.table::rbindlist(total_biom_df), ts_df)
# write.csv(ts_df, here("output", "TimeSeries_Summary.csv"))