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

# Paths
om_scenario_path <- here("output", "OM_Scenarios_Sensitivity") # path to OM folder
dir_out <- here("output", "Summary_Plots") # path to output folder
dir.create(dir_out)

# Collate Parameter Results -----------------------------------------------

# Get results from self tests
all_results <- get_results(om_scenario_path = om_scenario_path)
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
data.table::fwrite(ts_df, here("output", "OM_Scenarios_Sensitivity", "TimeSeries_Summary.csv"))
data.table::fwrite(aic_df, here("output", "OM_Scenarios_Sensitivity","AIC_Convergence_Summary.csv"))

# Parameter Results ------------------------------------------------------
# Create relative error and CV metrics for parameters
param_df <- param_df %>% mutate(RE = (mle_val - t) / t, # Relative Error
                                TE = mle_val - t, # Total Error
                                ARE = abs(mle_val - t/t)) %>%  # Absolute Relative Error
  dplyr::select(-X)

# Differentiate time components
param_df <- param_df %>% 
  mutate(
    EM_Scenario_Full = EM_Scenario, 
    time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int") ~ "Fleet Intersect", # Fleet Transition Intersects
  ),
  EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'), # remove preceeding letters
  time_comp = factor(time_comp, levels = c("Terminal", "Fleet Trans End",
                                           "Fleet Intersect"))) %>% 
  data.frame()

data.table::fwrite(param_df, here("output", "OM_Scenarios_Sensitivity", "Parameter_Summary.csv"))

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

data.table::fwrite(ts_re_df, here("output", "OM_Scenarios_Sensitivity", "TimeSeries_RE.csv"))
data.table::fwrite(ts_te_df, here("output", "OM_Scenarios_Sensitivity", "TimeSeries_TE.csv"))
data.table::fwrite(ts_are_df, here("output", "OM_Scenarios_Sensitivity", "TimeSeries_ARE.csv"))

# Get Selectivity estimates ---------------------------------------------------------
for(i in 1:length(unique_oms)) {
  
  # List files in folder
  files <- list.files(here(om_scenario_path, unique_oms[i]))
  # Load in OM data
  load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
  # Remove.RData and .pdf
  files <- files[str_detect(files, ".RData|.pdf|.csv") == FALSE]
  # Pre-allocate list to store all EMs
  em_slx_list <- list()
  
  for(j in 1:length(files)) {
    
    # Load in .RData from model runs
    load(here(om_scenario_path, unique_oms[i], files[j], paste(files[j], ".RData", sep = "")))
    
    # Pre-allocate list object here
    F_Slx_list <- list()
    
    # Get convergence statistics here
    convergence <- param_df[param_df$OM_Scenario == unique_oms[i] &
                            param_df$EM_Scenario_Full == files[j],]
    
    # Get total biomass estimate from models here
    for(m in 1:length(model_list)) {
      
      if(length(model_list[[m]])!=0) {
        # Get fishery selectivity from model
        F_Slx_df <- data.table::data.table(reshape2::melt(model_list[[m]]$model_fxn$rep$F_Slx))
        names(F_Slx_df) <- c("Year", "Age", "Fleet", "Sex", "Selex")
        
        # Construct our dataframe
        F_Slx_df_bind <- cbind(F_Slx_df, sim = m, 
                               conv = convergence$conv[m],
                               OM_Scenario = unique_oms[i],
                               EM_Scenario = files[j])
        
        # Now bind everything together
        F_Slx_list[[m]] <- F_Slx_df_bind
      }
    } # end m 
    
    # Turn from list to dataframe
    F_Slx_mod_df <- data.table::rbindlist(F_Slx_list)
    em_slx_list[[j]] <- F_Slx_mod_df # put biomass dataframes from EMs into list
    print(j)
  } # end j
  
  # Rbind out list
  EM_fish_Slx <- data.table::rbindlist(em_slx_list)
  
  # Do some residual munging and output into OM folders
  EM_fish_Slx <- EM_fish_Slx %>% 
    mutate(time_comp = case_when(
      str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
      str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
      str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
    ),
    time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")),
    Sex = case_when(
      str_detect(Sex, "1") ~ "Female",
      str_detect(Sex, "2") ~ "Male"), 
    EM_Scenario_Full = EM_Scenario, # Retaining full EM Scenario Name
    EM_Scenario = str_remove(EM_Scenario, "Term_|TrxE_|Int_"))
  
  data.table::fwrite(EM_fish_Slx, file = here('output', "OM_Scenarios_Sensitivity",
                                              unique_oms[i], "EM_Fish_Selex.csv"))
 
   # Get OM selectivity 
  OM_fish_Slx <- data.table::data.table(reshape2::melt(oms$Fish_selex_at_age[,,,,1]))
  names(OM_fish_Slx) <- c("Year", "Age", "Fleet", "Sex", "Selex")
  OM_fish_Slx$OM_Scenario <- unique_oms[i]
  
  # Now parse out the numbers from this
  OM_fish_Slx <- OM_fish_Slx %>% 
    mutate(Year = parse_number(paste(Year)),
           Age = parse_number(paste(Age)),
           Fleet = parse_number(paste(Fleet)),
           Sex = parse_number(paste(Sex)),
           Sex = case_when( # Differentiate sexes
             str_detect(Sex, "1") ~ "Female",
             str_detect(Sex, "2") ~ "Male"))
  
  # write out csv
  data.table::fwrite(OM_fish_Slx, file = here('output', "OM_Scenarios_Sensitivity",
                                              unique_oms[i], "OM_Fish_Selex.csv"))
  
  print(i)
} # end i

# Now read these all back in, bind them together and output as one single file
# Need to process huge outputs like this, or else R fucking crashes
om_list <- list()

for(i in 1:length(unique_oms)) {
  OM_Selex = data.table::fread(here('output', "OM_Scenarios_Sensitivity",
                                    unique_oms[i], "OM_Fish_Selex.csv"))
  om_list[[i]] <- OM_Selex
  print(i)
} # end i

# Write out csvs here
data.table::fwrite(data.table::rbindlist(om_list), here("output", "OM_Scenarios_Sensitivity", "OM_Fish_Selex.csv"), row.names = FALSE)

# Get Population Selex ----------------------------------------------------
time_comp <- c("First Year", "Fleet Intersect", "Fleet Trans End", "Terminal")

### OM ----------------------------------------------------------------------
pop_sel_om <- data.frame()
for(i in 1:length(unique_oms)) {
  
  # Load in a given OM
  load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
  # Assessment years to filter to
  if(str_detect(unique_oms[i], "Fast")) yr <- c(1, 27, 30, 50)
  if(str_detect(unique_oms[i], "Slow")) yr <- c(1, 36, 50, 70)
  
  for(j in 1:length(time_comp)) {
    
    # Get Fratio here
    fratio <- oms$fish_mort[yr[j],,1] / sum(oms$fish_mort[yr[j],,1])
    
    # Get weighted average pop selex here
    Female <- colSums(fratio * t(oms$Fish_selex_at_age[1,,,1,1])) # female
    Male <- colSums(fratio * t(oms$Fish_selex_at_age[1,,,2,1])) # male
    
    # Put this into a dataframe
    pop_sel_df <- data.frame(cbind(Female, Male), Age = 1:30, 
                             OM_Scenario = unique_oms[i],
                             time_comp = time_comp[j])
    
    pop_sel_om <- rbind(pop_sel_df, pop_sel_om)
    
  } # end j loop
} # end i loop

# Residual munging here
pop_sel_om <- pop_sel_om %>% 
  pivot_longer(names_to = "Sex", values_to = "Selex", 
               cols = c("Female", "Male"))

data.table::fwrite(pop_sel_om, here("output", "OM_Scenarios_Sensitivity", "Pop_Selex_OM.csv"))

### EM ----------------------------------------------------------------------
pop_sel_em <- data.frame()
twofleet_pop_sel_ems <- data.frame()
for(i in 1:length(unique_oms)) {
  
  # Assessment years to filter to
  if(str_detect(unique_oms[i], "Fast") & str_detect(unique_oms[i], "Ext")) {
    yr <- c(30, 99)
    time_comp <- "Terminal"
  }
  if(str_detect(unique_oms[i], "Fast") & !str_detect(unique_oms[i], "Ext")) {
    yr <- c(27, 30, 50)
    time_comp <- "Terminal"
  }
  if(str_detect(unique_oms[i], "Slow")) {
    yr <- c(36, 50, 70)
    time_comp <- "Terminal"
  }
  
  # Read in EM Fishery Selex
  EM_fish_Slx = data.table::fread(here('output', "OM_Scenarios_Sensitivity", unique_oms[i], "EM_Fish_Selex.csv"))
  EM_fish_Slx = EM_fish_Slx[EM_fish_Slx$conv == "Converged"]
  
  for(j in 1:length(time_comp)) {
    
    # Time component filtering
    tc_filter <- time_comp[j]

    # Subset to desired dataframe
    pop_sel_em_sub <- EM_fish_Slx %>% 
      filter(OM_Scenario == unique_oms[i], time_comp == tc_filter,
             Year %in% c((yr[j] - 4):yr[j]))
    
    # Subset to onefleet variants
    onefleet_pop_sel_em <- pop_sel_em_sub[str_detect(pop_sel_em_sub$EM_Scenario, "1Fl")]
    
    # Summarize one fleet variants
    onefleet_summary <- onefleet_pop_sel_em %>% 
      group_by(OM_Scenario, EM_Scenario, Sex, Age, sim, time_comp) %>% 
      summarize(Selex = mean(Selex)) %>% 
      ungroup() %>% 
      group_by(OM_Scenario, EM_Scenario, Sex, Age, time_comp) %>% 
      summarize(Median_Selex = median(Selex),
                Lwr_95 = quantile(Selex, 0.025),
                Upr_95 = quantile(Selex, 0.975),
                Lwr_75 = quantile(Selex, 0.125, na.rm = T),
                Upr_75 = quantile(Selex, 0.875, na.rm = T))
    
    # Subset to twofleet variants
    twofleet_pop_sel_em <- pop_sel_em_sub[str_detect(pop_sel_em_sub$EM_Scenario, "2Fl")]
    
    # Get unique 2 fleet ems 
    unique_models_full <- unique(twofleet_pop_sel_em$EM_Scenario_Full)
    unique_models <- unique(twofleet_pop_sel_em$EM_Scenario)
    
    # Load in model_lists from unique EM models
    for(k in 1:length(unique_models_full)) {
      
      # load in 2 fleet models
      load(here(om_scenario_path, unique_oms[i], unique_models_full[k], paste(unique_models_full[k], ".RData", sep = "")))
      
      # Loop through simulation runs to get population selex (fratio * selexAA for both males and females)
      for(m in 1:length(model_list)) {
        if(length(model_list[[m]]) > 0) { # only go thorugh this if the list isn't null
          # Get Fratio - munge into matrix format
          fratio <- matrix(exp(model_list[[m]]$sd_rep$par.fixed[names(model_list[[m]]$sd_rep$par.fixed) == "ln_Fy"]),
                           ncol = 2, nrow = length(1:yr[j]))
          fratio <- fratio[length(1:yr[j]),] / sum(fratio[length(1:yr[j]),])
          
          # Female Selectivity
          Female <- colSums(fratio * t(colMeans(model_list[[m]]$model_fxn$rep$F_Slx[(yr[j] - 4):yr[j],,,1])))
          # Male Selectivity
          Male <- colSums(fratio * t(colMeans(model_list[[m]]$model_fxn$rep$F_Slx[(yr[j] - 4):yr[j],,,2])))
          
          # Put this into a dataframe
          twofleet_pop_sel_df <- data.frame(cbind(Female, Male), Age = 1:30, 
                                            OM_Scenario = unique_oms[i],
                                            EM_Scenario = unique_models[k],
                                            time_comp = time_comp[j],
                                            sim = m)
          
          twofleet_pop_sel_ems <- rbind(twofleet_pop_sel_df, twofleet_pop_sel_ems)
        } # end if statement
      } # end m loop
    } # end k loop
    
    pop_sel_em <- rbind(pop_sel_em, onefleet_summary)
    
  } # end j loop
} # end i loop

# Residual munging for 2 fleet models
twofleet_pop_sel_ems <- twofleet_pop_sel_ems %>% 
  pivot_longer(names_to = "Sex", values_to = "Selex", cols = c("Female", "Male")) %>% 
  group_by(OM_Scenario, EM_Scenario, Sex, Age, sim, time_comp) %>% 
  summarize(Selex = mean(Selex)) %>% 
  ungroup() %>% 
  group_by(OM_Scenario, EM_Scenario, Sex, Age, time_comp) %>% 
  summarize(Median_Selex = median(Selex),
            Lwr_95 = quantile(Selex, 0.025),
            Upr_95 = quantile(Selex, 0.975),
            Lwr_75 = quantile(Selex, 0.125, na.rm = T),
            Upr_75 = quantile(Selex, 0.875, na.rm = T))

pop_sel_em <- rbind(twofleet_pop_sel_ems, pop_sel_em)
data.table::fwrite(pop_sel_em, here("output", "OM_Scenarios_Sensitivity", "Pop_Selex_EM.csv"))

# # Get Numbers at age ------------------------------------------------------
all_NAA_list <- list() # NAA estimates for both OM and EM
for(i in 1:length(unique_oms)) {
  
  # List files in folder
  files <- list.files(here(om_scenario_path, unique_oms[i]))
  # Load in OM data
  load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
  # Remove.RData and .pdf
  files <- files[str_detect(files, ".RData|.pdf|.csv") == FALSE]
  om_NAA_list <- list() # pre-allocate length
  
  # Get NAA for OMs
  om_NAA_df <- reshape2::melt(oms$N_at_age)
  names(om_NAA_df) <- c("Year", "Age", "Sex", "sim", "True_Numbers")
  
  # Parse numbers out
  om_NAA_df <- om_NAA_df %>% mutate(
    Year = parse_number(as.character(Year)),
    Age = parse_number(as.character(Age)),
    Sex = parse_number(as.character(Sex)),
    sim = parse_number(as.character(sim)))
  
  for(j in 1:length(files)) {
    
    # Load in .RData from model runs
    load(here(om_scenario_path, unique_oms[i],
              files[j], paste(files[j], ".RData", sep = "")))
    
    # Pre-allocate list object here
    em_NAA_list <- list()
    
    # Get convergence statistics here
    convergence <- param_df[param_df$OM_Scenario == unique_oms[i] &
                              param_df$EM_Scenario_Full == files[j],]
    
    # Get total biomass estimate from models here
    for(m in 1:length(model_list)) {
      
      if(length(model_list[[m]]) != 0) {
        # Get fishery selectivity from model
        em_NAA_df <- reshape2::melt(model_list[[m]]$model_fxn$rep$NAA)
        names(em_NAA_df) <- c("Year", "Age", "Sex", "Est_Numbers")
        
        # Subset to correct simulation
        om_sim_NAA = om_NAA_df[om_NAA_df$sim == m &
                                 om_NAA_df$Year %in% c(1:max(em_NAA_df$Year)),]$True_Numbers
        
        # Construct our dataframe
        em_NAA_df_bind <- cbind(em_NAA_df, sim = m,
                                True_Numbers = om_sim_NAA,
                                conv = convergence$conv[m],
                                OM_Scenario = unique_oms[i],
                                EM_Scenario = files[j])
        
        # Now bind everything together
        em_NAA_list[[m]] <- em_NAA_df_bind
      } # end if
    } # end m
    
    # Turn from list to dataframe
    em_NAA_df <- data.table::rbindlist(em_NAA_list)
    om_NAA_list[[j]] <- em_NAA_df # put dataframes from EMs into list
    print(j)
    
  } # end j
  
  # Rbind list and do some residual munging.
  om_NAA_df <- data.table::rbindlist(om_NAA_list)
  
  # Changing some names, etc.
  om_NAA_df <- om_NAA_df %>%
    mutate(time_comp = case_when(
      str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
      str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
      str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
    ),
    time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")),
    Sex = case_when(
      str_detect(Sex, "1") ~ "Female",
      str_detect(Sex, "2") ~ "Male"),
    EM_Scenario_Full = EM_Scenario, # Retaining full EM Scenario Name
    EM_Scenario = str_remove(EM_Scenario, "Term_|TrxE_|Int_"))
  
  # Now, output this
  data.table::fwrite(om_NAA_df, file = here('output', "OM_Scenarios_Sensitivity",
                                            unique_oms[i], "OM_EM_NAA.csv"))
  
  print(i)
} # end i


### Summarize NAA -----------------------------------------------------------
all_naa_sum <- data.frame()
for(i in 1:length(unique_oms)) {
  # read in data
  NAA_df_re = data.table::fread(here("output", "OM_Scenarios_Sensitivity", unique_oms[i], "OM_EM_NAA.csv"))
  
  # Do some residual munging
  NAA_df_re_plot <- NAA_df_re %>%
    filter(conv == "Converged") %>% 
    mutate(RE = (Est_Numbers - True_Numbers) / True_Numbers) %>%
    group_by(OM_Scenario, EM_Scenario, Year, Age, Sex, time_comp) %>%
    summarize(Median_RE = median(RE, na.rm = TRUE),
              Lwr_95 = quantile(RE, 0.025, na.rm = TRUE),
              Upr_95 = quantile(RE, 0.975, na.rm = TRUE))
  
  # Now relevel factor for organizing plot
  NAA_df_re_plot <- NAA_df_re_plot %>%
    mutate(time_comp = factor(time_comp,
                              levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")))
  all_naa_sum <- rbind(all_naa_sum, NAA_df_re_plot) # bind
  print(i)
}

write.csv(all_naa_sum, here("output", "OM_Scenarios_Sensitivity", "NAA_Summary.csv")) # write out csv

