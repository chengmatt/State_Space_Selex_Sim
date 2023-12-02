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
data.table::fwrite(aic_df, here("output", "AIC_Convergence_Summary.csv"))

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

data.table::fwrite(param_df, here("output", "Parameter_Summary.csv"))

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

data.table::fwrite(ts_re_df, here("output", "TimeSeries_RE.csv"))
data.table::fwrite(ts_te_df, here("output", "TimeSeries_TE.csv"))
data.table::fwrite(ts_are_df, here("output", "TimeSeries_ARE.csv"))

# Get Selectivity estimates ---------------------------------------------------------
for(i in 8:length(unique_oms)) {
  
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
  
  data.table::fwrite(EM_fish_Slx, file = here('output', "OM_Scenarios",
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
  data.table::fwrite(OM_fish_Slx, file = here('output', "OM_Scenarios",
                                              unique_oms[i], "OM_Fish_Selex.csv"))
  
  print(i)
} # end i

# Now read these all back in, bind them together and output as one single file
# Need to process huge outputs like this, or else R fucking crashes
om_list <- list()

for(i in 1:length(unique_oms)) {
  OM_Selex = data.table::fread(here('output', "OM_Scenarios",
                                    unique_oms[i], "OM_Fish_Selex.csv"))
  om_list[[i]] <- OM_Selex
  print(i)
} # end i

# Write out csvs here
data.table::fwrite(data.table::rbindlist(om_list), here("output", "OM_Fish_Selex.csv"), row.names = FALSE)

# Get Population Selex ----------------------------------------------------
time_comp <- c("Fleet Intersect", "Fleet Trans End", "Terminal")

### OM ----------------------------------------------------------------------
pop_sel_om <- data.frame()
for(i in 1:length(unique_oms)) {
  
  # Load in a given OM
  load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
  # Assessment years to filter to
  if(str_detect(unique_oms[i], "Fast")) yr <- c(27, 30, 50)
  if(str_detect(unique_oms[i], "Slow")) yr <- c(36, 50, 70)
  
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

data.table::fwrite(pop_sel_om, here("output", "Pop_Selex_OM.csv"))

### EM ----------------------------------------------------------------------
pop_sel_em <- data.frame()
twofleet_pop_sel_ems <- data.frame()
for(i in 1:length(unique_oms)) {
  
  # Assessment years to filter to
  if(str_detect(unique_oms[i], "Fast")) yr <- c(27, 30, 50)
  if(str_detect(unique_oms[i], "Slow")) yr <- c(36, 50, 70)
  
  # Read in EM Fishery Selex
  EM_fish_Slx = data.table::fread(here('output', "OM_Scenarios", unique_oms[i], "EM_Fish_Selex.csv"))
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
data.table::fwrite(pop_sel_em, here("output", "Pop_Selex_EM.csv"))

# Get SPR from population Selex -------------------------------------------

# Load in random OM to get WAA, MortAA, and MatAA
load(here(om_scenario_path, "Fast_LL_High", paste("Fast_LL_High", ".RData", sep = "")))

om_spr_df <- data.frame()
em_spr_df <- data.frame()
for(i in 1:length(unique_oms)) {
  
  # Filter out time-block stuff 
  em_df <- pop_sel_em %>% filter(OM_Scenario == unique_oms[i], Sex == "Female")
  # Get OM dataframe selex
  om_df <- pop_sel_om %>% filter(Sex == "Female", OM_Scenario == unique_oms[i])
  
  for(j in 1:length(unique(om_df$time_comp))) {
    
    # Get true OM SPR
    sub_om_df <- om_df %>% filter(time_comp == unique(om_df$time_comp)[j])
    
    # Get OM SPR Values
    om_spr_values <- get_trialF_spr(MortAA = oms$Mort_at_age[1,,1], 
                                    SelexAA = sub_om_df$Selex,
                                    MatAA = oms$mat_at_age[1,,1,1], 
                                    WAA = oms$wt_at_age[1,,1,1],
                                    trial_F = seq(0.001, 0.5, 0.001),
                                    F_x = 0.4) %>% 
      mutate(OM_Scenario = unique_oms[i],
             time_comp = unique(om_df$time_comp)[j])
    
    # Bind together OM SPR values
    om_spr_df <- rbind(om_spr_df, om_spr_values)
    
    # Get unique EMs
    unique_ems <- unique(em_df$EM_Scenario)
    
    for(e in 1:length(unique_ems)) {
      
      # Subset to time component
      sub_em_df <- em_df %>% filter(time_comp == unique(om_df$time_comp)[j],
                                    EM_Scenario == unique_ems[e])
      
      # Get EM SPR Values
      em_spr_values <- get_trialF_spr(MortAA = oms$Mort_at_age[1,,1], 
                                      SelexAA = sub_em_df$Median_Selex,
                                      MatAA = oms$mat_at_age[1,,1,1], 
                                      WAA = oms$wt_at_age[1,,1,1],
                                      trial_F = seq(0.001, 0.5, 0.001),
                                      F_x = 0.4) %>% 
        mutate(OM_Scenario = unique_oms[i],
               EM_Scenario = unique_ems[e],
               time_comp = unique(om_df$time_comp)[j])
      
      # Bind together OM SPR values
      em_spr_df <- rbind(em_spr_df, em_spr_values)
      
    } # end e loop
  } # end j loop
  print(i)
} # end i loop

data.table::fwrite(em_spr_df, here("output", "EM_SPR.csv"))
data.table::fwrite(om_spr_df, here("output", "OM_SPR.csv"))


# Sensitivity for SPR (Maturity and Selectivity Schedules) ----------------
om_spr_mat_df <- data.frame()
em_spr_mat_df <- data.frame()

# Maturity schedules
MatAA_OM = oms$mat_at_age[1,,1,1]
MatAA_10 = 1 / (1 + exp(0.4 * (10 - 1:30))) # Age at 50% maturity = 10
MatAA_4 = 1 / (1 + exp(1.25 * (4 - 1:30))) # Age at 50% maturity = 4

plot(MatAA_OM, type = "l", xlab = "Age", ylab = "Maturity")
lines(MatAA_10, type = "l", col = "blue")
lines(MatAA_4, type = "l", col = "red")
legend(x = 8, y = 0.5,
       legend = c("Maturity 3 years back",  "OM Maturity", "Maturity 3 years ahead"), 
       fill = c("red", "black", "blue"))

# Maturity options
mat_opt <- list(MatAA_OM = MatAA_OM, MatAA_10 = MatAA_10, MatAA_4 = MatAA_4)

for(mat in 1:length(mat_opt)) {
  for(i in 1:length(unique_oms)) {
    
    # Filter out time-block stuff 
    em_df <- pop_sel_em %>% filter(OM_Scenario == unique_oms[i], Sex == "Female")
    # Get OM dataframe selex
    om_df <- pop_sel_om %>% filter(Sex == "Female", OM_Scenario == unique_oms[i])
    
    for(j in 1:length(unique(om_df$time_comp))) {
      
      # Get true OM SPR
      sub_om_df <- om_df %>% filter(time_comp == unique(om_df$time_comp)[j])
      
      # Get OM SPR Values
      om_spr_values <- get_trialF_spr(MortAA = oms$Mort_at_age[1,,1], 
                                      SelexAA = sub_om_df$Selex,
                                      MatAA = mat_opt[[mat]], 
                                      WAA = oms$wt_at_age[1,,1,1],
                                      trial_F = seq(0.001, 0.5, 0.001),
                                      F_x = 0.4) %>% 
        mutate(OM_Scenario = unique_oms[i],
               time_comp = unique(om_df$time_comp)[j],
               Mat_Opt = names(mat_opt[mat]))
      
      # Bind together OM SPR values
      om_spr_mat_df <- rbind(om_spr_mat_df, om_spr_values)
      
      # Get unique EMs
      unique_ems <- unique(em_df$EM_Scenario)
      
      for(e in 1:length(unique_ems)) {
        
        # Subset to time component
        sub_em_df <- em_df %>% filter(time_comp == unique(om_df$time_comp)[j],
                                      EM_Scenario == unique_ems[e])
        
        # Get EM SPR Values
        em_spr_values <- get_trialF_spr(MortAA = oms$Mort_at_age[1,,1], 
                                        SelexAA = sub_em_df$Median_Selex,
                                        MatAA = mat_opt[[mat]], 
                                        WAA = oms$wt_at_age[1,,1,1],
                                        trial_F = seq(0.001, 0.5, 0.001),
                                        F_x = 0.4) %>% 
          mutate(OM_Scenario = unique_oms[i],
                 EM_Scenario = unique_ems[e],
                 time_comp = unique(om_df$time_comp)[j],
                 Mat_Opt = names(mat_opt[mat]))
        
        # Bind together OM SPR values
        em_spr_mat_df <- rbind(em_spr_mat_df, em_spr_values)
        
      } # end e loop
    } # end j loop
    print(i)
  } # end i loop
} # end maturity option loop

data.table::fwrite(em_spr_mat_df, here("output", "EM_SPR_MatSens.csv"))
data.table::fwrite(om_spr_mat_df, here("output", "OM_SPR_MatSens.csv"))
# 
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
  data.table::fwrite(om_NAA_df, file = here('output', "OM_Scenarios",
                                              unique_oms[i], "OM_EM_NAA.csv"))

  print(i)
} # end i

# # Get composition fits ----------------------------------------------------
# all_em_fish_comps_list <- list() # em pred composition list
# om_fish_comps_list <- list() # om composition list
# 
# for(i in 1:length(unique_oms)) {
#   
#   # List files in folder
#   files <- list.files(here(om_scenario_path, unique_oms[i]))
#   # Load in OM data
#   load(here(om_scenario_path, unique_oms[i], paste(unique_oms[i], ".RData", sep = "")))
#   # Remove.RData and .pdf
#   files <- files[str_detect(files, ".RData|.pdf|.csv") == FALSE]
#   # Pre-allocate list to store all EMs
#   em_fish_comps_list <- list()
#   
#   for(j in 1:length(files)) {
#     
#     # Load in .RData from model runs
#     load(here(om_scenario_path, unique_oms[i], files[j], paste(files[j], ".RData", sep = "")))
#     
#     # Pre-allocate list object here
#     Fish_comps_list <- list()
#     
#     # Get total biomass estimate from models here
#     for(m in 1:length(model_list)) {
#       
#       # Get fishery selectivity from model
#       Fish_Comps_df <- reshape2::melt(model_list[[m]]$model_fxn$rep$pred_fish_age_comps)
#       names(Fish_Comps_df) <- c("Year", "Age", "Fleet", "Sex", "Prop")
#       
#       # Check convergence here
#       convergence_status <- aic_df %>% 
#           filter(sim == m, OM_Scenario == unique_oms[i], EM_Scenario == files[j])
#         
#       # Construct our dataframe
#       F_Comps_df_bind <- cbind(Fish_Comps_df, sim = m, 
#                                conv = convergence_status$conv,
#                                OM_Scenario = unique_oms[i], EM_Scenario = files[j])
#       
#       
#       # Now bind everything together
#       Fish_comps_list[[m]] <- F_Comps_df_bind
#       
#     } # end m 
#     
#     # Turn from list to dataframe
#     Fish_Comps_mod_df <- data.table::rbindlist(Fish_comps_list)
#     em_fish_comps_list[[j]] <- Fish_Comps_mod_df # put biomass dataframes from EMs into list
#     print(j)
#     
#   } # end j
#   
#   # put selex dataframes from ems into all models list
#   all_em_fish_comps_list[[i]] <- data.table::rbindlist(em_fish_comps_list) 
#   
#   # OM Munging Comps and CAA ------------------------------------------------
#   
#   # Empty array to store stuff in for each OM 1 fleet combo
#   obs_fish_age_comps <- array(data = 0, dim = c(dim(oms$N_at_age)[1] - 1,  # years
#                                                 dim(oms$N_at_age)[2], # ages
#                                                 1, # fleets
#                                                 dim(oms$Fish_Age_Comps)[4], # sexes
#                                                 dim(oms$N_at_age)[4])) # simulations
#   
#   # Get fishery age comps for 1 fleet model
#   for(sim in 1:dim(oms$N_at_age)[4]) { # number of simulations loop 
#     # Observed catches 
#     obs_catches <- as.matrix(apply(as.matrix(oms$Catch_agg[1:50,,sim]), 1, FUN = sum), ncol = 1)
#     # Catch weighting here
#     catch_sum <- rowSums(as.matrix(oms$Catch_agg[1:50,,sim]))
#     catch_weight <- as.matrix(oms$Catch_agg[1:50,,sim]) / catch_sum 
#     
#     for(s in 1:dim(oms$Fish_Age_Comps)[4]) { # sex loop
#       for(f in 1:dim(oms$Fish_Age_Comps)[3]) { # fleet loop
#         # Filter to save as an object
#         fish_age_comps <-  oms$Fish_Age_Comps[1:(dim(oms$N_at_age)[1] - 1),,f,s,sim] * catch_weight[,f]
#         # Increment comps - fixing fleet index to 1 here
#         obs_fish_age_comps[,,1,s,sim] <- obs_fish_age_comps[,,1,s,sim] + fish_age_comps
#       } # end f loop
#     } # end s loop
#     
#     # Now, apply the proportion function over a single fleet
#     for(s in 1:dim(oms$Fish_Age_Comps)[4]) {
#       obs_fish_age_comps[,,,s, sim] <- t(apply(obs_fish_age_comps[,,,s,sim], MARGIN = 1, 
#                                                FUN=function(x) { x/sum(x) }))
#     } # s loop
#   } # end sim loop
#   
#   # one fleet comps
#   one_fleet_comps <- reshape2::melt(obs_fish_age_comps)
#   names(one_fleet_comps) <- c("Year", "Ages", "Fleets", "Sexes", "Sim", "Prop")
#   one_fleet_comps$fleet_type <- "One Fleet"
#   
#   # Empty array to store stuff in for each OM 2fleet combo
#   two_obs_fish_age_comps <- array(data = 0, dim = c(dim(oms$N_at_age)[1] - 1,  # years
#                                                     dim(oms$N_at_age)[2], # ages
#                                                     dim(oms$Fish_Age_Comps)[3], # fleets
#                                                     dim(oms$Fish_Age_Comps)[4], # sexes
#                                                     dim(oms$N_at_age)[4])) # simulations
#   
#   # Get two fleet comps
#   for(sim in 1:dim(oms$N_at_age)[4]) {
#     for(s in 1:dim(oms$Fish_Age_Comps)[4]) { # sex loop
#       for(f in 1:dim(oms$Fish_Age_Comps)[3]) { # fleet loop
#         two_obs_fish_age_comps[,,f,s, sim] <- t(
#           apply(oms$Fish_Age_Comps[1:(dim(oms$N_at_age)[1] - 1),,f,s,sim], MARGIN = 1, 
#                 FUN=function(x) { x/sum(x) })
#         )
#       } # end s loop
#     } # end f loop
#   } # end sim loop
#   
#   two_fleet_comps <- reshape2::melt(two_obs_fish_age_comps)
#   names(two_fleet_comps) <- c("Year", "Ages", "Fleets", "Sexes", "Sim", "Prop")
#   two_fleet_comps$fleet_type <- "Two Fleets"
#   
#   # Now, bind the 2 fleet and one fleet dataframe together
#   fleet_comps <- rbind(one_fleet_comps, two_fleet_comps)
#   fleet_comps$OM_Scenario <- unique_oms[i] # name OM
#   om_fish_comps_list[[i]] <- fleet_comps # input into our om list
#   print(i)
#   
# } # end i
# 
# # Rbind to list all ems 
# EM_fish_comps <- data.table::rbindlist(all_em_fish_comps_list)
# om_fish_comps_df <- data.table::rbindlist(om_fish_comps_list)
# 
# data.table::fwrite(EM_fish_comps, here("output", "EM_Fish_Comps.csv"), row.names = FALSE)
# data.table::fwrite(om_fish_comps_df, here("output", "OM_Fish_Comps.csv"), row.names = FALSE)
# 

# Get likelihood components (fishery comps) -----------------------------------------------
all_om_list = list()
for(i in 1:length(unique_oms)) {
  
  # List out EMs within each OM here
  ems = list.files(here(om_scenario_path, unique_oms[i]))
  ems <- ems[str_detect(ems, ".RData|.pdf|.csv") == FALSE]
  em_list = list()
  
  # Now load in unique EMs and extract values
  for(k in 1:length(unique(ems))) {
    
    # load in rdata
    load(here(om_scenario_path, unique_oms[i], ems[k], paste(ems[k], ".RData", sep = ""))) 
    # figure out number of fleets
    n_fleets = dim(model_list[[1]]$model_fxn$rep$fish_comp_nLL)[2]
    sim_list = list()
    
    for(sim in 1:length(model_list)) {
      fleet_df = data.frame() # to store fleet specific stuff
      for(f in 1:n_fleets) {
        females_cumsum = cumsum(model_list[[sim]]$model_fxn$rep$fish_comp_nLL[,f,1]) # females
        males_cumsum = cumsum(model_list[[sim]]$model_fxn$rep$fish_comp_nLL[,f,2]) # males
        store_df = data.frame(Females = females_cumsum, Males = males_cumsum, OM = unique_oms[i],
                   EM = ems[k], Year = 1:length(females_cumsum), sim = sim, Fleet = paste("Fleet", f))
        fleet_df = rbind(fleet_df, store_df)
      } # end f fleet loop
      sim_list[[sim]] = fleet_df # input into simulation list
    } # sim loop
    
    # Rbind list into dataframe for simulations of a given EM
    sim_df <- data.table::rbindlist(sim_list)
    em_list[[k]] = sim_df # input simulation for a given EM into em list
    print(k)
  } # end k
  
  # rbind list into om list
  em_df <- data.table::rbindlist(em_list)
  all_om_list[[i]] = em_df
} # end i

# rbind list into df and write it out
fish_comps_likes_om_df = data.table::rbindlist(all_om_list)
data.table::fwrite(fish_comps_likes_om_df, here("output", "FishComp_Likes.csv"))

