
# Purpose: To replot manuscript figures 
# Creator: Matthew LH. Cheng
# 11/29/23

# Set Up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggh4x)
library(data.table)
library(ggrepel)

# Paths
om_scenario_path <- here("output", "OM_Scenarios") # path to OM folder
dir_out <- here("output", "Summary_Plots") # path to output folder
dir.create(dir_out)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
files <- list.files(fxn_path) # Load in all functions from the functions folder
for(i in 1:length(files)) source(here(fxn_path, files[i]))

# read in csvs
# Parameter and Time Series Summaries
AIC_df <- data.table::fread(here("output", "AIC_Convergence_Summary.csv")) %>%   # time series total error (converged runs only)
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
  # time series total error (converged runs only)
param_df <- read.csv(here("output", "Parameter_Summary.csv")) %>% # parameters
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
ts_all_df <- data.table::fread(here("output", "TimeSeries_Summary.csv")) %>%  # all time series data
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
ts_re_df <- read.csv(here("output", "TimeSeries_RE.csv")) %>% # time series relative error (converged runs only)
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
ts_te_df <- read.csv(here("output", "TimeSeries_TE.csv")) %>%  # time series total error (converged runs only)
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
ts_are_df <- read.csv(here("output", "TimeSeries_ARE.csv")) %>%  # time series abs relative error (converged runs only)
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
naa_df <- read.csv(here("output", "NAA_Summary.csv")) %>%  # numbers at age 
  filter(!str_detect(OM_Scenario, "Rev|Ext"))

# Selectivity and Comps
om_slx_df <- data.table::fread(here("output", "OM_Fish_Selex.csv")) %>% # OM Selectivity Values
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
pop_sel_om <- data.table::fread(here("output", "Pop_Selex_OM.csv")) %>%  # OM Pop'n Selectivity Values
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
pop_sel_em <- data.table::fread(here("output", "Pop_Selex_EM.csv")) %>%  # EM Pop'n Selectivity Values
  filter(!str_detect(OM_Scenario, "Rev|Ext"))


# Unique oms and other components
unique_oms <- unique(param_df$OM_Scenario) # unique oms
ts_pars <- unique(ts_re_df$par_name) # unique parameter names

# Timeblock sensitivity models
blk_sens_mods <- "Blk\\b|Blk_1|Blk_3|Blk_5|Blk_-1|Blk_-3|Blk_-5"
all_models <- "2Fl_LL|2Fl_LGam|1Fl_L_TI|1Fl_Gam_TI|1Fl_LL_Blk\\b|1Fl_LGam_Blk\\b|1Fl_L_RW_1.25|1Fl_Gam_RW_2.0|1Fl_L_SP_0.75|1Fl_Gam_SP_0.75"
gen_om <- c("Fast", "Slow") # general OMs

# Get plot order for EM models
plot_order <- c("2Fleet_Logist_Logist", 
                 "2Fleet_Logist_Gamma",  
                "1Fleet_Logist_TimeInvar", 
                "1Fleet_Gamma_TimeInvar", 
                "1Fleet_Logist_Logist_TimeBlk", 
                "1Fleet_Logist_Gamma_TimeBlk", 
                "1Fleet_Logist_RandomWalk", 
                "1Fleet_Gamma_RandomWalk", 
                "1Fleet_Logist_SemiPar", 
                "1Fleet_Gamma_SemiPar")

# Fast and slow plot order
fast_om_plot_order <- c("Fast-Logistic", "Fast-Gamma-Old", "Fast-Gamma-Young")
slow_om_plot_order <- c("Slow-Logistic", "Slow-Gamma-Old", "Slow-Gamma-Young")

# Sensitivity block models order
blk_sens_order <- c("1Fl_LL_Blk_-5", "1Fl_LL_Blk_-4", "1Fl_LL_Blk_-3", "1Fl_LL_Blk_-2", "1Fl_LL_Blk_-1", "1Fl_LL_Blk", 
                    "1Fl_LL_Blk_1", "1Fl_LL_Blk_2", "1Fl_LL_Blk_3", "1Fl_LL_Blk_4", "1Fl_LL_Blk_5", 
                    "1Fl_LGam_Blk_-5", "1Fl_LGam_Blk_-4", "1Fl_LGam_Blk_-3", "1Fl_LGam_Blk_-2", "1Fl_LGam_Blk_-1", 
                    "1Fl_LGam_Blk", "1Fl_LGam_Blk_1",  "1Fl_LGam_Blk_2","1Fl_LGam_Blk_3", "1Fl_LGam_Blk_4", "1Fl_LGam_Blk_5")

blk_sens_labels <- c("1Fl_Logist_Logist_TimeBlk_-5", "1Fl_Logist_Logist_TimeBlk_-4", 
                     "1Fl_Logist_Logist_TimeBlk_-3", "1Fl_Logist_Logist_TimeBlk_-2",
                     "1Fl_Logist_Logist_TimeBlk_-1", "1Fl_Logist_Logist_TimeBlk",
                     "1Fl_Logist_Logist_TimeBlk_1", "1Fl_Logist_Logist_TimeBlk_2", 
                     "1Fl_Logist_Logist_TimeBlk_3", "1Fl_Logist_Logist_TimeBlk_4",
                     "1Fl_Logist_Logist_TimeBlk_5",
                     "1Fl_Logist_Gamma_TimeBlk_-5", "1Fl_Logist_Gamma_TimeBlk_-4", "1Fl_Logist_Gamma_TimeBlk_-3", 
                     "1Fl_Logist_Gamma_TimeBlk_-2", "1Fl_Logist_Gamma_TimeBlk_-1", "1Fl_Logist_Gamma_TimeBlk",
                     "1Fl_Logist_Gamma_TimeBlk_1", "1Fl_Logist_Gamma_TimeBlk_2",
                     "1Fl_Logist_Gamma_TimeBlk_3", "1Fl_Logist_Gamma_TimeBlk_4","1Fl_Logist_Gamma_TimeBlk_5")


# Residual Munging --------------------------------------------------------

### Convergence Stuff -------------------------------------------------------

# Convergence for all models we want to look at
conv_stat <- AIC_df %>% 
  filter(conv == "Converged") %>% 
  drop_na() %>% # drop nas for those that did not converge
  filter(str_detect(EM_Scenario, all_models)) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), 
  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
     EM_Scenario = str_replace(EM_Scenario, "Gam", "G")) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(converged = sum(conv == "Converged")/200) %>% # divide by 200 = number of sims run
  ungroup() %>% 
  mutate(OM_Scenario = factor(OM_Scenario, 
                                labels = c(fast_om_plot_order,
                                           slow_om_plot_order)),
           EM_Scenario = factor(EM_Scenario,
                                labels = c("1Fleet_Gamma_RandomWalk",
                                           "1Fleet_Gamma_SemiPar",
                                           "1Fleet_Gamma_TimeInvar",
                                           "1Fleet_Logist_RandomWalk",
                                           "1Fleet_Logist_SemiPar",
                                           "1Fleet_Logist_TimeInvar",
                                           "1Fleet_Logist_Gamma_TimeBlk",
                                           "1Fleet_Logist_Logist_TimeBlk",
                                           "2Fleet_Logist_Gamma",
                                           "2Fleet_Logist_Logist"))) %>% 
  mutate(func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
         abbrev = case_when(
           str_detect(EM_Scenario, "2Fleet") ~ "2Fleet",
           str_detect(EM_Scenario, "TimeInvar") ~ "1Fleet-TimeInvar",
           str_detect(EM_Scenario, "Blk") ~ "1Fleet-Block",
           str_detect(EM_Scenario, "Random") ~ "1Fleet-RandWlkPar",
           str_detect(EM_Scenario, "Semi") ~ "1Fleet-SemiPar"
         ),
         abbrev = factor(abbrev, levels = rev(c("2Fleet", "1Fleet-TimeInvar",
                                                "1Fleet-Block", "1Fleet-RandWlkPar",
                                                "1Fleet-SemiPar"))))
# Set order for plot
order <- vector()
for(o in 1:length(plot_order)) {
  order[o] <- which(grepl(plot_order[o], x = unique(conv_stat$EM_Scenario)))
} # end o loop

# relevel factors
conv_stat$EM_Scenario <- factor(conv_stat$EM_Scenario, levels = c(unique(conv_stat$EM_Scenario)[order]))

# Do some summart stats
conv_stat %>% 
  filter(Dat_Qual == "High") %>% 
  group_by(EM_Scenario) %>% 
  summarize(mean = mean(converged))
  
  
### AIC Stuff ---------------------------------------------------------------
AIC_df <- AIC_df %>% 
  filter(conv == "Converged") %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = case_when(
    OM_Scenario == "Fast_LL" ~ "Fast-Logistic",
    OM_Scenario == "Fast_LG_Y" ~ "Fast-Gamma-Young",
    OM_Scenario == "Fast_LG_O" ~ "Fast-Gamma-Old",
    OM_Scenario == "Slow_LL" ~ "Slow-Logistic",
    OM_Scenario == "Slow_LG_Y" ~ "Slow-Gamma-Young",
    OM_Scenario == "Slow_LG_O" ~ "Slow-Gamma-Old"
  ),
  OM_Scenario = factor(OM_Scenario,  levels = c(fast_om_plot_order, slow_om_plot_order))) 

# Two fleet AIC df
twofleet_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "2Fl_LL|2Fl_LGam")) %>% 
  group_by(sim, OM_Scenario, time_comp, Dat_Qual) %>% 
  mutate(min = min(AIC), min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "2 Fleet",
         EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
         EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
         EM_Scenario = case_when(
           EM_Scenario == "2Fl_LL" ~ "2Fleet_Logist_Logist",
           EM_Scenario == "2Fl_LG" ~ "2Fleet_Logist_Gamma",
         )) %>% 
  mutate(func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
         abbrev = case_when(
           str_detect(EM_Scenario, "2Fleet") ~ "2Fleet",
           str_detect(EM_Scenario, "TimeInvar") ~ "1Fleet-TimeInvar",
           str_detect(EM_Scenario, "Blk") ~ "1Fleet-Block",
           str_detect(EM_Scenario, "Random") ~ "1Fleet-RandWlkPar",
           str_detect(EM_Scenario, "Semi") ~ "1Fleet-SemiPar"
         ),
         abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet-TimeInvar",
                                            "1Fleet-Block", "1Fleet-RandWlkPar",
                                            "1Fleet-SemiPar")))

# Get mean values here
twofleet_aic_vals = AIC_df %>% 
  filter(str_detect(EM_Scenario, "2Fl_")) %>% 
  mutate(Selex_Groups = case_when(
    str_detect(OM_Scenario, "Logist_Logist") ~ "LL",
    str_detect(OM_Scenario, "Logist_Gamma") ~ "LG"
  )) %>% 
  group_by(Selex_Groups, EM_Scenario) %>% 
  summarize(mean = mean(AIC)) %>% 
  pivot_wider(names_from = EM_Scenario, values_from = mean) %>% 
  mutate(abs_diff = abs(`2Fl_LGam` - `2Fl_LL`)) 


### 1 Fleet Models AIC -------------------------------------------------------
# filter to 1 fleet models
onefleet_aic <- AIC_df %>% 
  drop_na() %>% 
  filter(str_detect(EM_Scenario, "1Fl_"),
         str_detect(EM_Scenario, all_models)) %>% 
  group_by(sim, OM_Scenario, time_comp, Dat_Qual) %>% 
  mutate(min = min(AIC, na.rm = TRUE),
         min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC, na.rm = TRUE) / 200) %>% 
  mutate(fleet = "1 Fleet",
         EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0|_0.75"),
         EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
         EM_Scenario = case_when(
         EM_Scenario == "1Fl_L_SP" ~ "1Fleet_Logist_SemiPar",
         EM_Scenario == "1Fl_G_SP" ~ "1Fleet_Gamma_SemiPar",
         EM_Scenario == "1Fl_L_RW" ~ "1Fleet_Logist_RandomWalk",
         EM_Scenario == "1Fl_G_RW" ~ "1Fleet_Gamma_RandomWalk",
         EM_Scenario == "1Fl_LG_Blk" ~ "1Fleet_Logist_Gamma_TimeBlk",
         EM_Scenario == "1Fl_LL_Blk" ~ "1Fleet_Logist_Logist_TimeBlk",
         EM_Scenario == "1Fl_L_TI" ~ "1Fleet_Logist_TimeInvar",
         EM_Scenario == "1Fl_G_TI" ~ "1Fleet_Gamma_TimeInvar"
         ), EM_Scenario = factor(EM_Scenario, levels = rev(plot_order))) %>% 
  mutate(func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
         abbrev = case_when(
           str_detect(EM_Scenario, "2Fleet") ~ "2Fleet",
           str_detect(EM_Scenario, "TimeInvar") ~ "1Fleet-TimeInvar",
           str_detect(EM_Scenario, "Blk") ~ "1Fleet-Block",
           str_detect(EM_Scenario, "Random") ~ "1Fleet-RandWlkPar",
           str_detect(EM_Scenario, "Semi") ~ "1Fleet-SemiPar"
         ),
         abbrev = factor(abbrev, levels = rev(c("2Fleet", "1Fleet-TimeInvar",
                                                "1Fleet-Block", "1Fleet-RandWlkPar",
                                                "1Fleet-SemiPar"))))


### Time Block AIC Stuff ----------------------------------------------------

# Filter to 1 fleet time block models
onefleet_blk_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "LL_Blk|LGam_Blk"),
         str_detect(OM_Scenario, "Fast")) %>% 
  mutate(
    selex_form = case_when(
      str_detect(EM_Scenario, "LGam") ~ "Logistic_Gamma",
      str_detect(EM_Scenario, "LL") ~ "Logistic_Logistic",
    )) %>% group_by(sim, OM_Scenario, time_comp, Dat_Qual) %>% 
  mutate(min = min(AIC), min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "1 Fleet",
         EM_Scenario = factor(EM_Scenario, levels = blk_sens_order, labels = blk_sens_labels)) %>% 
  mutate(func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
         blk_years = case_when(
           str_detect(EM_Scenario, "-1") ~ "-1",
           str_detect(EM_Scenario, "-2") ~ "-2",
           str_detect(EM_Scenario, "-3") ~ "-3",
           str_detect(EM_Scenario, "-4") ~ "-4",
           str_detect(EM_Scenario, "-5") ~ "-5",
           str_detect(EM_Scenario, "Blk\\b") ~ "0",
           str_detect(EM_Scenario, "_1") ~ "1",
           str_detect(EM_Scenario, "_2") ~ "2",
           str_detect(EM_Scenario, "_3") ~ "3",
           str_detect(EM_Scenario, "_4") ~ "4",
           str_detect(EM_Scenario, "_5") ~ "5",
         ),
         blk_years = factor(blk_years, levels = c("-5","-4", "-3", "-2", "-1",
                                                  "0", "1", "2", "3", "4", "5")))

# get the best time-block models based on aic
best_fleetblk_models_df = onefleet_blk_aic %>% 
  group_by(OM_Scenario, time_comp) %>% 
  filter(n_minAIC == max(n_minAIC)) %>% 
  mutate(Model = "Best")

# select dataframe components
best_fleetblk_models = best_fleetblk_models_df %>% 
  select(OM_Scenario, EM_Scenario, Model, time_comp, Dat_Qual) %>% 
  mutate(OM_Scenario = factor(OM_Scenario, levels = c("Fast_Logist_Logist", 
                                                      "Fast_Logist_Gamma_Old", 
                                                      "Fast_Logist_Gamma_Young"),
                              labels = c(fast_om_plot_order)))

### ABC Parameters ----------------------------------------------------------

# Filter to relevant components for parameters
om_scenario_params <- param_df %>% filter(str_detect(EM_Scenario, all_models), 
                                          type %in% c("F_0.4", "ABC")) %>%
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = factor(OM_Scenario, labels = c("Fast-Gamma-Old",
                                               "Fast-Gamma-Young",
                                               "Fast-Logistic",
                                               "Slow-Gamma-Old",
                                               "Slow-Gamma-Young",
                                               "Slow-Logistic")),
         EM_Scenario = factor(EM_Scenario,
                              labels = c("1Fleet_Gamma_RandomWalk",
                                         "1Fleet_Gamma_SemiPar",
                                         "1Fleet_Gamma_TimeInvar",
                                         "1Fleet_Logist_RandomWalk",
                                         "1Fleet_Logist_SemiPar",
                                         "1Fleet_Logist_TimeInvar",
                                         "1Fleet_Logist_Gamma_TimeBlk",
                                         "1Fleet_Logist_Logist_TimeBlk",
                                         "2Fleet_Logist_Gamma",
                                         "2Fleet_Logist_Logist")))

# Point ranges for relative error and total error
pt_rg_re <- om_scenario_params %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, type, Dat_Qual) %>% 
  summarize(median = median(RE), 
            sd = sd(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 =  quantile(RE, 0.975))

# Clarify names
pt_rg_re <- pt_rg_re %>% 
  mutate(type = case_when(
    type == "ABC" ~ "ABC",
    type == "F_0.4" ~ "F40%"),
    EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
    OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)),
    EM_Scenario = factor(EM_Scenario, levels = rev(plot_order))) %>% 
  mutate(func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
         abbrev = case_when(
           str_detect(EM_Scenario, "2Fleet") ~ "2Fleet",
           str_detect(EM_Scenario, "TimeInvar") ~ "1Fleet-TimeInvar",
           str_detect(EM_Scenario, "Blk") ~ "1Fleet-Block",
           str_detect(EM_Scenario, "Random") ~ "1Fleet-RandWlkPar",
           str_detect(EM_Scenario, "Semi") ~ "1Fleet-SemiPar"
         ),
         abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet-TimeInvar",
                                                "1Fleet-Block", "1Fleet-RandWlkPar",
                                                "1Fleet-SemiPar")))


# Selectivity Stuff -------------------------------------------------------

plot_df = pop_sel_em %>% 
  filter(str_detect(EM_Scenario, all_models)) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), 
  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0|_0.75"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
  OM_Scenario = factor(OM_Scenario, labels = c("Fast-Gamma-Old",
                                               "Fast-Gamma-Young",
                                               "Fast-Logistic",
                                               "Slow-Gamma-Old",
                                               "Slow-Gamma-Young",
                                               "Slow-Logistic")),
  EM_Scenario = factor(EM_Scenario,
                       labels = c("1Fleet_Gamma_RandomWalk",
                                  "1Fleet_Gamma_SemiPar",
                                  "1Fleet_Gamma_TimeInvar",
                                  "1Fleet_Logist_RandomWalk",
                                  "1Fleet_Logist_SemiPar",
                                  "1Fleet_Logist_TimeInvar",
                                  "1Fleet_Logist_Gamma_TimeBlk",
                                  "1Fleet_Logist_Logist_TimeBlk",
                                  "2Fleet_Logist_Gamma",
                                  "2Fleet_Logist_Logist")),
  OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)),
  EM_Scenario = factor(EM_Scenario, levels = plot_order)) %>% 
  mutate(func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
         abbrev = case_when(
           str_detect(EM_Scenario, "2Fleet") ~ "2Fleet",
           str_detect(EM_Scenario, "TimeInvar") ~ "1Fleet-TimeInvar",
           str_detect(EM_Scenario, "Blk") ~ "1Fleet-Block",
           str_detect(EM_Scenario, "Random") ~ "1Fleet-RandWlkPar",
           str_detect(EM_Scenario, "Semi") ~ "1Fleet-SemiPar"
         ),
         abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet-TimeInvar",
                                                "1Fleet-Block", "1Fleet-RandWlkPar",
                                                "1Fleet-SemiPar")))

plot_sel_om_df = pop_sel_om %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), 
  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = factor(OM_Scenario, labels = c("Fast-Gamma-Old",
                                               "Fast-Gamma-Young",
                                               "Fast-Logistic",
                                               "Slow-Gamma-Old",
                                               "Slow-Gamma-Young",
                                               "Slow-Logistic")),
  OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)))

### Time Series Stuff -------------------------------------------------------

# Relative error of time series
ts_re_om <- ts_re_df %>% filter(str_detect(EM_Scenario, all_models)) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
  OM_Scenario = factor(OM_Scenario,
                       labels = c("Fast-Gamma-Old",
                                  "Fast-Gamma-Young",
                                  "Fast-Logistic",
                                  "Slow-Gamma-Old",
                                  "Slow-Gamma-Young",
                                  "Slow-Logistic")),
  EM_Scenario = factor(EM_Scenario,
                       labels = c("1Fleet_Gamma_RandomWalk",
                                  "1Fleet_Gamma_SemiPar",
                                  "1Fleet_Gamma_TimeInvar",
                                  "1Fleet_Logist_RandomWalk",
                                  "1Fleet_Logist_SemiPar",
                                  "1Fleet_Logist_TimeInvar",
                                  "1Fleet_Logist_Gamma_TimeBlk",
                                  "1Fleet_Logist_Logist_TimeBlk",
                                  "2Fleet_Logist_Gamma",
                                  "2Fleet_Logist_Logist")))

# Set order for plot
order <- vector()
for(o in 1:length(plot_order)) {
  order[o] <- which(grepl(plot_order[o], x = unique(ts_re_om$EM_Scenario)))
} # end o loop

# Now relevel factor for organizing plot
ts_re_om <- ts_re_om %>% 
  mutate(EM_Scenario = factor(EM_Scenario, levels = unique(EM_Scenario)[order]),
         time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")))


### Time Block Time Series --------------------------------------------------
# Relative error of time series
tb_ts_re_om <- ts_re_df %>% filter(time_comp %in% c("Terminal", "Fleet Trans End"),
                                   str_detect(OM_Scenario, "Fast"),
                                   str_detect(EM_Scenario, "Blk"),
                                   par_name %in% c("Spawning Stock Biomass",
                                                   "Total Fishing Mortality"),
                                   EM_Scenario %in% blk_sens_order) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = factor(OM_Scenario, levels = c("Fast_LL", "Fast_LG_O", "Fast_LG_Y"),
                       labels = c("Fast-Logistic", "Fast-Gamma-Old",
                                  "Fast-Gamma-Young")),
  EM_Scenario = factor(EM_Scenario, levels = blk_sens_order,
                       labels = blk_sens_labels),
  BreakPoint = 
    case_when(
      str_detect(EM_Scenario, "-5") ~ 20,
      str_detect(EM_Scenario, "-3") ~ 22,
      str_detect(EM_Scenario, "-1") ~ 24,
      str_detect(EM_Scenario, "TimeBlk\\b") ~ 25,
      str_detect(EM_Scenario, "_1") ~ 26,
      str_detect(EM_Scenario, "_3") ~ 28,
      str_detect(EM_Scenario, "_5") ~ 30,
    )) %>% 
  mutate(func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
         blk_years = case_when(
           str_detect(EM_Scenario, "-1") ~ "-1",
           str_detect(EM_Scenario, "-2") ~ "-2",
           str_detect(EM_Scenario, "-3") ~ "-3",
           str_detect(EM_Scenario, "-4") ~ "-4",
           str_detect(EM_Scenario, "-5") ~ "-5",
           str_detect(EM_Scenario, "Blk\\b") ~ "0",
           str_detect(EM_Scenario, "_1") ~ "1",
           str_detect(EM_Scenario, "_2") ~ "2",
           str_detect(EM_Scenario, "_3") ~ "3",
           str_detect(EM_Scenario, "_4") ~ "4",
           str_detect(EM_Scenario, "_5") ~ "5",
         ),
         blk_years = factor(blk_years, levels = c("-5","-4", "-3", "-2", "-1",
                                                  "0", "1", "2", "3", "4", "5")))

# Left Join "Best AIC" Models
tb_ts_re_om = tb_ts_re_om %>% left_join(best_fleetblk_models,
                                        by = c("OM_Scenario", "EM_Scenario", "time_comp", "Dat_Qual")) %>% 
  mutate(Model = ifelse(is.na(Model), "Not Best", Model ))


### Minimax Solution --------------------------------------------------------

# Do some munging and compute maximum ARE and minimum ARE
minmax_df = ts_all_df %>% 
  filter(conv == "Converged") %>% 
  mutate(
    are = abs((mle_val - t) / t),
    time_comp = case_when(
      str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
      str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
      str_detect(EM_Scenario, "Int") ~ "Fleet Intersect", # Fleet Transition Intersects
    ),
    EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'), # remove preceeding letters
    time_comp = factor(time_comp, levels = c("Terminal", "Fleet Trans End",
                                             "Fleet Intersect")),
    Dat_Qual = case_when(
      str_detect(OM_Scenario, "High") ~ 'High',
      str_detect(OM_Scenario, "Low") ~ 'Low'),
    OM_Scenario = str_remove(OM_Scenario, "_High|_Low")) %>% 
  data.frame() %>% 
  filter(str_detect(EM_Scenario, all_models), type == "Spawning Stock Biomass") %>% 
  data.frame() %>% 
  group_by(OM_Scenario, EM_Scenario, Dat_Qual, time_comp) %>% 
  summarize(MARE = median(are)) %>% 
  ungroup() %>% 
  group_by(EM_Scenario, Dat_Qual, time_comp) %>% 
  mutate(max_median = max(MARE))

  
# Do some residual munging with names
minmax_df = minmax_df %>% 
  mutate(
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
  EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0|_0.75"),
  OM_Scenario = factor(OM_Scenario, 
                       levels = c("Fast_LL", "Fast_LG_O", "Fast_LG_Y",
                                  "Slow_LL", "Slow_LG_O", "Slow_LG_Y"),
                       labels = c("Fast_Logist_Logist",
                                  "Fast_Logist_Gamma_Old",
                                  "Fast_Logist_Gamma_Young",
                                  "Slow_Logist_Logist",
                                  "Slow_Logist_Gamma_Old",
                                  "Slow_Logist_Gamma_Young")),
  EM_Scenario = factor(EM_Scenario,
                       levels = c( "1Fl_G_SP" , "1Fl_L_SP" , 
                                   "1Fl_G_RW" , "1Fl_L_RW" ,
                                  "1Fl_LG_Blk" , "1Fl_LL_Blk" ,
                                  "1Fl_G_TI", "1Fl_L_TI" , 
                                  "2Fl_LG" , "2Fl_LL"), # relevel
                       labels = c("1Fleet_Gamma_SemiPar",
                                  "1Fleet_Logist_SemiPar",
                                  "1Fleet_Gamma_RandomWalk",
                                  "1Fleet_Logist_RandomWalk",
                                  "1Fleet_Logist_Gamma_TimeBlk",
                                  "1Fleet_Logist_Logist_TimeBlk",
                                  "1Fleet_Gamma_TimeInvar",
                                  "1Fleet_Logist_TimeInvar",
                                  "2Fleet_Logist_Gamma",
                                  "2Fleet_Logist_Logist")),
  time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End", "Terminal"))) %>% 
  ungroup()

# Find the minimum maximum median value for each time period
time_minmax_medians <- minmax_df %>%
  group_by(time_comp, Dat_Qual) %>%
  summarize(time_minmax_medians = min(max_median),
            time_maxmax_medians = max(max_median)) %>% 
  ungroup()

# Find the global minimuim maximum median value
global_minmax_medians <- minmax_df %>%
  group_by(Dat_Qual) %>% 
  summarize(global_minmax_medians = min(max_median),
            global_maxmax_medians = max(max_median)) 

# Now left join this
minmax_df = minmax_df %>% 
  left_join(time_minmax_medians, by = c("time_comp", "Dat_Qual")) %>% 
  left_join(global_minmax_medians, by = c("Dat_Qual"))

# NAA Summary -------------------------------------------------------------
# Relative error of NAA series
naa_df <- naa_df %>% filter(str_detect(EM_Scenario, all_models)) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
  OM_Scenario = factor(OM_Scenario,
                       labels = c("Fast-Gamma-Old",
                                  "Fast-Gamma-Young",
                                  "Fast-Logistic",
                                  "Slow-Gamma-Old",
                                  "Slow-Gamma-Young",
                                  "Slow-Logistic")),
  EM_Scenario = factor(EM_Scenario,
                       labels = c("1Fleet_Gamma_RandomWalk",
                                  "1Fleet_Gamma_SemiPar",
                                  "1Fleet_Gamma_TimeInvar",
                                  "1Fleet_Logist_RandomWalk",
                                  "1Fleet_Logist_SemiPar",
                                  "1Fleet_Logist_TimeInvar",
                                  "1Fleet_Logist_Gamma_TimeBlk",
                                  "1Fleet_Logist_Logist_TimeBlk",
                                  "2Fleet_Logist_Gamma",
                                  "2Fleet_Logist_Logist")))

# Set order for plot
order <- vector()
for(o in 1:length(plot_order)) {
  order[o] <- which(grepl(plot_order[o], x = unique(naa_df$EM_Scenario)))
} # end o loop

# Now relevel factor for organizing plot
naa_df <- naa_df %>% 
  mutate(EM_Scenario = factor(EM_Scenario, levels = unique(EM_Scenario)[order]),
         time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")))


# Plots -------------------------------------------------------------------
# Figure 2 (RE SSB Fast) --------------------------------------------------

# Filter to fast scenario, and differentiate models
ssb_fast_re_om = ts_re_om %>% 
  filter(par_name == "Spawning Stock Biomass", str_detect(OM_Scenario, "Fast")) %>% 
  mutate( OM_Scenario = factor(OM_Scenario, levels = fast_om_plot_order),
          model_type = case_when( # differenitate model types
            str_detect(EM_Scenario, 'Random') ~ "1Fleet Random Walk",
            str_detect(EM_Scenario, 'Semi') ~ "1Fleet Semi-Parametric",
            str_detect(EM_Scenario, 'TimeInvar') ~ "1Fleet TimeInvar",
            str_detect(EM_Scenario, 'TimeBlk') ~ "1Fleet TimeBlock",
            str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
          ),
          func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
          abbrev = case_when(
                str_detect(model_type, "2Fleet") ~ "2Fleet",
                str_detect(model_type, "TimeInvar") ~ "1Fleet-TimeInvar",
                str_detect(model_type, "Block") ~ "1Fleet-Block",
                str_detect(model_type, "Random") ~ "1Fleet-RandWlkPar",
                str_detect(model_type, "Semi") ~ "1Fleet-SemiPar"
              ),
          abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet-TimeInvar",
                                             "1Fleet-Block", "1Fleet-RandWlkPar",
                                             "1Fleet-SemiPar"))) %>%  # functional form
  filter(Dat_Qual == "High")


pdf(here("figs", "Manuscript_Figures_v2", "Fig2_SSB.pdf"), width = 19, height = 21)
ggplot() + 
  geom_hline(yintercept = 0, col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  geom_ribbon(ssb_fast_re_om %>% filter(time_comp == "Terminal"),
              mapping = aes(x = year, y = median, ymin = lwr_95, ymax = upr_95, fill = func), 
              color = NA, alpha = 0.3) +
  geom_line(ssb_fast_re_om, mapping = aes(x = year, y = median, size = time_comp,
                                          color = func, linetype = time_comp), alpha = 1) +
  facet_grid(abbrev ~ OM_Scenario) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  scale_linetype_manual(values = c("Fleet Intersect" = 3, "Fleet Trans End" = 2, "Terminal" = 1)) +
  scale_size_manual(values = c("Fleet Intersect" = 2.8, "Fleet Trans End" = 1.35, "Terminal" = 1.5)) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_matt() +
  theme(legend.key.width = unit(1.3,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in Spawning Stock Biomass", fill = "Functional Form",
       color = "Functional Form", lty = "Assessment Period", alpha = "Assessment Period",
       size = "Assessment Period", title = "Fast Transition")
dev.off()

# Figure 3 (RE SSB Slow) --------------------------------------------------

# Filter to fast scenario, and differentiate models
ssb_slow_re_om = ts_re_om %>% 
  filter(par_name == "Spawning Stock Biomass", str_detect(OM_Scenario, "Slow")) %>% 
  mutate( OM_Scenario = factor(OM_Scenario, levels = slow_om_plot_order),
          model_type = case_when( # differenitate model types
            str_detect(EM_Scenario, 'Random') ~ "1Fleet Random Walk",
            str_detect(EM_Scenario, 'Semi') ~ "1Fleet Semi-Parametric",
            str_detect(EM_Scenario, 'TimeInvar') ~ "1Fleet TimeInvar",
            str_detect(EM_Scenario, 'TimeBlk') ~ "1Fleet TimeBlock",
            str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
          ),
          func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
          abbrev = case_when(
            str_detect(model_type, "2Fleet") ~ "2Fleet",
            str_detect(model_type, "TimeInvar") ~ "1Fleet-TimeInvar",
            str_detect(model_type, "Block") ~ "1Fleet-Block",
            str_detect(model_type, "Random") ~ "1Fleet-RandWlkPar",
            str_detect(model_type, "Semi") ~ "1Fleet-SemiPar"
          ),
          abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet-TimeInvar",
                                             "1Fleet-Block", "1Fleet-RandWlkPar",
                                             "1Fleet-SemiPar"))) %>%  # functional form
  filter(Dat_Qual == "High")

pdf(here("figs", "Manuscript_Figures_v2", "Fig3_SSB.pdf"), width = 19, height = 21)
ggplot() + 
  geom_hline(yintercept = 0, col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  geom_ribbon(ssb_slow_re_om %>% filter(time_comp == "Terminal"),
              mapping = aes(x = year, y = median, ymin = lwr_95, ymax = upr_95, fill = func), 
              color = NA, alpha = 0.3) +
  geom_line(ssb_slow_re_om, mapping = aes(x = year, y = median, size = time_comp,
                                          color = func, linetype = time_comp), alpha = 1) +
  facet_grid(abbrev ~ OM_Scenario) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  scale_linetype_manual(values = c("Fleet Intersect" = 3, "Fleet Trans End" = 2, "Terminal" = 1)) +
  scale_size_manual(values = c("Fleet Intersect" = 2.8, "Fleet Trans End" = 1.35, "Terminal" = 1.5)) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_matt() +
  theme(legend.key.width = unit(1.3,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in Spawning Stock Biomass", fill = "Functional Form",
       color = "Functional Form", lty = "Assessment Period", alpha = "Assessment Period",
       size = "Assessment Period", title = "Slow Transition")
dev.off()


# Figure 4 (RE ABC Fast) --------------------------------------------------

pdf(here("figs", "Manuscript_Figures_v2", "Fig4_ABC.pdf"), width = 15, height = 15)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, "Fast"), type == "ABC"),
         aes(x = factor(abbrev), y = median, color = func, fill = func,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    facet_grid(time_comp~OM_Scenario, scales = "free_x") +
    scale_color_manual(values = c("#E69F00", "#0072B2")) +
    scale_fill_manual(values = c("#E69F00", "#0072B2")) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(x = "Estimation Models", y = "Relative Error in ABC", 
         fill = "Functional Form", color = "Functional Form", title = "Fast Transition") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text= element_text(size = 18),
          strip.text = element_text(size = 17),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20)) +
    coord_cartesian(ylim = c(-0.5, 0.5)))
dev.off()

# Figure 5 (RE ABC Slow) --------------------------------------------------
pdf(here("figs", "Manuscript_Figures_v2", "Fig5_ABC.pdf"), width = 15, height = 15)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, "Slow"), type == "ABC"),
         aes(x = factor(abbrev), y = median, color = func, fill = func,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    facet_grid(time_comp~OM_Scenario, scales = "free_x") +
    scale_color_manual(values = c("#E69F00", "#0072B2")) +
    scale_fill_manual(values = c("#E69F00", "#0072B2")) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(x = "Estimation Models", y = "Relative Error in ABC", 
         fill = "Functional Form", color = "Functional Form", title = "Slow Transition") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text= element_text(size = 18),
          strip.text = element_text(size = 17),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20)) +
    coord_cartesian(ylim = c(-0.5, 0.5)))
dev.off()


# Figure 6 (Quadratic Bias) -----------------------------------------------

# Get terminal NAA
sub_naa_df <- naa_df %>% filter(OM_Scenario == "Fast-Logistic",
                  str_detect(EM_Scenario, "TimeInvar"),
                  time_comp == "Terminal",
                  Sex == "Female") %>% 
  mutate(EM_Scenario = factor(EM_Scenario, labels = c("1Fleet-TimeInvar-Gamma",
                                                      "1Fleet-TimeInvar-Logist")))

# Filter to fast scenario, and differentiate models
f_fast_re_om = ts_re_om %>% 
  filter(par_name == "Total Fishing Mortality", str_detect(OM_Scenario, "Fast")) %>% 
  mutate( OM_Scenario = factor(OM_Scenario, levels = fast_om_plot_order),
          model_type = case_when( # differenitate model types
            str_detect(EM_Scenario, 'Random') ~ "1Fleet Random Walk",
            str_detect(EM_Scenario, 'Semi') ~ "1Fleet Semi-Parametric",
            str_detect(EM_Scenario, 'TimeInvar') ~ "1Fleet TimeInvar",
            str_detect(EM_Scenario, 'TimeBlk') ~ "1Fleet TimeBlock",
            str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
          ),
          func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
          func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
          abbrev = case_when(
            str_detect(model_type, "2Fleet") ~ "2Fleet",
            str_detect(model_type, "TimeInvar") ~ "1Fleet-TimeInvar",
            str_detect(model_type, "Block") ~ "1Fleet-Block",
            str_detect(model_type, "Random") ~ "1Fleet-RandWlkPar",
            str_detect(model_type, "Semi") ~ "1Fleet-SemiPar"
          ),
          abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet-TimeInvar",
                                             "1Fleet-Block", "1Fleet-RandWlkPar",
                                             "1Fleet-SemiPar"))) %>%  # functional form
  filter(Dat_Qual == "High")

# Get OM Age Structure
load(here(om_scenario_path, "Fast_LL_High", paste("Fast_LL_High", ".RData", sep = "")))

# Get Numbers at age
naa_om <- reshape2::melt(oms$N_at_age)
colnames(naa_om) <- c("Year", "Age", "Sex", "Sim", "Value")

# parse out numbers
naa_om <- naa_om %>% 
  mutate(Year = parse_number(as.character(Year)),
         Age = parse_number(as.character(Age)),
         Sex = parse_number(as.character(Sex)),
         Sim = parse_number(as.character(Sim)))

# get median here
med_naa_df <- naa_om %>% 
  group_by(Age, Year, Sex) %>% 
  summarize(median = median(Value))

### Numbers at age Females (Fast-Logistic) ------------------------------------------

naa_plot1 <- ggplot(sub_naa_df %>% filter(Age %in% c(6,26)) %>% 
                      mutate(Age = factor(Age, labels = c("Age 6 (Young)",
                                                          "Age 26 (Old)"))),
       aes(x = Year, y = Median_RE, fill = EM_Scenario)) +
  geom_rect(xmin = -Inf, xmax = 25, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.05) +
  geom_rect(xmin = 25, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "green3", alpha = 0.002) +
  geom_line(linewidth = 2, alpha = 1, aes(color = EM_Scenario)) +
  geom_ribbon(aes(ymin = Lwr_95, ymax = Upr_95), alpha = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black",
             lty = 2, linewidth = 0.5, alpha = 1) +
  facet_wrap(~Age, scales = "free_y", nrow = 2) +
  labs(x = "Year", y = "Relative Error in NAA (Females)", 
       fill = "Estimation Models", color = "Estimation Models") +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_matt() +
  theme(legend.position = "none", 
        title = element_text(size = 25),
        axis.title = element_text(size = 20),
        axis.text= element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 23), 
        panel.grid = element_blank(),
        axis.title.x = element_text(vjust = 5))

png(here("figs", "Manuscript_Figures_v2", "NAA1_EM.png"), width = 800, height = 400)
naa_plot1
dev.off()

selex_plot1 <- ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                               str_detect(OM_Scenario, "Fast-Logistic"),
                               str_detect(EM_Scenario, "TimeInvar")) %>% 
              mutate(EM_Scenario = factor(EM_Scenario, labels = c("1Fleet-TimeInvar-Logist",
                                                    "1Fleet-TimeInvar-Gamma"))),
            mapping = aes(x = Age, y = Median_Selex, color = EM_Scenario), lwd = 2, alpha = 1) +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast-Logistic"),
                                 str_detect(EM_Scenario, "TimeInvar")) %>% 
                mutate(EM_Scenario = factor(EM_Scenario, labels = c("1Fleet-TimeInvar-Logist",
                                                      "1Fleet-TimeInvar-Gamma"))),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95, fill = EM_Scenario),
              alpha = 0.25) +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", 
                                      time_comp %in% c( "First Year"),
                                      str_detect(OM_Scenario, "Fast-Logistic")),
            mapping = aes(x = Age, y = Selex), color = "black", 
            lwd = 1.5, alpha = 1, lty = 1) +
  scale_color_manual(values = c("#0072B2", "#E69F00")) +
  scale_fill_manual(values = c( "#0072B2", "#E69F00")) +
  labs(x = "Age", y = "Selectivity (Females)", lty = "Year", title = "Pre-Transition",
       color = "Estimation Model", fill = "Estimation Model") +
  theme_matt() +
  theme(legend.position = "top", 
        plot.title=element_text(hjust= 0.5),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(vjust = 5)) 


selex_plot2 <- ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                               str_detect(OM_Scenario, "Fast-Logistic"),
                               str_detect(EM_Scenario, "TimeInvar")),
            mapping = aes(x = Age, y = Median_Selex, color = func), lwd = 2, alpha = 1) +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast-Logistic"),
                                 str_detect(EM_Scenario, "TimeInvar")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95, fill = func),
              alpha = 0.25) +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", 
                                      time_comp %in% c( "Terminal"),
                                      str_detect(OM_Scenario, "Fast-Logistic")),
            mapping = aes(x = Age, y = Selex), color = "black", 
            lwd = 1.5, alpha = 1, lty = 2) +
  scale_color_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  labs(x = "Age", y = "Selectivity (Females)", lty = "Year", title = "Post-Transition") +
  theme_matt() +
  theme(legend.position = c(0.7, 0.15), 
        plot.title=element_text(hjust= 0.5),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(vjust = 5)) 

png(here("figs", "Manuscript_Figures_v2", "Selex_PrePost.png"), width = 850)
ggpubr::ggarrange(selex_plot1, selex_plot2, common.legend = TRUE)
dev.off()

# get f plot
# f_plot <- ggplot(f_fast_re_om %>% filter(time_comp == "Terminal",
#                                          str_detect(OM_Scenario, "Fast-Logistic"),
#                                          str_detect(EM_Scenario, "TimeInvar")),
#                  mapping = aes(x = year, y = median, ymin = lwr_95, ymax = upr_95, fill = func, color = func)) + 
#   geom_rect(xmin = -Inf, xmax = 25, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.05, color = NA) +
#   geom_rect(xmin = 25, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "green3", alpha = 0.002, color = NA) +
#   geom_hline(yintercept = 0, col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
#   geom_ribbon(color = NA, alpha = 0.3) +
#   geom_line(alpha = 1, lwd = 2) +
#   coord_cartesian(ylim = c(-0.5, 0.5)) +
#   scale_color_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
#   scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
#   theme_matt() +
#   theme(legend.key.width = unit(1.3,"cm")) +
#   theme(legend.position = "top",
#         axis.text = element_text(size = 20), panel.grid = element_blank(),
#         axis.title.x = element_text(vjust = 5),
#         axis.title = element_text(size = 20)) +
#   labs(x = "Year", y = "Relative Error in Total Fishing Mortality")
# 
# png(here("figs", "Manuscript_Figures_v2", "F_Quad.png"))
# f_plot
# dev.off()

# naa_om_plot <- ggplot(med_naa_df %>% filter(Sex == 1,
#                              Age %in% c(3, 6, 9, 12, 21, 24, 27, 30)), 
#        aes(x = Year, y = median, fill = factor(Age))) +
#   geom_col(color = "black", alpha = 0.75) +
#   theme_matt() +
#   ggthemes::scale_fill_colorblind() +
#   theme(legend.key.width = unit(1.3,"cm"),
#         axis.text.y = element_text(angle = 90)) +
#   guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
#   theme(legend.position = c(0.65, 0.8),
#         axis.text = element_text(size = 23)) +
#   labs(x = "Year", y = "True Median NAA (Females)", fill = "Age", title = "E)")

# png(here("figs", "Manuscript_Figures_v2", "NAA_OM.png"))
# naa_om_plot
# dev.off()

# ssb plot
# ssb_plot <- ggplot() +
#   geom_hline(yintercept = 0, col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
#   geom_ribbon(ssb_fast_re_om %>% filter(time_comp == "Terminal",
#                                         str_detect(EM_Scenario, "TimeInvar"),
#                                         str_detect(OM_Scenario, "Fast-Logistic")),
#               mapping = aes(x = year, y = median, ymin = lwr_95, ymax = upr_95, fill = EM_Scenario),
#               color = NA, alpha = 0.3) +
#   geom_line(ssb_fast_re_om %>% filter(time_comp == "Terminal",
#                                       str_detect(EM_Scenario, "TimeInvar"),
#                                       str_detect(OM_Scenario, "Fast-Logistic")),
#             mapping = aes(x = year, y = median, color = EM_Scenario), alpha = 1, size = 2) +
#   coord_cartesian(ylim = c(-0.5, 0.5)) +
#   scale_color_manual(values = c("#E69F00", "#0072B2")) +
#   scale_fill_manual(values = c("#E69F00", "#0072B2")) +
#   theme_matt() +
#   theme(legend.key.width = unit(1.3,"cm"),
#         axis.text.y = element_text(angle = 90),
#         axis.text = element_text(size = 23)) +
#   theme(legend.position = "none") +
#   labs(x = "Year", y = "Relative Error in Spawning Stock Biomass")
# 
# png(here("figs", "Manuscript_Figures_v2", "SSB_Quad.png"))
# ssb_plot
# dev.off()

# Combined SSB and F plot
f_sub <- f_fast_re_om %>% filter(time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast-Logistic"),
                                 str_detect(EM_Scenario, "TimeInvar")) %>% 
  mutate(Type = "F")

ssb_sub <- ssb_fast_re_om %>% filter(time_comp == "Terminal",
                                        str_detect(EM_Scenario, "TimeInvar"),
                                        str_detect(OM_Scenario, "Fast-Logistic")) %>% 
  mutate(Type = "SSB")

f_ssb_sub <- rbind(f_sub, ssb_sub)

f_ssb_plot <- ggplot(f_ssb_sub, mapping = aes(x = year, y = median, ymin = lwr_95, 
                                ymax = upr_95, color = EM_Scenario,
                                fill = EM_Scenario, lty = Type)) +
  geom_rect(xmin = -Inf, xmax = 25, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.05, color = NA) +
  geom_rect(xmin = 25, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "green3", alpha = 0.002, color = NA) +
  geom_hline(yintercept = 0, col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line(alpha = 1, lwd = 2) +
  scale_color_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  theme_matt() +
  theme(legend.position = c(0.875, 0.925),
        axis.text = element_text(size = 20), panel.grid = element_blank(),
        axis.title.x = element_text(vjust = 5),
        axis.title = element_text(size = 20),
        legend.key.width = unit(1.3,"cm"),
        legend.background = element_blank()) +
  labs(x = "Year", y = "Relative Error", lty = "")

png(here("figs", "Manuscript_Figures_v2", "F_SSB.png"), width = 800, height = 400)
f_ssb_plot
dev.off()


# f_ssb_plot <- ggpubr::ggarrange(ssb_plot, f_plot, nrow = 1)
# sel_naa_plot <- ggpubr::ggarrange(selex_plot, naa_om_plot)
# comb_f_ssb_sel_naa_plot <- ggpubr::ggarrange(f_ssb_plot, sel_naa_plot, nrow = 2)
# pdf(here("figs", "Manuscript_Figures_v2", "Fig6_Quad.pdf"), width = 17, height = 25)
# ggpubr::ggarrange(naa_plot, comb_f_ssb_sel_naa_plot, 
#                   ncol= 1, widths = c(0.5, 0.25, 0.25), common.legend = TRUE)
# dev.off()

# Figure S1 (Convergence) -------------------------------------------------
pdf(here("figs", "Manuscript_Figures_v2", "FigS1_Convergence.pdf"), width = 25, height = 15)
ggplot(conv_stat %>% 
         mutate(OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)),
                EM_Scenario = forcats::fct_rev(abbrev)) %>% 
                              filter(Dat_Qual == "High"),
       mapping = aes(x = abbrev, y = converged * 100, group = func, fill = func))  +
  geom_col(position = position_dodge2(width = 1), 
           alpha = 0.7, color = "black") +
  geom_hline(aes(yintercept = 90), col = "black", lty = 2, size = 0.7, alpha = 1) +
  geom_hline(aes(yintercept = 50), col = "black", lty = 2, size = 0.7, alpha = 1) +
  facet_grid(time_comp~OM_Scenario,  scales = "free_x") +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_matt() + 
  labs(fill = "Functional Form", y = "Convergence Rate (%)", x = "Estimation Models") +
  theme(legend.position = "top") +
  coord_flip()
dev.off()


# Figure S2 and S3 (AIC) --------------------------------------------------

pdf(here("figs", "Manuscript_Figures_v2", "FigS2_AIC_2Fl.pdf"), width = 17, height = 10)
ggplot(twofleet_aic %>% filter(Dat_Qual == "High") %>% 
         mutate(OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order))), 
       aes(x = abbrev, y = round(n_minAIC, 2), fill = func)) +
  geom_col(position = position_dodge2(width = 0.1), 
           alpha = 0.7, color = "black") +  
  facet_grid(time_comp~OM_Scenario) +
  geom_hline(aes(yintercept = 0.9), col = "black", lty = 2, size = 0.7, alpha = 1) +
  geom_hline(aes(yintercept = 0.5), col = "black", lty = 2, size = 0.7, alpha = 1) +
  coord_flip() +
  labs(x = "Operating Models", y = "Proportion of models with lowest AIC", fill = "Functional Form of 2Fleet Models") +
  theme_test() +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme(legend.position = "top",
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text= element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.key.width = unit(1, "cm"))
dev.off()

pdf(here("figs", "Manuscript_Figures_v2", "FigS3_AIC_1Fl.pdf"), width = 20, height = 13)
ggplot(onefleet_aic %>% filter(Dat_Qual == "High") %>% 
         mutate(OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order))), 
       aes(x = abbrev, y = round(n_minAIC, 2), fill = func)) +
  geom_col(position = position_dodge2(width = 0.1), 
           alpha = 0.7, color = "black") +  
  facet_grid(time_comp~OM_Scenario) +
  coord_flip() +
  geom_hline(aes(yintercept = 0.9), col = "black", lty = 2, size = 0.7, alpha = 1) +
  geom_hline(aes(yintercept = 0.5), col = "black", lty = 2, size = 0.7, alpha = 1) +
  labs(x = "Estimation Models", y = "Proportion of models with lowest AIC", fill = "Functional Form of 1Fleet Models") +
  theme_test() +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme(legend.position = "top",
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text= element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.key.width = unit(1, "cm"))
dev.off()

# Figure S4 (RE Total F Fast) --------------------------------------------------
# 
# # Filter to fast scenario, and differentiate models
# f_fast_re_om = ts_re_om %>% 
#   filter(par_name == "Total Fishing Mortality", str_detect(OM_Scenario, "Fast")) %>% 
#   mutate( OM_Scenario = factor(OM_Scenario, levels = fast_om_plot_order),
#           model_type = case_when( # differenitate model types
#             str_detect(EM_Scenario, 'Random') ~ "1Fleet Random Walk",
#             str_detect(EM_Scenario, 'Semi') ~ "1Fleet Semi-Parametric",
#             str_detect(EM_Scenario, 'TimeInvar') ~ "1Fleet TimeInvar",
#             str_detect(EM_Scenario, 'TimeBlk') ~ "1Fleet TimeBlock",
#             str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
#           ),
#           func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
#           func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
#           abbrev = case_when(
#             str_detect(model_type, "2Fleet") ~ "2Fleet",
#             str_detect(model_type, "TimeInvar") ~ "1Fleet-TimeInvar",
#             str_detect(model_type, "Block") ~ "1Fleet-Block",
#             str_detect(model_type, "Random") ~ "1Fleet-RandWlkPar",
#             str_detect(model_type, "Semi") ~ "1Fleet-SemiPar"
#           ),
#           abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet-TimeInvar",
#                                              "1Fleet-Block", "1Fleet-RandWlkPar",
#                                              "1Fleet-SemiPar"))) %>%  # functional form
#   filter(Dat_Qual == "High")
# 
# 
# pdf(here("figs", "Manuscript_Figures_v2", "FigS4_F.pdf"), width = 19, height = 21)
# ggplot() + 
#   geom_hline(yintercept = 0, col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
#   geom_ribbon(f_fast_re_om %>% filter(time_comp == "Terminal"),
#               mapping = aes(x = year, y = median, ymin = lwr_95, ymax = upr_95, fill = func), 
#               color = NA, alpha = 0.3) +
#   geom_line(f_fast_re_om, mapping = aes(x = year, y = median, size = time_comp,
#                                           color = func, linetype = time_comp), alpha = 1) +
#   facet_grid(abbrev ~ OM_Scenario) +
#   coord_cartesian(ylim = c(-0.5, 0.5)) +
#   scale_linetype_manual(values = c("Fleet Intersect" = 3, "Fleet Trans End" = 2, "Terminal" = 1)) +
#   scale_size_manual(values = c("Fleet Intersect" = 2.8, "Fleet Trans End" = 1.35, "Terminal" = 1.5)) +
#   scale_color_manual(values = c("#E69F00", "#0072B2")) +
#   scale_fill_manual(values = c("#E69F00", "#0072B2")) +
#   theme_matt() +
#   theme(legend.key.width = unit(1.3,"cm")) +
#   theme(legend.position = "top") +
#   labs(x = "Year", y = "Relative Error in Total Fishing Mortality", fill = "Functional Form",
#        color = "Functional Form", lty = "Assessment Period", alpha = "Assessment Period",
#        size = "Assessment Period", title = "Fast Transition")
# dev.off()
# 
# 
# # Figure S5 (RE Total F Slow) --------------------------------------------------
# 
# # Filter to fast scenario, and differentiate models
# f_slow_re_om = ts_re_om %>% 
#   filter(par_name == "Total Fishing Mortality", str_detect(OM_Scenario, "Slow")) %>% 
#   mutate( OM_Scenario = factor(OM_Scenario, levels = slow_om_plot_order),
#           model_type = case_when( # differenitate model types
#             str_detect(EM_Scenario, 'Random') ~ "1Fleet Random Walk",
#             str_detect(EM_Scenario, 'Semi') ~ "1Fleet Semi-Parametric",
#             str_detect(EM_Scenario, 'TimeInvar') ~ "1Fleet TimeInvar",
#             str_detect(EM_Scenario, 'TimeBlk') ~ "1Fleet TimeBlock",
#             str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
#           ),
#           func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
#           func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logistic"),
#           abbrev = case_when(
#             str_detect(model_type, "2Fleet") ~ "2Fleet",
#             str_detect(model_type, "TimeInvar") ~ "1Fleet-TimeInvar",
#             str_detect(model_type, "Block") ~ "1Fleet-Block",
#             str_detect(model_type, "Random") ~ "1Fleet-RandWlkPar",
#             str_detect(model_type, "Semi") ~ "1Fleet-SemiPar"
#           ),
#           abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet-TimeInvar",
#                                              "1Fleet-Block", "1Fleet-RandWlkPar",
#                                              "1Fleet-SemiPar"))) %>%  # functional form
#   filter(Dat_Qual == "High")
# 
# 
# pdf(here("figs", "Manuscript_Figures_v2", "FigS5_F.pdf"), width = 19, height = 21)
# ggplot() + 
#   geom_hline(yintercept = 0, col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
#   geom_ribbon(f_slow_re_om %>% filter(time_comp == "Terminal"),
#               mapping = aes(x = year, y = median, ymin = lwr_95, ymax = upr_95, fill = func), 
#               color = NA, alpha = 0.3) +
#   geom_line(f_slow_re_om, mapping = aes(x = year, y = median, size = time_comp,
#                                         color = func, linetype = time_comp), alpha = 1) +
#   facet_grid(abbrev ~ OM_Scenario) +
#   coord_cartesian(ylim = c(-0.5, 0.5)) +
#   scale_linetype_manual(values = c("Fleet Intersect" = 3, "Fleet Trans End" = 2, "Terminal" = 1)) +
#   scale_size_manual(values = c("Fleet Intersect" = 2.8, "Fleet Trans End" = 1.35, "Terminal" = 1.5)) +
#   scale_color_manual(values = c("#E69F00", "#0072B2")) +
#   scale_fill_manual(values = c("#E69F00", "#0072B2")) +
#   theme_matt() +
#   theme(legend.key.width = unit(1.3,"cm")) +
#   theme(legend.position = "top") +
#   labs(x = "Year", y = "Relative Error in Total Fishing Mortality", fill = "Functional Form",
#        color = "Functional Form", lty = "Assessment Period", alpha = "Assessment Period",
#        size = "Assessment Period", title = "Slow Transition")
# dev.off()

# Figure S6 (Time Block AIC) ----------------------------------------------

pdf(here("figs", "Manuscript_Figures_v2", "FigS6_BlkAIC.pdf"), width = 15, height = 10)
ggplot(onefleet_blk_aic %>% filter(time_comp %in% c("Terminal", "Fleet Trans End"), Dat_Qual == "High"),
       aes(x = blk_years, y = n_minAIC, fill = func)) +
  geom_col(position = position_dodge(width = 1), color = "black") +
  ggh4x::facet_nested(time_comp~OM_Scenario, scales = "free") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  ylim(0,1) +
  geom_hline(aes(yintercept = 0.9), col = "black", lty = 2, size = 0.7, alpha = 1) +
  geom_hline(aes(yintercept = 0.5), col = "black", lty = 2, size = 0.7, alpha = 1) +
  labs(x = "Years since fleet structure change (Time Block specification)",
       y = "Proportion of models with lowest AIC",
       fill = "Functional Form of 1Fleet Time Block Models") +
  theme_test() +
  theme(legend.position = "top", 
        title = element_text(size = 20),
        axis.title = element_text(size = 17),
        axis.text= element_text(size = 15),
        strip.text = element_text(size = 17),
        legend.text = element_text(size = 15)) 
dev.off()


# Figure S7 and S8 (Time Block SSB) ----------------------------------------------

(a = ggplot(tb_ts_re_om %>% filter(Dat_Qual == "High", par_name == "Spawning Stock Biomass", 
                              str_detect(EM_Scenario, "Logist_Logist"),
                              str_detect(EM_Scenario, "Blk\\b|_1|_3|_5|-1|-3|-5")), 
            aes(x = year, y = median, fill = EM_Scenario))  +
  geom_line(linewidth = 1.5, alpha = 0.85, aes(color = blk_years)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  coord_cartesian(ylim = c(-0.4,0.4)) +
  scale_linetype_manual(values = c(2,1)) +
  guides(linetype = "none") +
  facet_grid(time_comp~OM_Scenario) + 
  labs(x = "Year", y = "Relative Error in Spawning Biomass", 
       color = "Years since fleet structure change (Time Block specification)", title = "Fast Transition") +
  theme_matt() +
  theme(legend.position = "top", 
          title = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 17),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 17)))


(b = ggplot(tb_ts_re_om %>% filter(Dat_Qual == "High", par_name == "Spawning Stock Biomass", 
                                   str_detect(EM_Scenario, "Logist_Gamma"),
                                   str_detect(EM_Scenario, "Blk\\b|_1|_3|_5|-1|-3|-5")), 
            aes(x = year, y = median, fill = EM_Scenario))  +
    geom_line(linewidth = 1.5, alpha = 0.85, aes(color = blk_years)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    coord_cartesian(ylim = c(-0.4,0.4)) +
    scale_linetype_manual(values = c(2,1)) +
    guides(linetype = "none") +
    facet_grid(time_comp~OM_Scenario) + 
    labs(x = "Year", y = "Relative Error in Spawning Biomass", 
         color = "Years since fleet structure change (Time Block specification)", title = "Fast Transition") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 17),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 17)))
pdf(here("figs", "Manuscript_Figures_v2", "FigS7_BlkSSB_LL.pdf"), width = 15, height = 12)
a
dev.off()

pdf(here("figs", "Manuscript_Figures_v2", "FigS8_BlkSSB_LG.pdf"), width = 15, height = 12)
b
dev.off()


# Figure S9 (Selectivity Fast Females) --------------------------------------------

pdf(here("figs", "Manuscript_Figures_v2", "FigS9_FastSelex_F.pdf"), width = 12, height = 15)
ggplot() +
    geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast")),
              mapping = aes(x = Age, y = Median_Selex, color = func), lwd = 1, alpha = 1) +
    geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                   str_detect(OM_Scenario, "Fast")),
                mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95, fill = func), alpha = 0.25) +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", 
                                      time_comp %in% c("Terminal", "First Year"),
                                      str_detect(OM_Scenario, "Fast")),
            mapping = aes(x = Age, y = Selex,lty = time_comp), color = "black", 
            lwd = 1, alpha = 1) +
    facet_grid(abbrev ~ OM_Scenario) +
    scale_color_manual(values = c("#E69F00", "#0072B2")) +
    scale_fill_manual(values = c("#E69F00", "#0072B2")) +
    labs(x = "Age", y = "Selectivity (Females)", 
         fill = "Functional Form", color = "Functional Form",
         lty = "Year", title = "Fast Transition") +
    theme_matt() +
    coord_cartesian(ylim = c(0, 2.5))
dev.off()


# Figure S10 (Selectivity Fast Males) --------------------------------------------

pdf(here("figs", "Manuscript_Figures_v2", "FigS10_FastSelex_M.pdf"), width = 12, height = 15)
ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                               str_detect(OM_Scenario, "Fast")),
            mapping = aes(x = Age, y = Median_Selex, color = func), lwd = 1, alpha = 1) +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95, fill = func), alpha = 0.25) +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Male", 
                                      time_comp %in% c("Terminal", "First Year"),
                                      str_detect(OM_Scenario, "Fast")),
            mapping = aes(x = Age, y = Selex,lty = time_comp), color = "black", 
            lwd = 1, alpha = 1) +
  facet_grid(abbrev ~ OM_Scenario) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  labs(x = "Age", y = "Selectivity (Males)", 
       fill = "Functional Form", color = "Functional Form",
       lty = "Year", title = "Fast Transition") +
  theme_matt() +
  coord_cartesian(ylim = c(0, 2.5))
dev.off()



# Figure S11 (Selectivity Slow Females) --------------------------------------------


pdf(here("figs", "Manuscript_Figures_v2", "FigS11_SlowSelex_F.pdf"), width = 12, height = 15)
ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                               str_detect(OM_Scenario, "Slow")),
            mapping = aes(x = Age, y = Median_Selex, color = func), lwd = 1, alpha = 1) +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Slow")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95, fill = func), alpha = 0.25) +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", 
                                      time_comp %in% c("Terminal", "First Year"),
                                      str_detect(OM_Scenario, "Slow")),
            mapping = aes(x = Age, y = Selex,lty = time_comp), color = "black", 
            lwd = 1, alpha = 1) +
  facet_grid(abbrev ~ OM_Scenario) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  labs(x = "Age", y = "Selectivity (Females)", 
       fill = "Functional Form", color = "Functional Form",
       lty = "Year", title = "Slow Transition") +
  theme_matt() +
  coord_cartesian(ylim = c(0, 2.5))
dev.off()


# Figure S10 (Selectivity Slow Males) --------------------------------------------

pdf(here("figs", "Manuscript_Figures_v2", "FigS12_SlowSelex_M.pdf"), width = 12, height = 15)
ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                               str_detect(OM_Scenario, "Slow")),
            mapping = aes(x = Age, y = Median_Selex, color = func), lwd = 1, alpha = 1) +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Slow")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95, fill = func), alpha = 0.25) +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Male", 
                                      time_comp %in% c("Terminal", "First Year"),
                                      str_detect(OM_Scenario, "Slow")),
            mapping = aes(x = Age, y = Selex,lty = time_comp), color = "black", 
            lwd = 1, alpha = 1) +
  facet_grid(abbrev ~ OM_Scenario) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  labs(x = "Age", y = "Selectivity (Males)", 
       fill = "Functional Form", color = "Functional Form",
       lty = "Year", title = "Slow Transition") +
  theme_matt() +
  coord_cartesian(ylim = c(0, 2.5))
dev.off()

# Table (Minimax Solution) -------------------------------------------------
plot_minmax_df = minmax_df %>% filter(Dat_Qual == "High") # filter to high data quality
pdf(here("figs", "Manuscript_Figures_v2", "Fig_Minimax.pdf"), width = 20, height = 10)
print(
  ggplot(plot_minmax_df %>% filter(Dat_Qual == "High"), 
         aes(x = factor(OM_Scenario), y = factor(EM_Scenario), fill = MARE, label = round(MARE, 4))) +
    geom_tile(alpha = 0.25) +
    geom_text(color = ifelse(plot_minmax_df$MARE ==  plot_minmax_df$time_minmax_medians, "green", "black"), size = 5) +    
                facet_wrap(~time_comp, scales = "free_x") +
    scale_fill_distiller(palette = "Spectral", direction = -1) + 
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = "OM Scenario", y = "Estimation Models", fill = "SSB MARE") +
    theme_test() +
    theme(legend.position = "top",
          strip.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 17, color = "black"),
          legend.text = element_text(size = 17),
          legend.title = element_text(size = 20),
          legend.key.width = unit(1, "cm"))
)
dev.off()

# Output tables
# terminal year
terminal_minmax_df = minmax_df %>% 
  mutate(MARE = round(MARE, 4)) %>% 
  filter(Dat_Qual == "High", time_comp == "Terminal") %>% 
  select(OM_Scenario, EM_Scenario, MARE) %>% 
  pivot_wider(names_from = "OM_Scenario", values_from = "MARE")

# rearrange stuff
minmax_df <- minmax_df %>% 
  mutate(OM_Scenario = factor(OM_Scenario, 
                              levels = c("Fast_Logist_Logist",
                                         "Fast_Logist_Gamma_Old",
                                         "Fast_Logist_Gamma_Young",
                                         "Slow_Logist_Logist",
                                         "Slow_Logist_Gamma_Old",
                                         "Slow_Logist_Gamma_Young"),
                              labels = c(fast_om_plot_order, slow_om_plot_order)))

# fleet transition end
fleettrans_minmax_df = minmax_df %>% 
  mutate(MARE = round(MARE, 4)) %>% 
  filter(Dat_Qual == "High", time_comp == "Fleet Trans End") %>% 
  select(OM_Scenario, EM_Scenario, MARE) %>% 
  pivot_wider(names_from = "OM_Scenario", values_from = "MARE")

# fleet intersection
fleetint_minmax_df = minmax_df %>% 
  mutate(MARE = round(MARE, 4)) %>% 
  filter(Dat_Qual == "High", time_comp == "Fleet Intersect") %>% 
  select(OM_Scenario, EM_Scenario, MARE) %>% 
  pivot_wider(names_from = "OM_Scenario", values_from = "MARE")

# terminal
terminal_minmax_df = minmax_df %>% 
  mutate(MARE = round(MARE, 4)) %>% 
  filter(Dat_Qual == "High", time_comp == "Terminal") %>% 
  select(OM_Scenario, EM_Scenario, MARE) %>% 
  pivot_wider(names_from = "OM_Scenario", values_from = "MARE")

terminal_minmax_df = terminal_minmax_df[,c("EM_Scenario", fast_om_plot_order, slow_om_plot_order)] %>% 
  mutate(EM_Scenario = fct_relevel(EM_Scenario, plot_order)) %>% arrange(EM_Scenario)
fleettrans_minmax_df = fleettrans_minmax_df[,c("EM_Scenario", fast_om_plot_order, slow_om_plot_order)] %>% 
  mutate(EM_Scenario = fct_relevel(EM_Scenario, plot_order)) %>% arrange(EM_Scenario)
fleetint_minmax_df = fleetint_minmax_df[,c("EM_Scenario", fast_om_plot_order, slow_om_plot_order)] %>% 
  mutate(EM_Scenario = fct_relevel(EM_Scenario, plot_order)) %>% arrange(EM_Scenario)

write.table(terminal_minmax_df, here("output", "Minmax_Terminal.txt"), sep = ",", row.names = FALSE)
write.table(fleettrans_minmax_df, here("output", "Minmax_FleetTrans.txt"), sep = ",", row.names = FALSE)
write.table(fleetint_minmax_df, here("output", "Minmax_FleetInt.txt"), sep = ",", row.names = FALSE)


# Selectivity EM Plot -----------------------------------------------------

selex_path = here("output", "OM_Scenarios", "Fast_LL_High") # selectivity plot path
selex_em_df = data.table::fread(here(selex_path, "EM_Fish_Selex.csv")) 

plt_selex_em_df = selex_em_df %>% 
  filter(Sex == "Female", time_comp == "Terminal", str_detect(EM_Scenario, all_models)) %>%  # filter to terminal year and other stuff for plotting
  mutate(EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0|_0.05"),
         EM_Scenario = str_replace(EM_Scenario, "Gam", "G"), # filter out stuff and cleaning up stuff
         EM_Scenario = factor(EM_Scenario,
                                labels = c("1Fleet_RandWlkPar_Gamma",
                                           "1Fleet_SemiPar_Gamma",
                                           "1Fleet_TimeInvar_Gamma",
                                           "1Fleet_RandWlkPar_Logist",
                                           "1Fleet_SemiPar_Logist",
                                           "1Fleet_TimeInvar_Logist",
                                           "1Fleet_Block_Gamma",
                                           "1Fleet_Block_Logist",
                                           "2Fleet_Gamma",
                                           "2Fleet_Logist"))) %>% 
  group_by(Year, Age, EM_Scenario, Fleet) %>% 
  summarize(Selex = mean(Selex))


(pa = ggplot(plt_selex_em_df %>% filter(str_detect(EM_Scenario, "2Fl"), Year == 50), 
       aes(x = Age, y = Selex, color = factor(Fleet))) +
  geom_line(size = 1.3) +
  facet_wrap(~EM_Scenario) +
  theme_bw() +
  labs(x = "Age", y = "Selectivity", color = "Fleet") +
  theme(legend.position = c(0.85, 0.15),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.text= element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        title = element_text(size = 25),
        legend.key.width = unit(1.5, "cm")))

(pb = ggplot(plt_selex_em_df %>% filter(str_detect(EM_Scenario, "Invar"), Year == 50), 
       aes(x = Age, y = Selex, color = factor(Fleet))) +
  geom_line(size = 1.3) +
  facet_wrap(~EM_Scenario) +
  theme_bw() +
  labs(x = "Age", y = "", color = "Fleet") +
  theme(legend.position = c(0.85, 0.15),
          strip.text = element_text(size = 13),
          axis.title = element_text(size = 15),
          axis.text= element_text(size = 13, color = "black"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 13),
          title = element_text(size = 25),
          legend.key.width = unit(1.5, "cm")))

(pc = ggplot(plt_selex_em_df %>% filter(str_detect(EM_Scenario, "Block"), Year %in% c(1, 50)) %>% 
               mutate(Time_Block = ifelse(Year == 1, "Year 1 - 24",
                                          "Year 25 - 50")), 
       aes(x = Age, y = Selex, color = factor(Time_Block))) +
    geom_line(size = 1.3) +
    facet_wrap(~EM_Scenario) +
    theme_bw() +
    labs(x = "Age", y = "", color = "Block Period") +
    theme(legend.position = c(0.83, 0.15),
          strip.text = element_text(size = 13),
          axis.title = element_text(size = 15),
          axis.text= element_text(size = 13, color = "black"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 13),
          title = element_text(size = 25),
          legend.key.width = unit(1.5, "cm")))
  
age_year_levels <- as.factor(rev(seq(1, 50, 1))) 

(pd = ggplot(plt_selex_em_df %>% filter(str_detect(EM_Scenario, "Rand|Semi")),
             aes(Age, Year, height = Selex, group = Year, fill = Year)) + 
    ggridges::geom_density_ridges(stat = "identity", scale = 2, alpha = 0.3) +
    scale_fill_viridis_c() +
    facet_wrap(~EM_Scenario, nrow = 1) +
    scale_y_continuous(trans = "reverse") +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 13),
          axis.title = element_text(size = 15),
          axis.text= element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15),
          title = element_text(size = 25),
          legend.key.width = unit(1.5, "cm")))

a_b_c = ggpubr::ggarrange(pa, pb, pc, nrow = 1)

pdf(here("figs", "Manuscript_Figures_v2", "Fig_EM_Selex_A.pdf"), width = 20, height = 5)
a_b_c
dev.off()

pdf(here("figs", "Manuscript_Figures_v2", "Fig_EM_Selex_B.pdf"), width = 13, height = 15)
pd
dev.off()

pdf(here("figs", "Manuscript_Figures_v2", "Fig_EM_Selex.pdf"), width = 17.85, height = 15)
ggpubr::ggarrange(a_b_c, pd, nrow = 2, heights = c(0.4, 1), labels = c("A", "B"),
                  label.x = -0.003, label.y = 1.01,font.label = list(size = 28, color = "black"))
dev.off()

