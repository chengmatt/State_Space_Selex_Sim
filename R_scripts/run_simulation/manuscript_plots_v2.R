
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
AIC_df <- data.table::fread(here("output", "AIC_Convergence_Summary.csv")) # time series total error (converged runs only)
param_df <- read.csv(here("output", "Parameter_Summary.csv")) # parameters
ts_re_df <- read.csv(here("output", "TimeSeries_RE.csv")) # time series relative error (converged runs only)
ts_te_df <- read.csv(here("output", "TimeSeries_TE.csv")) # time series total error (converged runs only)
ts_are_df <- read.csv(here("output", "TimeSeries_ARE.csv")) # time series abs relative error (converged runs only)

# Selectivity and Comps
om_slx_df <- data.table::fread(here("output", "OM_Fish_Selex.csv")) # OM Selectivity Values
pop_sel_om <- data.table::fread(here("output", "Pop_Selex_OM.csv")) # OM Pop'n Selectivity Values
pop_sel_em <- data.table::fread(here("output", "Pop_Selex_EM.csv")) # EM Pop'n Selectivity Values

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
fast_om_plot_order <- c("Fast_Logist_Logist", "Fast_Logist_Gamma_Old", "Fast_Logist_Gamma_Young")
slow_om_plot_order <- c("Slow_Logist_Logist", "Slow_Logist_Gamma_Old", "Slow_Logist_Gamma_Young")

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
  ), OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
     EM_Scenario = str_replace(EM_Scenario, "Gam", "G")) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(converged = sum(conv == "Converged")/200) %>% # divide by 200 = number of sims run
  ungroup() %>% 
  mutate(OM_Scenario = factor(OM_Scenario, 
                                labels = c("Fast_Logist_Gamma_Old",
                                           "Fast_Logist_Gamma_Young",
                                           "Fast_Logist_Logist",
                                           "Slow_Logist_Gamma_Old",
                                           "Slow_Logist_Gamma_Young",
                                           "Slow_Logist_Logist")),
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
  order[o] <- which(grepl(plot_order[o], x = unique(conv_stat$EM_Scenario)))
} # end o loop

# relevel factors
conv_stat$EM_Scenario <- factor(conv_stat$EM_Scenario, levels = c(unique(conv_stat$EM_Scenario)[order]))

### AIC Stuff ---------------------------------------------------------------

AIC_df <- AIC_df %>% 
  filter(conv == "Converged") %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = case_when(
    OM_Scenario == "Fast_LL" ~ "Fast_Logist_Logist",
    OM_Scenario == "Fast_LG_Y" ~ "Fast_Logist_Gamma_Young",
    OM_Scenario == "Fast_LG_O" ~ "Fast_Logist_Gamma_Old",
    OM_Scenario == "Slow_LL" ~ "Slow_Logist_Logist",
    OM_Scenario == "Slow_LG_Y" ~ "Slow_Logist_Gamma_Young",
    OM_Scenario == "Slow_LG_O" ~ "Slow_Logist_Gamma_Old"
  ),
  OM_Scenario = factor(OM_Scenario,  levels = c(fast_om_plot_order, slow_om_plot_order)))

# Two fleet AIC df
twofleet_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "2Fl_")) %>% 
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
         ))

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
         ), EM_Scenario = factor(EM_Scenario, levels = rev(plot_order))) 


### Time Block AIC Stuff ----------------------------------------------------

# Filter to 1 fleet time block models
onefleet_blk_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "Blk"),
         str_detect(OM_Scenario, "Fast")) %>% 
  mutate(
    selex_form = case_when(
      str_detect(EM_Scenario, "LGam") ~ "Logistic_Gamma",
      str_detect(EM_Scenario, "LL") ~ "Logistic_Logistic",
    )) %>% group_by(sim, OM_Scenario, time_comp, selex_form, Dat_Qual) %>% 
  mutate(min = min(AIC), min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, selex_form, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "1 Fleet",
         EM_Scenario = factor(EM_Scenario, levels = blk_sens_order, labels = blk_sens_labels))

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
                              labels = c("Fast_Logist_Logist",
                                         "Fast_Logist_Gamma_Old",
                                         "Fast_Logist_Gamma_Young")))

### ABC Parameters ----------------------------------------------------------

# Filter to relevant components for parameters
om_scenario_params <- param_df %>% filter(str_detect(EM_Scenario, all_models), 
                                          type %in% c("F_0.4", "ABC")) %>%
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = factor(OM_Scenario, labels = c("Fast_Logist_Gamma_Old",
                                         "Fast_Logist_Gamma_Young",
                                         "Fast_Logist_Logist",
                                         "Slow_Logist_Gamma_Old",
                                         "Slow_Logist_Gamma_Young",
                                         "Slow_Logist_Logist")),
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
            lwr_95 = quantile(RE, 0.025),
            upr_95 =  quantile(RE, 0.975))

# Clarify names
pt_rg_re <- pt_rg_re %>% 
  mutate(type = case_when(
    type == "ABC" ~ "ABC",
    type == "F_0.4" ~ "F40%"),
    EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
    OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)),
    EM_Scenario = factor(EM_Scenario, levels = rev(plot_order)))


# Selectivity Stuff -------------------------------------------------------

plot_df = pop_sel_em %>% 
  filter(str_detect(EM_Scenario, all_models)) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), 
  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
  OM_Scenario = factor(OM_Scenario, labels = c("Fast_Logist_Gamma_Old",
                                               "Fast_Logist_Gamma_Young",
                                               "Fast_Logist_Logist",
                                               "Slow_Logist_Gamma_Old",
                                               "Slow_Logist_Gamma_Young",
                                               "Slow_Logist_Logist")),
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
  EM_Scenario = factor(EM_Scenario, levels = plot_order))

plot_sel_om_df = pop_sel_om %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), 
  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = factor(OM_Scenario, labels = c("Fast_Logist_Gamma_Old",
                                               "Fast_Logist_Gamma_Young",
                                               "Fast_Logist_Logist",
                                               "Slow_Logist_Gamma_Old",
                                               "Slow_Logist_Gamma_Young",
                                               "Slow_Logist_Logist")),
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
                       labels = c("Fast_Logist_Gamma_Old",
                                  "Fast_Logist_Gamma_Young",
                                  "Fast_Logist_Logist",
                                  "Slow_Logist_Gamma_Old",
                                  "Slow_Logist_Gamma_Young",
                                  "Slow_Logist_Logist")),
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
                       labels = c("Fast_Logist_Logist", "Fast_Logist_Gamma_Old",
                                  "Fast_Logist_Gamma_Young")),
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
    )) 

# Left Join "Best AIC" Models
tb_ts_re_om = tb_ts_re_om %>% left_join(best_fleetblk_models,
                                        by = c("OM_Scenario", "EM_Scenario", "time_comp", "Dat_Qual")) %>% 
  mutate(Model = ifelse(is.na(Model), "Not Best", Model ))

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
          func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logisitic")) # functional form

(a = ggplot(ssb_fast_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Fast_Logist_Logist"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
   facet_wrap(~EM_Scenario, nrow = 2) +
   scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
   scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
   coord_cartesian(ylim = c(-0.4, 0.4)) +
   geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
   theme_matt() +
   labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Fast_Logist_Logist"))

(b = ggplot(ssb_fast_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Fast_Logist_Gamma_Old"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.4, 0.4)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "Relative Error in Spawning Biomass", fill = "Time", color = "Time", title = "Fast_Logist_Gamma_Old"))

(c = ggplot(ssb_fast_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Fast_Logist_Gamma_Young"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.4, 0.4)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Fast_Logist_Gamma_Young"))

pdf(here("figs", "Manuscript_Figures_v2", "Fig2_SSB.pdf"), width = 20, height = 23)
ggpubr::ggarrange(a,b,c, common.legend = TRUE, ncol = 1, labels = c("A", "B", "C"), widths = c(0.975,1,1),
                  label.x = 0.02, label.y = 1.01,font.label = list(size = 34, color = "black"))
dev.off()


# Figure 3 (RE ABC Fast) --------------------------------------------------

pdf(here("figs", "Manuscript_Figures_v2", "Fig3_ABC.pdf"), width = 15, height = 10)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, "Fast_"),
                             type == "ABC"),
         aes(x = factor(EM_Scenario), y = median, color = time_comp, fill = time_comp,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), 
                    size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    geom_vline(xintercept = c(seq(1.5, 7.5, 1)), lwd = 0.25, alpha = 0.75) +
    facet_grid(~OM_Scenario, scales = "free_x") +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_x_discrete(guide = guide_axis(angle = 0)) +
    labs(x = "Estimation Models", y = "Relative Error in ABC",
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          axis.title = element_text(size = 17),
          axis.text= element_text(size = 15),
          strip.text = element_text(size = 17),
          axis.text.y = element_text(angle = 90),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17)) +
    coord_flip(ylim = c(-0.4, 0.4))
)
dev.off()


# Figure 4 (RE SSB Slow) --------------------------------------------------

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
          func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logisitic")) # functional form

(a = ggplot(ssb_slow_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Slow_Logist_Logist"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.4, 0.4)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Slow_Logist_Logist"))

(b = ggplot(ssb_slow_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Slow_Logist_Gamma_Old"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.4, 0.4)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "Relative Error in Spawning Biomass", fill = "Time", color = "Time", title = "Slow_Logist_Gamma_Old"))

(c = ggplot(ssb_slow_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Slow_Logist_Gamma_Young"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.4, 0.4)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Slow_Logist_Gamma_Young"))

pdf(here("figs", "Manuscript_Figures_v2", "Fig4_SSB.pdf"), width = 20, height = 23)
ggpubr::ggarrange(a,b,c, common.legend = TRUE, ncol = 1, labels = c("A", "B", "C"), widths = c(0.975,1,1),
                  label.x = 0.02, label.y = 1.01,font.label = list(size = 34, color = "black"))
dev.off()

# Figure 5 (RE ABC Slow) --------------------------------------------------
pdf(here("figs", "Manuscript_Figures_v2", "Fig5_ABC.pdf"), width = 15, height = 10)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, "Slow_"),
                             type == "ABC"),
         aes(x = factor(EM_Scenario), y = median, color = time_comp, fill = time_comp,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), 
                    size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    geom_vline(xintercept = c(seq(1.5, 7.5, 1)), lwd = 0.25, alpha = 0.75) +
    facet_grid(~OM_Scenario, scales = "free_x") +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_x_discrete(guide = guide_axis(angle = 0)) +
    labs(x = "Estimation Models", y = "Relative Error in ABC",
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          axis.title = element_text(size = 17),
          axis.text= element_text(size = 15),
          strip.text = element_text(size = 17),
          axis.text.y = element_text(angle = 90),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17)) +
    coord_flip(ylim = c(-0.4, 0.4))
)
dev.off()

# Figure S1 (Convergence) -------------------------------------------------
pdf(here("figs", "Manuscript_Figures_v2", "FigS1_Convergence.pdf"), width = 25, height = 13)
ggplot(conv_stat %>% 
         mutate(OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)),
                EM_Scenario = forcats::fct_rev(EM_Scenario)) %>% 
                              filter(Dat_Qual == "High"),
       mapping = aes(x = EM_Scenario, y = converged * 100, group = time_comp, fill = time_comp))  +
  geom_col(position = position_dodge2(width = 1), 
           alpha = 0.7, color = "black") +
  geom_hline(aes(yintercept = 90), col = "black", lty = 2, size = 0.7, alpha = 1) +
  geom_hline(aes(yintercept = 50), col = "black", lty = 2, size = 0.7, alpha = 1) +
  facet_grid(time_comp~OM_Scenario,  scales = "free_x") +
  theme_matt() + 
  labs(fill = "Assessment Period", y = "Convergence Rate (%)", x = "Estimation Models") +
  theme(legend.position = "none") +
  coord_flip() 
dev.off()


# Figure S2 and S3 (AIC) --------------------------------------------------

pdf(here("figs", "Manuscript_Figures_v2", "FigS2_AIC_2Fl.pdf"), width = 15, height = 5)
ggplot(twofleet_aic %>% filter(Dat_Qual == "High"), 
       aes(x = OM_Scenario, y = EM_Scenario, label = round(n_minAIC, 2))) +
  geom_tile(alpha = 0.35, color = "black") +
  geom_text(size = 5) +
  facet_wrap(~time_comp) +
  coord_flip() +
  labs(x = "Operating Models", y = "Estimation Models", fill = "Proportion of models with lowest AIC") +
  theme_test() +
  theme(legend.position = "top",
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text= element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.key.width = unit(1, "cm"))
dev.off()

pdf(here("figs", "Manuscript_Figures_v2", "FigS3_AIC_1Fl.pdf"), width = 20, height = 8)
ggplot(onefleet_aic %>% filter(Dat_Qual == "High"), 
       aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC, label = round(n_minAIC, 2))) +
  geom_tile(alpha = 0.65) +
  geom_text(size = 6) +
  facet_grid(~time_comp, scales = "free_y") +
  scale_y_discrete(guide = guide_axis(angle = 90)) + 
  coord_flip() +
  scale_fill_distiller(palette = "Spectral", direction = 1) + 
  labs(x = "Operating Models", y = "Estimation Models", fill = "Proportion of models with lowest AIC") +
  theme_test() +
  theme(legend.position = "top",
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text= element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.key.width = unit(1, "cm"))
dev.off()

# Figure S4 (RE Total F Fast) --------------------------------------------------

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
          func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logisitic")) # functional form

(a = ggplot(f_fast_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Fast_Logist_Logist"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Fast_Logist_Logist"))

(b = ggplot(f_fast_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Fast_Logist_Gamma_Old"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Fast_Logist_Gamma_Old"))

(c = ggplot(f_fast_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Fast_Logist_Gamma_Young"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Fast_Logist_Gamma_Young"))

pdf(here("figs", "Manuscript_Figures_v2", "FigS4_F.pdf"), width = 20, height = 23)
ggpubr::ggarrange(a,b,c, common.legend = TRUE, ncol = 1, labels = c("A", "B", "C"), widths = c(0.975,1,1),
                  label.x = 0.02, label.y = 1.01,font.label = list(size = 34, color = "black"))
dev.off()


# Figure S5 (RE Total F Slow) --------------------------------------------------

# Filter to fast scenario, and differentiate models
f_slow_re_om = ts_re_om %>% 
  filter(par_name == "Total Fishing Mortality", str_detect(OM_Scenario, "Slow")) %>% 
  mutate( OM_Scenario = factor(OM_Scenario, levels = slow_om_plot_order),
          model_type = case_when( # differenitate model types
            str_detect(EM_Scenario, 'Random') ~ "1Fleet Random Walk",
            str_detect(EM_Scenario, 'Semi') ~ "1Fleet Semi-Parametric",
            str_detect(EM_Scenario, 'TimeInvar') ~ "1Fleet TimeInvar",
            str_detect(EM_Scenario, 'TimeBlk') ~ "1Fleet TimeBlock",
            str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
          ),
          func = ifelse(str_detect(EM_Scenario, "Gamma"), "Gamma", "Logisitic")) # functional form

(a = ggplot(f_slow_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Slow_Logist_Logist"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Slow_Logist_Logist"))

(b = ggplot(f_slow_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Slow_Logist_Gamma_Old"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Slow_Logist_Gamma_Old"))

(c = ggplot(f_slow_re_om %>% filter(Dat_Qual == "High", OM_Scenario == "Slow_Logist_Gamma_Young"),
            aes(x = year, y = median, color=time_comp)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill=time_comp),
                alpha = 0.35, color = NA) +
    geom_line(show.legend = F, size = 2, alpha = 1) +
    facet_wrap(~EM_Scenario, nrow = 2) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    theme_matt() +
    labs(x = "Year", y = "", fill = "Time", color = "Time", title = "Slow_Logist_Gamma_Young"))

pdf(here("figs", "Manuscript_Figures_v2", "FigS5_F.pdf"), width = 20, height = 23)
ggpubr::ggarrange(a,b,c, common.legend = TRUE, ncol = 1, labels = c("A", "B", "C"), widths = c(0.975,1,1),
                  label.x = 0.02, label.y = 1.01,font.label = list(size = 34, color = "black"))
dev.off()


# Figure S6 (Time Block AIC) ----------------------------------------------

pdf(here("figs", "Manuscript_Figures_v2", "FigS6_BlkAIC.pdf"), width = 15, height = 10)
ggplot(onefleet_blk_aic %>% filter(time_comp %in% c("Terminal", "Fleet Trans End"), Dat_Qual == "High"),
       aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC, label = round(n_minAIC, 2))) +
  geom_tile(alpha = 0.65) +
  geom_text(size = 5) +
  ggh4x::facet_nested(selex_form~time_comp, scales = "free") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_distiller(palette = "Spectral", direction = 1) + 
  guides(fill = guide_colourbar(barwidth = 10)) +
  labs(x = "OM Scenario", y = "Estimation Models", fill = "Proportion of models with lowest AIC") +
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
  geom_line(linewidth = 1.5, alpha = 0.8, aes(color = EM_Scenario)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  coord_cartesian(ylim = c(-0.4,0.4)) +
  scale_linetype_manual(values = c(2,1)) +
  guides(linetype = "none") +
  facet_grid(OM_Scenario~time_comp) + 
  labs(x = "Year", y = "Relative Error in Spawning Biomass", 
       fill = "Estimation Models", color = "Estimation Models") +
  theme_matt() +
  theme(legend.position = "right", 
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
    geom_line(linewidth = 1.5, alpha = 0.8, aes(color = EM_Scenario)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    coord_cartesian(ylim = c(-0.4,0.4)) +
    scale_linetype_manual(values = c(2,1)) +
    guides(linetype = "none") +
    facet_grid(OM_Scenario~time_comp) + 
    labs(x = "Year", y = "Relative Error in Spawning Biomass", 
         fill = "Estimation Models", color = "Estimation Models") +
    theme_matt() +
    theme(legend.position = "right", 
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

(a = ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                               str_detect(OM_Scenario, "Fast_Logist_Logist")),
            mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast_Logist_Logist")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
              fill = "red") +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal", 
                                  str_detect(OM_Scenario, "Fast_Logist_Logist")),
            mapping = aes(x = Age, y = Selex), color = "black",
            position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
  facet_wrap( ~ EM_Scenario, nrow = 2) +
  labs(x = "Age", y = "", fill = "Assessment Period",
       title = "Fast_Logist_Logist") +
  theme_matt() +
   coord_cartesian(ylim = c(0, 2.5)))

(b = ggplot() +
    geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast_Logist_Gamma_Old")),
              mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
    geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                   str_detect(OM_Scenario, "Fast_Logist_Gamma_Old")),
                mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
                fill = "red") +
    geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal", 
                                        str_detect(OM_Scenario, "Fast_Logist_Gamma_Old")),
              mapping = aes(x = Age, y = Selex), color = "black",
              position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
    facet_wrap( ~ EM_Scenario, nrow = 2) +
    labs(x = "Age", y = "Population Selectivity (Terminal Period)", fill = "Assessment Period",
         title = "Fast_Logist_Gamma_Old") +
    theme_matt() +
    coord_cartesian(ylim = c(0, 2.5)) )

(c = ggplot() +
    geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast_Logist_Gamma_Young")),
              mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
    geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                   str_detect(OM_Scenario, "Fast_Logist_Gamma_Young")),
                mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
                fill = "red") +
    geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal", 
                                        str_detect(OM_Scenario, "Fast_Logist_Gamma_Young")),
              mapping = aes(x = Age, y = Selex), color = "black",
              position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
    facet_wrap( ~ EM_Scenario, nrow = 2) +
    labs(x = "Age", y = "", fill = "Assessment Period",
         title = "Fast_Logist_Gamma_Young") +
    theme_matt() +
    coord_cartesian(ylim = c(0, 2.5)))


pdf(here("figs", "Manuscript_Figures_v2", "FigS9_FastSelex_F.pdf"), width = 20, height = 25)
ggpubr::ggarrange(a,b,c, common.legend = TRUE, ncol = 1, labels = c("A", "B", "C"), widths = c(0.975,1,1),
                  label.x = 0.01, label.y = 1.01,font.label = list(size = 34, color = "black"))
dev.off()


# Figure S10 (Selectivity Fast Males) --------------------------------------------

(a = ggplot() +
   geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                str_detect(OM_Scenario, "Fast_Logist_Logist")),
             mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "blue") +
   geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                  str_detect(OM_Scenario, "Fast_Logist_Logist")),
               mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
               fill = "blue") +
   geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal", 
                                       str_detect(OM_Scenario, "Fast_Logist_Logist")),
             mapping = aes(x = Age, y = Selex), color = "black",
             position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
   facet_wrap( ~ EM_Scenario, nrow = 2) +
   labs(x = "Age", y = "", fill = "Assessment Period",
        title = "Fast_Logist_Logist") +
   theme_matt() +
   coord_cartesian(ylim = c(0, 2.5)))

(b = ggplot() +
    geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast_Logist_Gamma_Old")),
              mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "blue") +
    geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                   str_detect(OM_Scenario, "Fast_Logist_Gamma_Old")),
                mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
                fill = "blue") +
    geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal", 
                                        str_detect(OM_Scenario, "Fast_Logist_Gamma_Old")),
              mapping = aes(x = Age, y = Selex), color = "black",
              position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
    facet_wrap( ~ EM_Scenario, nrow = 2) +
    labs(x = "Age", y = "Population Selectivity (Terminal Period)", fill = "Assessment Period",
         title = "Fast_Logist_Gamma_Old") +
    theme_matt() +
    coord_cartesian(ylim = c(0, 2.5)) )

(c = ggplot() +
    geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast_Logist_Gamma_Young")),
              mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "blue") +
    geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                   str_detect(OM_Scenario, "Fast_Logist_Gamma_Young")),
                mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
                fill = "blue") +
    geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal", 
                                        str_detect(OM_Scenario, "Fast_Logist_Gamma_Young")),
              mapping = aes(x = Age, y = Selex), color = "black",
              position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
    facet_wrap( ~ EM_Scenario, nrow = 2) +
    labs(x = "Age", y = "", fill = "Assessment Period",
         title = "Fast_Logist_Gamma_Young") +
    theme_matt() +
    coord_cartesian(ylim = c(0, 2.5)))


pdf(here("figs", "Manuscript_Figures_v2", "FigS10_FastSelex_M.pdf"), width = 20, height = 25)
ggpubr::ggarrange(a,b,c, common.legend = TRUE, ncol = 1, labels = c("A", "B", "C"), widths = c(0.975,1,1),
                  label.x = 0.01, label.y = 1.01,font.label = list(size = 34, color = "black"))
dev.off()


# Figure S11 (Selectivity Slow Females) --------------------------------------------

(a = ggplot() +
   geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                str_detect(OM_Scenario, "Slow_Logist_Logist")),
             mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
   geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                  str_detect(OM_Scenario, "Slow_Logist_Logist")),
               mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
               fill = "red") +
   geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal", 
                                       str_detect(OM_Scenario, "Slow_Logist_Logist")),
             mapping = aes(x = Age, y = Selex), color = "black",
             position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
   facet_wrap( ~ EM_Scenario, nrow = 2) +
   labs(x = "Age", y = "", fill = "Assessment Period",
        title = "Slow_Logist_Logist") +
   theme_matt() +
   coord_cartesian(ylim = c(0, 2.5)))

(b = ggplot() +
    geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Slow_Logist_Gamma_Old")),
              mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
    geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                   str_detect(OM_Scenario, "Slow_Logist_Gamma_Old")),
                mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
                fill = "red") +
    geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal", 
                                        str_detect(OM_Scenario, "Slow_Logist_Gamma_Old")),
              mapping = aes(x = Age, y = Selex), color = "black",
              position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
    facet_wrap( ~ EM_Scenario, nrow = 2) +
    labs(x = "Age", y = "Population Selectivity (Terminal Period)", fill = "Assessment Period",
         title = "Slow_Logist_Gamma_Old") +
    theme_matt() +
    coord_cartesian(ylim = c(0, 2.5)) )

(c = ggplot() +
    geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Slow_Logist_Gamma_Young")),
              mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
    geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                   str_detect(OM_Scenario, "Slow_Logist_Gamma_Young")),
                mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
                fill = "red") +
    geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal", 
                                        str_detect(OM_Scenario, "Slow_Logist_Gamma_Young")),
              mapping = aes(x = Age, y = Selex), color = "black",
              position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
    facet_wrap( ~ EM_Scenario, nrow = 2) +
    labs(x = "Age", y = "", fill = "Assessment Period",
         title = "Slow_Logist_Gamma_Young") +
    theme_matt() +
    coord_cartesian(ylim = c(0, 2.5)))


pdf(here("figs", "Manuscript_Figures_v2", "FigS11_SlowSelex_F.pdf"), width = 20, height = 25)
ggpubr::ggarrange(a,b,c, common.legend = TRUE, ncol = 1, labels = c("A", "B", "C"), widths = c(0.975,1,1),
                  label.x = 0.01, label.y = 1.01,font.label = list(size = 34, color = "black"))
dev.off()


# Figure S10 (Selectivity Slow Males) --------------------------------------------

(a = ggplot() +
   geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                str_detect(OM_Scenario, "Slow_Logist_Logist")),
             mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "blue") +
   geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                  str_detect(OM_Scenario, "Slow_Logist_Logist")),
               mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
               fill = "blue") +
   geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal", 
                                       str_detect(OM_Scenario, "Slow_Logist_Logist")),
             mapping = aes(x = Age, y = Selex), color = "black",
             position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
   facet_wrap( ~ EM_Scenario, nrow = 2) +
   labs(x = "Age", y = "", fill = "Assessment Period",
        title = "Slow_Logist_Logist") +
   theme_matt() +
   coord_cartesian(ylim = c(0, 2.5)))

(b = ggplot() +
    geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Slow_Logist_Gamma_Old")),
              mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "blue") +
    geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                   str_detect(OM_Scenario, "Slow_Logist_Gamma_Old")),
                mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
                fill = "blue") +
    geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal", 
                                        str_detect(OM_Scenario, "Slow_Logist_Gamma_Old")),
              mapping = aes(x = Age, y = Selex), color = "black",
              position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
    facet_wrap( ~ EM_Scenario, nrow = 2) +
    labs(x = "Age", y = "Population Selectivity (Terminal Period)", fill = "Assessment Period",
         title = "Slow_Logist_Gamma_Old") +
    theme_matt() +
    coord_cartesian(ylim = c(0, 2.5)) )

(c = ggplot() +
    geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Slow_Logist_Gamma_Young")),
              mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "blue") +
    geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal",
                                   str_detect(OM_Scenario, "Slow_Logist_Gamma_Young")),
                mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
                fill = "blue") +
    geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Male", time_comp == "Terminal", 
                                        str_detect(OM_Scenario, "Slow_Logist_Gamma_Young")),
              mapping = aes(x = Age, y = Selex), color = "black",
              position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
    facet_wrap( ~ EM_Scenario, nrow = 2) +
    labs(x = "Age", y = "", fill = "Assessment Period",
         title = "Slow_Logist_Gamma_Young") +
    theme_matt() +
    coord_cartesian(ylim = c(0, 2.5)))


pdf(here("figs", "Manuscript_Figures_v2", "FigS12_SlowSelex_M.pdf"), width = 20, height = 25)
ggpubr::ggarrange(a,b,c, common.legend = TRUE, ncol = 1, labels = c("A", "B", "C"), widths = c(0.975,1,1),
                  label.x = 0.01, label.y = 1.01,font.label = list(size = 34, color = "black"))
dev.off()

