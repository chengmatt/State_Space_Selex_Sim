# Purpose: To replot sensitivity figures 
# Creator: Matthew LH. Cheng
# 11/29/23

# Set Up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggh4x)
library(data.table)
library(ggrepel)

# Paths
om_scenario_path <- here("output", "OM_Scenarios_Sensitivity") # path to OM folder
dir_out <- here("output", "Summary_Plots") # path to output folder
dir.create(dir_out)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
files <- list.files(fxn_path) # Load in all functions from the functions folder
for(i in 1:length(files)) source(here(fxn_path, files[i]))

# read in csvs
# Parameter and Time Series Summaries
AIC_df <- data.table::fread(here("output", "OM_Scenarios_Sensitivity", "AIC_Convergence_Summary.csv")) %>%   # time series total error (converged runs only)
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
  # time series total error (converged runs only)
param_df <- read.csv(here("output", "OM_Scenarios_Sensitivity","Parameter_Summary.csv")) %>% # parameters
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
ts_all_df <- data.table::fread(here("output", "OM_Scenarios_Sensitivity", "TimeSeries_Summary.csv")) %>%  # all time series data
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
ts_re_df <- read.csv(here("output", "OM_Scenarios_Sensitivity", "TimeSeries_RE.csv")) %>% # time series relative error (converged runs only)
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
ts_te_df <- read.csv(here("output", "OM_Scenarios_Sensitivity", "TimeSeries_TE.csv")) %>%  # time series total error (converged runs only)
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
ts_are_df <- read.csv(here("output", "OM_Scenarios_Sensitivity", "TimeSeries_ARE.csv")) %>%  # time series abs relative error (converged runs only)
  filter(!str_detect(OM_Scenario, "Rev|Ext"))

# Selectivity and Comps
om_slx_df <- data.table::fread(here("output", "OM_Scenarios_Sensitivity", "OM_Fish_Selex.csv")) %>% # OM Selectivity Values
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
pop_sel_om <- data.table::fread(here("output", "OM_Scenarios_Sensitivity", "Pop_Selex_OM.csv")) %>%  # OM Pop'n Selectivity Values
  filter(!str_detect(OM_Scenario, "Rev|Ext"))
pop_sel_em <- data.table::fread(here("output", "OM_Scenarios_Sensitivity", "Pop_Selex_EM.csv")) %>%  # EM Pop'n Selectivity Values
  filter(!str_detect(OM_Scenario, "Rev|Ext"))


# Unique oms and other components
unique_oms <- unique(param_df$OM_Scenario) # unique oms
ts_pars <- unique(ts_re_df$par_name) # unique parameter names

# Timeblock sensitivity models
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


ssb_fast <- ggplot() + 
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
       size = "Assessment Period", title = "Fast Transition (Survey Sensitivity Test)")

pdf(here("figs", "Sensitivity_Manuscript_Survey", "Fig2_SSB.pdf"), width = 19, height = 21)
ssb_fast
dev.off()

pdf(here("figs", "Manuscript_Figures_v2", "FigS7_SSB.pdf"), width = 19, height = 21)
ssb_fast
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

ssb_slow <- ggplot() + 
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
       size = "Assessment Period", title = "Slow Transition (Survey Sensitivity Test)")

pdf(here("figs", "Sensitivity_Manuscript_Survey", "Fig3_SSB.pdf"), width = 19, height = 21)
ssb_slow
dev.off()

pdf(here("figs", "Manuscript_Figures_v2", "FigS8_SSB.pdf"), width = 19, height = 21)
ssb_slow
dev.off()


# Figure 4 (RE ABC Fast) --------------------------------------------------

pdf(here("figs", "Sensitivity_Manuscript_Survey", "Fig4_ABC.pdf"), width = 15, height = 15)
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
pdf(here("figs", "Sensitivity_Manuscript_Survey", "Fig5_ABC.pdf"), width = 15, height = 15)
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


# Figure S1 (Convergence) -------------------------------------------------
pdf(here("figs", "Sensitivity_Manuscript_Survey", "FigS1_Convergence.pdf"), width = 25, height = 15)
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

pdf(here("figs", "Sensitivity_Manuscript_Survey", "FigS2_AIC_2Fl.pdf"), width = 17, height = 10)
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

pdf(here("figs", "Sensitivity_Manuscript_Survey", "FigS3_AIC_1Fl.pdf"), width = 20, height = 13)
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


pdf(here("figs", "Sensitivity_Manuscript_Survey", "FigS4_F.pdf"), width = 19, height = 21)
ggplot() + 
  geom_hline(yintercept = 0, col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  geom_ribbon(f_fast_re_om %>% filter(time_comp == "Terminal"),
              mapping = aes(x = year, y = median, ymin = lwr_95, ymax = upr_95, fill = func), 
              color = NA, alpha = 0.3) +
  geom_line(f_fast_re_om, mapping = aes(x = year, y = median, size = time_comp,
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
  labs(x = "Year", y = "Relative Error in Total Fishing Mortality", fill = "Functional Form",
       color = "Functional Form", lty = "Assessment Period", alpha = "Assessment Period",
       size = "Assessment Period", title = "Fast Transition")
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


pdf(here("figs", "Sensitivity_Manuscript_Survey", "FigS5_F.pdf"), width = 19, height = 21)
ggplot() + 
  geom_hline(yintercept = 0, col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  geom_ribbon(f_slow_re_om %>% filter(time_comp == "Terminal"),
              mapping = aes(x = year, y = median, ymin = lwr_95, ymax = upr_95, fill = func), 
              color = NA, alpha = 0.3) +
  geom_line(f_slow_re_om, mapping = aes(x = year, y = median, size = time_comp,
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
  labs(x = "Year", y = "Relative Error in Total Fishing Mortality", fill = "Functional Form",
       color = "Functional Form", lty = "Assessment Period", alpha = "Assessment Period",
       size = "Assessment Period", title = "Slow Transition")
dev.off()

# Figure S9 (Selectivity Fast Females) --------------------------------------------

pdf(here("figs", "Sensitivity_Manuscript_Survey", "FigS9_FastSelex_F.pdf"), width = 12, height = 15)
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

pdf(here("figs", "Sensitivity_Manuscript_Survey", "FigS10_FastSelex_M.pdf"), width = 12, height = 15)
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


pdf(here("figs", "Sensitivity_Manuscript_Survey", "FigS11_SlowSelex_F.pdf"), width = 12, height = 15)
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

pdf(here("figs", "Sensitivity_Manuscript_Survey", "FigS12_SlowSelex_M.pdf"), width = 12, height = 15)
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

