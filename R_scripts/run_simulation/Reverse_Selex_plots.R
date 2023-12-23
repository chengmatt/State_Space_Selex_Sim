# Purpose: To plot selectivity scenarios, but reversed from the orignial OM
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

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
files <- list.files(fxn_path) # Load in all functions from the functions folder
for(i in 1:length(files)) source(here(fxn_path, files[i]))

# read in csvs
# Parameter and Time Series Summaries
AIC_df <- data.table::fread(here("output", "AIC_Convergence_Summary.csv")) # time series total error (converged runs only)
param_df <- read.csv(here("output", "Parameter_Summary.csv")) # parameters
ts_all_df <- data.table::fread(here("output", "TimeSeries_Summary.csv")) # all time series data
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
all_models <- "2Fl_LL|2Fl_GamL|1Fl_L_TI|1Fl_Gam_TI|1Fl_LL_Blk\\b|1Fl_GamL_Blk\\b|1Fl_L_RW_1.25|1Fl_Gam_RW_2.0|1Fl_L_SP_0.75|1Fl_Gam_SP_0.75"
gen_om <- c("Fast", "Slow") # general OMs

# Get plot order for EM models
plot_order <- c("2Fleet_Logist_Logist", 
                "2Fleet_Gamma_Logist",  
                "1Fleet_Logist_TimeInvar", 
                "1Fleet_Gamma_TimeInvar", 
                "1Fleet_Logist_Logist_TimeBlk", 
                "1Fleet_Gamma_Logist_TimeBlk", 
                "1Fleet_Logist_RandomWalk", 
                "1Fleet_Gamma_RandomWalk", 
                "1Fleet_Logist_SemiPar", 
                "1Fleet_Gamma_SemiPar")

# Fast and slow plot order
fast_om_plot_order <- c("Fast_Logist_Logist_Rev", "Fast_Gamma_Logist_Old_Rev", "Fast_Gamma_Logist_Young_Rev")
slow_om_plot_order <- c("Slow_Logist_Logist_Rev", "Slow_Gamma_Logist_Old_Rev", "Slow_Gamma_Logist_Young_Rev")

# Residual Munging --------------------------------------------------------

### ABC Parameters ----------------------------------------------------------

# Filter to relevant components for parameters
om_scenario_params <- param_df %>% filter(str_detect(EM_Scenario, all_models), 
                                          type %in% c("F_0.4", "ABC")) %>%
  filter(str_detect(OM_Scenario, "Rev")) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = factor(OM_Scenario, 
                              labels = c("Fast_Gamma_Logist_Old_Rev",
                                         "Fast_Gamma_Logist_Young_Rev",
                                         "Fast_Logist_Logist_Rev",
                                         "Slow_Gamma_Logist_Old_Rev",
                                         "Slow_Gamma_Logist_Young_Rev",
                                         "Slow_Logist_Logist_Rev")))

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
    OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order))) %>% 
  mutate(func = ifelse(str_detect(EM_Scenario, "G"), "Gamma", "Logistic"),
         abbrev = case_when(
           str_detect(EM_Scenario, "2Fl") ~ "2Fleet",
           str_detect(EM_Scenario, "TI") ~ "1Fleet TimeInvar",
           str_detect(EM_Scenario, "Blk") ~ "1Fleet Block",
           str_detect(EM_Scenario, "RW") ~ "1Fleet RandWlkPar",
           str_detect(EM_Scenario, "SP") ~ "1Fleet SemiPar"
         ),
         abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet TimeInvar",
                                            "1Fleet Block", "1Fleet RandWlkPar",
                                            "1Fleet SemiPar")))

### Time Series Stuff -------------------------------------------------------

# Relative error of time series
ts_re_om <- ts_re_df %>% filter(str_detect(EM_Scenario, all_models)) %>% 
  filter(str_detect(OM_Scenario, "Rev")) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
  OM_Scenario = factor(OM_Scenario, 
                       labels = c("Fast_Gamma_Logist_Old_Rev",
                                  "Fast_Gamma_Logist_Young_Rev",
                                  "Fast_Logist_Logist_Rev",
                                  "Slow_Gamma_Logist_Old_Rev",
                                  "Slow_Gamma_Logist_Young_Rev",
                                  "Slow_Logist_Logist_Rev")))

# Now relevel factor for organizing plot
ts_re_om <- ts_re_om %>% 
  mutate(time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")))

# Plots -------------------------------------------------------------------
# Figure 2 (RE SSB Fast) --------------------------------------------------

# Filter to fast scenario, and differentiate models
ssb_fast_re_om = ts_re_om %>% 
  filter(par_name == "Spawning Stock Biomass", str_detect(OM_Scenario, "Fast")) %>% 
  mutate( OM_Scenario = factor(OM_Scenario, levels = fast_om_plot_order),
          model_type = case_when( # differenitate model types
            str_detect(EM_Scenario, 'RW') ~ "1Fleet Random Walk",
            str_detect(EM_Scenario, 'SP') ~ "1Fleet Semi-Parametric",
            str_detect(EM_Scenario, 'TI') ~ "1Fleet TimeInvar",
            str_detect(EM_Scenario, 'Blk') ~ "1Fleet TimeBlock",
            str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
          ),
          func = ifelse(str_detect(EM_Scenario, "G"), "Gamma", "Logistic"),
          abbrev = case_when(
            str_detect(model_type, "2Fleet") ~ "2Fleet",
            str_detect(model_type, "TimeInvar") ~ "1Fleet TimeInvar",
            str_detect(model_type, "Block") ~ "1Fleet Block",
            str_detect(model_type, "Random") ~ "1Fleet RandWlkPar",
            str_detect(model_type, "Semi") ~ "1Fleet SemiPar"
          ),
          abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet TimeInvar",
                                             "1Fleet Block", "1Fleet RandWlkPar",
                                             "1Fleet SemiPar"))) # functional form


pdf(here("figs", "Reverse_Selex_Figs", "Fig2_SSB.pdf"), width = 17, height = 19)
ggplot(ssb_fast_re_om %>% filter(time_comp %in% c("Fleet Trans End", "Terminal"), Dat_Qual == "High"),
       aes(x = year, y = median, color = func, group = interaction(func, time_comp), alpha = time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = func, alpha = time_comp), color = NA) +
  geom_line(aes(color = func, linetype = time_comp, size = time_comp), alpha = 1) +
  facet_grid(abbrev ~ OM_Scenario) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  scale_linetype_manual(values = c("Fleet Trans End" = 2, "Terminal" = 1)) +
  scale_alpha_manual(values = c("Fleet Trans End" = 0.25, "Terminal" = 0.35)) +
  scale_size_manual(values = c("Fleet Trans End" = 1.25, "Terminal" = 1.5)) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in Spawning Stock Biomass", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
dev.off()


# Figure 3 (RE ABC Fast) --------------------------------------------------

pdf(here("figs", "Reverse_Selex_Figs", "Fig3_ABC.pdf"), width = 25, height = 13)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, "Fast_"),
                             type == "ABC", time_comp != "Fleet Intersect"),
         aes(x = factor(abbrev), y = median, color = func, fill = func,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    facet_grid(time_comp~OM_Scenario, scales = "free_x") +
    scale_color_manual(values = c("#E69F00", "#0072B2")) +
    scale_fill_manual(values = c("#E69F00", "#0072B2")) +
    scale_x_discrete(guide = guide_axis(angle = 0)) +
    labs(x = "Estimation Models", y = "Relative Error in ABC", 
         fill = "Functional Form", color = "Functional Form") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          axis.title = element_text(size = 17),
          axis.text= element_text(size = 15),
          strip.text = element_text(size = 17),
          axis.text.y = element_text(angle = 90),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17)) +
    coord_cartesian(ylim = c(-0.5, 0.5)))
dev.off()

# Figure 4 (RE SSB Slow) --------------------------------------------------

# Filter to fast scenario, and differentiate models
ssb_slow_re_om = ts_re_om %>% 
  filter(par_name == "Spawning Stock Biomass", str_detect(OM_Scenario, "Slow")) %>% 
  mutate( OM_Scenario = factor(OM_Scenario, levels = slow_om_plot_order),
          model_type = case_when( # differenitate model types
            str_detect(EM_Scenario, 'RW') ~ "1Fleet Random Walk",
            str_detect(EM_Scenario, 'SP') ~ "1Fleet Semi-Parametric",
            str_detect(EM_Scenario, 'TI') ~ "1Fleet TimeInvar",
            str_detect(EM_Scenario, 'Blk') ~ "1Fleet TimeBlock",
            str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
          ),
          func = ifelse(str_detect(EM_Scenario, "G"), "Gamma", "Logistic"),
          abbrev = case_when(
            str_detect(model_type, "2Fleet") ~ "2Fleet",
            str_detect(model_type, "TimeInvar") ~ "1Fleet TimeInvar",
            str_detect(model_type, "Block") ~ "1Fleet Block",
            str_detect(model_type, "Random") ~ "1Fleet RandWlkPar",
            str_detect(model_type, "Semi") ~ "1Fleet SemiPar"
          ),
          abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet TimeInvar",
                                             "1Fleet Block", "1Fleet RandWlkPar",
                                             "1Fleet SemiPar"))) # functional form

pdf(here("figs", "Reverse_Selex_Figs", "Fig4_SSB.pdf"), width = 17, height = 19)
ggplot(ssb_slow_re_om %>% filter(time_comp %in% c("Fleet Trans End", "Terminal"), Dat_Qual == "High"),
       aes(x = year, y = median, color = func, group = interaction(func, time_comp), alpha = time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = func, alpha = time_comp), color = NA) +
  geom_line(aes(color = func, linetype = time_comp, size = time_comp), alpha = 1) +
  facet_grid(abbrev ~ OM_Scenario) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  scale_linetype_manual(values = c("Fleet Trans End" = 2, "Terminal" = 1)) +
  scale_alpha_manual(values = c("Fleet Trans End" = 0.25, "Terminal" = 0.35)) +
  scale_size_manual(values = c("Fleet Trans End" = 1.25, "Terminal" = 1.5)) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in Spawning Stock Biomass", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
dev.off()

# Figure 5 (RE ABC Slow) --------------------------------------------------
pdf(here("figs", "Reverse_Selex_Figs", "Fig5_ABC.pdf"), width = 25, height = 13)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, "Slow_"),
                             type == "ABC", time_comp != "Fleet Intersect"),
         aes(x = factor(abbrev), y = median, color = func, fill = func,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    facet_grid(time_comp~OM_Scenario, scales = "free_x") +
    scale_color_manual(values = c("#E69F00", "#0072B2")) +
    scale_fill_manual(values = c("#E69F00", "#0072B2")) +
    scale_x_discrete(guide = guide_axis(angle = 0)) +
    labs(x = "Estimation Models", y = "Relative Error in ABC", 
         fill = "Functional Form", color = "Functional Form") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          axis.title = element_text(size = 17),
          axis.text= element_text(size = 15),
          strip.text = element_text(size = 17),
          axis.text.y = element_text(angle = 90),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17)) +
    coord_cartesian(ylim = c(-0.5, 0.5)))
dev.off()

# Figure S4 (RE Total F Fast) --------------------------------------------------

# Filter to fast scenario, and differentiate models
f_fast_re_om = ts_re_om %>% 
  filter(par_name == "Total Fishing Mortality", str_detect(OM_Scenario, "Fast")) %>% 
  mutate( OM_Scenario = factor(OM_Scenario, levels = fast_om_plot_order),
          model_type = case_when( # differenitate model types
            str_detect(EM_Scenario, 'RW') ~ "1Fleet Random Walk",
            str_detect(EM_Scenario, 'SP') ~ "1Fleet Semi-Parametric",
            str_detect(EM_Scenario, 'TI') ~ "1Fleet TimeInvar",
            str_detect(EM_Scenario, 'Blk') ~ "1Fleet TimeBlock",
            str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
          ),
          func = ifelse(str_detect(EM_Scenario, "G"), "Gamma", "Logistic"),
          abbrev = case_when(
            str_detect(model_type, "2Fleet") ~ "2Fleet",
            str_detect(model_type, "TimeInvar") ~ "1Fleet TimeInvar",
            str_detect(model_type, "Block") ~ "1Fleet Block",
            str_detect(model_type, "Random") ~ "1Fleet RandWlkPar",
            str_detect(model_type, "Semi") ~ "1Fleet SemiPar"
          ),
          abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet TimeInvar",
                                             "1Fleet Block", "1Fleet RandWlkPar",
                                             "1Fleet SemiPar"))) # functional form


pdf(here("figs", "Reverse_Selex_Figs", "FigS4_F.pdf"), width = 17, height = 19)
ggplot(f_fast_re_om %>% filter(time_comp %in% c("Fleet Trans End", "Terminal"), Dat_Qual == "High"),
       aes(x = year, y = median, color = func, group = interaction(func, time_comp), alpha = time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = func, alpha = time_comp), color = NA) +
  geom_line(aes(color = func, linetype = time_comp, size = time_comp), alpha = 1) +
  facet_grid(abbrev ~ OM_Scenario) +
  coord_cartesian(ylim = c(-1, 1)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  scale_linetype_manual(values = c("Fleet Trans End" = 2, "Terminal" = 1)) +
  scale_alpha_manual(values = c("Fleet Trans End" = 0.25, "Terminal" = 0.35)) +
  scale_size_manual(values = c("Fleet Trans End" = 1.25, "Terminal" = 1.5)) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in Total Fishing Mortality", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
dev.off()


# Figure S5 (RE Total F Slow) --------------------------------------------------

# Filter to fast scenario, and differentiate models
f_slow_re_om = ts_re_om %>% 
  filter(par_name == "Total Fishing Mortality", str_detect(OM_Scenario, "Slow")) %>% 
  mutate( OM_Scenario = factor(OM_Scenario, levels = slow_om_plot_order),
          model_type = case_when( # differenitate model types
            str_detect(EM_Scenario, 'RW') ~ "1Fleet Random Walk",
            str_detect(EM_Scenario, 'SP') ~ "1Fleet Semi-Parametric",
            str_detect(EM_Scenario, 'TI') ~ "1Fleet TimeInvar",
            str_detect(EM_Scenario, 'Blk') ~ "1Fleet TimeBlock",
            str_detect(EM_Scenario, '2Fl') ~ "2Fleet TimeInvar"
          ),
          func = ifelse(str_detect(EM_Scenario, "G"), "Gamma", "Logistic"),
          abbrev = case_when(
            str_detect(model_type, "2Fleet") ~ "2Fleet",
            str_detect(model_type, "TimeInvar") ~ "1Fleet TimeInvar",
            str_detect(model_type, "Block") ~ "1Fleet Block",
            str_detect(model_type, "Random") ~ "1Fleet RandWlkPar",
            str_detect(model_type, "Semi") ~ "1Fleet SemiPar"
          ),
          abbrev = factor(abbrev, levels = c("2Fleet", "1Fleet TimeInvar",
                                             "1Fleet Block", "1Fleet RandWlkPar",
                                             "1Fleet SemiPar"))) # functional form

pdf(here("figs", "Reverse_Selex_Figs", "FigS5_F.pdf"), width = 17, height = 19)
ggplot(f_slow_re_om %>% filter( Dat_Qual == "High"),
       aes(x = year, y = median, color = func, group = interaction(func, time_comp), alpha = time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = func, alpha = time_comp), color = NA) +
  geom_line(aes(color = func, linetype = time_comp, size = time_comp), alpha = 1) +
  facet_grid(abbrev ~ OM_Scenario) +
  coord_cartesian(ylim = c(-1, 1)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  scale_linetype_manual(values = c("Fleet Intersect" = 3, "Fleet Trans End" = 2, "Terminal" = 1)) +
  scale_alpha_manual(values = c("Fleet Intersect" = 0.1, "Fleet Trans End" = 0.25, "Terminal" = 0.35)) +
  scale_size_manual(values = c("Fleet Trans End" = 1.25, "Terminal" = 1.5)) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in Total Fishing Mortality", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
dev.off()

