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
all_models <- "2Fl_LL|2Fl_GamL|2Fl_GamGam|2Fl_LGam|1Fl_L_TI|1Fl_Gam_TI|1Fl_LL_Blk\\b|1Fl_GamL_Blk\\b|1Fl_GamGam_Blk\\b|1Fl_LGam_Blk\\b|1Fl_L_RW_1.25|1Fl_Gam_RW_2.0|1Fl_L_SP_0.75|1Fl_Gam_SP_0.75"

# Residual Munging --------------------------------------------------------

### Selex Stuff -------------------------------------------------------------
plot_df = pop_sel_em %>% 
  filter(str_detect(EM_Scenario, all_models)) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), 
  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0|_0.75"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"))

plot_sel_om_df = pop_sel_om %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), 
  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"))

### ABC Parameters ----------------------------------------------------------

# Filter to relevant components for parameters
om_scenario_params <- param_df %>% filter(str_detect(EM_Scenario, all_models), 
                                          type %in% c("F_0.4", "ABC")) %>%
  filter(str_detect(OM_Scenario, "Ext")) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), OM_Scenario = str_remove(OM_Scenario, "_High|_Low")
  )

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
    # OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order))
    ) %>% 
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
  filter(str_detect(OM_Scenario, "Ext")) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G")
  )

# Now relevel factor for organizing plot
ts_re_om <- ts_re_om %>% 
  mutate(time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")))

# Plots -------------------------------------------------------------------
# Figure 2 (RE SSB Fast) --------------------------------------------------

# Filter to fast scenario, and differentiate models
ssb_fast_re_om = ts_re_om %>% 
  filter(par_name == "Spawning Stock Biomass", str_detect(OM_Scenario, "Fast"),
         time_comp == "Terminal") 


pdf(here("figs", "Ext_Figs", "Fig2_SSB.pdf"), width = 15, height = 15)
ggplot(ssb_fast_re_om %>% filter(str_detect(EM_Scenario, "Blk")),
       aes(x = year, y = median, color=time_comp,fill=time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), color = NA, alpha = 0.35) +
  geom_line(alpha = 1) +
  facet_grid(EM_Scenario ~ OM_Scenario) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in Spawning Stock Biomass", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
ggplot(ssb_fast_re_om %>% filter(str_detect(EM_Scenario, "2Fl")),
       aes(x = year, y = median,color=time_comp,fill=time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), color = NA, alpha = 0.35) +
  geom_line(alpha = 1) +
  facet_grid(EM_Scenario ~ OM_Scenario) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in Spawning Stock Biomass", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
ggplot(ssb_fast_re_om %>% filter(str_detect(EM_Scenario, "RW|SP")),
       aes(x = year, y = median,color=time_comp,fill=time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), color = NA, alpha = 0.35) +
  geom_line(alpha = 1) +
  facet_grid(EM_Scenario ~ OM_Scenario) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in Spawning Stock Biomass", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
dev.off()


# Figure 3 (RE ABC Fast) --------------------------------------------------

pdf(here("figs", "Ext_Figs", "Fig3_ABC.pdf"), width = 15, height = 15)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, "Fast_"),
                             str_detect(EM_Scenario, "1Fl"), str_detect(EM_Scenario, "Gam"),
                             time_comp == "Terminal"),
         aes(x = OM_Scenario, y = median, color=time_comp,fill=time_comp,  ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge(width = 1)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    facet_grid(type~EM_Scenario, scales = "free") +
    coord_cartesian(ylim = c(-1, 1)) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
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
          legend.title = element_text(size = 17))
  )
dev.off()

# Figure S4 (RE Total F Fast) --------------------------------------------------

# Filter to fast scenario, and differentiate models
f_fast_re_om = ts_re_om %>% 
  filter(par_name == "Total Fishing Mortality", str_detect(OM_Scenario, "Fast"),
         time_comp == "Terminal") 

pdf(here("figs", "Ext_Figs", "FigS4_F.pdf"), width = 15, height = 15)
ggplot(f_fast_re_om %>% filter(str_detect(EM_Scenario, "Blk")),
       aes(x = year, y = median,color=time_comp,fill=time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), color = NA, alpha = 0.35) +
  geom_line(alpha = 1) +
  facet_grid(EM_Scenario ~ OM_Scenario) +
  coord_cartesian(ylim = c(-1, 1)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in F", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
ggplot(f_fast_re_om %>% filter(str_detect(EM_Scenario, "2Fl")),
       aes(x = year, y = median,color=time_comp,fill=time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), color = NA, alpha = 0.35) +
  geom_line(alpha = 1) +
  facet_grid(EM_Scenario ~ OM_Scenario) +
  coord_cartesian(ylim = c(-1, 1)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in F", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
ggplot(f_fast_re_om %>% filter(str_detect(EM_Scenario, "RW|SP")),
       aes(x = year, y = median,color=time_comp,fill=time_comp)) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), color = NA, alpha = 0.35) +
  geom_line(alpha = 1) +
  facet_grid(EM_Scenario ~ OM_Scenario) +
  coord_cartesian(ylim = c(-1, 1)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  theme_matt() +
  theme(legend.key.width = unit(1,"cm")) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Relative Error in F", fill = "Functional Form",
       color = "Functional Form", lty = "Time Component", alpha = "Time Component",
       size = "Time Component")
dev.off()



# Selectivity -------------------------------------------------------------
pdf(here("figs", "Ext_Figs", "Selex_Females.pdf"), width = 15, height = 15)

ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female",
                               str_detect(OM_Scenario, "Ext"), time_comp == "Terminal",
                               str_detect(EM_Scenario, "Blk")),
            mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1) +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Ext"), str_detect(EM_Scenario, "Blk")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25) +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                      str_detect(OM_Scenario, "Ext")),
            mapping = aes(x = Age, y = Selex), color = "black",
            position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
  facet_grid(EM_Scenario ~ OM_Scenario) +
  labs(x = "Age", y = "Selectivity (Females)", fill = "Functional Form", color = "Functional Form") +
  theme_matt() 

ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female",
                               str_detect(OM_Scenario, "Ext"), time_comp == "Terminal",
                               str_detect(EM_Scenario, "2Fl")),
            mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1) +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Ext"), str_detect(EM_Scenario, "2Fl")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25) +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                      str_detect(OM_Scenario, "Ext")),
            mapping = aes(x = Age, y = Selex), color = "black",
            position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
  facet_grid(EM_Scenario ~ OM_Scenario) +
  labs(x = "Age", y = "Selectivity (Females)", fill = "Functional Form", color = "Functional Form") +
  theme_matt() 
ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female",
                               str_detect(OM_Scenario, "Ext"), time_comp == "Terminal",
                               str_detect(EM_Scenario, "SP|RW")),
            mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1) +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Ext"), 
                                 str_detect(EM_Scenario, "SP|RW")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25) +
  geom_line(plot_sel_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                      str_detect(OM_Scenario, "Ext")),
            mapping = aes(x = Age, y = Selex), color = "black",
            position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
  facet_grid(EM_Scenario ~ OM_Scenario) +
  labs(x = "Age", y = "Selectivity (Females)", fill = "Functional Form", color = "Functional Form") +
  theme_matt() 
dev.off()
