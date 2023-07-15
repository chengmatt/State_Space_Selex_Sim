# Purpose: To plot simulation runs for manuscript purposes
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 6/8/23

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
ts_df <- data.table::fread(here("output", "TimeSeries_Summary.csv")) # time series dataframe
param_df <- data.table::fread(here("output", "Parameter_Summary.csv")) # parameters
ts_re_df <- data.table::fread(here("output", "TimeSeries_RE.csv")) # time series relative error (converged runs only)
ts_te_df <- data.table::fread(here("output", "TimeSeries_TE.csv")) # time series total error (converged runs only)
ts_are_df <- data.table::fread(here("output", "TimeSeries_ARE.csv")) # time series abs relative error (converged runs only)

# Selectivity and Comps
om_slx_df <- data.table::fread(here("output", "OM_Fish_Selex.csv")) # OM Selectivity Values
pop_sel_om <- data.table::fread(here("output", "Pop_Selex_OM.csv")) # OM Pop'n Selectivity Values
pop_sel_em <- data.table::fread(here("output", "Pop_Selex_EM.csv")) # EM Pop'n Selectivity Values

# SPR from population selex
om_spr_df <- data.table::fread(here("output", "OM_SPR.csv")) # SPR from Pop Selex OM
em_spr_df <- data.table::fread(here("output", "EM_SPR.csv")) # SPR from Pop Selex EM

# SPR Maturity Sensitivities
om_spr_mat_df <- data.table::fread(here("output", "OM_SPR_MatSens.csv")) # SPR from Pop Selex OM
em_spr_mat_df <- data.table::fread(here("output", "EM_SPR_MatSens.csv")) # SPR from Pop Selex EM

# Unique oms and other components
unique_oms <- unique(param_df$OM_Scenario) # unique oms
ts_pars <- unique(ts_re_df$par_name) # unique parameter names

# Timeblock sensitivity models
blk_sens_mods <- "Blk\\b|Blk_1|Blk_3|Blk_5|Blk_-1|Blk_-3|Blk_-5"
all_models <- "2Fl_LL|2Fl_LGam|1Fl_L_TI|1Fl_Gam_TI|1Fl_LL_Blk\\b|1Fl_LGam_Blk\\b|1Fl_L_RW_1.25|1Fl_Gam_RW_2.0"
gen_om <- c("Fast", "Slow") # general OMs

# Get plot order for EM models
plot_order <- c("2Fl_LL", "2Fl_LG",  "1Fl_L_TI", "1Fl_G_TI", 
                "1Fl_LL_Blk", "1Fl_LG_Blk", "1Fl_L_RW", "1Fl_G_RW")

# Fast and slow plot order
fast_om_plot_order <- c("Fast_LL", "Fast_LG_O", "Fast_LG_Y")
slow_om_plot_order <- c("Slow_LL", "Slow_LG_O", "Slow_LG_Y")

# Sensitivity block models order
blk_sens_order <- c("1Fl_LL_Blk_-5", "1Fl_LL_Blk_-3", "1Fl_LL_Blk_-1", "1Fl_LL_Blk",
                    "1Fl_LL_Blk_1", "1Fl_LL_Blk_3", "1Fl_LL_Blk_5",
                    "1Fl_LGam_Blk_-5", "1Fl_LGam_Blk_-3", "1Fl_LGam_Blk_-1", "1Fl_LGam_Blk",
                    "1Fl_LGam_Blk_1", "1Fl_LGam_Blk_3", "1Fl_LGam_Blk_5")


# Main Figures ------------------------------------------------------------

### Figure 2 (Fast SSB Plot) ----------------------------------------------------------------

# Relative error of time series
ts_re_om <- ts_re_df %>% filter(str_detect(EM_Scenario, all_models)) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
  EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"))

# Set order for plot
order <- vector()
for(o in 1:length(plot_order)) {
  order[o] <- which(grepl(plot_order[o], x = unique(ts_re_om$EM_Scenario)))
} # end o loop

# Now relevel factor for organizing plot
ts_re_om <- ts_re_om %>% 
  mutate(EM_Scenario = factor(EM_Scenario, levels = unique(EM_Scenario)[order]),
         time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")))

pdf(here("figs", "Manuscript_Figures", "Fig2_FastSSB_High.pdf"), width = 15, height = 8)
print(
  ggplot(ts_re_om %>% filter(Dat_Qual == "High",
                             par_name == "Spawning Stock Biomass",
                             str_detect(OM_Scenario, "Fast")) %>% 
           mutate(OM_Scenario = factor(OM_Scenario, levels = fast_om_plot_order)),
         aes(x = year, y = median))  +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = time_comp, group = time_comp), alpha = 0.35) +
    geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    facet_grid(OM_Scenario~EM_Scenario) +
    coord_cartesian(ylim = c(-0.5,0.5)) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    labs(x = "Year", y = "Relative Error in SSB", 
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top",
          title = element_text(size = 20),
          axis.text = element_text(size = 13), 
          strip.text = element_text(size = 13)) 
)

dev.off()

### Figure 3 (Slow SSB Plot) ----------------------------------------------------------------

pdf(here("figs", "Manuscript_Figures", "Fig3_SlowSSB_High.pdf"), width = 15, height = 8)
print(
  ggplot(ts_re_om %>% filter(Dat_Qual == "High",
                             par_name == "Spawning Stock Biomass",
                             str_detect(OM_Scenario, "Slow")) %>% 
           mutate(OM_Scenario = factor(OM_Scenario, levels = slow_om_plot_order)),
         aes(x = year, y = median))  +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = time_comp, group = time_comp), alpha = 0.35) +
    geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    facet_grid(OM_Scenario~EM_Scenario) +
    coord_cartesian(ylim = c(-0.5,0.5)) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    labs(x = "Year", y = "Relative Error in SSB", 
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top",
          title = element_text(size = 20),
          axis.text = element_text(size = 13), 
          strip.text = element_text(size = 13)) 
)
dev.off()

### Figure 4 (Fast Parameter Summary Plot) ----------------------------------------------------------------

# Filter to relevant components for parameters
om_scenario_params <- param_df %>% filter(str_detect(EM_Scenario, all_models), 
                                          type %in% c("F_0.4", "ABC")) %>%
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), OM_Scenario = str_remove(OM_Scenario, "_High|_Low"))

# Get Terminal SSB and Biomass
term_biom_ssb <- ts_re_df %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low")) %>% 
  filter(
    (year == c(50) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Terminal") |
      (year == c(70) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Terminal") |
      (year == c(30) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Fleet Trans End") |
      (year == c(50) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Fleet Trans End") |
      (year == c(27) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Fleet Intersect") |
      (year == c(40) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Fleet Intersect"),
    par_name %in% c("Spawning Stock Biomass"),
    str_detect(EM_Scenario, all_models)) %>% 
  dplyr::select(OM_Scenario, EM_Scenario, 
                time_comp, par_name, Dat_Qual,
                median, lwr_95, upr_95) %>% 
  rename(type = par_name)

# Point ranges for relative error and total error
pt_rg_re <- om_scenario_params %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, type, Dat_Qual) %>% 
  summarize(median = median(RE), 
            lwr_95 = quantile(RE, 0.025),
            upr_95 =  quantile(RE, 0.975))

# Get terminal SSB and biomass
pt_rg_re <- rbind(pt_rg_re, term_biom_ssb)

# Clarify names
pt_rg_re <- pt_rg_re %>% 
  mutate(type = case_when(
    type == "ABC" ~ "ABC",
    type == "F_0.4" ~ "F40%",
    type == "Spawning Stock Biomass" ~ 'Terminal SSB',
  ),
  EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
  OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)),
  EM_Scenario = factor(EM_Scenario, levels = plot_order))


pdf(here("figs", "Manuscript_Figures", "Fig4_FastParam_High.pdf"), width = 15, height = 8)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, "Fast")),
         aes(x = factor(EM_Scenario), y = median, color = time_comp, fill = time_comp,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), 
                    size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    geom_vline(xintercept = c(seq(1.5, 7.5, 1)), lwd = 0.25, alpha = 0.75) +
    facet_grid(type~OM_Scenario, scales = "free_x") +
    coord_cartesian(ylim = c(-0.5,0.5)) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = "EMs", y = "Relative Error", 
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          strip.text = element_text(size = 15))
)
dev.off()  

### Figure 5 (Slow Parameter Summary Plot) ----------------------------------------------------------------

pdf(here("figs", "Manuscript_Figures", "Fig5_SlowParam_High.pdf"), width = 15, height = 8)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, "Slow")),
         aes(x = factor(EM_Scenario), y = median, color = time_comp, fill = time_comp,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), 
                    size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    geom_vline(xintercept = c(seq(1.5, 7.5, 1)), lwd = 0.25, alpha = 0.75) +
    facet_grid(type~OM_Scenario, scales = "free_x") +
    coord_cartesian(ylim = c(-0.5,0.5)) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = "EMs", y = "Relative Error", 
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          strip.text = element_text(size = 15))
)
dev.off() 

### Figure 6 (Fast Selectivity Plot) ----------------------------------------------------------------

plot_df = pop_sel_em %>% 
  filter(str_detect(EM_Scenario, all_models)) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), 
  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
  EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
  OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)),
  EM_Scenario = factor(EM_Scenario, levels = plot_order))

plot_om_df = pop_sel_om %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ), 
  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)))

pdf(here("figs", "Manuscript_Figures", "Fig6_SelexFast_High.pdf"), width = 13, height = 5)
ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                str_detect(OM_Scenario, "Fast")),
             mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
              fill = "red") +
  geom_line(plot_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal", 
                                  str_detect(OM_Scenario, "Fast")),
            mapping = aes(x = Age, y = Selex), color = "black",
            position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
  facet_grid(OM_Scenario ~ EM_Scenario) +
  labs(x = "Age", y = "Normalized Population Selectivity (Terminal Year)", fill = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "none", 
        title = element_text(size = 20),
        strip.text = element_text(size = 15))
dev.off()

### Figure 7 (Slow Selectivity Plot) ----------------------------------------------------------------

pdf(here("figs", "Manuscript_Figures", "Fig7_SelexSlow_High.pdf"), width = 13, height = 5)
ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                               str_detect(OM_Scenario, "Slow")),
            mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Slow")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
              fill = "red") +
  geom_line(plot_om_df %>% filter(Dat_Qual == "High", Sex == "Female", time_comp == "Terminal", 
                                  str_detect(OM_Scenario, "Slow")),
            mapping = aes(x = Age, y = Selex), color = "black",
            position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
  facet_grid(OM_Scenario ~ EM_Scenario) +
  labs(x = "Age", y = "Normalized Population Selectivity (Terminal Year)", fill = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "none", 
        title = element_text(size = 20),
        strip.text = element_text(size = 15))
dev.off()

### Figure 8 (Variant of Minimax Solution) ----------------------------------------------------------------

plot_ts_df = ts_df %>% 
  filter(str_detect(EM_Scenario, all_models),
         type == "Spawning Stock Biomass") %>% 
  mutate(
    Dat_Qual = case_when(
      str_detect(OM_Scenario, "High") ~ 'Data Quality: High',
      str_detect(OM_Scenario, "Low") ~ 'Data Quality: Low'
    ), 
    OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
    EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
    EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
    OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order)),
    ARE = abs(mle_val - t) / t,
    time_comp = case_when(
      str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
      str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
      str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
    ), 
    EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'),
    time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End",
                                             "Terminal")),
    EM_Scenario = factor(EM_Scenario, levels = plot_order)
  )

# Get median values
median_vals = plot_ts_df %>% 
  group_by(EM_Scenario, time_comp) %>%
  summarize(median = median(ARE))

pdf(here("figs", "Manuscript_Figures", "Fig8_RobustSoln.pdf"), width = 19, height = 10)
# Plot
ggplot(plot_ts_df, aes(x = EM_Scenario, y = ARE, fill = EM_Scenario)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.75) +
  facet_wrap(~time_comp) +
  geom_text(data = median_vals, aes(y = median, label = round(median, 3)), size = 5, vjust = -0.5) +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_manual(values = c("2Fl_LG" = "green", "white")) +
  labs(x = "EMs", y = "Absolute Relative Error in SSB") +
  theme_matt() +
  theme(legend.position = "none")
dev.off()


### Figure 8 (Minimax Solution) ----------------------------------------------

minmax_df <- ts_are_df %>% 
  # Selective filtering down here to terminal years
  filter(str_detect(EM_Scenario, all_models),
         # Get terminal year of peels/EMs
         (year == c(50) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Terminal") |
           (year == c(70) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Terminal") |
           (year == c(30) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Fleet Trans End") |
           (year == c(50) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Fleet Trans End") |
           (year == c(27) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Fleet Intersect") |
           (year == c(40) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Fleet Intersect"),
         par_name == "Spawning Stock Biomass") %>% 
  data.frame() %>%
  group_by(EM_Scenario, par_name) %>%
  mutate(max_median = max(median)) %>% # find the maximum median MARE
  ungroup() %>% 
  mutate(EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
         EM_Scenario = str_replace(EM_Scenario, "Gam", "G"),
         EM_Scenario = factor(EM_Scenario, levels = rev(plot_order)))

# Find the minimum maximum median value
min_max_medians <- minmax_df %>%
  group_by(par_name) %>%
  summarize(min_max_medians = min(max_median))

# Now left join this
minmax_df = minmax_df %>% 
  left_join(min_max_medians, by = c("par_name"))

# Plot!
pdf(file = here("figs", "Manuscript_Figures", "MinMax_Summary.pdf"), width = 25, height = 8)
  print(
    ggplot(minmax_df,
           aes(x = factor(OM_Scenario), y = factor(EM_Scenario), 
               fill = median, label = round(median, 3))) +
      geom_tile(alpha = 0.35) +
      facet_wrap(~time_comp, scales = "free_x") +
      geom_text(color = ifelse(minmax_df$median ==  minmax_df$max_median &
                                 minmax_df$median != minmax_df$min_max_medians, "red", 
                               ifelse(minmax_df$median ==  minmax_df$min_max_medians, "green2", "black")),
                size = 5.5) +
      scale_fill_distiller(palette = "Spectral", direction = -1) + 
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      labs(x = "OM Scenario", y = "EM Scenario", fill = "Median MARE") +
      theme_test() +
      theme(legend.position = "top",
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            axis.text= element_text(size = 13, color = "black"),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 15),
            legend.key.width = unit(1, "cm"))
  )
dev.off()

# Supplemental Figures ----------------------------------------------------
### Supp Figure 1 and 2 (AIC) -----------------------------------------------------------
AIC_df <- AIC_df %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'Data Quality: High',
    str_detect(OM_Scenario, "Low") ~ 'Data Quality: Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order))
  )

# Two fleet AIC df
twofleet_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "2Fl_")) %>% 
  group_by(sim, OM_Scenario, time_comp, Dat_Qual) %>% 
  mutate(min = min(AIC), min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "2 Fleet")

pdf(here("figs", "Manuscript_Figures", "S1_AIC_2Fl.pdf"), width = 13, height = 5)
ggplot(twofleet_aic, aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC,
                             label = round(n_minAIC, 2))) +
    geom_tile(alpha = 0.65) +
    geom_text(size = 5) +
    facet_grid(Dat_Qual~time_comp, scales = "free") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_fill_distiller(palette = "Spectral", direction = 1) + 
    labs(x = "OM Scenario", y = "EMs", fill = "Proportion of models with lowest AIC") +
    theme_test() +
    theme(legend.position = "top",
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text= element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15),
          legend.key.width = unit(1, "cm"))
dev.off()

# filter to 1 fleet models
onefleet_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "1Fl_"),
         str_detect(EM_Scenario, all_models),
         !str_detect(EM_Scenario, blk_sens_mods)) %>% 
  group_by(sim, OM_Scenario, time_comp, Dat_Qual) %>% 
  mutate(min = min(AIC),
         min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "1 Fleet",
         EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
         EM_Scenario = str_replace(EM_Scenario, "Gam", "G")) 

# One fleet plot
pdf(here("figs", "Manuscript_Figures", "S2_AIC_1Fl.pdf"), width = 13, height = 7)
ggplot(onefleet_aic, aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC,
                             label = round(n_minAIC, 2))) +
    geom_tile(alpha = 0.65) +
    geom_text(size = 6) +
    facet_grid(Dat_Qual~time_comp, scales = "free_y") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_fill_distiller(palette = "Spectral", direction = 1) + 
    labs(x = "OM Scenario", y = "EMs", fill = "Proportion of models with lowest AIC") +
    theme_test() +
    theme(legend.position = "top",
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text= element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15),
          legend.key.width = unit(1, "cm"))
dev.off()

### Supp Figure 3 (Convergence) -----------------------------------------------------------

# Convergence for all models we want to look at
conv_stat <- AIC_df %>% 
  filter(str_detect(EM_Scenario, all_models)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(converged = sum(conv == "Converged")/200) %>% 
  mutate(EM_Scenario = str_remove(EM_Scenario, "_1.25|_2.0"),
         EM_Scenario = str_replace(EM_Scenario, "Gam", "G")) 

# Set order for plot
order <- vector()
for(o in 1:length(plot_order)) {
  order[o] <- which(grepl(plot_order[o], x = unique(conv_stat$EM_Scenario)))
} # end o loop

# relevel factors
conv_stat$EM_Scenario <- factor(conv_stat$EM_Scenario, levels = c(unique(conv_stat$EM_Scenario)[order]))

pdf(here("figs", "Manuscript_Figures", "S3_Convergence.pdf"), width = 25, height = 10)
print(
  ggplot(conv_stat %>%
           mutate(OM_Scenario = factor(OM_Scenario, levels = c(fast_om_plot_order, slow_om_plot_order))),
         mapping = aes(x = EM_Scenario, y = converged * 100, group = time_comp, fill = time_comp))  +
    geom_col(position = position_dodge2(width = 1), 
             alpha = 0.7, color = "black") +
    geom_hline(aes(yintercept = 50), col = "black", lty = 2, size = 0.7, alpha = 1) +
    facet_grid(Dat_Qual~OM_Scenario,  scales = "free_x") +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) + 
    theme_matt() +
    labs(fill = "Assessment Period", y = "Convergence Rate",
         x = "EMs") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          strip.text = element_text(size = 15))
)
dev.off()

### Supp Figure 4 (Fast SSB Plot - Low Data) -----------------------------------------------------------
pdf(here("figs", "Manuscript_Figures", "S4_FastSSB_Low.pdf"), width = 15, height = 8)
print(
  ggplot(ts_re_om %>% filter(Dat_Qual == "Low",
                             par_name == "Spawning Stock Biomass",
                             str_detect(OM_Scenario, "Fast")) %>% 
           mutate(OM_Scenario = factor(OM_Scenario, levels = fast_om_plot_order)), aes(x = year, y = median))  +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = time_comp, group = time_comp), alpha = 0.35) +
    geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    facet_grid(OM_Scenario~EM_Scenario) +
    coord_cartesian(ylim = c(-0.8,0.8)) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    labs(x = "Year", y = "Relative Error in SSB", 
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top",
          title = element_text(size = 20),
          axis.text = element_text(size = 13), 
          strip.text = element_text(size = 13)) 
)
dev.off()

### Supp Figure 5 (Slow SSB Plot - Low Data) -----------------------------------------------------------
pdf(here("figs", "Manuscript_Figures", "S5_SlowSSB_Low.pdf"), width = 15, height = 8)
print(
  ggplot(ts_re_om %>% filter(Dat_Qual == "Low",
                             par_name == "Spawning Stock Biomass",
                             str_detect(OM_Scenario, "Slow")) %>% 
           mutate(OM_Scenario = factor(OM_Scenario, levels = slow_om_plot_order)), aes(x = year, y = median))  +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = time_comp, group = time_comp), alpha = 0.35) +
    geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
    facet_grid(OM_Scenario~EM_Scenario) +
    coord_cartesian(ylim = c(-0.8,0.8)) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    labs(x = "Year", y = "Relative Error in SSB", 
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top",
          title = element_text(size = 20),
          axis.text = element_text(size = 13), 
          strip.text = element_text(size = 13)) 
)
dev.off()



### Supp Figure 6 Fast Total F Plot - High Data -------------------------------------------

pdf(here("figs", "Manuscript_Figures", "S6_TotalF_Fast_High.pdf"), width = 23.5, height = 8)
ggplot(ts_re_om %>% filter(Dat_Qual == "High",
                           par_name == "Total Fishing Mortality",
                           str_detect(OM_Scenario, "Fast")), 
       aes(x = year, y = median, fill = time_comp))  +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.35) +
  geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  facet_grid(OM_Scenario~EM_Scenario) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
  scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
  labs(x = "Year", y = "Relative Error in Fishing Mortality Multiplier", fill = "Assessment Period", color = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "top", title = element_text(size = 20),
        axis.text = element_text(size = 13), strip.text = element_text(size = 13)) 
dev.off()

### Supp Figure 7 Fast Total F Plot - Low Data -------------------------------------------

pdf(here("figs", "Manuscript_Figures", "S7_TotalF_Fast_Low.pdf"), width = 23.5, height = 8)
ggplot(ts_re_om %>% filter(Dat_Qual == "Low",
                           par_name == "Total Fishing Mortality",
                           str_detect(OM_Scenario, "Fast")), 
       aes(x = year, y = median, fill = time_comp))  +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.35) +
  geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  facet_grid(OM_Scenario~EM_Scenario) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
  scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
  labs(x = "Year", y = "Relative Error in Fishing Mortality Multiplier", fill = "Assessment Period", color = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "top", title = element_text(size = 20),
        axis.text = element_text(size = 13), strip.text = element_text(size = 13)) 
dev.off()

### Supp Figure 8 Fast Total F Plot - Low Data -------------------------------------------

pdf(here("figs", "Manuscript_Figures", "S8_TotalF_Slow_High.pdf"), width = 23.5, height = 8)
ggplot(ts_re_om %>% filter(Dat_Qual == "High",
                           par_name == "Total Fishing Mortality",
                           str_detect(OM_Scenario, "Slow")), 
       aes(x = year, y = median, fill = time_comp))  +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.35) +
  geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  facet_grid(OM_Scenario~EM_Scenario) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
  scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
  labs(x = "Year", y = "Relative Error in Fishing Mortality Multiplier", fill = "Assessment Period", color = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "top", title = element_text(size = 20),
        axis.text = element_text(size = 13), strip.text = element_text(size = 13)) 
dev.off()

### Supp Figure 9 Slow Total F Plot - Low Data -------------------------------------------

pdf(here("figs", "Manuscript_Figures", "S8_TotalF_Slow_High.pdf"), width = 23.5, height = 8)
ggplot(ts_re_om %>% filter(Dat_Qual == "Low",
                           par_name == "Total Fishing Mortality",
                           str_detect(OM_Scenario, "Slow")), 
       aes(x = year, y = median, fill = time_comp))  +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.35) +
  geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  facet_grid(OM_Scenario~EM_Scenario) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
  scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
  labs(x = "Year", y = "Relative Error in Fishing Mortality Multiplier", fill = "Assessment Period", color = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "top", title = element_text(size = 20),
        axis.text = element_text(size = 13), strip.text = element_text(size = 13)) 
dev.off()


### Supp Figure 10 (Fast Param Plot - Low Data) -----------------------------------------------------------

pdf(here("figs", "Manuscript_Figures", "S10_FastParam_Low.pdf"), width = 15, height = 8)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "Low", str_detect(OM_Scenario, "Fast")),
         aes(x = factor(EM_Scenario), y = median, color = time_comp, fill = time_comp,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), 
                    size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    geom_vline(xintercept = c(seq(1.5, 7.5, 1)), lwd = 0.25, alpha = 0.75) +
    facet_grid(type~OM_Scenario, scales = "free_x") +
    coord_cartesian(ylim = c(-0.5,0.5)) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = "EMs", y = "Relative Error", 
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          strip.text = element_text(size = 15))
)
dev.off()  

### Supp Figure 11 (Fast Param Plot - Low Data) -----------------------------------------------------------

pdf(here("figs", "Manuscript_Figures", "S11_SlowParam_Low.pdf"), width = 15, height = 8)
print(
  ggplot(pt_rg_re %>% filter(Dat_Qual == "Low", str_detect(OM_Scenario, "Slow")),
         aes(x = factor(EM_Scenario), y = median, color = time_comp, fill = time_comp,
             ymin = lwr_95, ymax = upr_95)) +
    geom_pointrange(position = position_dodge2(width = 0.65), 
                    size = 1, linewidth = 1) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
    geom_vline(xintercept = c(seq(1.5, 7.5, 1)), lwd = 0.25, alpha = 0.75) +
    facet_grid(type~OM_Scenario, scales = "free_x") +
    coord_cartesian(ylim = c(-0.5,0.5)) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 40)]) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = "EMs", y = "Relative Error", 
         fill = "Assessment Period", color = "Assessment Period") +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20),
          strip.text = element_text(size = 15))
)
dev.off()  



### Supp Figure 12 (Fast Selectivity Plot - Low Data) ------------------------

pdf(here("figs", "Manuscript_Figures", "S12_SelexFast_Low.pdf"), width = 13, height = 5)
ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "Low", Sex == "Female", time_comp == "Terminal",
                               str_detect(OM_Scenario, "Fast")),
            mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "Low", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Fast")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
              fill = "red") +
  geom_line(plot_om_df %>% filter(Dat_Qual == "Low", Sex == "Female", time_comp == "Terminal", 
                                  str_detect(OM_Scenario, "Fast")),
            mapping = aes(x = Age, y = Selex), color = "black",
            position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
  facet_grid(OM_Scenario ~ EM_Scenario) +
  labs(x = "Age", y = "Normalized Population Selectivity (Terminal Year)", fill = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "none", 
        title = element_text(size = 20),
        strip.text = element_text(size = 15))
dev.off()

### Supp Figure 13 (Slow Selectivity Plot - Low Data) ------------------------

pdf(here("figs", "Manuscript_Figures", "S13_SelexSlow_Low.pdf"), width = 13, height = 5)
ggplot() +
  geom_line(plot_df %>% filter(Dat_Qual == "Low", Sex == "Female", time_comp == "Terminal",
                               str_detect(OM_Scenario, "Slow")),
            mapping = aes(x = Age, y = Median_Selex), lwd = 1, alpha = 1, color = "red") +
  geom_ribbon(plot_df %>% filter(Dat_Qual == "Low", Sex == "Female", time_comp == "Terminal",
                                 str_detect(OM_Scenario, "Slow")),
              mapping = aes(x = Age, y = Median_Selex, ymin = Lwr_95, ymax = Upr_95), alpha = 0.25,
              fill = "red") +
  geom_line(plot_om_df %>% filter(Dat_Qual == "Low", Sex == "Female", time_comp == "Terminal", 
                                  str_detect(OM_Scenario, "Slow")),
            mapping = aes(x = Age, y = Selex), color = "black",
            position = position_jitter(width = 0.3), lwd = 1, alpha = 1, lty = 2) +
  facet_grid(OM_Scenario ~ EM_Scenario) +
  labs(x = "Age", y = "Normalized Population Selectivity (Terminal Year)", fill = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "none", 
        title = element_text(size = 20),
        strip.text = element_text(size = 15))
dev.off()

### Supp Figure 14 (Time-Block AIC) -----------------------------------------

# Filter to 1 fleet time block models
onefleet_blk_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "Blk"),
         str_detect(OM_Scenario, "Fast")) %>% 
  mutate(
    selex_form = case_when(
      str_detect(EM_Scenario, "LGam") ~ "Logistic-Gamma",
      str_detect(EM_Scenario, "LL") ~ "Logistic-Logistic",
    )) %>% group_by(sim, OM_Scenario, time_comp, selex_form, Dat_Qual) %>% 
  mutate(min = min(AIC),
         min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, selex_form, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "1 Fleet")

pdf(here("figs", "Manuscript_Figures", "S14_AIC_Blk.pdf"), width = 15, height = 8.5)
ggplot(onefleet_blk_aic %>% filter(time_comp == "Terminal"),
       aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC, label = round(n_minAIC, 2))) +
    geom_tile(alpha = 0.65) +
    geom_text(size = 6) +
    facet_grid(selex_form~Dat_Qual, scales = "free") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_fill_distiller(palette = "Spectral", direction = 1) + 
    labs(x = "OM Scenario", y = "EMs", fill = "Proportion of models with lowest AIC") +
    theme_test() +
    theme(legend.position = "top",
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text= element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15),
          legend.key.width = unit(1, "cm"))
dev.off()

### Supp Figure 15 (Time-Block Sensitivity SSB - High) ---------------------------------
# Relative error of time series
ts_re_om <- ts_re_df %>% filter(time_comp %in% c("Terminal", "Fleet Trans End"),
                                str_detect(OM_Scenario, "Fast"),
                                str_detect(EM_Scenario, blk_sens_mods),
                                par_name %in% c("Spawning Stock Biomass",
                                                "Total Fishing Mortality")) %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'Data Quality: High',
    str_detect(OM_Scenario, "Low") ~ 'Data Quality: Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"),
  OM_Scenario = factor(OM_Scenario, levels = c("Fast_LL", "Fast_LG_O", "Fast_LG_Y")),
  EM_Scenario = factor(EM_Scenario, levels = blk_sens_order))

pdf(here("figs", "Manuscript_Figures", "S15_BlkSSB_High.pdf"), width = 23.5, height = 8)
ggplot(ts_re_om %>% filter(Dat_Qual == "Data Quality: High",
                           par_name == "Spawning Stock Biomass"), 
           aes(x = year, y = median, fill = time_comp))  +
      geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.35) +
      geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
      geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
      facet_grid(OM_Scenario~EM_Scenario) +
      coord_cartesian(ylim = c(-0.5,0.5)) +
      scale_color_manual(values = viridis::viridis(n = 50)[c(20, 43)]) +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(20, 43)]) +
      labs(x = "Year", y = "Relative Error in SSB", fill = "Assessment Period", color = "Assessment Period") +
      theme_matt() +
      theme(legend.position = "top", title = element_text(size = 20),
            axis.text = element_text(size = 13), strip.text = element_text(size = 13)) 
dev.off()

### Supp Figure 16 (Time-Block Sensitivity SSB - Low) ---------------------------------

pdf(here("figs", "Manuscript_Figures", "S16_BlkSSB_Low.pdf"), width = 23.5, height = 8)
ggplot(ts_re_om %>% filter(Dat_Qual == "Data Quality: Low",
                           par_name == "Spawning Stock Biomass"), 
       aes(x = year, y = median, fill = time_comp))  +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.35) +
  geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  facet_grid(OM_Scenario~EM_Scenario) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  scale_color_manual(values = viridis::viridis(n = 50)[c(20, 43)]) +
  scale_fill_manual(values = viridis::viridis(n = 50)[c(20, 43)]) +
  labs(x = "Year", y = "Relative Error in SSB", fill = "Assessment Period", color = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "top", title = element_text(size = 20),
        axis.text = element_text(size = 13), strip.text = element_text(size = 13)) 
dev.off()

### Supp Figure 17 (Time-Block Sensitivity Total F - High) ---------------------------------

pdf(here("figs", "Manuscript_Figures", "S17_BlkTotalF_High.pdf"), width = 23.5, height = 8)
ggplot(ts_re_om %>% filter(Dat_Qual == "Data Quality: High",
                           par_name == "Total Fishing Mortality"), 
       aes(x = year, y = median, fill = time_comp))  +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.35) +
  geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  facet_grid(OM_Scenario~EM_Scenario) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  scale_color_manual(values = viridis::viridis(n = 50)[c(20, 43)]) +
  scale_fill_manual(values = viridis::viridis(n = 50)[c(20, 43)]) +
  labs(x = "Year", y = "Relative Error in Fishing Mortality Multiplier", fill = "Assessment Period", color = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "top", title = element_text(size = 20),
        axis.text = element_text(size = 13), strip.text = element_text(size = 13)) 
dev.off()

### Supp Figure 18 (Time-Block Sensitivity Total F - Low) ---------------------------------

pdf(here("figs", "Manuscript_Figures", "S18_BlkTotalF_Low.pdf"), width = 23.5, height = 8)
ggplot(ts_re_om %>% filter(Dat_Qual == "Data Quality: Low",
                           par_name == "Total Fishing Mortality"), 
       aes(x = year, y = median, fill = time_comp))  +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.35) +
  geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
  facet_grid(OM_Scenario~EM_Scenario) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  scale_color_manual(values = viridis::viridis(n = 50)[c(20, 43)]) +
  scale_fill_manual(values = viridis::viridis(n = 50)[c(20, 43)]) +
  labs(x = "Year", y = "Relative Error in Fishing Mortality Multiplier", fill = "Assessment Period", color = "Assessment Period") +
  theme_matt() +
  theme(legend.position = "top", title = element_text(size = 20),
        axis.text = element_text(size = 13), strip.text = element_text(size = 13)) 
dev.off()

