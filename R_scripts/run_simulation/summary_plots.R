# Purpose: To make summary plots of simulation runs
# Creator: Matthew LH. Cheng
# Date 3/26/23

# Set up -----------------------------------------------------------------

library(here)
library(tidyverse)
library(ggh4x)
library(data.table)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

# Paths
om_scenario_path <- here("output", "OM_Scenarios") # path to OM folder
dir_out <- here("output", "Summary_Plots") # path to output folder
dir.create(dir_out)

# read in csvs
# Parameter and Time Series Summaries
AIC_df <- data.table::fread(here("output", "AIC_Convergence_Summary.csv")) # time series total error (converged runs only)
param_df <- data.table::fread(here("output", "Parameter_Summary.csv")) # parameters
# ts_df <- data.table::fread(here("output", "TimeSeries_Summary.csv")) # time series values 
ts_re_df <- data.table::fread(here("output", "TimeSeries_RE.csv")) # time series relative error (converged runs only)
ts_te_df <- data.table::fread(here("output", "TimeSeries_TE.csv")) # time series total error (converged runs only)
ts_are_df <- data.table::fread(here("output", "TimeSeries_ARE.csv")) # time series abs relative error (converged runs only)
NAA_df <- data.table::fread(here("output", "NAA_Summary.csv")) # numbers at age 

# Selectivity and Comps
om_slx_df <- data.table::fread(here("output", "OM_Fish_Selex.csv")) # OM Selectivity Values
pop_sel_om <- data.table::fread(here("output", "Pop_Selex_OM.csv")) # OM Pop'n Selectivity Values
pop_sel_em <- data.table::fread(here("output", "Pop_Selex_EM.csv")) # EM Pop'n Selectivity Values
om_spr_df <- data.table::fread(here("output", "OM_SPR.csv")) # SPR from Pop Selex OM
em_spr_df <- data.table::fread(here("output", "EM_SPR.csv")) # SPR from Pop Selex EM
om_spr_mat_df <- data.table::fread(here("output", "OM_SPR_MatSens.csv")) # SPR from Pop Selex OM
em_spr_mat_df <- data.table::fread(here("output", "EM_SPR_MatSens.csv")) # SPR from Pop Selex EM
om_comps_df <- data.table::fread(here("output", "OM_Fish_Comps.csv")) # OM Composition values
em_comps_df <- data.table::fread(here("output", "EM_Fish_Comps.csv")) # EM Composition values

# Unique oms and other components
unique_oms <- unique(param_df$OM_Scenario) # unique oms
ts_pars <- unique(ts_re_df$par_name) # unique parameter names

# Timeblock sensitivity models
blk_sens_mods <- "Blk_1|Blk_2|Blk_3|Blk_4|Blk_5|Blk_-1|Blk_-2|Blk_-3|Blk_-4|Blk_-5"
all_models <- "2Fl_LL|2Fl_LGam|2Fl_LExpL|1Fl_L_TI|1Fl_Gam_TI|1Fl_ExpL_TI|1Fl_LL_Blk\\b|1Fl_LGam_Blk\\b|1Fl_LExpL_Blk\\b|1Fl_L_RW_1.25|1Fl_Gam_RW_2.0|1Fl_ExpL_RW_2.0"
gen_om <- c("Fast_LL", "Slow_LL", "Fast_LG", "Slow_LG") # general OMs

# Get plot order for EM models
plot_order <- c("2Fl_LL", "2Fl_LGam", "2Fl_LExpL", 
                "1Fl_L_TI", "1Fl_Gam_TI", "1Fl_ExpL_TI",
                "1Fl_LL_Blk", "1Fl_LGam_Blk", "1Fl_LExpL_Blk",
                "1Fl_L_RW_1.25", "1Fl_Gam_RW_2.0", "1Fl_ExpL_RW_2.0")


# AIC Summary -------------------------------------------------------------

AIC_df <- AIC_df %>% 
  mutate(Dat_Qual = case_when(
    str_detect(OM_Scenario, "High") ~ 'High',
    str_detect(OM_Scenario, "Low") ~ 'Low'
  ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"))

# To determine which models to plot out, etc.
# AIC plot
pdf(file = here(dir_out, "AIC_SummaryPlot.pdf"), width = 17.5, height = 7)

# Two fleet AIC df
twofleet_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "2Fl_")) %>% 
  group_by(sim, OM_Scenario, time_comp, Dat_Qual) %>% 
  mutate(min = min(AIC), min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "2 Fleet")

# Two fleet plot
(twofleet_plot <- ggplot(twofleet_aic, 
                         aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC,
                             label = round(n_minAIC, 2))) +
    geom_tile(alpha = 0.85) +
    geom_text(size = 5) +
    facet_grid(Dat_Qual~time_comp, scales = "free") +
    scale_fill_viridis_c() +
    labs(x = "OM Scenario", y = "EMs", fill = "Proportion of lowest AIC") +
    theme_test() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90))) 

# filter to 1 fleet random walk models
onefleet_rw_aic <- AIC_df %>% 
  mutate(
    selex_form = case_when(
      str_detect(EM_Scenario, "ExpL") ~ "Exponential Logistic",
      str_detect(EM_Scenario, "1Fl_L_RW") ~ "Logistic",
      str_detect(EM_Scenario, "1Fl_Gam_") ~ "Gamma",
    )) %>%  filter(str_detect(EM_Scenario, "1Fl_"),
                   str_detect(EM_Scenario, "RW_")) %>% 
  group_by(sim, OM_Scenario, time_comp, selex_form, Dat_Qual) %>% 
  mutate(min = min(AIC),
         min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, selex_form, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "1 Fleet")

# one fleet RW models 
# For random walk models, Logistic = 1.25, Gamma= 2.0, and ExpL = 2.0 = best AIC
(onefleet_rw_low_plot <- ggplot(onefleet_rw_aic %>% filter(Dat_Qual == "Low"), 
                            aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC,
                                label = round(n_minAIC, 2))) +
    geom_tile(alpha = 0.85) +
    geom_text(size = 6) +
    facet_grid(selex_form~time_comp, scales = "free") +
    scale_fill_viridis_c() +
    theme_test() +
    labs(x = "OM Scenario", y = "EMs", 
         fill = "Proportion of lowest AIC",
         title = "Low Data Quality") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90))) 

(onefleet_rw_high_plot <- ggplot(onefleet_rw_aic %>% filter(Dat_Qual == "High"), 
                                aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC,
                                    label = round(n_minAIC, 2))) +
    geom_tile(alpha = 0.85) +
    geom_text(size = 6) +
    facet_grid(selex_form~time_comp, scales = "free") +
    scale_fill_viridis_c() +
    theme_test() +
    labs(x = "OM Scenario", y = "EMs", 
         fill = "Proportion of lowest AIC",
         title = "High Data Quality") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90))) 

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
  mutate(fleet = "1 Fleet")

# One fleet plot
(onefleet_plot <- ggplot(onefleet_aic, 
                         aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC,
                             label = round(n_minAIC, 2))) +
    geom_tile(alpha = 0.85) +
    geom_text(size = 6) +
    facet_grid(Dat_Qual~time_comp, scales = "free_y") +
    scale_fill_viridis_c() +
    theme_test() +
    labs(x = "OM Scenario", y = "EMs", fill = "Proportion of lowest AIC") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90))) 

# Filter to 1 fleet time block models
onefleet_blk_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "Blk"),
         time_comp == "Terminal",
         str_detect(OM_Scenario, "Fast")) %>% 
  mutate(
    selex_form = case_when(
      str_detect(EM_Scenario, "LExpL") ~ "Logistic-Exponential Logistic",
      str_detect(EM_Scenario, "LGam") ~ "Logistic-Gamma",
      str_detect(EM_Scenario, "LL") ~ "Logistic-Logistic",
    )) %>% group_by(sim, OM_Scenario, time_comp, selex_form, Dat_Qual) %>% 
  mutate(min = min(AIC),
         min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, selex_form, Dat_Qual) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "1 Fleet")

# One fleet time block models (Seems like time-block 3 models are the best here?)
(onefleet_blk_plot <- ggplot(onefleet_blk_aic,
                             aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC,
                                 label = round(n_minAIC, 2))) +
    geom_tile(alpha = 0.85) +
    geom_text(size = 4) +
    facet_grid(selex_form~Dat_Qual, scales = "free") +
    scale_fill_viridis_c() +
    theme_test() +
    labs(x = "OM Scenario", y = "EMs", 
         fill = "Proportion of lowest AIC") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90))) 

dev.off()


# Parameter Summary (And Terminal SSB + Biomass) plot -------------------------------------------------------------

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
  ))

# Set order for plot
order <- vector()
for(o in 1:length(plot_order)) {
  order[o] <- which(grepl(plot_order[o], x = unique(pt_rg_re$EM_Scenario)))
} # end o loop

# Now relevel factor for organizing plot
pt_rg_re$EM_Scenario <- with(pt_rg_re,factor(EM_Scenario, levels = unique(EM_Scenario)[order]))

pdf(here(dir_out, "Param_Sum.pdf"), width = 20, height = 8)

for(i in 1:length(gen_om)) {
  # plot now!
  print(
    ggplot(pt_rg_re %>% filter(Dat_Qual == "Low", str_detect(OM_Scenario, gen_om[i])), 
           aes(x = factor(EM_Scenario), y = median, ymin = lwr_95, ymax = upr_95,
               color = time_comp, fill = time_comp, 
               label = round(median, 2))) +
      geom_pointrange(position = position_dodge2(width = 0.65), 
                      size = 1, linewidth = 1) +
      geom_vline(xintercept = c(seq(1.5, 13.5, 1)), lwd = 0.5) +
      geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
      facet_grid(type~OM_Scenario, scales = "free_x") +
      coord_cartesian(ylim = c(-0.4,0.4)) +
      scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      labs(x = "EM Scenarios", y = "Relative Error", fill = "Time Component",
           color = "Time Component", title = "Low Data Quality") +
      theme_matt() +
      theme(legend.position = "top", title = element_text(size = 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90)) 
  )
  
  print(
    ggplot(pt_rg_re %>% filter(Dat_Qual == "High", str_detect(OM_Scenario, gen_om[i])), 
           aes(x = factor(EM_Scenario), y = median, ymin = lwr_95, ymax = upr_95,
               color = time_comp, fill = time_comp, 
               label = round(median, 2))) +
      geom_pointrange(position = position_dodge2(width = 0.65), 
                      size = 1, linewidth = 1) +
      geom_vline(xintercept = c(seq(1.5, 13.5, 1)), lwd = 0.5) +
      geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
      facet_grid(type~OM_Scenario, scales = "free_x") +
      coord_cartesian(ylim = c(-0.4,0.4)) +
      scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      labs(x = "EM Scenarios", y = "Relative Error", fill = "Time Component",
           color = "Time Component", title = "High Data Quality") +
      theme_matt() +
      theme(legend.position = "top", title = element_text(size = 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90)) 
  )
}

dev.off()


# Parameter Summary (Time-block models) ------------------------------

# Filter to relevant components for parameters
om_scenario_params <- param_df %>% filter(time_comp == "Terminal",
                                          str_detect(OM_Scenario, "Fast"),
                                          str_detect(EM_Scenario, "Blk"), 
                                          type %in% c("F_0.4", "ABC"))

# Get terminal biomass as well
term_biom_ssb_blk <- ts_re_df %>% 
  filter(
    (year == c(50) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Terminal") |
      (year == c(70) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Terminal") |
      (year == c(30) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Fleet Trans End") |
      (year == c(50) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Fleet Trans End") |
      (year == c(27) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Fleet Intersect") |
      (year == c(40) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Fleet Intersect"),
    par_name %in% c("Spawning Stock Biomass"),
    str_detect(EM_Scenario, "Blk"),
    time_comp == "Terminal",
    str_detect(OM_Scenario, "Fast")) %>% 
  dplyr::select(OM_Scenario, EM_Scenario, time_comp, par_name,
                median, lwr_95, upr_95) %>% 
  rename(type = par_name)

# Point ranges for relative error and total error
pt_rg_re_blk <- om_scenario_params %>% 
 group_by(EM_Scenario, time_comp, type, OM_Scenario) %>% 
 summarize(median = median(RE), 
           lwr_95 = quantile(RE, 0.025),
           upr_95 =  quantile(RE, 0.975))

pt_rg_re_blk <- rbind(pt_rg_re_blk, term_biom_ssb_blk)

pdf(here(dir_out, "Param_Sum_BlkModels.pdf"), width = 35, height = 15)
  
  # plot low data quality
  print(
    ggplot(pt_rg_re_blk %>% 
             filter(str_detect(OM_Scenario, "Low")), aes(x = factor(EM_Scenario), y = median, ymin = lwr_95, ymax = upr_95)) +
      geom_pointrange(position = position_dodge2(width = 0.5), size = 1, linewidth = 1,
                      alpha = 0.5) +
      geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
      coord_cartesian(ylim = c(-0.5,0.5)) +
      facet_grid(type~OM_Scenario) +
      labs(x = "EM Scenarios", y = "Relative Error", fill = "Time Component",
           color = "Time Component", title = "Low Data Quality") +
      theme_matt() +
      theme(legend.position = "top", title = element_text(size = 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90),
            axis.ticks.x = element_blank()) 
  )
  
  # plot high data quality
  print(
    ggplot(pt_rg_re_blk %>% 
             filter(str_detect(OM_Scenario, "High")), aes(x = factor(EM_Scenario), y = median, ymin = lwr_95, ymax = upr_95)) +
      geom_pointrange(position = position_dodge2(width = 0.5), size = 1, linewidth = 1,
                      alpha = 0.5) +
      geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
      coord_cartesian(ylim = c(-0.5,0.5)) +
      facet_grid(type~OM_Scenario) +
      labs(x = "EM Scenarios", y = "Relative Error", fill = "Time Component",
           color = "Time Component", title = "High Data Quality") +
      theme_matt() +
      theme(legend.position = "top", title = element_text(size = 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90),
            axis.ticks.x = element_blank()) 
  )
  
dev.off()

# Time Series plot -----------------------------------------------------
  
pdf(here(dir_out, "TS_Sum.pdf"), width = 25, height = 6)
for(i in 1:length(gen_om)) {
  
    # Relative error of time series
    ts_re_om <- ts_re_df %>% filter(str_detect(OM_Scenario, gen_om[i]),
                                    str_detect(EM_Scenario, all_models)) %>% 
      mutate(Dat_Qual = case_when(
        str_detect(OM_Scenario, "High") ~ 'High',
        str_detect(OM_Scenario, "Low") ~ 'Low'
      ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"))
    
    # Set order for plot
      order <- vector()
      for(o in 1:length(plot_order)) {
        order[o] <- which(grepl(plot_order[o], x = unique(ts_re_om$EM_Scenario)))
      } # end o loop

    # Now relevel factor for organizing plot
    ts_re_om <- ts_re_om %>% 
      mutate(EM_Scenario = factor(EM_Scenario, levels = unique(EM_Scenario)[order]),
             time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End",
                                                      "Terminal")))
    
    for(j in 1:length(ts_pars)) {
      # Now loop through each time series component and print the plot out
      print(
        ggplot(ts_re_om %>% filter(Dat_Qual == "Low",
                                   par_name == ts_pars[j]), aes(x = year, y = median))  +
          geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = time_comp, group = time_comp), alpha = 0.3) +
          geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
          geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
          facet_grid(OM_Scenario~EM_Scenario) +
          coord_cartesian(ylim = c(-0.5,0.5)) +
          scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
          scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
          labs(x = "Year", y = paste(ts_pars[j], "RE"), 
               fill = "Time", color = "Time", title = "Low Data Quality") +
          theme_matt() +
          # theme(aspect.ratio=1)
          theme(legend.position = "top",
                title = element_text(size = 20),
                axis.text = element_text(size = 13), 
                strip.text = element_text(size = 13)) 
      )
      
      print(
        ggplot(ts_re_om %>% filter(Dat_Qual == "High",
                                   par_name == ts_pars[j]), aes(x = year, y = median))  +
          geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = time_comp, group = time_comp), alpha = 0.3) +
          geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
          geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
          facet_grid(OM_Scenario~EM_Scenario) +
          coord_cartesian(ylim = c(-0.5,0.5)) +
          scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
          scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
          labs(x = "Year", y = paste(ts_pars[j], "RE"), 
               fill = "Time", color = "Time", title = "High Data Quality") +
          theme_matt() +
          # theme(aspect.ratio=1)
          theme(legend.position = "top", 
                title = element_text(size = 20),
                axis.text = element_text(size = 13),
                strip.text = element_text(size = 13)) 
      )
    } # end j

} # end i

dev.off()


# Time Series (Time Block Models) -----------------------------------------

pdf(here(dir_out, "TS_SumBlkModels.pdf"), width = 45, height = 10)

for(i in 1:length(ts_pars)) {
  
  # Relative error of time series
  ts_re_om <- ts_re_df %>% filter(time_comp == "Terminal", 
                                  str_detect(OM_Scenario, "Fast"),
                                  str_detect(EM_Scenario, "Blk")) %>% 
    mutate(Dat_Qual = case_when(
      str_detect(OM_Scenario, "High") ~ 'High',
      str_detect(OM_Scenario, "Low") ~ 'Low'
    ),  OM_Scenario = str_remove(OM_Scenario, "_High|_Low"))
  
    # Now loop through each time series component and print the plot out
    print(
      ggplot(ts_re_om %>% filter(Dat_Qual == "Low",
                                 par_name == ts_pars[i]), aes(x = year, y = median))  +
        geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.3) +
        geom_line(linewidth = 2, alpha = 1) +
        geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
        facet_grid(OM_Scenario~EM_Scenario) +
        coord_cartesian(ylim = c(-0.4,0.4)) +
        scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
        scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
        labs(x = "Year", y = paste(ts_pars[i], "RE"), fill = "Time", color = "Time",
             title = "Low Data Quality") +
        theme_matt() +
        # theme(aspect.ratio=1)
        theme(legend.position = "top", title = element_text(size = 20),
              axis.text = element_text(size = 13), strip.text = element_text(size = 11)) 
    )
    
    print(
      ggplot(ts_re_om %>% filter(Dat_Qual == "High", par_name == ts_pars[i]), 
             aes(x = year, y = median))  +
        geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.3) +
        geom_line(linewidth = 2, alpha = 1) +
        geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
        facet_grid(OM_Scenario~EM_Scenario) +
        coord_cartesian(ylim = c(-0.4,0.4)) +
        scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
        scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
        labs(x = "Year", y = paste(ts_pars[i], "RE"), 
             fill = "Time", color = "Time", title = "High Data Quality") +
        theme_matt() +
        # theme(aspect.ratio=1)
        theme(legend.position = "top", title = element_text(size = 20),
              axis.text = element_text(size = 13), strip.text = element_text(size = 11)) 
    )
} # end i

dev.off()


# Numbers At Age Plots ----------------------------------------------------

# Do some residual munging
NAA_df_re <- NAA_df %>% 
  mutate(RE = (Est_Numbers - True_Numbers) / True_Numbers) %>% 
  group_by(OM_Scenario, EM_Scenario, Year, Age, Sex, time_comp) %>% 
  summarize(Median_RE = median(RE, na.rm = TRUE),
            Lwr_95 = quantile(RE, 0.025, na.rm = TRUE),
            Upr_95 = quantile(RE, 0.975, na.rm = TRUE))

pdf(here(dir_out, "NAA_Summary.pdf"), width = 35, height = 25)
for(i in 1:length(unique_oms)) {
  
  # Relative error of time series
  NAA_df_re_plot <- NAA_df_re %>% filter(OM_Scenario == unique_oms[i],
                                  str_detect(EM_Scenario, all_models))
  
  # Set order for plot
  order <- vector()
  for(o in 1:length(plot_order)) {
    order[o] <- which(grepl(plot_order[o], x = unique(NAA_df_re_plot$EM_Scenario)))
  } # end o loop
  
  # Now relevel factor for organizing plot
  NAA_df_re_plot <- NAA_df_re_plot %>% 
    mutate(time_comp = factor(time_comp, 
                              levels = c("Fleet Intersect", "Fleet Trans End", "Terminal")))
  # relvel factors
  NAA_df_re_plot$EM_Scenario <- factor(NAA_df_re_plot$EM_Scenario, 
                                  levels = c(unique(NAA_df_re_plot$EM_Scenario)[order]))
  
  print(
    ggplot(NAA_df_re_plot %>% 
             filter(Sex == "Female",
                    Age %in% c(seq(3, 30, 3))), 
           aes(x = Year, y = Median_RE, fill = time_comp)) +
      geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
      geom_ribbon(aes(ymin = Lwr_95, ymax = Upr_95), alpha = 0.3) +
      geom_hline(aes(yintercept = 0), col = "black", 
                 lty = 2, linewidth = 0.5, alpha = 1) +
      coord_cartesian(ylim = c(-0.4,0.4)) +
      facet_grid(EM_Scenario~Age) +
      labs(x = "Year", y = "RE", color = "Time Component",
           fill = "Time Component", title = paste("Female", unique_oms[i])) +
      scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      theme_bw() +
      theme(legend.position = "top", 
            title = element_text(size = 20),
            axis.text = element_text(size = 13), 
            strip.text = element_text(size = 13)) 
  )
  
  print(
    ggplot(NAA_df_re_plot %>% 
             filter(Sex == "Male",
                    Age %in% c(seq(3, 30, 3))), 
           aes(x = Year, y = Median_RE, fill = time_comp)) +
      geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
      geom_ribbon(aes(ymin = Lwr_95, ymax = Upr_95), alpha = 0.3) +
      geom_hline(aes(yintercept = 0), col = "black", 
                 lty = 2, linewidth = 0.5, alpha = 1) +
      coord_cartesian(ylim = c(-0.4,0.4)) +
      facet_grid(EM_Scenario~Age) +
      labs(x = "Year", y = "RE", color = "Time Component",
           fill = "Time Component", title = paste("Male", unique_oms[i])) +
      scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      theme_bw() +
      theme(legend.position = "top", 
            title = element_text(size = 20),
            axis.text = element_text(size = 13), 
            strip.text = element_text(size = 13)) 
  )
} # end i loop
dev.off()

# Convergence Summary plot ------------------------------------------------
# Convergence for all models we want to look at
conv_stat <- AIC_df %>% 
  filter(str_detect(EM_Scenario, all_models)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(converged = sum(conv == "Converged")/200)

# Random walk convergence statistics
RW_convergence <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "RW_")) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp, Dat_Qual) %>% 
  summarize(converged = sum(conv == "Converged")/200)

# Set order for plot
order <- vector()
for(o in 1:length(plot_order)) {
  order[o] <- which(grepl(plot_order[o], x = unique(conv_stat$EM_Scenario)))
} # end o loop

# relvel factors
conv_stat$EM_Scenario <- factor(conv_stat$EM_Scenario, 
                                levels = c(unique(conv_stat$EM_Scenario)[order]))

pdf(here(dir_out, "Convg_Sum.pdf"), width = 25, height = 10)

  print(
    ggplot(conv_stat, mapping = aes(x = EM_Scenario, y = converged * 100,
                         group = time_comp, fill = time_comp))  +
      geom_col(position = "dodge", alpha = 0.85) +
      geom_hline(aes(yintercept = 50), col = "black", lty = 2, size = 0.5, alpha = 1) +
      facet_grid(Dat_Qual~OM_Scenario,  scales = "free_x") +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      theme_matt() +
      labs(fill = "Time Component", y = "Convergence Rate",
           x = "EM Scenario") +
      # theme(aspect.ratio=1) +
      theme(legend.position = "top", title = element_text(size = 20),
            strip.text = element_text(size = 15)) +
      coord_flip()
  )
  
  # Random walk convergence statistics
  print(
    ggplot(RW_convergence, mapping = aes(x = EM_Scenario, y = converged * 100,
                                    group = time_comp, fill = time_comp))  +
      geom_col(position = "dodge", alpha = 0.85) +
      geom_hline(aes(yintercept = 50), col = "black", lty = 2, size = 0.5, alpha = 1) +
      facet_grid(Dat_Qual~OM_Scenario,  scales = "free_x") +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      theme_matt() +
      labs(fill = "Time Component", y = "Convergence Rate",
           x = "EM Scenario") +
      # theme(aspect.ratio=1) +
      theme(legend.position = "top", title = element_text(size = 20),
            strip.text = element_text(size = 15)) +
      coord_flip()
      
  )
  
dev.off()


# Min Max Solution Plot ---------------------------------------------------

minmax_df <- ts_are_df %>% 
  # Selective filtering down here to terminal years
  filter(str_detect(EM_Scenario, all_models),
           # Get terminal year of peels/EMs
           (year == c(50) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Terminal") |
           (year == c(70) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Terminal") |
           (year == c(30) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Fleet Trans End") |
           (year == c(50) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Fleet Trans End") |
           (year == c(27) & str_detect(OM_Scenario, "Fast_") ) & str_detect(time_comp, "Fleet Intersect") |
           (year == c(40) & str_detect(OM_Scenario, "Slow_") ) & str_detect(time_comp, "Fleet Intersect")) %>% 
  data.frame() %>%
  group_by(EM_Scenario, time_comp, par_name) %>%
  mutate(max_median = max(median)) %>% # find the maximum median MARE
  ungroup()

# Find the minimum maximum median value
min_max_medians <- minmax_df %>%
  group_by(par_name, time_comp) %>%
  summarize(min_max_medians = min(max_median))

# Now left join this
minmax_df = minmax_df %>% 
  left_join(min_max_medians, by = c("par_name", "time_comp")) 

# Plot!
pdf(file = here(dir_out, "MinMax_Summary.pdf"), width = 35, height = 10)
for(i in 1:length(ts_pars)) {
  print(
    ggplot(minmax_df %>% filter(par_name == ts_pars[i]), 
           aes(x = factor(OM_Scenario), y = factor(EM_Scenario), 
               fill = median, label = round(median, 3))) +
      geom_tile(alpha = 0.85) +
      facet_wrap(~time_comp, scales = "free_x") +
      geom_text(color = ifelse(minmax_df$median[minmax_df$par_name == ts_pars[i]] == 
                               minmax_df$max_median[minmax_df$par_name == ts_pars[i]] &
                               minmax_df$median[minmax_df$par_name == ts_pars[i]]
                               != minmax_df$min_max_medians[minmax_df$par_name == ts_pars[i]], "red", 
                               ifelse(minmax_df$median[minmax_df$par_name == ts_pars[i]] == 
                               minmax_df$min_max_medians[minmax_df$par_name == ts_pars[i]], "green", "black")),
                size = 5.5) +
      scale_fill_viridis_c() +
      labs(x = "OM Scenario", y = "EM Scenario", fill = "Median MARE",
           title = ts_pars[i]) +
      theme_test() +
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 90)) 
  )
} # end i
dev.off()

# Selectivity plots (OM Selex vs. EM Comparison) -------------------------------------------------------
# Subset to first year for OMs
om_slx_df_sub <- om_slx_df %>% filter(Year == 1)

pdf(here(dir_out, "Selex_Sum.pdf"), width = 22, height = 5)
for(i in 1:length(unique_oms)) {
  
  # Read in selectivity dataframe for a given OM
  em_slx_df = data.table::fread(here('output', "OM_Scenarios", unique_oms[i], "EM_Fish_Selex.csv"))
  
  # Filter only to terminal years of terminal, 10, 5; also do some residual munging on the dataframe
  trunc_em <- em_slx_df %>% 
    filter((Year %in% c(50) & str_detect(OM_Scenario, "Fast_") ) & time_comp == "Terminal" |
           (Year %in% c(70) & str_detect(OM_Scenario, "Slow_") ) & time_comp == "Terminal"  |
           (Year %in% c(30) & str_detect(OM_Scenario, "Fast_") ) & time_comp == "Fleet Trans End"  |
           (Year %in% c(50) & str_detect(OM_Scenario, "Slow_") ) & time_comp == "Fleet Trans End" |
           (Year %in% c(27) & str_detect(OM_Scenario, "Fast_") ) & time_comp == "Fleet Intersect" |
           (Year %in% c(40) & str_detect(OM_Scenario, "Slow_") ) & time_comp == "Fleet Intersect")
  
  # Subset plot here for plotting purposes
  em_plot_df <- trunc_em %>% 
    filter(OM_Scenario == unique_oms[i],
           str_detect(EM_Scenario, all_models)) %>% 
    group_by(EM_Scenario, Age, Fleet, Sex, time_comp) %>% 
    summarize(mean = mean(Selex, na.rm = T)) 
  
  # Subset om dataframe plot to a specific om
  om_plot_df <- om_slx_df_sub %>% 
    filter(OM_Scenario == unique_oms[i])
    
  # Set order for plot
  order <- vector()
  for(o in 1:length(plot_order)) {
    order[o] <- which(grepl(plot_order[o], x = unique(em_plot_df$EM_Scenario)))
  } # end o loop

  # Now relevel factor for organizing plot
  em_plot_df$EM_Scenario <- factor(em_plot_df$EM_Scenario, levels = unique(em_plot_df$EM_Scenario)[order])
  
  print(
    ggplot() +
      geom_line(om_plot_df, mapping = aes(x = Age, y = Selex, 
                                          linetype = as.factor(Fleet)), size = 1, alpha = 0.5)  +
      geom_point(om_plot_df, mapping = aes(x = Age, y = Selex, 
                                           shape = as.factor(Fleet)),  size = 2, alpha = 0.75)  +
      geom_line(em_plot_df, mapping = aes(x = Age, y = mean, color = as.factor(time_comp), 
                                          lty = as.factor(Fleet)), size = 1.3, alpha = 1) +
      scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      facet_grid(Sex~EM_Scenario) +
      theme_matt() +
      labs(x = "Age", y = "Mean Proportion Selected", color = "EM Time Component", 
           shape = "Fleet", linetype = "Fleet", title = unique_oms[i]) +
      theme(title = element_text(size = 20),
            legend.box = "vertical", strip.text = element_text(size = 12)) 
  )
  
} # end i

dev.off()

# Get Population Selection Curves -----------------------------------------
### Models of Interest ------------------------------------------------------

# Now plot this out!
pdf(file = here(dir_out, "PopSelex.pdf"), width = 25, height = 15)
for(i in 1:length(unique_oms)) {
  
    # Filter out time-block stuff 
    plot_df <- pop_sel_em %>% 
      filter(str_detect(EM_Scenario, all_models), 
             OM_Scenario == unique_oms[i])
    
    # Set order for plot
    order <- vector()
    for(o in 1:length(plot_order)) {
      order[o] <- which(grepl(plot_order[o], x = unique(plot_df$EM_Scenario)))
    } # end o loop
    
    # Now relevel factor for organizing plot
    plot_df$EM_Scenario <- factor(plot_df$EM_Scenario, levels = unique(plot_df$EM_Scenario)[order])
    
    # Females
    female_plot <- ggplot() +
    # EM Median Population Selex
    geom_line(plot_df %>% filter(Sex == "Female"), 
              mapping = aes(x = Age, y = Median_Selex),  lwd = 1.5,
              color = "red", alpha = 0.5) +
    geom_ribbon(plot_df %>% filter(Sex == "Female"),
                mapping = aes(x = Age, y = Median_Selex, 
                              ymin = Lwr_95, ymax = Upr_95), alpha = 0.3,
                fill = "red") +
      # OM Median Population Selex
      geom_line(pop_sel_om %>% filter(Sex == "Female", OM_Scenario == unique_oms[i]), 
                mapping = aes(x = Age, y = Selex),
                color = "black", lty = 2, lwd = 1.5) +
    facet_grid(time_comp~EM_Scenario) +
    labs(x = "Age", y = "Population Selectivity",
         color = "Time Component", fill = "Time Component",
         title = paste(unique_oms[i], "Female")) +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20))
    
    # Males
    male_plot <- ggplot() +
      # EM Median Population Selex
      geom_line(plot_df %>% filter(Sex == "Male"), 
                mapping = aes(x = Age, y = Median_Selex), lwd = 1.5,color = 'blue',
                alpha = 0.5) +
      geom_ribbon(plot_df %>% filter(Sex == "Male"),
                  mapping = aes(x = Age, y = Median_Selex, 
                                ymin = Lwr_95, ymax = Upr_95), alpha = 0.3, fill = "blue") +
      # OM Median Population Selex
      geom_line(pop_sel_om %>% filter(Sex == "Male",OM_Scenario == unique_oms[i]), 
                mapping = aes(x = Age, y = Selex),
                color = "black", lty = 2, lwd = 1.5) +
      facet_grid(time_comp~EM_Scenario) +
      labs(x = "Age", y = "Population Selectivity",
           color = "Time Component", fill = "Time Component",
           title = paste(unique_oms[i], "Male")) +
      theme_matt() +
      theme(legend.position = "top", 
            title = element_text(size = 20))
    
      print(cowplot::plot_grid(female_plot, male_plot, ncol = 1))
      
      # Plots of relative error
      re_slx_plot <- plot_df %>% 
        left_join(pop_sel_om %>% 
                  rename(True = Selex), by = c("OM_Scenario", "Age", "Sex", "time_comp")) %>% 
        mutate(re = (Median_Selex - True) / True,
               lwr_95_re = (Lwr_95 - True) / True,
               upr_95_re = (Upr_95 - True) / True)
      
      print(
        ggplot(re_slx_plot, aes(x = Age, y = re, ymin = lwr_95_re, ymax = upr_95_re,
                                fill = Sex)) +
          geom_hline(yintercept = 0, lty = 2, lwd = 1.5) +
          geom_ribbon(alpha = 0.3) +
          geom_line(lwd = 1.5, alpha = 0.5, aes(color = Sex)) +
          coord_cartesian(ylim = c(-1.5, 1.5)) +
          scale_color_manual(values = c("red", "blue")) +
          scale_fill_manual(values = c("red", "blue")) +
          facet_grid(time_comp~EM_Scenario) +
          labs(x = "Age", y = "Relative Error",
               color = "Sex", fill = "Sex",
               title = paste(unique_oms[i])) +
          theme_matt() +
          theme(legend.position = "top", 
                title = element_text(size = 20))
      )

} # end i loop
dev.off()


### Block Models ------------------------------------------------------------

# Now plot this out!
pdf(file = here(dir_out, "PopSelex_BlkModels.pdf"), width = 20, height = 25)
fast_oms <- unique_oms[str_detect(unique_oms, "Fast")]
for(i in 1:length(fast_oms)) {
  
  # Filter out time-block stuff 
  plot_df <- pop_sel_em %>% 
    filter(str_detect(EM_Scenario, "Blk"), OM_Scenario == fast_oms[i])
  
  # Females
  female_plot <- ggplot() +
    # EM Median Population Selex
    geom_line(plot_df %>% filter(Sex == "Female",
                                 time_comp == "Terminal"), 
              mapping = aes(x = Age, y = Median_Selex), lwd = 1.5,
              color = "red", alpha = 0.5) +
    geom_ribbon(plot_df %>% filter(Sex == "Female",
                                   time_comp == "Terminal"),
                mapping = aes(x = Age, y = Median_Selex, 
                              ymin = Lwr_95, ymax = Upr_95), alpha = 0.3,
                fill = "red") +
    # OM Median Population Selex
    geom_line(pop_sel_om %>% filter(Sex == "Female", 
                                    OM_Scenario == unique_oms[i],
                                    time_comp == "Terminal"), 
              mapping = aes(x = Age, y = Selex),
              color = "black", lty = 2, lwd = 1.5) +
    facet_wrap(~EM_Scenario) +
    labs(x = "Age", y = "Population Selectivity",
         color = "Time Component", fill = "Time Component",
         title = paste(unique_oms[i], "Female")) +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20))
  
  # Males
  male_plot <- ggplot() +
    # EM Median Population Selex
    geom_line(plot_df %>% filter(Sex == "Male",
                                 time_comp == "Terminal"), 
              mapping = aes(x = Age, y = Median_Selex), lwd = 1.5,color = 'blue',
              alpha = 0.5) +
    geom_ribbon(plot_df %>% filter(Sex == "Male",
                                   time_comp == "Terminal"),
                mapping = aes(x = Age, y = Median_Selex, 
                              ymin = Lwr_95, ymax = Upr_95), alpha = 0.3, fill = "blue") +
    # OM Median Population Selex
    geom_line(pop_sel_om %>% filter(Sex == "Male",
                                    OM_Scenario == unique_oms[i],
                                    time_comp == "Terminal"), 
              mapping = aes(x = Age, y = Selex),
              color = "black", lty = 2, lwd = 1.5) +
    facet_wrap(~EM_Scenario) +
    labs(x = "Age", y = "Population Selectivity",
         color = "Time Component", fill = "Time Component",
         title = paste(unique_oms[i], "Male")) +
    theme_matt() +
    theme(legend.position = "top", 
          title = element_text(size = 20))
  
  print(cowplot::plot_grid(female_plot, male_plot, ncol = 1))
  
  # Plots of relative error
  re_slx_plot <- plot_df %>% 
    left_join(pop_sel_om %>% 
                rename(True = Selex), by = c("OM_Scenario", "Age", "Sex", "time_comp")) %>% 
    mutate(re = (Median_Selex - True) / True,
           lwr_95_re = (Lwr_95 - True) / True,
           upr_95_re = (Upr_95 - True) / True)
  
  print(
    ggplot(re_slx_plot %>% 
             filter(time_comp == "Terminal"), aes(x = Age, y = re, ymin = lwr_95_re, ymax = upr_95_re,
                                                  fill = Sex)) +
      geom_hline(yintercept = 0, lty = 2, lwd = 1.5) +
      geom_ribbon(alpha = 0.3) +
      geom_line(lwd = 1.5, alpha = 0.5, aes(color = Sex)) +
      coord_cartesian(ylim = c(-1.5, 1.5)) +
      scale_color_manual(values = c("red", "blue")) +
      scale_fill_manual(values = c("red", "blue")) +
      facet_wrap(~EM_Scenario) +
      labs(x = "Age", y = "Relative Error",
           color = "Sex", fill = "Sex",
           title = paste(unique_oms[i])) +
      theme_matt() +
      theme(legend.position = "top", 
            title = element_text(size = 20))
  )
  
} # end i loop
dev.off()

# Get SPR Differences -----------------------------------------------------

pdf(file = here(dir_out, "SPR_PopSelex.pdf"), width = 27, height = 13)
# Now plot this out!
for(i in 1:length(unique_oms)) {
  
  # Filter out time-block stuff 
  em_plot_df <- em_spr_df %>% 
    filter(str_detect(EM_Scenario, all_models), 
           OM_Scenario == unique_oms[i])
  
  # Get OM SPR
  om_plot_df <- om_spr_df %>% filter(OM_Scenario == unique_oms[i])
  
  # Get F40 lines for EM
  em_f40 <- em_plot_df %>% 
    group_by(EM_Scenario, time_comp) %>% 
    summarize(f40_em = trial_F[which.min(Diff)])
  
  # Get F40 lines for oM
  om_f40 <- om_plot_df %>% 
    group_by(time_comp) %>% 
    summarize(f40_om = trial_F[which.min(Diff)])
  
  # Calculate relative error
  em_f40 <- em_f40 %>% 
    left_join(om_f40, by = c("time_comp")) %>% 
    mutate(re = (f40_em - f40_om) / f40_om)
  
  # Set order for plot
  em_f40_order <- vector()
  em_scenario_order <- vector()
  for(o in 1:length(plot_order)) {
    em_f40_order[o] <- which(grepl(plot_order[o], x = unique(em_f40$EM_Scenario)))
    em_scenario_order[o] <- which(grepl(plot_order[o], x = unique(em_plot_df$EM_Scenario)))
  } # end o loop
  
  # Now relevel factor for organizing plot
  em_f40$EM_Scenario <- factor(em_f40$EM_Scenario, levels = unique(em_f40$EM_Scenario)[em_f40_order])
  em_plot_df$EM_Scenario <- factor(em_plot_df$EM_Scenario, levels = unique(em_plot_df$EM_Scenario)[em_scenario_order])
  
  print(
    ggplot() +
      # SPR lines - EM
      geom_line(em_plot_df, mapping = aes(x = trial_F, y = SPR, color = "Estimated"),
                lwd = 1, alpha = 0.75) +
      # SPR lines - OM
      geom_line(om_plot_df, mapping = aes(x = trial_F, y = SPR, color = "True"), 
                lwd = 1, alpha = 0.75, lty = 2) +
      # F40 lines - EM
      geom_vline(em_f40, mapping = aes(xintercept = f40_em, color = "Estimated"),
                 lwd = 1, lty = 1, alpha = 0.75) +
      # F40 lines - OM
      geom_vline(om_f40, mapping = aes(xintercept = f40_om, color = "True"),
                 lwd = 1, lty = 2, alpha = 0.75) +
      # Geom text for relative error
      geom_text(em_f40, mapping = aes(x = 0.4, y = 0.85, 
                                      label = paste("RE:", round(re, 3) * 100)))+
      
      scale_color_manual(values = c("blue", "black"), label = c("Estimated", "True")) +
      facet_grid(time_comp~EM_Scenario) +
      labs(x = "Exploitation Rate", y = "SPR", color = "", title = unique_oms[i]) +
      theme_matt() +
      theme(legend.position = "top", 
            title = element_text(size = 20))
  )
  
} # end i 
dev.off()


# Get SPR differences (Maturity Sensitivity) ------------------------------

pdf(file = here(dir_out, "SPR_PopSelex_MatSense.pdf"), width = 27, height = 13)
# Now plot this out!
for(i in 1:length(unique_oms)) {
  
  # Filter out time-block stuff 
  em_plot_df <- em_spr_mat_df %>% 
    filter(str_detect(EM_Scenario, all_models), 
           OM_Scenario == unique_oms[i])
  
  # Get OM SPR
  om_plot_df <- om_spr_mat_df %>% filter(OM_Scenario == unique_oms[i])
  
  # Get F40 lines for EM
  em_f40 <- em_plot_df %>% 
    group_by(EM_Scenario, time_comp, Mat_Opt) %>% 
    summarize(f40_em = trial_F[which.min(Diff)])
  
  # Get F40 lines for oM
  om_f40 <- om_plot_df %>% 
    group_by(time_comp, Mat_Opt) %>% 
    summarize(f40_om = trial_F[which.min(Diff)])
  
  # Calculate relative error
  em_f40 <- em_f40 %>% 
    left_join(om_f40, by = c("time_comp", "Mat_Opt")) %>% 
    mutate(re = (f40_em - f40_om) / f40_om)
  
  # Set order for plot
  em_f40_order <- vector()
  em_scenario_order <- vector()
  for(o in 1:length(plot_order)) {
    em_f40_order[o] <- which(grepl(plot_order[o], x = unique(em_f40$EM_Scenario)))
    em_scenario_order[o] <- which(grepl(plot_order[o], x = unique(em_plot_df$EM_Scenario)))
  } # end o loop
  
  # Now relevel factor for organizing plot
  em_f40$EM_Scenario <- factor(em_f40$EM_Scenario, levels = unique(em_f40$EM_Scenario)[em_f40_order])
  em_plot_df$EM_Scenario <- factor(em_plot_df$EM_Scenario, levels = unique(em_plot_df$EM_Scenario)[em_scenario_order])
  
  print(
    ggplot() +
      # SPR lines - EM
      geom_line(em_plot_df %>% filter(time_comp == "Terminal"),
                mapping = aes(x = trial_F, y = SPR, color = "Estimated"),
                lwd = 1, alpha = 0.75) +
      # SPR lines - OM
      geom_line(om_plot_df %>% filter(time_comp == "Terminal"), 
                mapping = aes(x = trial_F, y = SPR, color = "True"), 
                lwd = 1, alpha = 0.75, lty = 2) +
      # F40 lines - EM
      geom_vline(em_f40 %>% filter(time_comp == "Terminal"), 
                 mapping = aes(xintercept = f40_em, color = "Estimated"),
                 lwd = 1, lty = 1, alpha = 0.75) +
      # F40 lines - OM
      geom_vline(om_f40 %>% filter(time_comp == "Terminal"),
                 mapping = aes(xintercept = f40_om, color = "True"),
                 lwd = 1, lty = 2, alpha = 0.75) +
      # Geom text for relative error
      geom_text(em_f40 %>% filter(time_comp == "Terminal"),
                mapping = aes(x = 0.4, y = 0.85, label = paste("RE:", round(re, 3) * 100)))+
      scale_color_manual(values = c("blue", "black"), label = c("Estimated", "True")) +
      facet_grid(Mat_Opt~EM_Scenario) +
      labs(x = "Exploitation Rate", y = "SPR", color = "", title = unique_oms[i]) +
      theme_matt() +
      theme(legend.position = "top", 
            title = element_text(size = 20))
  )
  
} # end i 
dev.off()

# 
# # Composition plots --------------------------------------------------------
# 
# # Get manual colors
# colors <- viridis::viridis(n = 10)[c(1, 4, 10)]
# 
# # Filter comps to terminal year
# om_prop_term <- om_comps_df %>% 
#   group_by(OM_Scenario, Ages, Sexes, fleet_type) %>% 
#   summarize(Prop = mean(Prop)) %>% 
#   mutate(time_comp = "Terminal",
#          Sex = case_when(
#            str_detect(Sexes, "1") ~ "Female",
#            str_detect(Sexes, "2") ~ "Male"),
#          Fleet = case_when(
#            str_detect(fleet_type, "One Fleet") ~ 1,
#            str_detect(fleet_type, "Two Fleets") ~ 2
#          )) %>% 
#   filter(fleet_type == "One Fleet") %>% 
#   dplyr::select(-Sexes)
# 
# # EM proportion munging
# em_prop <- em_comps_df %>% 
#   filter(!str_detect(EM_Scenario, "0.5|1.0|2.0"),
#          !str_detect(EM_Scenario, "5_|10_"),
#          conv == "Converged") %>%  # remove 2 fleet model
#   mutate(EM_Scenario = str_remove(EM_Scenario, "Term_|10_|5_"),
#          Sex = case_when(
#            str_detect(Sex, "1") ~ "Female",
#            str_detect(Sex, "2") ~ "Male")) %>% 
#   group_by(Age, Fleet, Sex, OM_Scenario, EM_Scenario) %>% 
#   summarize(Pred = mean(Prop)) 
# 
# 
# # Set order for plot
# order <- vector()
# no2fl_plot_order <- plot_order[-1] # remove 2 fleets
# for(o in 1:length(no2fl_plot_order)) {
#   order[o] <- which(grepl(no2fl_plot_order[o], x = unique(em_prop$EM_Scenario)))
# } # end o loop
# 
# # Now relevel factor for organizing plot
# em_prop$EM_Scenario <- factor(em_prop$EM_Scenario, levels = unique(em_prop$EM_Scenario)[order])
#   
# # Left join EM and OM
# em_om_df <- em_prop %>% 
#   left_join(om_prop_term, 
#             by = c("Sex", "Age" = "Ages", "Fleet", "OM_Scenario")) %>% 
#   drop_na() %>% 
#   mutate(Diff = Prop - Pred)
# 
# # now plot!
# pdf(here(dir_out, "Comp_Fits_Sum.pdf"), width = 22, height = 8)
# for(i in 1:length(unique_oms)) {
#   
#   # Filter to a given OM
#   plot_df <- em_om_df %>% 
#     filter(OM_Scenario == unique_oms[i])
#   
#   # Comps and predictions
#   print(
#     ggplot() +
#       geom_col(plot_df, mapping = aes(x = Age, y = Prop, fill = Diff), alpha = 1,
#                position = "identity", color = "grey", size = 0.5) +
#       geom_line(plot_df, mapping = aes(x = Age, y = Pred), size = 0.85,
#                 color = "black", lty = 2) +
#       facet_grid(Sex~EM_Scenario) +
#       scale_fill_gradient2(midpoint = 0, oob = scales::squish) + 
#       guides (fill = guide_colourbar(barwidth = 16, barheight = 1,
#                                      title.position = "top")) +
#       theme_matt() +
#       labs(x = "Age", y = "Mean Proportion", title = unique_oms[i],
#            fill = "Difference (Observed - Predicted)") +
#       theme(title = element_text(size = 20), 
#             aspect.ratio = 0.8)
#   )
# 
# }
# dev.off()
