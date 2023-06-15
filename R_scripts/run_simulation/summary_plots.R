# Purpose: To make summary plots of simulation runs
# Notes: For TVCont - picking sd of 1.5 for all runs since it seemed to work well
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
param_df <- data.table::fread(here("output", "Parameter_Summary.csv")) # parameters
ts_df <- data.table::fread(here("output", "TimeSeries_Summary.csv")) # time series values 
ts_re_df <- data.table::fread(here("output", "TimeSeries_RE.csv")) # time series relative error (converged runs only)
ts_te_df <- data.table::fread(here("output", "TimeSeries_TE.csv")) # time series total error (converged runs only)
ts_are_df <- data.table::fread(here("output", "TimeSeries_ARE.csv")) # time series abs relative error (converged runs only)
AIC_df <- data.table::fread(here("output", "AIC_Convergence_Summary.csv")) # time series total error (converged runs only)

# Selectivity and Comps
om_slx_df <- data.table::fread(here("output", "OM_Fish_Selex.csv")) # OM Selectivity Values
em_slx_df <- data.table::fread(here("output", "EM_Fish_Selex.csv")) # EM Selectivity Values
om_comps_df <- data.table::fread(here("output", "OM_Fish_Comps.csv")) # OM Composition values
em_comps_df <- data.table::fread(here("output", "EM_Fish_Comps.csv")) # EM Composition values

# Unique oms and other components
unique_oms <- unique(param_df$OM_Scenario) # unique oms
ts_pars <- unique(ts_re_df$par_name) # unique parameter names
 
# Get plot order for EM models
plot_order <- c("2Fl_LL", "2Fl_LGam", "2Fl_LExpL", "1Fl_L_TI", "1Fl_Gam_TI", "1Fl_ExpL_TI",
                "1Fl_LL_Blk", "1Fl_LGam_Blk", "1Fl_LExpL_Blk",
                "1Fl_L_RW", "1Fl_Gam_RW", "1Fl_ExpL_RW")

# Parameter Summary plot -------------------------------------------------------------
param_df <-  param_df %>% 
  mutate(
  EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'),
  time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End",
                                           "Terminal"))) %>% 
  filter(conv == "Converged") # only converged runs

pdf(here(dir_out, "Param_Sum.pdf"), width = 25, height = 10)

for(i in 1:length(unique_oms)) {
  
  # Filter to relevant components for parameters
  om_scenario_params <- param_df %>% filter(OM_Scenario == unique_oms[i],
                                            !str_detect(EM_Scenario, "0.5|1.0|2.0"), 
                                            type %in% c("F_0.4", "ABC"))
  
  # Set order for plot
  if(length(plot_order) == length(unique(om_scenario_params$EM_Scenario))) {
    order <- vector()
    for(o in 1:length(plot_order)) {
      order[o] <- which(grepl(plot_order[o], x = unique(om_scenario_params$EM_Scenario)))
    } # end o loop
  } # end if
  
  # Now relevel factor for organizing plot
  om_scenario_params$EM_Scenario <- with(om_scenario_params,
                                         factor(EM_Scenario, 
                                         levels = unique(EM_Scenario)[order]))
  
  # Point ranges for relative error and total error
  pt_rg_re <- om_scenario_params %>% 
    group_by(EM_Scenario, time_comp, type) %>% 
    summarize(median = median(RE), 
              lwr_95 = quantile(RE, 0.025),
              upr_95 =  quantile(RE, 0.975))
  
  # plot now!
  print(
    ggplot(pt_rg_re, aes(x = factor(EM_Scenario), y = median, ymin = lwr_95, ymax = upr_95,
                                   color = time_comp, fill = time_comp)) +
      geom_pointrange(position = position_dodge2(width = 1), size = 1.5, linewidth = 1) +
      geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
      facet_grid(type~EM_Scenario, scales = "free_x") +
      coord_cartesian(ylim = c(-0.4,0.4)) +
      scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      labs(x = "EM Scenarios", y = "Relative Error", fill = "Time Component",
           color = "Time Component", title = unique_oms[i]) +
      theme_matt() +
      theme(legend.position = "top", title = element_text(size = 20),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) 
  )
  
} #end i

dev.off()



# Time Series plot -----------------------------------------------------

# Quick munging here
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
  
pdf(here(dir_out, "TS_Sum.pdf"), width = 25, height = 10)

for(i in 1:length(unique_oms)) {
  
    # Relative error of time series
    ts_re_om <- ts_re_df %>% filter(OM_Scenario == unique_oms[i],
                                    !str_detect(EM_Scenario, "0.5|1.0|2.0"),
                                    par_name %in% c("Spawning Stock Biomass",
                                                    "Total Biomass",
                                                    "Total Recruitment"))
    
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
    
    # Now loop through each time series component and print the plot out
    print(
      ggplot(ts_re_om %>% 
               filter(par_name != "Depletion"), aes(x = year, y = median))  +
        geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = time_comp, group = time_comp), alpha = 0.3) +
        geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
        geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
        facet_grid(par_name~EM_Scenario) +
        coord_cartesian(ylim = c(-0.4,0.4)) +
        scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
        scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
        labs(x = "Year", y = paste("RE"), fill = "Time", color = "Time",
             title = unique(ts_re_om$OM_Scenario)) +
        theme_matt() +
        # theme(aspect.ratio=1)
        theme(legend.position = "top", title = element_text(size = 20),
              axis.text = element_text(size = 13), strip.text = element_text(size = 13)) 
    )

} # end i

dev.off()


# Convergence Summary plot ------------------------------------------------

pdf(here(dir_out, "Convg_Sum.pdf"), width = 25, height = 10)

for(i in 1:length(unique_oms)) {
  
  conv_stat <- AIC_df %>% 
    mutate(time_comp = case_when(
      str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
      str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
      str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
    ),
      EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'),
      time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End",
                                                    "Terminal"))) %>% 
    filter(!str_detect(EM_Scenario, "0.5|1.0|2.0"),
           OM_Scenario == unique_oms[i]) %>% 
    group_by(EM_Scenario, time_comp) %>% 
    summarize(converged = sum(conv == "Converged")/200)

  # Set order for plot
    order <- vector()
    for(o in 1:length(plot_order)) {
      order[o] <- which(grepl(plot_order[o], x = unique(conv_stat$EM_Scenario)))
    } # end o loop

  # relvel factors
  conv_stat$EM_Scenario <- factor(conv_stat$EM_Scenario, 
                                  levels = c(unique(conv_stat$EM_Scenario)[order]))
  
  print(
    ggplot(conv_stat %>% 
             filter(EM_Scenario %in% c(unique(conv_stat$EM_Scenario))),
           mapping = aes(x = EM_Scenario, y = converged * 100,
                         group = time_comp, fill = time_comp))  +
      geom_col(position = "dodge", alpha = 0.85) +
      geom_hline(aes(yintercept = 50), col = "black", lty = 2, size = 0.5, alpha = 1) +
      facet_wrap(~EM_Scenario, nrow = 1, scales = "free_x") +
      scale_fill_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)]) +
      theme_matt() +
      labs(fill = "Time Component", y = "Convg Rate", 
           title = unique_oms[i]) +
      # theme(aspect.ratio=1) +
      theme(legend.position = "top", title = element_text(size = 20),
            strip.text = element_text(size = 15), 
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())  
  )
  
}

dev.off()


# AIC Plot ----------------------------------------------------------------

AIC_df <- AIC_df %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
  ),
  EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_'),
  time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End",
                                           "Terminal")))

# filter to 2 fleet models (different data)
twofleet_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "2Fl_")) %>% 
  group_by(sim, OM_Scenario, time_comp) %>% 
  mutate(min = min(AIC),
         min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "2 Fleet")

# filter to 1 fleet models
onefleet_aic <- AIC_df %>% 
  filter(str_detect(EM_Scenario, "1Fl_"),
         !str_detect(EM_Scenario, "0.5|1.0|2.0")) %>% 
  group_by(sim, OM_Scenario, time_comp) %>% 
  mutate(min = min(AIC),
         min_AIC = ifelse(AIC == min, 1, 0)) %>% 
  group_by(OM_Scenario, EM_Scenario, time_comp) %>% 
  summarize(n_minAIC = sum(min_AIC) / 200) %>% 
  mutate(fleet = "1 Fleet")

# AIC plot
pdf(file = here(dir_out, "AIC_SummaryPlot.pdf"), width = 10, height = 7)
# Two fleet plot
(twofleet_plot <- ggplot(twofleet_aic, 
                        aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC,
                        label = round(n_minAIC, 2))) +
  geom_tile(alpha = 0.85) +
  geom_text(size = 5) +
  facet_grid(~time_comp, scales = "free") +
  scale_fill_viridis_c() +
  theme(legend.position = "top")) 

# One fleet plot
(onefleet_plot <- ggplot(onefleet_aic, 
                        aes(x = OM_Scenario, y = EM_Scenario, fill = n_minAIC,
                       label = round(n_minAIC, 2))) +
  geom_tile(alpha = 0.85) +
  geom_text(size = 6) +
  facet_grid(~time_comp, scales = "free") +
  scale_fill_viridis_c() +
  theme_test() +
  theme(legend.position = "top")) 
  
dev.off()


# Min Max Solution Plot ---------------------------------------------------

minmax_df <- ts_are_df %>% 
  filter(!str_detect(EM_Scenario, "0.5|1.0|2.0"),
         # Get terminal year of peels/EMs
           (year %in% c(50) & str_detect(OM_Scenario, "Fast_") ) & str_detect(EM_Scenario, "Term_") |
           (year %in% c(70) & str_detect(OM_Scenario, "Slow_") ) & str_detect(EM_Scenario, "Term_") |
           (year %in% c(30) & str_detect(OM_Scenario, "Fast_") ) & str_detect(EM_Scenario, "TrxE_") |
           (year %in% c(50) & str_detect(OM_Scenario, "Slow_") ) & str_detect(EM_Scenario, "TrxE_") |
           (year %in% c(28) & str_detect(OM_Scenario, "Fast_") ) & str_detect(EM_Scenario, "Int_") |
           (year %in% c(39) & str_detect(OM_Scenario, "Slow_") ) & str_detect(EM_Scenario, "Int_"),
         par_name %in% c("Total Recruitment",
                         "Spawning Stock Biomass",
                         "Total Biomass")) %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE_") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int_") ~ "Fleet Intersect" # Fleet Transition Intersects
  ),
  time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End",
                                           "Terminal")),
  EM_Scenario = str_remove(EM_Scenario, 'Term_|TrxE_|Int_')) %>% 
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
pdf(file = here(dir_out, "MinMax_Summary.pdf"), width = 15, height = 10)
ggplot(minmax_df, aes(x = factor(OM_Scenario), y = factor(EM_Scenario), 
              fill = median, label = round(median, 2))) +
  geom_tile(alpha = 0.85) +
  geom_text(color = ifelse(minmax_df$median == minmax_df$max_median &
                             minmax_df$median != minmax_df$min_max_medians, "red", 
                           ifelse(a$median == minmax_df$min_max_medians, "green", "black")),
            size = 5.5) +
  facet_grid(par_name~time_comp, scales = "free_y") +
  scale_fill_viridis_c() +
  theme_test()
dev.off()

# Selectivity plots -------------------------------------------------------

# Subset to first year for OMs
om_slx_df_sub <- om_slx_df %>% 
  mutate(Sex = case_when(
    str_detect(Sex, "1") ~ "Female",
    str_detect(Sex, "2") ~ "Male")) %>% 
    filter(Year == 1)

# Filter only to terminal years of terminal, 10, 5
# Also do some residual munging on the dataframe
trunc_em <- em_slx_df %>% 
  filter(!str_detect(EM_Scenario, "0.5|1.0|2.0"),
           (Year %in% c(50) & str_detect(OM_Scenario, "Fast_") ) & str_detect(EM_Scenario, "Term_") |
           (Year %in% c(70) & str_detect(OM_Scenario, "Slow_") ) & str_detect(EM_Scenario, "Term_") |
           (Year %in% c(30) & str_detect(OM_Scenario, "Fast_") ) & str_detect(EM_Scenario, "TrxE_") |
           (Year %in% c(50) & str_detect(OM_Scenario, "Slow_") ) & str_detect(EM_Scenario, "TrxE_") |
           (Year %in% c(28) & str_detect(OM_Scenario, "Fast_") ) & str_detect(EM_Scenario, "Int_") |
           (Year %in% c(39) & str_detect(OM_Scenario, "Slow_") ) & str_detect(EM_Scenario, "Int_"),
            conv == "Converged") %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal", # Terminal Year
    str_detect(EM_Scenario, "TrxE") ~ "Fleet Trans End", # Fleet Transition End
    str_detect(EM_Scenario, "Int") ~ "Fleet Intersect" # Fleet Transition Intersects
  ),
  time_comp = factor(time_comp, levels = c("Fleet Intersect", "Fleet Trans End",
                               "Terminal")),
  Sex = case_when(
    str_detect(Sex, "1") ~ "Female",
    str_detect(Sex, "2") ~ "Male"), 
  EM_Scenario = str_remove(EM_Scenario, "Term_|TrxE_|Int_"))

pdf(here(dir_out, "Selex_Sum.pdf"), width = 22, height = 5)

for(i in 1:length(unique_oms)) {
  
  # Subset plot here for plotting purposes
  em_plot_df <- trunc_em %>% 
    filter(OM_Scenario == unique_oms[i]) %>% 
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
            legend.box = "vertical", strip.text = element_text(size = 12),
            aspect.ratio = 0.8) 
  )
  
} # end i

dev.off()


# Composition plots --------------------------------------------------------

# Get manual colors
colors <- viridis::viridis(n = 10)[c(1, 4, 10)]

# Filter comps to terminal year
om_prop_term <- om_comps_df %>% 
  group_by(OM_Scenario, Ages, Sexes, fleet_type) %>% 
  summarize(Prop = mean(Prop)) %>% 
  mutate(time_comp = "Terminal",
         Sex = case_when(
           str_detect(Sexes, "1") ~ "Female",
           str_detect(Sexes, "2") ~ "Male"),
         Fleet = case_when(
           str_detect(fleet_type, "One Fleet") ~ 1,
           str_detect(fleet_type, "Two Fleets") ~ 2
         )) %>% 
  filter(fleet_type == "One Fleet") %>% 
  dplyr::select(-Sexes)

# EM proportion munging
em_prop <- em_comps_df %>% 
  filter(!str_detect(EM_Scenario, "0.5|1.0|2.0"),
         !str_detect(EM_Scenario, "5_|10_"),
         conv == "Converged") %>%  # remove 2 fleet model
  mutate(EM_Scenario = str_remove(EM_Scenario, "Term_|10_|5_"),
         Sex = case_when(
           str_detect(Sex, "1") ~ "Female",
           str_detect(Sex, "2") ~ "Male")) %>% 
  group_by(Age, Fleet, Sex, OM_Scenario, EM_Scenario) %>% 
  summarize(Pred = mean(Prop)) 


# Set order for plot
order <- vector()
no2fl_plot_order <- plot_order[-1] # remove 2 fleets
for(o in 1:length(no2fl_plot_order)) {
  order[o] <- which(grepl(no2fl_plot_order[o], x = unique(em_prop$EM_Scenario)))
} # end o loop

# Now relevel factor for organizing plot
em_prop$EM_Scenario <- factor(em_prop$EM_Scenario, levels = unique(em_prop$EM_Scenario)[order])
  
# Left join EM and OM
em_om_df <- em_prop %>% 
  left_join(om_prop_term, 
            by = c("Sex", "Age" = "Ages", "Fleet", "OM_Scenario")) %>% 
  drop_na() %>% 
  mutate(Diff = Prop - Pred)

# now plot!
pdf(here(dir_out, "Comp_Fits_Sum.pdf"), width = 22, height = 8)
for(i in 1:length(unique_oms)) {
  
  # Filter to a given OM
  plot_df <- em_om_df %>% 
    filter(OM_Scenario == unique_oms[i])
  
  # Comps and predictions
  print(
    ggplot() +
      geom_col(plot_df, mapping = aes(x = Age, y = Prop, fill = Diff), alpha = 1,
               position = "identity", color = "grey", size = 0.5) +
      geom_line(plot_df, mapping = aes(x = Age, y = Pred), size = 0.85,
                color = "black", lty = 2) +
      facet_grid(Sex~EM_Scenario) +
      scale_fill_gradient2(midpoint = 0, oob = scales::squish) + 
      guides (fill = guide_colourbar(barwidth = 16, barheight = 1,
                                     title.position = "top")) +
      theme_matt() +
      labs(x = "Age", y = "Mean Proportion", title = unique_oms[i],
           fill = "Difference (Observed - Predicted)") +
      theme(title = element_text(size = 20), 
            aspect.ratio = 0.8)
  )

}
dev.off()
