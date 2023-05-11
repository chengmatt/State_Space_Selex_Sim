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
param_df <- data.table::fread(here("output", "Parameter_Summary.csv")) # parameters
ts_df <- data.table::fread(here("output", "TimeSeries_Summary.csv")) # time series values 
ts_re_df <- data.table::fread(here("output", "TimeSeries_RE.csv")) # time series relative error (converged runs only)
om_slx_df <- data.table::fread(here("output", "OM_Fish_Selex.csv")) # OM Selectivity Values
em_slx_df <- data.table::fread(here("output", "EM_Fish_Selex.csv")) # EM Selectivity Values
om_comps_df <- data.table::fread(here("output", "OM_Fish_Comps.csv")) # OM Composition values
em_comps_df <- data.table::fread(here("output", "EM_Fish_Comps.csv")) # EM Composition values


# Unique oms and other components
unique_oms <- unique(param_df$OM_Scenario) # unique oms
ts_pars <- unique(ts_re_df$par_name) # unique parameter names
 
# Get plot order for EM models
plot_order <- c("2Fl", "1Fl_L_TI", "1Fl_Gam_TI", "1Fl_ExpL_TI",
                "1Fl_L_TVBlk", "1Fl_Gam_TVBlk", "1Fl_ExpL_TVBlk",
                "1Fl_L_TVCont", "1Fl_Gam_TVCont", "1Fl_ExpL_TVCont")

# Parameter Summary plot -------------------------------------------------------------
param_df <-  param_df %>% 
  mutate(
time_comp = case_when(
      str_detect(EM_Scenario, "Term_") ~ "Terminal",
      str_detect(EM_Scenario, "10_") ~ "10 Years Post Fl Chg",
      str_detect(EM_Scenario, "5_") ~ "5 Years Post Fl Chg"
    ),
  EM_Scenario = str_remove(EM_Scenario, "Term_|10_|5_"),
  mle_val = as.numeric(mle_val),
  t = as.numeric(t),
  mle_sd = as.numeric(mle_sd),
  RE = (mle_val - t) / t,
  CV = mle_sd / mle_val,
  time_comp = factor(time_comp, levels = c("5 Years Post Fl Chg", 
                                           "10 Years Post Fl Chg",
                                           "Terminal"))) %>% 
  filter(conv == "Converged") # only converged runs

pdf(here(dir_out, "Param_Sum.pdf"), width = 25, height = 10)

for(i in 1:length(unique_oms)) {
  
  # Filter to relevant components for parameters
  om_scenario_params <- param_df %>% filter(str_detect(OM_Scenario, unique_oms[i]),
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
  
  pt_rg <- om_scenario_params %>% 
    group_by(EM_Scenario, time_comp, type) %>% 
    summarize(median = median(RE), 
              lwr_95 = quantile(RE, 0.025),
              upr_95 =  quantile(RE, 0.975))
    
  
  # plot now!
  print(
    ggplot(pt_rg, aes(x = factor(EM_Scenario), y = median, ymin = lwr_95, ymax = upr_95,
                                   color = time_comp, fill = time_comp)) +
      geom_pointrange(position = position_dodge2(width = 1), size = 1.5, linewidth = 1) +
      geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 0.5, alpha = 1) +
      facet_grid(type~EM_Scenario, scales = "free") +
      coord_cartesian(ylim = c(-0.75,0.75)) +
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

ts_re_df <-  ts_re_df %>% 
  mutate(time_comp = case_when(
    str_detect(EM_Scenario, "Term_") ~ "Terminal",
    str_detect(EM_Scenario, "10_") ~ "10 Years Post Fl Chg",
    str_detect(EM_Scenario, "5_") ~ "5 Years Post Fl Chg"
  ),
  EM_Scenario = str_remove(EM_Scenario, "Term_|10_|5_"))
  
pdf(here(dir_out, "TS_Sum.pdf"), width = 25, height = 10)

for(i in 1:length(unique_oms)) {
  
    # Relative error of time series
    ts_re_om <- ts_re_df %>% filter(OM_Scenario == unique_oms[i],
                                    !str_detect(EM_Scenario, "0.5|1.0|2.0"))
    
    # Set order for plot
    if(length(plot_order) == length(unique(ts_re_om$EM_Scenario))) {
      order <- vector()
      for(o in 1:length(plot_order)) {
        order[o] <- which(grepl(plot_order[o], x = unique(ts_re_om$EM_Scenario)))
      } # end o loop
    } # end if
    
    # Now relevel factor for organizing plot
    ts_re_om <- ts_re_om %>% 
      mutate(EM_Scenario = factor(EM_Scenario, levels = unique(EM_Scenario)[order]),
             time_comp = factor(time_comp, levels = c("5 Years Post Fl Chg",
                                           "10 Years Post Fl Chg",
                                           "Terminal")))
    
    # Now loop through each time series component and print the plot out
    print(
      ggplot(ts_re_om %>% 
               filter(par_name != "Depletion"), aes(x = year, y = median))  +
        geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = time_comp, group = time_comp), alpha = 0.3) +
        geom_line(linewidth = 2, alpha = 1, aes(color = time_comp)) +
        geom_hline(aes(yintercept = 0), col = "black", lty = 2, linewidth = 0.5, alpha = 1) +
        facet_grid(par_name~EM_Scenario) +
        coord_cartesian(ylim = c(-0.75,0.75)) +
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
  
  # Get convergence statistics here
  conv_stat <- ts_df %>% 
    mutate(time_comp = case_when(
      str_detect(EM_Scenario, "Term_") ~ "Terminal",
      str_detect(EM_Scenario, "10_") ~ "10 Years Post Fl Chg",
      str_detect(EM_Scenario, "5_") ~ "5 Years Post Fl Chg"
    ),
    EM_Scenario = str_remove(EM_Scenario, "Term_|10_|5_")) %>% 
    filter(year == 1, type == "Spawning Stock Biomass",
           OM_Scenario == unique_oms[i],
           !str_detect(EM_Scenario, "0.5|1.0|2.0")) %>% 
    group_by(EM_Scenario, time_comp) %>% 
    summarize(converged = sum(conv == "Converged") / 200) %>% 
    mutate(time_comp = factor(time_comp, levels = c("5 Years Post Fl Chg",
                                                    "10 Years Post Fl Chg",
                                                    "Terminal")))
  
  # Set order for plot
  if(length(plot_order) == length(unique(conv_stat$EM_Scenario))) {
    order <- vector()
    for(o in 1:length(plot_order)) {
      order[o] <- which(grepl(plot_order[o], x = unique(conv_stat$EM_Scenario)))
    } # end o loop
  } # end if
  
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
  filter(Year %in% c(30, 35, 50),
         !str_detect(EM_Scenario, "0.5|1.0|2.0"),
         conv == "Converged") %>% 
  mutate(time_comp = case_when(
    Year == 50 ~ "Terminal",
    Year == 35 ~ "10 Years Post Fl Chg",
    Year == 30 ~ "5 Years Post Fl Chg"
  ),
  time_comp = factor(time_comp, levels = c('5 Years Post Fl Chg', "10 Years Post Fl Chg",
                                           "Terminal")),
  Sex = case_when(
    str_detect(Sex, "1") ~ "Female",
    str_detect(Sex, "2") ~ "Male"), 
  EM_Scenario = str_remove(EM_Scenario, "Term_|10_|5_"))

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
  if(length(plot_order) == length(unique(em_plot_df$EM_Scenario))) {
  order <- vector()
  for(o in 1:length(plot_order)) {
    order[o] <- which(grepl(plot_order[o], x = unique(em_plot_df$EM_Scenario)))
  } # end o loop
} # end if
  
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
  
  # Residuals of comps
  # print(
  #   ggplot(plot_df, mapping = aes(x = Age, y = Diff)) +
  #     geom_line(color = "black", size = 3.3) +
  #     geom_line(mapping = aes(color = Diff), size = 3) +
  #     geom_point(pch = 21, color = "black", mapping = aes(fill = Diff),
  #                size = 5) +
  #     geom_hline(aes(yintercept = 0), lty = 2, size = 1.5) +
  #     facet_grid(Sex~EM_Scenario) +
  #     scale_color_gradient2(midpoint = 0, oob = scales::squish) +
  #     scale_fill_gradient2(midpoint = 0, oob = scales::squish) +
  #     guides (fill = guide_colourbar(barwidth = 16, barheight = 1,
  #                                    title.position = "top"),
  #             color = guide_colourbar(barwidth = 16, barheight = 1,
  #                                     title.position = "top")) +
  #     theme_matt() +
  #     labs(x = "Age", title = unique_oms[i],
  #          y = "Difference (Observed - Predicted)",
  #          fill = "Difference (Observed - Predicted)",
  #          color = "Difference (Observed - Predicted)") 
  # )
  
}
dev.off()
