# Purpose: To make plots of self tests
# Creator: Matthew LH. Cheng
# Date 3/26/23

theme_matt <- function() {
  theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text= element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15)) 
}


# Set up -----------------------------------------------------------------

library(here)
library(tidyverse)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

# Paths
om_scenario_path <- here("output", "OM_Scenarios") # path to OM folder
dir_out <- here("output", "Self_Tests") # path to output folder

# Get results from self tests
self_test_res <- get_results(om_scenario_path = om_scenario_path)
param_df <- self_test_res$Parameter_Sum # parameter dataframe
ts_df <- self_test_res$TimeSeries_Sum # time series dataframe

# Create relative error and CV metrics for parameters
param_df <- param_df %>% mutate(RE = (mle_val - t) / t,
                                CV = (mle_var / mle_val) * 100) %>% 
  dplyr::select(-X)


# Plots -------------------------------------------------------------------

# Relative Error of Parameters --------------------------------------------

pdf(here(dir_out, "Parameter_RE_SelfTests.pdf"), height = 8, width = 15)
ggplot(param_df %>% filter(!is.na(t),
                           type != "SSB0"), aes(x = EM_Scenario, y = RE,
                                           fill = type)) +
  geom_boxplot(alpha = 0.65) +
  ggthemes::scale_fill_colorblind() +
  geom_hline(yintercept = 0, col = "blue", lty = 2, lwd = 1.3) +
  facet_wrap(~OM_Scenario, scales = "free", ncol = 4) +
  labs(x = "EM Scenarios", y = "Relative Error", fill = "Parameter Type") +
  theme_matt() +
  ylim(-0.3, 0.3) 
dev.off() 

# Relative Error of Time Series -------------------------------------------

# Time series parameters
ts_pars <- unique(ts_df$type)
ts_re_df <- data.frame() # empty dataframe to store

# Loop through to extract time series components
for(i in 1:length(ts_pars)) {
  # Get time series components here
  ts_comp <- get_RE_precentiles(df = ts_all %>% filter(type == ts_pars[i]),
                                est_val_col = 2, true_val_col = 6,
                                par_name = ts_pars[i], group_vars = c("year", 
                                                              "OM_Scenario",
                                                              "EM_Scenario"))
  
  ts_re_df <- rbind(ts_re_df, ts_comp)

} # end i loop


# Time Series Relative Error (Logistic Logistic) ---------------------------------------------------------

pdf(here(dir_out, "RE_Time_Series_Fl_LL.pdf"), height = 8, width = 15)
for(i in 1:length(ts_pars)) {
  print(
    ggplot(ts_re_df %>% filter(str_detect(EM_Scenario, "LL"), par_name == ts_pars[i]), 
           aes(x = year, y = median) ) +
      geom_ribbon(aes(ymin = lwr_80, ymax = upr_80), alpha = 0.8, fill = "grey4") +
      geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.5, fill = "grey2") +
      geom_line(col = "white", size = 3) + 
      geom_hline(aes(yintercept = 0), col = "red", lty = 2, size = 1.3, alpha = 1) +
      facet_grid(EM_Scenario~OM_Scenario) +
      labs(x = "Year", y = "Relative Error", scales = "free_y", title = ts_pars[i]) +
      theme_matt()
  )
} # end i

dev.off()

# Time Series Relative Error (Logistic Gamma) ---------------------------------------------------------

pdf(here(dir_out, "RE_Time_Series_Fl_LG.pdf"), height = 8, width = 15)
for(i in 1:length(ts_pars)) {
  print(
    ggplot(ts_re_df %>% filter(str_detect(EM_Scenario, "LG"), par_name == ts_pars[i]), 
           aes(x = year, y = median) ) +
      geom_ribbon(aes(ymin = lwr_80, ymax = upr_80), alpha = 0.8, fill = "grey4") +
      geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.5, fill = "grey2") +
      geom_line(col = "white", size = 3) + 
      # geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.8, alpha = 0.85) +
      geom_hline(aes(yintercept = 0), col = "red", lty = 2, size = 1.3, alpha = 1) +
      facet_grid(EM_Scenario~OM_Scenario) +
      labs(x = "Year", y = "Relative Error", scales = "free_y", title = ts_pars[i]) +
      theme_matt()
  )
} # end i

dev.off()






