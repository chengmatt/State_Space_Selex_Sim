# Purpose: To test TMB EM biases, debug, troubleshoot, and develop fxns for EM model
# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(TMB)
library(cowplot)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
source(here(fxn_path, "simulate_data.R"))
source(here(fxn_path, "Utility_fxns.R"))
source(here(fxn_path, "prepare_EM_input.R"))

compile_tmb(wd = here("src"), cpp = "EM.cpp")

# Path to general input biological parameters
# spreadsheet_path <- here("input", "EBS_Pollock_Inputs.xlsx")
spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")

# simulate data
simulate_data(fxn_path = fxn_path, 
              check_equil = FALSE,
              spreadsheet_path = spreadsheet_path, 
              rec_type = "BH",
              Start_F = c(0.108, 0.01 * 0.108), 
              Fish_Start_yr = c(1, 1), 
              Surv_Start_yr = c(1),  
              max_rel_F_M = c(0.01, 0.75), 
              desc_rel_F_M = c(NULL, NULL), 
              F_type = c("Const_Ramp_Const", "Const_Ramp_Const"),
              yr_chng = c(25, 25), 
              yr_chng_end = c(30, 30),
              fish_likelihood = c("multinomial", "multinomial"),
              srv_likelihood = "multinomial",
              Input_Fish_N_Max = c(1e3, 1e3), 
              Input_Srv_N_Max = c(200),
              fish_CV = c(0.1, 0.1),
              srv_CV = c(0.01), 
              catch_CV = c(0.0, 0.0), 
              Input_N_Fish_Time = "F_Vary", 
              Input_N_Fish_Fixed = c(1e3, 1e3),
              Mort_Time = "Constant", 
              q_Mean_Fish = c(0.05, 0.05), 
              q_Mean_Surv = 0.01, 
              fish_selex = c("logistic", "logistic"), 
              srv_selex = c("logistic"), 
              fish_pars = list(Fleet_1_L = matrix(data = c(3, 0.85, 5, 0.6), 
                                                  nrow = 2, byrow = TRUE),
                               Fleet_2_L = matrix(data = c(6, 0.85, 10, 0.4), 
                                                  nrow = 2, byrow = TRUE)), # fish fleet 2
              srv_pars = list(Fleet_3_SL = matrix(data = c(2, 0.85, 5, 0.3), 
                                                  nrow = 2, byrow = TRUE)), # survey fleet 1
              f_ratio = 0.5, m_ratio = 0.5)

# if switching to a single sex, be sure to change the nrow to the number of sexes,
# and to make sure the selex parameters for the fleets align n_pars * n_sexes
# e.g., (7, 0.8, 4, 0.3) for a logistic with two sexes, nrow = 2

plot_OM(path = here("figs", "Base_OM_Figs"), file_name = "OM_Check.pdf")

# Load in data ------------------------------------------------------------

par_all <- data.frame()
ts_all <- data.frame()

# Specify years here
years <- 1:(n_years - 1)

start.time <- Sys.time()
compile_tmb(wd = here("src"), cpp = "EM.cpp")

for(sim in 1:n_sims){

  # Prepare inputs here
  input <- prepare_EM_input(years = years,
                            n_fleets = 2, 
                            catch_cv = c(0.01, 0.01),
                            F_Slx_Blocks_Input = matrix(c(rep(0)),
                                                        nrow = length(years),
                                                        ncol = 2), # fishery blocks
                            S_Slx_Blocks_Input = matrix(c(0), # selectivity blocks
                                                        nrow = length(years), 
                                                        ncol = 1),
                            use_fish_index = FALSE,
                            Fish_Comp_Like_Model = c("multinomial", "multinomial"),
                            Srv_Comp_Like_Model = c("multinomial"),
                            rec_model = "BH", 
                            F_Slx_Model_Input = c("logistic", "logistic"),
                            S_Slx_Model_Input = c("logistic"), 
                            time_selex = "None",
                            n_time_selex_pars = NULL,
                            fix_pars = c(
                                         "ln_SigmaRec",
                                         "ln_q_fish", 
                                         "ln_h"),
                            sim = sim)


  # Run EM model here and get sdrep
  tryCatch(expr = model <- run_EM(data = input$data, parameters = input$parameters, 
                                  map = rlist::list.append(input$map), 
                                  n.newton = 3, 
                                  # random = c("ln_fish_selpars_re"),
                                  silent = F, getsdrep = TRUE), error = function(e){e})
  
  # Check model convergence
  convergence_status <- check_model_convergence(mle_optim = model$mle_optim,
                                                mod_rep = model$model_fxn,
                                                sd_rep = model$sd_rep,
                                                min_grad = 0.01)
  conv[sim] <- convergence_status$Convergence
  max_par[sim] <- convergence_status$Max_Grad_Par
  
  quants_df <- get_quants(sd_rep = model$sd_rep,
                 model_fxn = model$model_fxn, 
                 sim = sim, 
                 conv = conv[sim],
                 n_years = length(years),
                 ntrue_fish_fleets = dim(Fish_Age_Comps)[3], 
                 nmod_fish_fleets = input$data$n_fleets,
                 ages = ages, 
                 F_x = 0.4)
  
  ts_all <- rbind(quants_df$TS_df, ts_all)
  par_all <- rbind(quants_df$Par_df, par_all)
  print(paste("done w/  sim = ", sim))
  print(conv[sim])
  
}


# Summary Checks ----------------------------------------------------------

# Get percentiles
f_sum <- get_RE_precentiles(df = f_all %>% filter(conv == "Converged"), 
                            est_val_col = 1, true_val_col = 5,  
                            par_name = "Total Fishing Mortality", group_vars = "year")

ssb_sum <- get_RE_precentiles(df = ssb_all %>% filter(conv == "Converged"), 
                              est_val_col = 1, true_val_col = 5, 
                              par_name = "Spawning Stock Biomass", group_vars = "year")

rec_sum <- get_RE_precentiles(df = rec_all %>% filter(conv == "Converged"), 
                              est_val_col = 1, true_val_col = 5, 
                              par_name = "Total Recruitment", group_vars = "year")

fish_mu_age_sum <- get_RE_precentiles(df = fish_mu_age %>% filter(conv == "Converged"),
                                      est_val_col = 4, true_val_col = 6,
                                      par_name = "Mean Predicted Fishery Age",
                                      group_vars = c("year","fleet", "sex"))

srv_mu_age_sum <- get_RE_precentiles(df = srv_mu_age %>% filter(conv == "Converged"),
                                     est_val_col = 4, true_val_col = 6,
                                     par_name = "Mean Predicted Survey Age",
                                     group_vars = c("year","fleet", "sex"))

depletion_sum <-  get_RE_precentiles(df = depletion_all %>% filter(conv == "Converged"),
                                     est_val_col = 1, true_val_col = 5,
                                     par_name = "Depletion (SSB / SSB0)", group_vars = "year")

all <- rbind(rec_sum, ssb_sum, f_sum,fish_mu_age_sum, srv_mu_age_sum,
             depletion_sum) 


# Get relative error time series
(est_plot <- plot_RE_ts_ggplot(data = all, x = year, y = median, 
                               lwr_1 = lwr_80, upr_1 = upr_80,
                               lwr_2 = lwr_95, upr_2 = upr_95, 
                               facet_name = par_name))


# par(mfrow = c(3, 1))
plot_RE_ts_base(data = all, par_name = "Spawning Stock Biomass", ylim = c(-0.9, 0.9))
plot_RE_ts_base(data = all, par_name = "Total Fishing Mortality", ylim = c(-0.9, 0.9))
plot_RE_ts_base(data = all, par_name = "Total Recruitment", ylim = c(-0.9, 0.9))

# Parameter estimates
par_df <- par_all %>% mutate(RE = (mle_val - t ) / t) %>% 
  filter(conv == "Converged") %>% 
  mutate(CV = abs(mle_sd / trans_mle_val))
ssb0_all_df <- ssb0_all %>% mutate(RE = (mle_val - t ) / t) %>% 
  filter(conv == "Converged")

# hist((par_df$mle_val[par_df$type == "DMF1"] - 1.5) / 1.5)

# par_df %>% group_by(type) %>% summarize(mle_val = median(mle_val))

# Quick summary stats
par_sum <- get_RE_precentiles(df = par_all %>% filter(conv == "Converged"), 
                              est_val_col = 2, true_val_col = 7, 
                              par_name = "", group_vars = "type")

ggplot(par_df, aes(x = type, y = RE, fill = type)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(aes(yintercept = 0), size = 0.85, lty = 2, col = "black") +
  theme_bw() +
  # coord_cartesian(ylim = c(-1.5, 1.5)) +
  facet_wrap(~type, scale = "free") +
 ggthemes::scale_fill_colorblind() +
  theme(legend.position = "none")


(par_plot <- ggplot(par_df %>% filter(RE < 1), aes(x = RE, fill = type)) +
    geom_density(alpha = 0.2) +
    facet_wrap(~type, scales = "free") +
    geom_errorbarh(inherit.aes = FALSE, data = par_sum, 
                   aes(xmin = lwr_95, xmax = upr_95, y = 0, color = type),
                   height = 0, size = 2.5, alpha = 0.55, linetype = 1) +
    geom_point(inherit.aes = FALSE, data = par_sum, 
               aes(x= median, y = 0, color = type), size = 5, alpha = 0.95) +
    geom_vline(aes(xintercept = 0), linetype = 2,
               size = 0.85, col = "black", alpha = 1) +
    # coord_cartesian(xlim = c(-0.3, 0.3)) +
    ggsci::scale_color_jco() +
    ggsci::scale_fill_jco() +
    labs(x = "Relative Error", y = "Probability Density", linetype = "", color = "") +
    theme_bw() + 
    theme(strip.text = element_text(size = 13),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13, color = "black"),
          legend.position = "none", legend.text = element_text(size = 15)))

(ssb0_all_plot <- ggplot(ssb0_all_df %>% filter, aes(x = RE)) +
    geom_density(alpha = 0.2) +
    geom_vline(aes(xintercept = median(RE)), linetype = 2,
               size = 0.85, col = "red", alpha = 1) +
    geom_vline(aes(xintercept = 0), linetype = 2,
               size = 0.85, col = "black", alpha = 1) +
    # coord_cartesian(xlim = c(-0.3, 0.3)) +
    ggsci::scale_color_jco() +
    ggsci::scale_fill_jco() +
    labs(x = "Relative Error", y = "Probability Density", linetype = "", color = "") +
    theme_bw() + 
    theme(strip.text = element_text(size = 13),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13, color = "black"),
          legend.position = "none", legend.text = element_text(size = 15)))

plot_grid(par_plot, est_plot, ncol = 2, align = "hv", axis = "bl",
          rel_widths  = c(0.75, 1))

f_all %>% 
  filter(conv == "Converged") %>% 
  ggplot(aes(x = year, y = mle_val, group = sim)) +
  geom_line(color = "grey", size = 1) +
  geom_line(aes(y = t), color = "red", lty = 2, size = 1) +
  theme_bw()



# Look at Comps by Year ---------------------------------------------------

pred_fish_comps <- reshape2::melt(model$model_fxn$rep$pred_fish_age_comps) %>% 
  mutate(type = "pred")
names(pred_fish_comps) <- c("year", "age", "fleet", "sex", "props", "type")

true_fish_comps <- reshape2::melt(input$data$obs_fish_age_comps) %>% 
  mutate(type = "pred")
names(true_fish_comps) <- c("year", "age", "fleet", "sex", "props", "type")

ggplot() +
  geom_col(pred_fish_comps, mapping = aes(x = age, y = props), alpha = 0.5) +
  geom_line(true_fish_comps, mapping = aes(x = age, y = props,
                                           col = "pred"), size = 1) +
  # geom_line(caa, mapping = aes(x = age, y = props,
  #                              col = "true_caa"), size = 1, lty = 2) +
  scale_color_manual(name = c("pred", "true_caa"), 
                     values = c("black", "red")) +
  facet_wrap(~year) +
  theme_bw()

# Sum to 1 check
a = pred_fish_comps %>% group_by(year, fleet, sex) %>% 
  summarize(sum = sum(props))

b = true_fish_comps %>% group_by(year, fleet, sex) %>% 
  summarize(sum = sum(props))
