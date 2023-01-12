# Purpose: To test TMB EM biases, debug, troubleshoot, and develop fxns and EM model
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
  spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")
  
  # simulate data
  simulate_data(fxn_path = fxn_path, spreadsheet_path = spreadsheet_path, 
                check_equil = FALSE, rec_type = "mean_rec",
                n_years = 101, Start_F = c(0.01, 0.01), 
                Fish_Start_yr = c(70, 70), Surv_Start_yr = c(70), 
                max_rel_F_M = c(1.5, 1.5), desc_rel_F_M = c(0.15, NULL), 
                F_type = c("Contrast", "Const_Inc"), yr_chng = c(96), 
                fish_Neff_max = c(150, 150), srv_Neff_max = c(10), fish_CV = c(0.1, 0.1),
                srv_CV = c(0.1), catch_CV = c(0, 0), Neff_Fish_Time = "F_Vary", fixed_Neff = c(5, 5),
                Mort_Time = "Constant", q_Mean_Fish = c(0.05, 0.08), q_Mean_Surv = 0.01, 
                Rec_Dev_Type = "iid", rho_rec = NA, 
                fish_selex = c("logistic", "logistic"), srv_selex = c("logistic"), 
                fish_pars = list(Fleet_1_L = matrix(data = c(6, 0.8), nrow = 1, byrow = TRUE),
                                 Fleet_1_L = matrix(data = c(7, 0.8), nrow = 1, byrow = TRUE)),
                srv_pars = list(Fleet_3_SL = matrix(data = c(4,0.8), nrow = 1, byrow = TRUE)), 
                f_ratio = 1, m_ratio = 0)
  
  plot_OM(path = here("figs", "Base_OM_Figs"), file_name = "OM_Check.pdf")

  # Load in data ------------------------------------------------------------
  
  ssb_all <- data.frame()
  f_all <- data.frame()
  rec_all <- data.frame()
  biom_all <- data.frame()
  par_all <- data.frame()
  max_par <- vector()
  fish_mu_age <- data.frame()
  srv_mu_age <- data.frame()
  conv <- vector()
  depletion_all <- data.frame()
  
  for(sim in 1:n_sims){
  
  # Prepare inputs here
  input <- prepare_EM_input(years = Fish_Start_yr[1]:(n_years - 1),
                   n_fleets = 1, 
                   catch_cv = c(0.03, 0.03),
                   F_Slx_Blocks_Input = matrix(c(rep(0, 31)),
                                               nrow = length(years), ncol = 1), # fishery blocks
                   S_Slx_Blocks_Input = matrix(c(0), # selectivity blocks
                                               nrow = length(years), ncol = 1),
                   use_catch = TRUE,
                   use_fish_index = FALSE,
                   use_srv_index = TRUE,
                   use_fish_comps = TRUE,
                   use_srv_comps = TRUE,
                   rec_model = 0, 
                   F_Slx_Model_Input = c("logistic", "logistic"),
                   S_Slx_Model_Input = c("logistic"), 
                   Sex_Ratio = as.vector(1),
                   sim = sim)
    
    # Map to fix parameters
    map <- list(ln_SigmaRec = factor(NA),
    # ln_fish_selpars = factor(c(c(1,2,3), NA, c(4,5,6), NA)),    
    ln_q_fish = factor(rep(NA, 1)))

  # Run EM model here and get sdrep
  model <- run_EM(data = input$data, parameters = input$parameters, 
         map = map, n.newton = 5, silent = TRUE, getsdrep = TRUE)

  # Check model convergence
  convergence_status <- check_model_convergence(mle_optim = model$mle_optim, 
                                                mod_rep = model$model_fxn,
                                                sd_rep = model$sd_rep, min_grad = 0.001)
  conv[sim] <- convergence_status$Convergence
  max_par[sim] <- convergence_status$Max_Grad_Par

  # Get parameter estimates
  M_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_M", log = TRUE) %>% 
    mutate(t = mean(Mort_at_age), type = "mortality", sim = sim, conv = conv[sim])
  q_srv_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_q_srv", log = TRUE) %>% 
    mutate(t = mean(q_Surv), type = "q_surv", sim = sim, conv = conv[sim])
  srv_sel_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_srv_selpars", log = TRUE) %>% 
    mutate(t = c(4, 0.8), type = c("a50_srv", "k_srv"), sim = sim, conv = conv[sim])
  meanrec_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_MeanRec", log = TRUE) %>% 
    mutate(t = exp(2.75), type = "meanrec", sim = sim, conv = conv[sim])
  
  # Bind parameter estimates
  par_all <- rbind(M_df, q_srv_df, srv_sel_df, meanrec_df, par_all)

  # Recruitment
  rec_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Total_Rec") %>% 
    mutate(t = rec_total[70:(n_years-1),sim], sim = sim, conv = conv[sim],
           year = 70:(n_years-1))
  rec_all <- rbind(rec_df, rec_all)
  
  # Check SSB
  ssb_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "SSB") %>% 
    mutate(t = SSB[Fish_Start_yr[1]:(n_years-1), sim], sim = sim, conv = conv[sim],
           year = 70:(n_years-1))
  ssb_all <- rbind(ssb_df, ssb_all)
  
  # ssb_all %>%
  #   filter(sim %in% c(1:10)) %>%
  #   ggplot(aes(x = year, y = mle_val, ymin = lwr_95, ymax = upr_95)) +
  #   geom_line() +
  #   facet_wrap(~sim) +
  #   geom_line(aes(y = t), color = "red") +
  #   geom_ribbon(alpha = 0.3)

  # Check F
  f_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Total_Fy") %>% 
    mutate(t = rowSums(fish_mort[Fish_Start_yr[1]:(n_years-1),,sim]), sim = sim, conv = conv[sim],
           year = 70:(n_years-1))
  f_all <- rbind(f_all, f_df)
  
  # Check depletion rates
  depletion_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Depletion") %>%
    mutate(t = (SSB[Fish_Start_yr[1]:(n_years-1),sim]/SSB[Fish_Start_yr[1],sim]),
           sim = sim, conv = conv[sim], year = Fish_Start_yr[1]:(n_years-1))
  depletion_all <- rbind(depletion_df, depletion_all)
  # 
  # # Check survey mean age
  srv_mean_ages <- extract_mean_age_vals(mod_rep = model$model_fxn, comp_name = "pred_srv_age_comps",
                        bins = ages, comp_start_yr = Fish_Start_yr[1], sim = sim,
                        n_fish_true_fleets = NULL) %>% mutate(conv = conv[sim])
  
  srv_mu_age <- rbind(srv_mu_age, srv_mean_ages)
  
  # # Check fishery mean age
  fish_mean_ages <- extract_mean_age_vals(mod_rep = model$model_fxn, comp_name = "pred_fish_age_comps",
                                         bins = ages, comp_start_yr = Fish_Start_yr[1], sim = sim,
                                         n_fish_true_fleets = 2) %>% mutate(conv = conv[sim])
  fish_mu_age <- rbind(fish_mu_age, fish_mean_ages)

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

# biom_sum <- get_RE_precentiles(df = biom_all %>% filter(conv == "Converged"), 
#                                est_val_col = 1, true_val_col = 8, 
#                                par_name = "Total Biomass", group_vars = "year")

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

all <- rbind(rec_sum, ssb_sum, f_sum, depletion_sum, srv_mu_age_sum, fish_mu_age_sum) 

# Parameter estimates
par_df <- par_all %>% mutate(RE = (mle_val - t ) / t) 

# Quick summary stats
par_sum <- get_RE_precentiles(df = par_all %>% filter(conv == "Converged"), 
                              est_val_col = 2, true_val_col = 7, 
                              par_name = "", group_vars = "type")
  
(est_plot <- ggplot(all, aes(x = year, y = median)) +
  # geom_ribbon(aes(ymin = lwr_70, ymax = upr_70), alpha = 0.7, fill = "grey") +
  geom_ribbon(aes(ymin = lwr_80, ymax = upr_80), alpha = 0.5, fill = "grey4") +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.3, fill = "grey2") +
  # geom_ribbon(aes(ymin = lwr_100, ymax = upr_100), alpha = 0.4, fill = "grey") +
  geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.8, alpha = 0.85) +
  # geom_line( color = "white", size = 1,alpha = 1) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 1, alpha = 1) +
  facet_wrap(~par_name, scales = "free", ncol = 2) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  labs(x = "Year", y = "Relative Error") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13, color = "black")))

(par_plot <- ggplot(par_df, aes(x = RE, fill = type)) +
geom_density(alpha = 0.2) +
facet_wrap(~type, scales = "free", nrow = 2) +
geom_errorbarh(inherit.aes = FALSE, data = par_sum, 
               aes(xmin = lwr_95, xmax = upr_95, y = 0, color = type),
               height = 0, size = 2.5, alpha = 0.55, linetype = 1) +
# geom_errorbarh(inherit.aes = FALSE, data = par_sum, 
#                  aes(xmin = lwr_70, xmax = upr_70, y = 0, color = type),
#                  height = 0, size = 2.5, alpha = 1, linetype = 1) +
geom_point(inherit.aes = FALSE, data = par_sum, 
           aes(x= median, y = 0, color = type), size = 5, alpha = 0.95) +
  geom_vline(aes(xintercept = 0), linetype = 2,
             size = 0.85, col = "black", alpha = 1) +
    # coord_cartesian(xlim = c(-1, 1)) +
# ggsci::scale_color_jco() +
# ggsci::scale_fill_jco() +
  labs(x = "Relative Error", y = "Probability Density", linetype = "", color = "") +
theme_bw() + 
theme(strip.text = element_text(size = 13),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13, color = "black"),
      legend.position = "none", legend.text = element_text(size = 15)))

plot_grid(par_plot, est_plot, ncol = 1, align = "hv", axis = "bl",
          rel_heights = c(0.70, 1))

