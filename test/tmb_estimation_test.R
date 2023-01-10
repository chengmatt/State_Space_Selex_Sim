# Purpose: To test TMB EM biases and to debug + troubleshoot
  # Set up ------------------------------------------------------------------
  
  library(here)
  library(tidyverse)
  library(TMB)
  library(cowplot)
  
  # Load in all functions into the environment
  fxn_path <- here("R_scripts", "functions")
  source(here(fxn_path, "simulate_data.R"))
  source(here("R_scripts", "functions", "Utility_fxns.R"))
  
  compile_tmb(wd = here("src"), cpp = "EM.cpp")
  
  # Path to general input biological parameters
  spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")
  
  # simulate data
  simulate_data(fxn_path = fxn_path, spreadsheet_path = spreadsheet_path, 
                check_equil = FALSE, rec_type = "mean_rec",
                n_years = 101, Start_F = c(0.01, 0.01), 
                Fish_Start_yr = c(70, 70), Surv_Start_yr = c(70), 
                max_rel_F_M = c(1.5, 1), desc_rel_F_M = c(0.15, NULL), 
                F_type = c("Contrast", "Const_Inc"), yr_chng = c(86), 
                fish_Neff = c(150, 150), srv_Neff = c(150), fish_CV = c(0.1, 0.1),
                srv_CV = c(0.1), catch_CV = c(0, 0), Neff_Fish_Time = "F_Vary", fixed_Neff = c(30, 30),
                Mort_Time = "Constant", q_Mean_Fish = c(0.05, 0.08), q_Mean_Surv = 0.01, 
                Rec_Dev_Type = "iid", rho_rec = NA, 
                fish_selex = c("logistic", "logistic"), srv_selex = c("logistic"), 
                fish_pars = list(Fleet_1_L = matrix(data = c(6, 0.8), nrow = 1, byrow = TRUE),
                                 Fleet_1_L = matrix(data = c(4, 0.8), nrow = 1, byrow = TRUE)),
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
  
  ages <- 1:30
  years <- Fish_Start_yr[1]:(n_years - 1)
  n_sexes <- n_sex
  n_fleets <- n_fish_fleets
  n_fish_indices <- 2
  n_srv_indices <- 1
  
  for(sim in 1:n_sims){
  
  # Get observed catches
  obs_catches <- matrix(Catch_agg[(Fish_Start_yr[1]:(n_years-1)),,sim], nrow = length(years), ncol = n_fish_fleets)

  # Observed fishery age comps
  obs_fish_age_comps <- array(data = NA, dim = c(length(years), length(ages), n_fleets))
  
  for(f in 1:n_fish_fleets) { # needs to loop thorugh transpose and apply for ez retnetion of array dimensions
    obs_fish_age_comps[,,f] <- t(
      apply(Fish_Age_Comps[Fish_Start_yr[1]:(n_years - 1),,f,,sim], MARGIN = 1, 
            FUN=function(x) { x/sum(x) })
    )
  }
  
  obs_fish_age_Neff <- matrix(fish_Neff[Fish_Start_yr[1]:(n_years - 1),], nrow = length(years), ncol = n_fish_fleets)
  
  # Observed survey age comps
  obs_srv_age_comps <- t(apply(Survey_Age_Comps[Fish_Start_yr[1]:(n_years - 1),,,,sim], MARGIN = 1, 
                               FUN=function(x) { x/sum(x) }))
  
  obs_srv_age_Neff <- matrix(srv_Neff[Fish_Start_yr[1]:(n_years - 1),], nrow = length(years), ncol = n_srv_fleets)
  
  # Observed fishery indices
  obs_fish_indices <-as.matrix( Fishery_Index_Agg[Fish_Start_yr[1]:(n_years - 1),,sim], nrow = length(years), ncol = n_fish_fleets)
  obs_srv_indices <- as.matrix(Survey_Index_Agg[Fish_Start_yr[1]:(n_years - 1),,sim], nrow = length(years), ncol = n_srv_fleets)
  
  WAA <- wt_at_age[Fish_Start_yr[1]:(n_years),,1,sim]
  MatAA <- mat_at_age[Fish_Start_yr[1]:(n_years),,1,sim]
  Sex_Ratio <- c(1,0)
  
  biom_df <- melt(Biom_at_age)
  names(biom_df) <- c("Year", "Age", "Sex", "Sim", "Biomass")
  
  # Biomass aggregated
  biom_df <- biom_df %>% 
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim))) %>% 
    filter(Sim == sim,
           Year >= Fish_Start_yr[1] & Year < 101) %>% 
    group_by(Year) %>% 
    summarize(Biomass = sum(Biomass, na.rm = TRUE))
  
  fish_cv <- c(0.1, 0.1)
  srv_cv <- 0.1
  catch_cv <- c(0.01, 0.01)
  
  n_fish_comps = 2
  n_srv_comps = 1
  
  F_Slx_Blocks <- matrix(c(0), nrow = length(years), ncol = n_fleets)
  S_Slx_Blocks <- matrix(c(0), nrow = length(years), ncol = n_srv_fleets)
  
  # set up indicators for whether or not we want to use certain data sources
  use_catch <- matrix(1, nrow = n_years, ncol = n_fish_fleets)
  use_fish_index <- matrix(1, nrow = n_years, ncol = n_fish_indices)
  use_srv_index <- matrix(1, nrow = n_years, ncol = n_srv_indices)
  use_fish_comps <- array(1, dim = c(n_years, n_fish_comps, n_sex))
  use_srv_comps <- array(1, dim = c(n_years, n_srv_comps, n_sex))
  
  # TMB Section -------------------------------------------------------------
  
  
  # Fill in list for data
  data <- list( ages = ages, years = years,
                n_sexes = n_sexes, n_fleets = 2,
                n_fish_indices = n_fish_indices, n_srv_indices = n_srv_indices,
                obs_catches = obs_catches, 
                obs_fish_age_comps = array(obs_fish_age_comps, dim = c(31, 30, 2, 1)),
                obs_fish_age_Neff = obs_fish_age_Neff, 
                obs_srv_age_comps = array(obs_srv_age_comps, dim = c(31, 30, 1, 1)),
                obs_srv_age_Neff = obs_srv_age_Neff, obs_fish_indices =  obs_fish_indices,
                obs_srv_indices = obs_srv_indices, WAA = array(WAA, dim = c(32, 30, 1)), 
                MatAA = array(MatAA, dim = c(32, 30, 1)), F_Slx_Blocks = F_Slx_Blocks,
                S_Slx_Blocks = S_Slx_Blocks, 
                # Init_N = as.vector(N_at_age[70,,,sim]),
                Sex_Ratio = as.vector(c(1)),  rec_model = 0, 
                fish_cv = fish_cv, srv_cv = srv_cv, catch_cv = catch_cv,
                F_Slx_model = as.vector(c(0, 0)), n_fish_comps = 2, n_srv_comps = 1,
                S_Slx_model = as.vector(0),
                use_catch = use_catch, use_fish_index= use_fish_index,
                use_srv_index= use_srv_index, use_fish_comps = use_fish_comps,
                use_srv_comps = use_srv_comps
  )
  
  # Define parameter inits here
  parameters <- list(ln_SigmaRec = 0.6, ln_MeanRec = 2.75,
                     ln_M = log(0.1),  
                     ln_fish_selpars = log(array(c(6, 4, 0.8, 0.8), dim = c(2, 1, 1, 2))),
                     ln_srv_selpars = array(5, dim = c(1, 1, 1, 2)),
                     ln_N1_Devs = log(rnorm(length(ages)-2,5, 1)),
                     ln_Fy = log(as.matrix(fish_mort[Fish_Start_yr[1]:((n_years) -1),,sim])),
                     ln_q_fish = as.matrix(rep(log(0.1), n_fish_fleets)), 
                     ln_q_srv = as.matrix(rep(log(0.01), n_srv_fleets)),
                     ln_RecDevs = rec_devs[Fish_Start_yr[1]:((n_years) -1),sim])
  
  map <- list(ln_SigmaRec = factor(NA))
              # ln_fish_selpars = factor(rep(NA, 4)),
              # ln_M = factor(NA),
              # ln_N1_Devs = factor(rep(NA, length(ages)-2)),
              # ln_MeanRec = factor(NA),
              # ln_q_fish = factor(rep(NA, 2)),
              # ln_q_srv = factor(NA),
              # ln_Fy = factor(rep(NA, 62)),
              # ln_RecDevs = factor(rep(NA, 31)))
  
  compile_tmb(wd = here("src"), cpp = "EM.cpp")
  
  # sum_FAA
  # (sum_FAA = fish_mort[70,1,sim] * Fish_selex_at_age[1,26,1,,sim] +
  #   fish_mort[70,2,sim] * Fish_selex_at_age[1,26,1,,sim])
  # # 
  # # # ZAA
  # (ZAA = (fish_mort[70,1,sim] * Fish_selex_at_age[1,26,1,,sim]) +
  # (fish_mort[70,2,sim] * Fish_selex_at_age[1,26,2,,sim]) + 0.1)
  # # 
  # # # NAA
  # (NAA = N_at_age[70, 26,,sim])
  # my_model$rep$NAA[4,26,1] * my_model$rep$FAA[4,26,1,1] * (1 - exp( -my_model$rep$ZAA[4,26,1])) / my_model$rep$ZAA[4,26,1]
  # # 
  # # # CAA
  # (CAA = Catch_at_age[70,26 , 1,, sim])
  # # 
  # # # Agg catch
  # Catch_agg[70, 1, sim]
  # # 
  # FISH_MORT = (fish_mort[70,1,sim] * Fish_selex_at_age[1,26,1,,sim] )
  # NAA * FISH_MORT * (1 - exp(-ZAA)) / ZAA
  # 
  # FISH_MORT * NAA * (1-exp(-ZAA)) / ZAA
  
  # plot(my_model$rep$NAA[4,,], type = "l")
  # lines(N_at_age[73,,,sim], col = "red")
  
  # Make ADFun
  my_model <- MakeADFun(data, parameters, map, DLL="EM", silent = T)
  mle_optim <- stats::nlminb(my_model$par, my_model$fn, my_model$gr, 
                             control = list(iter.max = 1e5, eval.max = 1e5))
  
  # Additional newton steps to take
  add_newton(n.newton = 5, ad_model = my_model, mle_optim = mle_optim)
  
  my_model$rep <- my_model$report(par = mle_optim$par)
  sd_rep <- TMB::sdreport(my_model)
  
  # Check model convergence
  convergence_status <- check_model_convergence(mle_optim = mle_optim, mod_rep = my_model,
                                                sd_rep = sd_rep, min_grad = 0.001)
  conv[sim] <- convergence_status$Convergence
  max_par[sim] <- convergence_status$Max_Grad_Par
  
  # plot(my_model$rep$pred_catches[,1], type = "l")
  # lines(Catch_agg[Fish_Start_yr[1]:(n_years-1),1,sim], col = "red")
  # plot(my_model$rep$pred_catches[,2], type = "l")
  # lines(Catch_agg[Fish_Start_yr[1]:(n_years-1),2,sim], col = "red")
  # plot(sd_rep$value[names(sd_rep$value)=="Total_Fy"], type = "l")
  # lines(rowSums(fish_mort[Fish_Start_yr[1]:(n_years-1),,sim]), col = "red")
  # plot(my_model$rep$NAA[1,,], type = "l")
  # lines(N_at_age[70,,,sim], col = "red")
  # plot(my_model$rep$F_Slx[1,,1,], type = "l")
  # lines(Fish_selex_at_age[1,,1,1,1], type = "l", col = "red")
  # plot(my_model$rep$F_Slx[1,,2,], type = "l")
  # lines(Fish_selex_at_age[1,,2,1,1], type = "l", col = "red")
  # plot(my_model$rep$S_Slx[1,,1,], type = "l")
  # lines(Surv_selex_at_age[1,,1,1,1], type = "l", col = "red")
  # plot(sd_rep$value[names(sd_rep$value)=="SSB"], type = "l")
  # lines(SSB[Fish_Start_yr[1]:(n_years-1),sim], col = "red")

  # Get parameter estimates
  M_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_M", log = TRUE) %>% 
    mutate(t = mean(Mort_at_age), type = "mortality", sim = sim, conv = conv[sim])
  q_fish_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_q_fish", log = TRUE) %>%
    mutate(t = q_Fish[1,,sim], type = "q_fish", sim = sim, conv = conv[sim])
  q_srv_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_q_srv", log = TRUE) %>% 
    mutate(t = mean(q_Surv), type = "q_surv", sim = sim, conv = conv[sim])
  fish_sel_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_fish_selpars", log = TRUE) %>%
    mutate(t = c(6,6,0.8,0.8), type = c("a50_f1", "a50_f2", "d1", "d2"), sim = sim, conv = conv[sim])
  srv_sel_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_srv_selpars", log = TRUE) %>% 
    mutate(t = c(4, 0.8), type = c("a50_srv", "k_srv"), sim = sim, conv = conv[sim])
  meanrec_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_MeanRec", log = TRUE) %>% 
    mutate(t = exp(2.70), type = "meanrec", sim = sim, conv = conv[sim])
  
  # Bind parameter estimates
  par_all <- rbind(M_df, q_srv_df, srv_sel_df, meanrec_df, par_all)

  # Recruitment
  rec_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "Total_Rec") %>% 
    mutate(t = rec_total[70:(n_years-1),sim], sim = sim, conv = conv[sim],
           year = 70:(n_years-1))
  rec_all <- rbind(rec_df, rec_all)
  
  # Check SSB
  ssb_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "SSB") %>% 
    mutate(t = SSB[Fish_Start_yr[1]:(n_years-1), sim], sim = sim, conv = conv[sim],
           year = 70:(n_years-1))
  ssb_all <- rbind(ssb_df, ssb_all)
  
  # ssb_df %>% ggplot(aes(x = year, y = mle_val, ymin = lwr_95, ymax = upr_95)) +
  #   geom_line() +
  #   geom_line(aes(y = t), color = "red") +
  #   geom_ribbon(alpha = 0.3)

  # Check F
  f_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "Total_Fy") %>% 
    mutate(t = rowSums(fish_mort[Fish_Start_yr[1]:(n_years-1),,sim]), sim = sim, conv = conv[sim],
           year = 70:(n_years-1))
  f_all <- rbind(f_all, f_df)
  
  # Check total biomass
  t_biom <- extract_ADREP_vals(sd_rep = sd_rep, par = "Total_Biom") %>% 
    mutate(sim = sim, conv = conv[sim], year = 70:(n_years-1), t = biom_df$Biomass)
  biom_all <- rbind(t_biom, biom_all)
  
  # Check depletion rates
  # depletion_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "Depletion") %>% 
  #   mutate(t = (SSB[Fish_Start_yr:(n_years-1),sim]/SSB[Fish_Start_yr,sim]), 
  #          sim = sim, conv = conv[sim], year = 70:(n_years-1))
  # depletion_all <- rbind(depletion_df, depletion_all)
  # 
  # # Check survey mean age
  # srv_mean_ages <- extract_mean_age_vals(mod_rep = my_model, comp_name = "pred_srv_age_comps", 
  #                       bins = ages, comp_start_yr = Fish_Start_yr, sim = sim, 
  #                       n_fish_true_fleets = NULL) %>% mutate(conv = conv[sim])
  # srv_mu_age <- rbind(srv_mu_age, srv_mean_ages)
  # 
  # # Check fishery mean age
  # fish_mean_ages <- extract_mean_age_vals(mod_rep = my_model, comp_name = "pred_fish_age_comps", 
  #                                        bins = ages, comp_start_yr = Fish_Start_yr, sim = sim, 
  #                                        n_fish_true_fleets = 1) %>% mutate(conv = conv[sim])
  # fish_mu_age <- rbind(fish_mu_age, fish_mean_ages)

  print(paste("done w/  sim = ", sim))
  print(conv[sim])
  
}
  
  
# Summary Checks ----------------------------------------------------------

# Get percentiles
f_sum <- get_RE_precentiles(df = f_all %>% filter(conv == "Converged"), 
                     est_val_col = 1, true_val_col = 5, 
                     par_name = "Total Fishing Mortality", year)
  
ssb_sum <- get_RE_precentiles(df = ssb_all %>% filter(conv == "Converged"), 
                              est_val_col = 1, true_val_col = 5, 
                              par_name = "Spawning Stock Biomass", year)

rec_sum <- get_RE_precentiles(df = rec_all %>% filter(conv == "Converged"), 
                              est_val_col = 1, true_val_col = 5, 
                              par_name = "Total Recruitment", year)

biom_sum <- get_RE_precentiles(df = biom_all %>% filter(conv == "Converged"), 
                               est_val_col = 1, true_val_col = 8, 
                               par_name = "Total Biomass", year)

fish_mu_age_sum <- get_RE_precentiles(df = fish_mu_age %>% filter(conv == "Converged"), 
                                      est_val_col = 5, true_val_col = 6, 
                                      par_name = "Mean Predicted Fishery Age", year)

srv_mu_age_sum <- get_RE_precentiles(df = srv_mu_age %>% filter(conv == "Converged"), 
                                     est_val_col = 5, true_val_col = 6, 
                                     par_name = "Mean Predicted Survey Age", year)

depletion_sum <-  get_RE_precentiles(df = depletion_all %>% filter(conv == "Converged"), 
                                     est_val_col = 1, true_val_col = 5, 
                                     par_name = "Depletion (SSB / SSB0)", year)

all <- rbind(rec_sum, ssb_sum, f_sum, biom_sum) 

# Parameter estimates
par_df <- par_all %>% mutate(RE = (mle_val - t ) / t) 

# Quick summary stats
par_sum <- get_RE_precentiles(df = par_all %>% filter(conv == "Converged"), 
                              est_val_col = 2, true_val_col = 7, 
                              par_name = NULL, type)
  
(est_plot <- ggplot(all %>% filter(par_name != "Depletion (SSB / SSB0)"), aes(x = year, y = median)) +
  # geom_ribbon(aes(ymin = lwr_70, ymax = upr_70), alpha = 0.7, fill = "grey") +
  geom_ribbon(aes(ymin = lwr_80, ymax = upr_80), alpha = 0.5, fill = "grey4") +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.3, fill = "grey2") +
  # geom_ribbon(aes(ymin = lwr_100, ymax = upr_100), alpha = 0.4, fill = "grey") +
  geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.8, alpha = 0.85) +
  # geom_line( color = "white", size = 1,alpha = 1) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 1, alpha = 1) +
  facet_wrap(~par_name, scales = "free") +
  # coord_cartesian(ylim = c(-0.4, 0.4)) +
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
ggsci::scale_color_jco() +
ggsci::scale_fill_jco() +
  labs(x = "Relative Error", y = "Probability Density", linetype = "", color = "") +
theme_bw() + 
theme(strip.text = element_text(size = 13),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13, color = "black"),
      legend.position = "none", legend.text = element_text(size = 15)))

plot_grid(par_plot, est_plot, ncol = 1, align = "hv", axis = "bl",
          rel_heights = c(0.70, 1))

