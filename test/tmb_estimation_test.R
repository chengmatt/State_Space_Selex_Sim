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
                n_years = 101, Start_F = c(0.01), 
                Fish_Start_yr = c(70), Surv_Start_yr = c(70), 
                max_rel_F_M = c(1.5), desc_rel_F_M = c(0.15), 
                F_type = c("Contrast"), yr_chng = c(86), 
                fish_Neff = c(150), srv_Neff = c(150), fish_CV = c(0.1),
                srv_CV = c(0.1), Neff_Fish_Time = "F_Vary", fixed_Neff = 30,
                Mort_Time = "Constant", q_Mean_Fish = 0.08, q_Mean_Surv = 0.01, 
                Rec_Dev_Type = "iid", rho_rec = NA, 
                fish_selex = c("double_logistic"), srv_selex = c("logistic"), 
                fish_pars = list(Fleet_1_L = matrix(data = c(0.3, 0.5, 0.5, 16 ), nrow = 1, byrow = TRUE)),
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
  years <- Fish_Start_yr:(n_years - 1)
  n_sexes <- n_sex
  n_fleets <- n_fish_fleets
  n_fish_indices <- 1
  n_srv_indices <- 1
  
  for(sim in 1:n_sims){
  
  # Calculate catches 
  obs_catches <- melt(Catch_at_age) %>%
    drop_na() %>% # drop nas in the last year
    rename(Year = Var1, Age = Var2, # Rename varialbes
           Fleet = Var3, Sex = Var4, Sim = Var5, Catch = value) %>%
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim))) %>%  # Parse number for year and simulation
    filter(Year >= Fish_Start_yr[1],
           Sim == sim)
  
  obs_catches <- matrix(with(obs_catches, tapply(Catch, list(Year), FUN = sum)))

  # Observed fishery age comps
  obs_fish_age_comps <- t(apply(Fish_Age_Comps[Fish_Start_yr:(n_years - 1),,,,sim], MARGIN = 1, 
                                FUN=function(x) { x/sum(x) }))
  obs_fish_age_Neff <- matrix(fish_Neff[Fish_Start_yr:(n_years - 1),])
  
  # Observed survey age comps
  obs_srv_age_comps <- t(apply(Survey_Age_Comps[Fish_Start_yr:(n_years - 1),,,,sim], MARGIN = 1, 
                               FUN=function(x) { x/sum(x) }))
  obs_srv_age_Neff <- matrix(srv_Neff[Fish_Start_yr:(n_years - 1),])
  
  # Observed fishery indices
  obs_fish_indices <-as.matrix( Fishery_Index_Agg[Fish_Start_yr:(n_years - 1),,sim])
  obs_srv_indices <- as.matrix(Survey_Index_Agg[Fish_Start_yr:(n_years - 1),,sim])
  
  WAA <- wt_at_age[Fish_Start_yr:(n_years),,1,sim]
  MatAA <- mat_at_age[Fish_Start_yr:(n_years),,1,sim]
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
  
  fish_cv <- 0.1
  srv_cv <- 0.1
  catch_cv <- 0.05
  
  F_Slx_Blocks <- matrix(c(0), nrow = length(years), ncol = n_fleets)
  S_Slx_Blocks <- matrix(c(0), nrow = length(years), ncol = n_srv_fleets)
  
  # TMB Section -------------------------------------------------------------
  
  
  # Fill in list for data
  data <- list( ages = ages, years = years,
                n_sexes = n_sexes, n_fleets = 1,
                n_fish_indices = n_fish_indices, n_srv_indices = n_srv_indices,
                obs_catches = obs_catches, 
                obs_fish_age_comps = array(obs_fish_age_comps, dim = c(31, 30, 1, 1)),
                obs_fish_age_Neff = obs_fish_age_Neff, 
                obs_srv_age_comps = array(obs_srv_age_comps, dim = c(31, 30, 1, 1)),
                obs_srv_age_Neff = obs_srv_age_Neff, obs_fish_indices =  obs_fish_indices,
                obs_srv_indices = obs_srv_indices, WAA = array(WAA, dim = c(32, 30, 1)), 
                MatAA = array(MatAA, dim = c(32, 30, 1)), F_Slx_Blocks = F_Slx_Blocks,
                S_Slx_Blocks = S_Slx_Blocks,
                Sex_Ratio = as.vector(c(1)),  rec_model = 0, 
                fish_cv = fish_cv, srv_cv = srv_cv, catch_cv = catch_cv,
                Init_N = as.vector(N_at_age[(Fish_Start_yr),,,sim]),
                F_Slx_model = as.vector(2), n_fish_comps = 1, n_srv_comps = 1,
                S_Slx_model = as.vector(0)
  )
  
  # Define parameter inits here
  parameters <- list(ln_SigmaRec = 0.6, ln_MeanRec = 2.75,
                     ln_M = log(0.1),  
                     ln_fish_selpars = log(array(c(0.3, 0.5, 0.5, 16 ), dim = c(1, 1, 1, 4))),
                     ln_srv_selpars = array(5, dim = c(1, 1, 1, 2)),
                     ln_N1_Devs = log(rnorm(length(ages)-2,5, 1)),
                     ln_Fy = log(as.matrix(fish_mort[Fish_Start_yr:((n_years) -1),,sim])),
                     ln_q_fish = as.matrix(rep(log(0.08), n_fish_fleets)), 
                     ln_q_srv = as.matrix(rep(log(0.01), n_srv_fleets)),
                     ln_RecDevs = rec_devs[Fish_Start_yr:((n_years) -1),sim])
  
  
  map <- list(ln_SigmaRec = factor(NA))
              # ln_fish_selpars = factor(rep(NA, 4)))
              # ln_M = factor(NA),
              # ln_MeanF = factor(NA)
              # ln_MeanRec = factor(NA),
              # ln_q_fish = factor(NA),
              # ln_q_srv = factor(NA),
              # ln_Fy = factor(rep(NA, 31)))
              # ln_RecDevs = factor(rep(NA, 31)))
  
  compile_tmb(wd = here("src"), cpp = "EM.cpp")
  
  # Make ADFun
  my_model <- MakeADFun(data, parameters, map, DLL="EM", silent = T)
  mle_optim <- stats::nlminb(my_model$par, my_model$fn, my_model$gr,
                             control = list(iter.max = 1e5, eval.max = 1e5))
  
  # Additional newton steps to take
  add_newton(n.newton = 3, ad_model = my_model, mle_optim = mle_optim)
  
  my_model$rep <- my_model$report(my_model$env$last.par.best)
  sd_rep <- TMB::sdreport(my_model)
  
  # Check model convergence
  convergence_status <- check_model_convergence(mle_optim = mle_optim, mod_rep = my_model,
                                                sd_rep = sd_rep, min_grad = 0.01)
  conv[sim] <- convergence_status$Convergence
  max_par[sim] <- convergence_status$Max_Grad_Par
  
  # Get parameter estimates
  M_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_M", log = TRUE) %>% 
    mutate(t = mean(Mort_at_age), type = "mortality", sim = sim, conv = conv[sim])
  q_fish_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_q_fish", log = TRUE) %>% 
    mutate(t = mean(q_Fish), type = "q_fish", sim = sim, conv = conv[sim])
  q_srv_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_q_srv", log = TRUE) %>% 
    mutate(t = mean(q_Surv), type = "q_surv", sim = sim, conv = conv[sim])
  # fish_sel_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_fish_selpars", log = TRUE) %>%
    # mutate(t = c(4, 6), type = c("delta_fish", "amax_fish"), sim = sim, conv = conv[sim])
  srv_sel_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_srv_selpars", log = TRUE) %>% 
    mutate(t = c(4, 0.8), type = c("a50_srv", "k_srv"), sim = sim, conv = conv[sim])
  meanrec_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_MeanRec", log = TRUE) %>% 
    mutate(t = exp(2.75), type = "meanrec", sim = sim, conv = conv[sim])
  
  # Bind parameter estimates
  par_all <- rbind(M_df, q_fish_df, q_srv_df, srv_sel_df, meanrec_df, par_all)

  # Recruitment
  rec_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "Total_Rec") %>% 
    mutate(t = rec_total[70:(n_years-1),sim], sim = sim, conv = conv[sim],
           year = 70:(n_years-1))
  rec_all <- rbind(rec_df, rec_all)
  
  # Check SSB
  ssb_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "SSB") %>% 
    mutate(t = SSB[Fish_Start_yr:(n_years-1), sim], sim = sim, conv = conv[sim],
           year = 70:(n_years-1))
  ssb_all <- rbind(ssb_df, ssb_all)

  # Check F
  f_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "Total_Fy") %>% 
    mutate(t = fish_mort[Fish_Start_yr:(n_years-1),,sim], sim = sim, conv = conv[sim],
           year = 70:(n_years-1))
  f_all <- rbind(f_all, f_df)
  
  # Check total biomass
  t_biom <- extract_ADREP_vals(sd_rep = sd_rep, par = "Total_Biom") %>% 
    mutate(sim = sim, conv = conv[sim], year = 70:(n_years-1), t = biom_df$Biomass)
  biom_all <- rbind(t_biom, biom_all)
  
  # Check depletion rates
  depletion_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "Depletion") %>% 
    mutate(t = (SSB[Fish_Start_yr:(n_years-1),sim]/SSB[Fish_Start_yr,sim]), 
           sim = sim, conv = conv[sim], year = 70:(n_years-1))
  depletion_all <- rbind(depletion_df, depletion_all)
  
  # Check survey mean age
  srv_mean_ages <- extract_mean_age_vals(mod_rep = my_model, comp_name = "pred_srv_age_comps", 
                        bins = ages, comp_start_yr = Fish_Start_yr, sim = sim, 
                        n_fish_true_fleets = NULL) %>% mutate(conv = conv[sim])
  srv_mu_age <- rbind(srv_mu_age, srv_mean_ages)
  
  # Check fishery mean age
  fish_mean_ages <- extract_mean_age_vals(mod_rep = my_model, comp_name = "pred_fish_age_comps", 
                                         bins = ages, comp_start_yr = Fish_Start_yr, sim = sim, 
                                         n_fish_true_fleets = 1) %>% mutate(conv = conv[sim])
  fish_mu_age <- rbind(fish_mu_age, fish_mean_ages)

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

all <- rbind(rec_sum, ssb_sum, f_sum, biom_sum, depletion_sum, srv_mu_age_sum, fish_mu_age_sum) 

# Parameter estimates
par_df <- par_all %>% mutate(RE = (mle_val - t ) / t) 

# Quick summary stats
par_sum <- get_RE_precentiles(df = par_all %>% filter(conv == "Converged"), 
                              est_val_col = 2, true_val_col = 7, 
                              par_name = NULL, type)
  
(est_plot <- ggplot(all %>% filter(par_name != "Depletion (SSB / SSB0)"), aes(x = year, y = median)) +
  # geom_ribbon(aes(ymin = lwr_75, ymax = upr_75), alpha = 0.7, fill = "grey") +
  geom_ribbon(aes(ymin = lwr_80, ymax = upr_80), alpha = 0.5, fill = "grey4") +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.3, fill = "grey2") +
  # geom_ribbon(aes(ymin = lwr_100, ymax = upr_100), alpha = 0.4, fill = "grey") +
  geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.8, alpha = 0.85) +
  # geom_line( color = "white", size = 1,alpha = 1) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 1, alpha = 1) +
  facet_wrap(~par_name, scales = "free") +
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
geom_errorbarh(inherit.aes = FALSE, data = par_sum, 
                 aes(xmin = lwr_75, xmax = upr_75, y = 0, color = type),
                 height = 0, size = 2.5, alpha = 1, linetype = 1) +
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
          rel_heights = c(0.75, 1))

