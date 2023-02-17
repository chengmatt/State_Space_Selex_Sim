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
      spreadsheet_path <- here("input", "EBS_Pollock_Inputs.xlsx")
      # spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")
      
      # simulate data
      simulate_data(fxn_path = fxn_path, 
                    check_equil = FALSE,
                    spreadsheet_path = spreadsheet_path, 
                    rec_type = "BH",
                    Start_F = c(0.01), 
                    Fish_Start_yr = c(70), 
                    Surv_Start_yr = c(70), 
                    max_rel_F_M = c(1), 
                    desc_rel_F_M = c(0.01), 
                    F_type = c("Contrast"),
                    yr_chng = 85, 
                    yr_chng_end = 100,
                    fish_Neff_max = c(200), 
                    srv_Neff_max = c(200),
                    fish_CV = c(0.1, 0.1),
                    srv_CV = c(0.1), 
                    catch_CV = c(0, 0), 
                    Neff_Fish_Time = "Constant", 
                    fixed_Neff = c(30),
                    Mort_Time = "Constant", 
                    q_Mean_Fish = c(0.05), 
                    q_Mean_Surv = 0.01, 
                    fish_selex = c("logistic"), 
                    srv_selex = c("logistic"), 
                    # if switching to a single sex, be sure to change the nrow to the number of sexes,
                    # and to make sure the selex parameters for the fleets align n_pars * n_sexes
                    # e.g., (7, 0.8, 4, 0.3) for a logistic with two sexes, nrow = 2
                    fish_pars = list(Fleet_1_L = matrix(data = c(4, 0.8 ), 
                                                        nrow = 1, byrow = TRUE)), # fish fleet 2
                    srv_pars = list(Fleet_3_SL = matrix(data = c(4,0.8 ), 
                                                        nrow = 1, byrow = TRUE)), # survey fleet 1
                    f_ratio = 0.5, m_ratio = 0.5)
      
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
      all_harv_rates <- data.frame()
      conv <- vector()
      depletion_all <- data.frame()
      ssb0_all <- data.frame()
      
      # Specify years here
      years <- Fish_Start_yr[1]:(n_years - 1)
      
      start.time <- Sys.time()
      
    for(sim in 1:n_sims){
    
    # Prepare inputs here
    input <- prepare_EM_input(years = years,
                              n_fleets = 1, 
                              catch_cv = c(0.01),
                              F_Slx_Blocks_Input = matrix(c(rep(0)),
                                                          nrow = length(years),
                                                          ncol = 1), # fishery blocks
                              S_Slx_Blocks_Input = matrix(c(0), # selectivity blocks
                                                          nrow = length(years), 
                                                          ncol = 1),
                              use_fish_index = FALSE,
                              rec_model = "BH", 
                              F_Slx_Model_Input = c("logistic"),
                              S_Slx_Model_Input = c("logistic"), 
                              time_selex = "None",
                              n_time_selex_pars = NULL,
                              fix_pars = c("ln_SigmaRec", "logit_q_fish", "ln_h"),
                              sim = sim)
    
    # input$parameters$fixed_sel_re_fish[] <- c(0.35, 0.3)
    
    # input$parameters$ln_M <- log(0.125)
    
    compile_tmb(wd = here("src"), cpp = "EM.cpp")
    
    # Run EM model here and get sdrep
    tryCatch(expr = model <- run_EM(data = input$data, parameters = input$parameters, 
                                    map = input$map, n.newton = 1, 
                                    # random = "ln_fish_selpars_re",
                                    silent = T, getsdrep = TRUE), error = function(e){e})
    
    # plot(model$model_fxn$rep$SBPR_SSB0, type = "l")
    # lines(SPR_SSB0, col = "red")
    
    # Checking fixed effects parameter length
    # sum(names(model$sd_rep$par.fixed) == "ln_fish_selpars_re")
    # model$sd_rep
    
    # f_pars <- model$sd_rep$par.fixed[names(model$sd_rep$par.fixed) == "ln_fish_selpars"]
    # f_repars <- model$sd_rep$par.random[names(model$sd_rep$par.random) == "ln_fish_selpars_re"]
    # 
    # sel_res <- model$model_fxn$rep$ln_fish_selpars_re
    # plot(sel_res[,,1,1], type = "l")
    # 
    # # Checking fixed effects
    # 
    # for(i in 1:31) {
    #   if(i == 1)  plot(model$model_fxn$rep$F_Slx[i,,,1], type = "l", ylim = c(0,1))
    #   else lines(model$model_fxn$rep$F_Slx[i,,,1])
    # }
    # 
    # 
    # lines(model$model_fxn$rep$F_Slx[31,,1,1], lwd= 3, col = "purple")
    # lines(Fish_selex_at_age[1,,1,1,1], col = "red", lwd= 3)
    # lines(Fish_selex_at_age[1,,2,1,1], col = "blue", lwd= 3)
    
    # Matrix::image(model$model_fxn$env$spHess(random=TRUE))
    
    # Check model convergence
    convergence_status <- check_model_convergence(mle_optim = model$mle_optim, 
                                                  mod_rep = model$model_fxn,
                                                  sd_rep = model$sd_rep, 
                                                  min_grad = 0.001)
    conv[sim] <- convergence_status$Convergence
    max_par[sim] <- convergence_status$Max_Grad_Par
    
    # Get parameter estimates
    M_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_M", trans = "log") %>% 
      mutate(t = mean(Mort_at_age), type = "mortality", sim = sim, conv = conv[sim])
    q_srv_df <- extract_parameter_vals(sd_rep = model$sd_rep,   par = "logit_q_srv", trans = "logit",
                                       logit_bounds = c(0, 1)) %>% 
      mutate(t = mean(q_Surv), type = "q_surv", sim = sim, conv = conv[sim])
    meanrec_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_RecPars", trans = "log") %>% 
      mutate(t = r0, type = "r0/meanrec", sim = sim, conv = conv[sim])
    # fish_sel_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_fish_selpars", trans = "log") %>% 
    #   mutate(t = c(7, 4, 0.8, 0.8), type = c("f1", "f2", "f1d", "f2d"), sim = sim, conv = conv[sim])
    ssb0_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "ssb0") %>%
    mutate(t = ssb0, type = "ssb0", sim = sim, conv = conv[sim])
    ssb0_all <- rbind(ssb0_df, ssb0_all)
    # Bind parameter estimates
    par_all <- rbind(M_df, q_srv_df, meanrec_df, par_all)
    

    if(sim > 3) {
      rec_stuff <- par_all %>% filter(type == "r0/meanrec",  conv == "Converged")
      hist((rec_stuff$mle_val - rec_stuff$t) / rec_stuff$mle_val, xlim = c(-1, 0.5))
      abline(v = 0, col = "red", lwd = 5)
      median_rec <- median((rec_stuff$mle_val - rec_stuff$t) / rec_stuff$mle_val)
      abline(v = median_rec, lwd = 5)
    }
    
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
    
    # Check F
    f_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Total_Fy") %>% 
      mutate(t = fish_mort[Fish_Start_yr[1]:(n_years-1),,sim], sim = sim, conv = conv[sim],
             year = 70:(n_years-1))
    f_all <- rbind(f_all, f_df)
    
    # Check depletion rates
    # depletion_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Depletion") %>%
    # mutate(t = (SSB[Fish_Start_yr[1]:(n_years-1),sim]/ssb0),
    # sim = sim, conv = conv[sim], year = Fish_Start_yr[1]:(n_years-1))
    # depletion_all <- rbind(depletion_df, depletion_all)
    
    # # Check survey mean age
    srv_mean_ages <- extract_mean_age_vals(mod_rep = model$model_fxn, comp_name = "pred_srv_age_comps",
                                           bins = ages, comp_start_yr = Fish_Start_yr[1], sim = sim,
                                           n_fish_true_fleets = NULL) %>% mutate(conv = conv[sim])
    
    srv_mu_age <- rbind(srv_mu_age, srv_mean_ages)
    
    # Check fishery mean age
    fish_mean_ages <- extract_mean_age_vals(mod_rep = model$model_fxn, comp_name = "pred_fish_age_comps",
                                            bins = ages, comp_start_yr = Fish_Start_yr[1], sim = sim,
                                            n_fish_true_fleets = 1) %>% mutate(conv = conv[sim])
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

fish_mu_age_sum <- get_RE_precentiles(df = fish_mu_age %>% filter(conv == "Converged"),
                                      est_val_col = 4, true_val_col = 6,
                                      par_name = "Mean Predicted Fishery Age", 
                                      group_vars = c("year","fleet", "sex"))

srv_mu_age_sum <- get_RE_precentiles(df = srv_mu_age %>% filter(conv == "Converged"),
                                     est_val_col = 4, true_val_col = 6,
                                     par_name = "Mean Predicted Survey Age", 
                                     group_vars = c("year","fleet", "sex"))

# depletion_sum <-  get_RE_precentiles(df = depletion_all %>% filter(conv == "Converged"), 
#                                      est_val_col = 1, true_val_col = 5, 
#                                      par_name = "Depletion (SSB / SSB0)", group_vars = "year")

all <- rbind(rec_sum, ssb_sum, f_sum, srv_mu_age_sum, fish_mu_age_sum) 
  

# Get relative error time series
(est_plot <- plot_RE_ts(data = all, x = year, y = median, 
           lwr_1 = lwr_80, upr_1 = upr_80,
           lwr_2 = lwr_95, upr_2 = upr_95, 
           facet_name = par_name))


# Parameter estimates
par_df <- par_all %>% mutate(RE = (mle_val - t ) / t) %>% 
  filter(conv == "Converged")

# Quick summary stats
par_sum <- get_RE_precentiles(df = par_all %>% filter(conv == "Converged"), 
                              est_val_col = 2, true_val_col = 7, 
                              par_name = "", group_vars = "type")

(par_plot <- ggplot(par_df, aes(x = RE, fill = type)) +
geom_density(alpha = 0.2) +
facet_wrap(~type, scales = "free", nrow = 2) +
geom_errorbarh(inherit.aes = FALSE, data = par_sum, 
               aes(xmin = lwr_95, xmax = upr_95, y = 0, color = type),
               height = 0, size = 2.5, alpha = 0.55, linetype = 1) +
geom_point(inherit.aes = FALSE, data = par_sum, 
           aes(x= median, y = 0, color = type), size = 5, alpha = 0.95) +
  geom_vline(aes(xintercept = 0), linetype = 2,
             size = 0.85, col = "black", alpha = 1) +
coord_cartesian(xlim = c(-1, 0.5)) +
ggsci::scale_color_jco() +
ggsci::scale_fill_jco() +
  labs(x = "Relative Error", y = "Probability Density", linetype = "", color = "") +
theme_bw() + 
theme(strip.text = element_text(size = 13),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13, color = "black"),
      legend.position = "none", legend.text = element_text(size = 15)))

# ssb0 par sums
ssb0_par_sum <- get_RE_precentiles(df = ssb0_all %>% filter(conv == "Converged"), 
                              est_val_col = 1, true_val_col = 5, 
                              par_name = "", group_vars = "type")

ssb0_all <- ssb0_all %>% mutate(RE = (mle_val - t ) / t) %>% 
  filter(conv == "Converged")

ggplot(ssb0_all, aes(x = RE, fill = type)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~type, scales = "free", nrow = 2) +
  geom_errorbarh(inherit.aes = FALSE, data = ssb0_par_sum, 
                 aes(xmin = lwr_95, xmax = upr_95, y = 0, color = type),
                 height = 0, size = 2.5, alpha = 0.55, linetype = 1) +
  geom_point(inherit.aes = FALSE, data = ssb0_par_sum, 
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
        legend.position = "none", legend.text = element_text(size = 15))


plot_grid(par_plot, est_plot, ncol = 1, align = "hv", axis = "bl",
          rel_heights = c(0.70, 1))

plot(rowSums(fish_mort[Fish_Start_yr[1]:(n_years-1),,sim]), type = "l")
lines(exp(model$sd_rep$par.fixed[names(model$sd_rep$par.fixed) == "ln_Fy"]),
      col = "red")

f_all %>% 
  ggplot(aes(x = year, y = mle_val, group = sim)) +
  geom_line(color = "grey", size = 1) +
  geom_line(aes(y = t), color = "red") +
  theme_bw()

