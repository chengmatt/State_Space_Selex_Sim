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
spreadsheet_path <- here("input", "EBS_Pollock_Inputs.xlsx")
# spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")

# simulate data
simulate_data(fxn_path = fxn_path, 
              check_equil = FALSE,
              spreadsheet_path = spreadsheet_path, 
              rec_type = "BH",
              Start_F = c(0.01, 0.01), 
              Fish_Start_yr = c(1, 1), 
              Surv_Start_yr = c(1),  
              max_rel_F_M = c(1, 1), 
              desc_rel_F_M = c(0.05, 0.05), 
              F_type = c("Contrast", "Const_Inc"),
              yr_chng = c(25, 25), 
              yr_chng_end = c(30, 30),
              fish_likelihood = c("multinomial", "multinomial"),
              srv_likelihood = "multinomial",
              Input_Fish_N_Max = c(400, 400), 
              Input_Srv_N_Max = c(100),
              fish_CV = c(0.1, 0.1),
              srv_CV = c(0.1), 
              catch_CV = c(0.01, 0.01), 
              Input_N_Fish_Time = c("Constant", "Constant"), 
              Input_N_Fish_Fixed = c(100),
              Mort_Time = "Constant", 
              q_Mean_Fish = c(0.05, 0.05), 
              q_Mean_Surv = 0.01, 
              fish_selex = c("logistic", "logistic"), 
              srv_selex = c("logistic"), 
              # if switching to a single sex, be sure to change the nrow to the number of sexes,
              # and to make sure the selex parameters for the fleets align n_pars * n_sexes
              # e.g., (7, 0.8, 4, 0.3) for a logistic with two sexes, nrow = 2
              fish_pars = list(Fleet_1_L = matrix(data = c(5, 0.85), 
                                                  nrow = 1, byrow = TRUE),
                               Fleet_2_L = matrix(data = c(2, 0.85), 
                                                  nrow = 1, byrow = TRUE)), # fish fleet 2
              srv_pars = list(Fleet_3_SL = matrix(data = c(2, 0.85), 
                                                  nrow = 1, byrow = TRUE)), # survey fleet 1
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
all_harv_rates <- data.frame()
conv <- vector()
depletion_all <- data.frame()
ssb0_all <- data.frame()
fish_neff_all <- data.frame()
f40 <- vector()

# Specify years here
years <- Fish_Start_yr[1]:(n_years - 1)

start.time <- Sys.time()
compile_tmb(wd = here("src"), cpp = "EM.cpp")

for(sim in sim:n_sims){

  # Prepare inputs here
  input <- prepare_EM_input(years = years,
                            n_fleets = 2, 
                            catch_cv = c(0.01, 0.01),
                            F_Slx_Blocks_Input = matrix(c(rep(0)),
                                                        nrow = length(years),
                                                        ncol = 2), # fishery blocks
                            S_Slx_Blocks_Input = matrix(c(0), # selectivity blocks
                                                        nrow = length(years), 
                                                        ncol = 2),
                            # use_srv_comps = FALSE,
                            # use_srv_index = FALSE,
                            use_fish_index = FALSE,
                            # use_catch = FALSE,
                            Fish_Comp_Like_Model = c("multinomial" ),
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
                                         # "fixed_sel_re_fish",
                                         # "ln_N1Devs",
                                         # "ln_Fy",
                                         # "ln_RecDevs",
                                         # "ln_fish_selpars",
                                         # "ln_srv_selpars",
                                         # "ln_q_srv",
                                         # "ln_RecPars",
                                         # "ln_M"),
                                         # "ln_DM_Fish_Param"),
                            sim = sim)
  
  
  # input$parameters$ln_srv_selpars[] <- log(c(2, 3, 0.5, 0.4))
  # input$parameters$ln_fish_selpars[] <- log(c(3, 5, 0.8, 0.8))
  input$data$N1_Sex_Test <- log(matrix(N_at_age[1,,,sim], ncol = 2, nrow = 30))
  input$parameters$ln_N1Devs <- log(init_age_devs[,sim])
  # input$data$obs_fish_age_Input_N[] <- (Input_N_Fish[Fish_Start_yr[1]:(n_years-1),] + 
  #                                         (DM_Fish_Param * Input_N_Fish[Fish_Start_yr[1]:(n_years-1),])) /
  #   (Input_N_Fish[Fish_Start_yr[1]:(n_years-1),] + DM_Fish_Param)
  # input$parameters$ln_q_srv <- log(mean(q_Surv))
  # input$parameters$ln_q_fish <- log(mean(q_Fish))
  # 
  # input$parameters$ln_fish_selpars[] <- log(c(5, 0.85))
  # input$parameters$ln_srv_selpars[] <- log(c(2, 0.85))
  # input$parameters$ln_DM_Fish_Param[] <- log(2)
  # input$parameters$ln_DM_Srv_Param[] <- log(1)
  input$data$ssb0_dat <- ssb0
  # input$data$obs_fish_age_Input_N[] <- input$data$obs_fish_age_Input_N[] * 1
  # input$data$obs_fish_age_Input_N[] <- 1/(1+DM_Fish_Param) + Input_N_Fish[Fish_Start_yr[1]:(n_years-1),]*
                                      # (DM_Fish_Param/(1+DM_Fish_Param)) 
  # input$parameters$ln_SigmaRec <- 0.01
  # input$parameters$ln_SigmaRec <- log(input$parameters$ln_SigmaRec )
  # input$parameters$ln_Fy[] <- log(0.01)
  # input$parameters$ln_N1Devs[] <- 0.1
  # input$data$srv_cv <- 0.05
  # input$parameters$fixed_sel_re_fish[] <- log(2)

  # Run EM model here and get sdrep
  tryCatch(expr = model <- run_EM(data = input$data, parameters = input$parameters, 
                                  map = rlist::list.append(input$map), 
                                  n.newton = 3, 
                                  # random = "ln_RecDevs",
                                  silent = TRUE, getsdrep = TRUE), error = function(e){e})
  
  # Check model convergence
  convergence_status <- check_model_convergence(mle_optim = model$mle_optim,
                                                mod_rep = model$model_fxn,
                                                sd_rep = model$sd_rep,
                                                min_grad = 0.01)
  conv[sim] <- convergence_status$Convergence
  max_par[sim] <- convergence_status$Max_Grad_Par
  
  # if(conv[sim] == "Converged") {
  #   par(mfrow = c(1,2))
  #   year <- 1
  #   plot(model$model_fxn$rep$NAA[year,,1], type = "l", xlab ="Age",
  #        ylab = "Numbers")
  #   # lines(fish_mort[1+year-1,,sim] * Fish_selex_at_age[year,,,1,sim])
  #   lines(N_at_age[1+year-1,,1,sim], type = "l", col = "blue")
  #   plot(model$model_fxn$rep$NAA[year+1,,1], type = "l", xlab ="Age",
  #        ylab = "Numbers")
  #   # lines(fish_mort[1+year-1,,sim] * Fish_selex_at_age[year,,,1,sim])
  #   lines(N_at_age[1+year-1+1,,1,sim], type = "l", col = "blue")
  #   # for(i in 1:length(years)) {
  #   #   if(i == 1) plot(model$model_fxn$rep$F_Slx[1,,1,1], type = 'l')
  #   #   else lines(model$model_fxn$rep$F_Slx[i,,1,1])
  #   # }
  # }
  
  # Get parameter estimates
  M_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_M", trans = "log") %>%
    mutate(t = mean(Mort_at_age), type = "mortality", sim = sim, conv = conv[sim])
  q_srv_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_q_srv", trans = "log") %>%
    mutate(t = mean(q_Surv), type = "q_surv", sim = sim, conv = conv[sim])
  meanrec_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_RecPars", trans = "log") %>%
    mutate(t = c(r0), type = c("r0"), sim = sim, conv = conv[sim])
  # fish_sel_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_fish_selpars", trans = "log") %>%
  # mutate(t = c(5,0.85), type = c("f1", "f2"), sim = sim, conv = conv[sim])
  # srv_sel_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_srv_selpars", trans = "log") %>%
  #   mutate(t = c(2, 0.85), type = c("s1", "s2"), sim = sim, conv = conv[sim])
  fish_dm_theta_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_DM_Fish_Param", trans = "log") %>%
    mutate(t = c(DM_Fish_Param), type = c("DMF1"), sim = sim, conv = conv[sim])
  # srv_dm_theta_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_DM_Srv_Param", trans = "log") %>%
    # mutate(t = c(1), type = c("DMS1"), sim = sim, conv = conv[sim])
  # recsig_theta_df <- extract_parameter_vals(sd_rep = model$sd_rep, par = "ln_SigmaRec", trans = "log") %>%
    # mutate(t = c(sigma_rec), type = c("sigmarec"), sim = sim, conv = conv[sim])
  # Bind parameter estimates
  par_all <- rbind( fish_dm_theta_df, par_all)
                    # fish_sel_df, srv_sel_df)

  # Recruitment
  rec_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Total_Rec") %>% 
    mutate(t = rec_total[1:(n_years-1),sim], sim = sim, conv = conv[sim],
           year = 1:(n_years-1))
  rec_all <- rbind(rec_df, rec_all)
  
  # Check SSB
  ssb_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "SSB") %>% 
    mutate(t = SSB[Fish_Start_yr[1]:(n_years-1), sim], sim = sim, conv = conv[sim],
           year = 1:(n_years-1))
  ssb_all <- rbind(ssb_df, ssb_all)
  
  # Check F
  f_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Total_Fy") %>% 
    mutate(t = rowSums(fish_mort[Fish_Start_yr[1]:(n_years-1),,sim]), sim = sim, conv = conv[sim],
           year = 1:(n_years-1))
  f_all <- rbind(f_all, f_df)
  
  # Check depletion rates
  depletion_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Depletion") %>%
    mutate(t = (SSB[Fish_Start_yr[1]:(n_years-1),sim]/ssb0),
           sim = sim, conv = conv[sim], year = Fish_Start_yr[1]:(n_years-1))
  depletion_all <- rbind(depletion_df, depletion_all)
  
  # # Check survey mean age
  srv_mean_ages <- extract_mean_age_vals(mod_rep = model$model_fxn, comp_name = "pred_srv_age_comps",
                                         bins = ages, comp_start_yr = Fish_Start_yr[1], sim = sim,
                                         n_fish_true_fleets = NULL) %>% mutate(conv = conv[sim])
  
  srv_mu_age <- rbind(srv_mu_age, srv_mean_ages)
  
  # Check fishery mean age
  fish_mean_ages <- extract_mean_age_vals(mod_rep = model$model_fxn, comp_name = "pred_fish_age_comps",
                                          bins = ages, comp_start_yr = Fish_Start_yr[1], sim = sim,
                                          n_fish_true_fleets = 2) %>% mutate(conv = conv[sim])
  fish_mu_age <- rbind(fish_mu_age, fish_mean_ages)
  
  # Get SSB0
  ssb0_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "ssb0") %>% 
    mutate(t = ssb0, sim = sim, conv = conv[sim])
  ssb0_all <- rbind(ssb0_all, ssb0_df)
  
  # Get Fishery effective sample size
  # fish_neff_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Fish_Neff") %>%
  # mutate(t = Input_N_Fish[Fish_Start_yr[1]:(n_years-1),] 
  # , sim = sim, conv = conv[sim], year = 1:(n_years-1))
  
  # fish_neff_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Fish_Neff") %>%
  #   mutate(t =   (Input_N_Fish[Fish_Start_yr[1]:(n_years-1),] + (Input_N_Fish[Fish_Start_yr[1]:(n_years-1),] * DM_Fish_Param)) / (Input_N_Fish[Fish_Start_yr[1]:(n_years-1),] +  DM_Fish_Param), sim = sim, conv = conv[sim], year = 1:(n_years-1))
  
  # fish_neff_df <- extract_ADREP_vals(sd_rep = model$sd_rep, par = "Fish_Neff") %>%
  # mutate(t =  1/(1+DM_Fish_Param) + Input_N_Fish[Fish_Start_yr[1]:(n_years-1),]*(DM_Fish_Param/(1+DM_Fish_Param)), sim = sim, conv = conv[sim], year = 1:(n_years-1))
  # 
  # 
  # fish_neff_all <- rbind(fish_neff_all, fish_neff_df)
  # (Input_N_Fish[-31] / DM_Fish_Param^2)
  # 1/(1+DM_Fish_Param) + Input_N_Fish[Fish_Start_yr[1]:(n_years-1),]*(DM_Fish_Param/(1+DM_Fish_Param))
  # 1/(1+0.04496562) + (Input_N_Fish[Fish_Start_yr[1]:(n_years-1),] * 10) *(0.04496562/(1+0.04496562))
  
  est_M = exp(model$sd_rep$par.fixed[names(model$sd_rep$par.fixed) == "ln_M"])
  est_TermF = exp(model$sd_rep$par.fixed[names(model$sd_rep$par.fixed) == "ln_Fy"])[c(30, 60)]
  est_Selex = model$model_fxn$rep$F_Slx[30,,,]
  
  est_F40 <- get_Fx_refpt(ages = ages,
                 MortAA = rep(est_M, length = length(ages)), 
                 SelexAA = t(est_Selex), 
                 MatAA = mat_at_age[30,,1,sim], 
                 WAA = wt_at_age[30,,1,sim], 
                 Terminal_F = est_TermF, 
                 F_x = 0.4)
  
  f40[sim] <- est_F40
  
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

# fish_neff_sum <- get_RE_precentiles(df = fish_neff_all ,
#                                  est_val_col = 1, true_val_col = 5,
#                                  par_name = "Fish Neff", group_vars = "year")

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

dm <- par_df %>% filter(type == "q_surv")
rec0 <- par_df %>% filter(type == "DMF1")
plot(rec0$RE, dm$RE)

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

caa <- reshape2::melt(Catch_at_age[-31,,,,sim]/rowSums(Catch_at_age[-31,,,,sim])) %>% 
  mutate(type = "True_CAA",
         fleet = 1, sex = 1,
         Var1 = parse_number(paste(Var1)), Var2 = parse_number(paste(Var2)))
names(caa) <- c("year", "age", "props", "type", "fleet", "sex")



ggplot() +
  geom_col(pred_fish_comps, mapping = aes(x = age, y = props), alpha = 0.5) +
  geom_line(true_fish_comps, mapping = aes(x = age, y = props,
                                           col = "pred"), size = 1) +
  geom_line(caa, mapping = aes(x = age, y = props,
                               col = "true_caa"), size = 1, lty = 2) +
  scale_color_manual(name = c("pred", "true_caa"), 
                     values = c("black", "red")) +
  facet_wrap(~year) +
  theme_bw()

# Sum to 1 check
a = pred_fish_comps %>% group_by(year, fleet, sex) %>% 
  summarize(sum = sum(props))

b = true_fish_comps %>% group_by(year, fleet, sex) %>% 
  summarize(sum = sum(props))
