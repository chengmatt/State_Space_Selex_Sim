  
  # mean recruitment is positviely biased..
  # M may be causing pos bias in mean rec
  # rec devs are also causing positive bias in ssb
  # Set up ------------------------------------------------------------------
  
  library(here)
  library(tidyverse)
  library(TMB)
  
  source(here("R_scripts", "functions", "TMB_Utils.R"))
  
  compile_tmb(wd = here("src"), cpp = "EM.cpp")
  
  # Load in data ------------------------------------------------------------
  
  ssb_all <- data.frame()
  f_all <- data.frame()
  rec_all <- data.frame()
  biom_all <- data.frame()
  par_all <- data.frame()
  max_par <- vector()
  
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
                MatAA = array(MatAA, dim = c(32, 30, 1)),
                Sex_Ratio = as.vector(c(1)),  rec_model = 0, 
                fish_cv = fish_cv, srv_cv = srv_cv, catch_cv = catch_cv,
                Init_N = as.vector(N_at_age[(Fish_Start_yr),,,sim]),
                F_Slx_model = array(0, dim=c(1, 31, 1)), n_fish_comps = 1, n_srv_comps = 1,
                S_Slx_model = array(0, dim=c(1, 31, 1))
  )
  
  # Define parameter inits here
  parameters <- list(ln_SigmaRec = 0.6, ln_MeanRec = 2.75,
                     ln_M = log(0.1),  ln_a50_f = log(6), ln_k_f = log(0.8), 
                     ln_a50_s = log(4), ln_k_s = log(0.8),
                     ln_N1_Devs = log(rnorm(length(ages)-2,5, 1)),
                     ln_Fy = log(as.matrix(fish_mort[Fish_Start_yr:((n_years) -1),,sim])),
                     ln_q_fish = as.matrix(rep(log(0.08), n_fish_fleets)), 
                     ln_q_srv = as.matrix(rep(log(0.01), n_srv_fleets)),
                     ln_RecDevs = rec_devs[Fish_Start_yr:((n_years) -1),sim])
  
  
  map <- list(ln_SigmaRec = factor(NA))
              # ln_M = factor(NA),
              # ln_MeanF = factor(NA)
              # ln_MeanRec = factor(NA),
              # ln_a50_s = factor(NA), ln_k_s = factor(NA),
              # ln_a50_f = factor(NA), ln_k_f = factor(NA),
              # ln_q_fish = factor(NA),
              # ln_q_srv = factor(NA),
              # ln_Fy = factor(rep(NA, 31)))
              # ln_RecDevs = factor(rep(NA, 31)))
  
  compile_tmb(wd = here("src"), cpp = "EM.cpp")
  
  # Make ADFun
  my_model <- MakeADFun(data, parameters, map, DLL="EM", silent = T)
  mle_optim <- stats::nlminb(my_model$par, my_model$fn, my_model$gr, upper = 15, lower = -15,
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
  a50_fish_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_a50_f", log = TRUE) %>% 
    mutate(t = 6, type = "a50_fish", sim = sim, conv = conv[sim])
  a50_srv_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_a50_s", log = TRUE) %>% 
    mutate(t = 4, type = "a50_srv", sim = sim, conv = conv[sim])
  k_fish_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_k_f", log = TRUE) %>% 
    mutate(t = 0.8, type = "k_fish", sim = sim, conv = conv[sim])
  k_srv_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_k_s", log = TRUE) %>% 
    mutate(t = 0.8, type = "k_srv", sim = sim, conv = conv[sim])
  meanrec_df <- extract_parameter_vals(sd_rep = sd_rep, par = "ln_MeanRec", log = TRUE) %>% 
    mutate(t = exp(2.75), type = "meanrec", sim = sim, conv = conv[sim])
  
  # Bind parameter estimates
  par_all <- rbind(M_df, q_fish_df, q_srv_df, a50_fish_df, a50_srv_df,
                   k_fish_df, k_srv_df, meanrec_df, par_all)

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
    mutate(sim = sim, conv = conv[sim], year = 70:(n_years-1), t =biom_df$Biomass)
  biom_all <- rbind(t_biom, biom_all)

  print(paste("done w/  sim = ", sim))
  print(conv[sim])
}
  
# Quick checks
f_sum <- f_all %>% 
  filter(conv == "Converged") %>%
  mutate(type = "F",
         RE = (mle_val - t) /  t )%>% 
  group_by(year, type) %>% 
  summarize(median = median(RE), 
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975),
            lwr_80 = quantile(RE, 0.1),
            upr_80 = quantile(RE, 0.9))

ssb_sum <- ssb_all %>% 
  filter(conv == "Converged") %>%
  mutate(type = "ssb", RE = (mle_val - t) /  t) %>% 
  group_by(year, type) %>% 
  summarize(median = median(RE), 
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975),
            lwr_80 = quantile(RE, 0.1),
            upr_80 = quantile(RE, 0.9))

rec_sum <- rec_all %>% 
  filter(conv == "Converged") %>%
  mutate(type = "rec", RE = (mle_val - t) /  t) %>% 
  group_by(year, type) %>% 
  summarize(median = median(RE), 
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975),
            lwr_80 = quantile(RE, 0.1),
            upr_80 = quantile(RE, 0.9))  

biom_sum <- biom_all %>% 
  filter(conv == "Converged") %>%
  mutate(type = "biom", RE = (mle_val - t) /  t) %>% 
  group_by(year, type) %>% 
  summarize(median = median(RE), 
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975),
            lwr_80 = quantile(RE, 0.1),
            upr_80 = quantile(RE, 0.9)) 
  
all <- rbind(rec_sum, ssb_sum, f_sum, biom_sum) 

ggplot(all,aes(x = year, y = median)) +
  geom_ribbon(aes(ymin = lwr_80, ymax = upr_80), alpha = 0.6, fill = "grey4") +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.4, fill = "grey4") +
  geom_line( color = "white", size = 1,alpha = 1) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 3.8, stroke = 1, alpha = 1) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 1, alpha = 0.85) +
  facet_wrap(~type, scales = "free") +
  theme_bw()

# Parameter estimates
par_all %>% 
  mutate(RE = (mle_val - t ) / t) %>% 
  ggplot(aes(x = RE, fill = type)) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = 0), lty = 2, size = 1, col = "blue") +
  facet_wrap(~type, scales = "free") +
  ggthemes::scale_fill_colorblind() +
  theme_bw() +
  theme(legend.position = "none")
