  
  # mean recruitment is positviely biased..
  # M may be causing pos bias in mean rec
  # rec devs are also causing positive bias in ssb
  # Set up ------------------------------------------------------------------
  
  library(here)
  library(tidyverse)
  library(TMB)
  
  # Need to set wd in source folder
  setwd(here("src"))
  
  compile("EM.cpp") # Compile .cpp file
  dyn.load(dynlib("EM")) # Load in .cpp
  
  # Load in data ------------------------------------------------------------
  
  ssb_all <- data.frame()
  f_all <- data.frame()
  rec_all <- data.frame()
  biom_all <- data.frame()
  M_all <- vector()
  q_fish <- vector()
  q_srv <- vector()
  a50_f<- vector()
  k_f<- vector()
  a50_s<- vector()
  k_s<- vector()
  meanrec <- vector()
  conv <- vector()
  
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
                fish_cv = fish_cv, srv_cv = srv_cv,
                Init_N = as.vector(N_at_age[(Fish_Start_yr),,,sim]),
                F_Slx_model = array(0, dim=c(1, 31, 1)), n_fish_comps = 1, n_srv_comps = 1,
                S_Slx_model = array(0, dim=c(1, 31, 1))
  )
  
  # Define parameter inits here
  parameters <- list(ln_SigmaRec = 0.8, ln_MeanRec = 2.75,
                     ln_M = log(0.1),  ln_a50_f = log(6), ln_k_f = log(0.8), 
                     ln_a50_s = log(4), ln_k_s = log(0.8),
                     ln_N1_Devs = log(rnorm(length(ages)-2,5, 1)),
                     # ln_MeanF = log(0.1),
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
  
  compile("EM.cpp") # Compile .cpp file
  dyn.load(dynlib("EM")) # Load in .cpp
  
  # Make ADFun
  my_model <- MakeADFun(data, parameters, map, DLL="EM", silent = T)
  mle_optim <- stats::nlminb(my_model$par, my_model$fn, my_model$gr, upper = 5, lower = -5,
                             control = list(iter.max = 1e5, eval.max = 1e5))
  
  # Additional newton steps to take
  n.newton <- 5
  try_improve <- tryCatch(expr =
                            for(i in 1:n.newton) {
                              g = as.numeric(my_model$gr(mle_optim$par))
                              h = optimHess(mle_optim$par, fn = my_model$fn, 
                                            gr = my_model$gr)
                              mle_optim$par = mle_optim$par - solve(h,g)
                              mle_optim$objective = my_model$fn(mle_optim$par)
                            }, error = function(e){e})
  
  my_model$rep <- my_model$report(my_model$env$last.par.best)
  sd_rep <- TMB::sdreport(my_model)
  
  if(mle_optim$convergence == 0 &
     sd_rep$pdHess == TRUE &
     max(abs(sd_rep$gradient.fixed)) < 0.001) conv[sim] = "Converged"
  else conv[sim] = "Not Converged"
  
  # Get naturla mortality
  # M_all[sim] <- (exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_M"]) - 0.1 )/ 0.1
  # q_fish[sim] <- (exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_q_fish"]) - 0.08) / 0.08
  # q_srv[sim] <-  (exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_q_srv"]) - 0.01) / 0.01
  # a50_f[sim] <- (exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_a50_f"]) - 6) / 6
  # a50_s[sim] <- (exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_a50_s"]) - 4) / 4
  # k_f[sim] <- (exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_k_f"]) - 0.8 )/ 0.8
  # k_s[sim] <- (exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_k_s"]) - 0.8)/0.8
  # meanrec[sim] <- (sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_MeanRec"] - 2.75)/2.75
  # 
  # Recruitment
  my_model$report()$rec_nLL
  rec = sd_rep$value[str_detect(names(sd_rep$value), "Total_Rec")]
  sd = sd_rep$sd[str_detect(names(sd_rep$value), "Total_Rec")]
  
  # Years 69 - 99 because we don;t have info for the last year projected SSB yet...
  # Recruitment enters after the fact
  rec_df <- data.frame(year = 70:100, rec = rec, sd = sd,
                       upr = rec + (1.96 * sd), downr = rec - (1.96 * sd),
                       t = rec_total[70:(n_years-1),sim], sim = sim, conv = conv[sim])
  
  
  # ggplot(rec_df, aes(x = year, y = rec, ymin = downr, ymax = upr)) +
  #   geom_line() +
  #   geom_line(aes(y = t), col = "red") +
  #   geom_ribbon(alpha = 0.3)
  
  rec_all <- rbind(rec_df, rec_all)
  
  Init_N_at_age <- N_at_age[(Fish_Start_yr),,,sim]
  plot(my_model$rep$NAA[1,,1], type = "l") # Model
  lines(Init_N_at_age, col = "red")
  
  # Check SSB
  ssb_df <- data.frame(year = 70:100, ssb = sd_rep$value[1:31], 
                       upr = sd_rep$value[1:31] + (1.96 * sd_rep$sd[1:31] ),
                       downr = sd_rep$value[1:31] - (1.96 * sd_rep$sd[1:31] ),
                       t = SSB[Fish_Start_yr:(n_years-1), sim], sim = sim, conv = conv[sim])
  
  ssb_all <- rbind(ssb_df, ssb_all)
  
  # ggplot(ssb_df, aes(x = year, y = ssb, ymin = downr, ymax = upr)) +
  #   geom_line() +
  #   geom_line(aes(y = t), col = "red") +
  #   geom_ribbon(alpha = 0.3)

  # Check F
  f_df <- data.frame(year = 70:100, f = sd_rep$value[33:63], 
                     upr = sd_rep$value[33:63] + (1.96 * sd_rep$sd[33:63] ),
                     downr = sd_rep$value[33:63] - (1.96 * sd_rep$sd[33:63]),
                     t = fish_mort[Fish_Start_yr:(n_years-1),,sim], sim = sim, conv = conv[sim])
  
  f_all <- rbind(f_all, f_df)
  
  t_biom <- data.frame(year = 70:100, f = sd_rep$value[str_detect(names(sd_rep$value), "Total_Biom")], 
                       t = biom_df$Biomass, sim = sim, conv = conv[sim])
  
  biom_all <- rbind(t_biom, biom_all)

  # ggplot(f_df, aes(x = year, y = f, ymin = downr, ymax = upr)) +
  #   geom_line() +
  #   geom_line(aes(y = t), col = "red") +
  #   geom_ribbon(alpha = 0.3)
  
  print(paste("done w/  sim = ", sim))
  print(conv[sim])
  }
  
  # Mean Rec, M are good
  # When rec devs get turned on, total biomass and ssb reuslt in negative bias?
  # When init rec devs get turned on recruitment results in postivie bias,
  # and biomass is still negative bias at terminal year (although less of this is in ssb)
  # overall total biomass has a small degree of postivie bias
  # Whensurvey selex gets turned on, everything else is similar although recruitmetn
  # increases in bias by a bit
  # When survey q gets turned on - rec changes a bit but in general, same thing
  # When fishery q gets turned on - rec changes a bit but in general, same things
  # When fishery selex gets turned on - rec reduces bias, biom increases bias, but in general, same thing
  # SSB and biom neg bias increases
  # When fishing mortality gets turned on, everything looks the same, but the negative bias
  # increases in total biomass, and the variability increases
  

# Quick checks
f_sum <- f_all %>% 
  # filter(conv == "Converged") %>%
  mutate(type = "F",
         RE = (f - t) /  t )%>% 
  group_by(year, type) %>% 
  summarize(median = median(RE), 
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975),
            lwr_80 = quantile(RE, 0.1),
            upr_80 = quantile(RE, 0.9))


ssb_sum <- ssb_all %>% 
  # filter(conv == "Converged") %>%
  mutate(type = "ssb", RE = (ssb - t) /  t) %>% 
  group_by(year, type) %>% 
  summarize(median = median(RE), 
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975),
            lwr_80 = quantile(RE, 0.1),
            upr_80 = quantile(RE, 0.9))

rec_sum <- rec_all %>% 
  # filter(conv == "Converged") %>%
  mutate(type = "rec", RE = (rec - t) /  t) %>% 
  group_by(year, type) %>% 
  summarize(median = median(RE), 
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975),
            lwr_80 = quantile(RE, 0.1),
            upr_80 = quantile(RE, 0.9))  

biom_sum <- biom_all %>% 
  # filter(conv == "Converged") %>%
  mutate(type = "biom", RE = (f - t) /  t) %>% 
  group_by(year, type) %>% 
  summarize(median = median(RE), 
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975),
            lwr_80 = quantile(RE, 0.1),
            upr_80 = quantile(RE, 0.9)) 
  
all <- rbind(rec_sum, ssb_sum, f_sum, biom_sum) 

ggplot(all, aes(x = year, y = median)) +
  geom_ribbon(aes(ymin = lwr_80, ymax = upr_80), alpha = 0.6, fill = "grey4") +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.4, fill = "grey4") +
  geom_line( color = "white", size = 1,alpha = 1) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 3.8, stroke = 1, alpha = 1) +
  geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 1, alpha = 0.85) +
  facet_wrap(~type, scales = "free") +
  theme_bw() +
  ylim(-0.3, 0.3)

# Check Selex
plot(my_model$rep$F_Slx[1,,,])
lines(Fish_selex_at_age[1,,,,1], col = "red")

plot(my_model$rep$S_Slx[1,,,])
lines(Surv_selex_at_age[1,,,,1], col = "red")

# Check catch
plot(sd_rep$value[str_detect(names(sd_rep$value), "pred_catches")])
# plot(my_model$rep$pred_catches)
lines(obs_catches, col = "red")

# Check indices
plot(sd_rep$value[str_detect(names(sd_rep$value), "pred_srv_indices")])
# plot(my_model$rep$pred_srv_indices)
lines(obs_srv_indices, col = "red")

plot(sd_rep$value[str_detect(names(sd_rep$value), "pred_fish_indices")])
# plot(my_model$rep$pred_fish_indices)
lines(obs_fish_indices, col = "red")

