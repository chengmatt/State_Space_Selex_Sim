# Purpose: A customizable operating model that generates an age-structured population
# Options are to develop a fishery and survey index, as well as fishery and survey comps
# Creator: Matthew LH. Cheng
# Date: 10/30/22

# Set up ------------------------------------------------------------------

library(here)
library(reshape2)
library(tidyverse)
library(wham)
library(TMB)
library(crayon)

# load in functions here 
source(here("R_scripts", "functions", "create_OM_objects.R")) # Create OM objects to hold stuff in
source(here("R_scripts", "functions", "beverton_holt_SR.R")) # Beverton Holt stock recruit
source(here("R_scripts", "functions", "plot_OM_objects.R")) # Plot base operating model objects
source(here("R_scripts", "functions", "read_input_params.R")) # Read in input parameters from an excel sheet
source(here("R_scripts", "functions", "specify_data_scenarios.R")) # Specify data scenarios
source(here("R_scripts", "functions", "specify_selex.R")) # Specify selectivity scenarios
source(here("R_scripts", "functions", "specify_nat_mort.R")) # Specify natural mortality scenarios
source(here("R_scripts", "functions", "specify_q.R")) # Specify catchability scenarios
source(here("R_scripts", "functions", "specify_F_pattern.R")) # Specify fishing mortality scenario
source(here("R_scripts", "functions", "specify_rec_devs.R")) # Specify recruitment deviates here
source(here("R_scripts", "functions", "specify_sex.R")) # Specify recruitment deviates here
source(here("R_scripts", "functions", "sample_index.R")) # Generate index data
source(here("R_scripts", "functions", "sample_comps.R")) # Generate comps data
source(here("R_scripts", "functions", "make_input.R")) # Put stuff into WHAM format
source(here("R_scripts", "functions", "get_results.R")) # Get WHAM data from model

# Path to general input biological parameters
spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")
# spreadsheet_path <- here("input", "PHalibut_Inputs.xlsx") # single-sex and multi-fleet

# Set up OM --------------------------------------------------------

# Read in excel sheet parameters and create OM objects to store values in
read_params_create_OM_objects(spreadsheet_path = spreadsheet_path)

# Specify fishing mortality pattern
get_Fs(Start_F = c(0.01, 0.01), 
       Fish_Start_yr = c(70, 70), 
       F_type = c("Contrast", "Const_Inc"),
       n_years = n_years, 
       max_rel_F_M = c(1.5, 1.5), 
       desc_rel_F_M = c(0.15,NULL), 
       mean_nat_mort = Mean_M,
       yr_chng = 86)

# Specify data scenarios here
fish_surv_data_scenarios(Fish_Start_yr = c(70, 70), 
                         Surv_Start_yr = c(70),
                         fish_Neff = c(100, 100),
                         srv_Neff = c(200), 
                         fish_CV = c(0.1, 0.1), 
                         srv_CV = c(0.1), 
                         Neff_Fish_Time = "F_Vary", 
                         fish_mort = fish_mort,
                         fixed_Neff = 50)

# Specify Natural Mortality
specify_nat_mort(Mort_Time = "Constant", 
                 Mean_M = Mean_M)

# Specify q for the fishery and survey
specify_q(q_Mean_Fish = c(0.08, 0.05),
          q_Mean_Surv = c(0.05))

# Specify recruitment deviates here - loops though each simulation
specify_rec_devs(Rec_Dev_Type = "iid", 
                 rho_rec = NA) 

# Specify selectivity parameterizations here
specify_selex(fish_selex = c("logistic", "logistic"), srv_selex = c("logistic"), 
# Fishery parameters
fish_pars = list(Fleet_1_L = matrix(data = c(9,0.3), nrow = n_sex, byrow = TRUE),
                 Fleet_2_L = matrix(data = c(3,2), nrow = n_sex, byrow = TRUE)),
# Survey parameters
srv_pars = list(Fleet_3_SL = matrix(data = c(5,0.95), nrow = n_sex, byrow = TRUE)), 
bins = ages)

# Specify sex ratios
specify_sex(f_ratio = 1, 
            m_ratio = 0) 

check_equil <- FALSE
rec_type <- "BH"

# Simulation Loop ---------------------------------------------------------

for(sim in 1:n_sims) {
  
  # Check equilibrium -------------------------------------------------------

  if(check_equil == TRUE) {
    
    # Turn fishing off
    rec_devs[,] <- 0 # Turn rec devs off
    fish_mort[,,] <- 0
    
    print("### Checking whether equilibrium conditions have been met ###")
    
  } # checking equilibrium conditions

  # Print simulation iteration
  print(paste("### Simulation",sim,"out of", n_sims, "###"))
  
  # Years loop  -------------------------------------------------------------

  for(y in 1:n_years) {
    
    if(y == 1) { 
      
      for(s in 1:n_sex) {
        
        # Initialize the population here first. We are going to seed the population with a starting number and the sex ratio
        N_at_age[y,1,s,sim] <- N_1 * sex_ratio[y,s] # Put that into our N_at_age array

        # Update Biomass at age after sex ratios have been assigned
        Biom_at_age[y,1,s,sim] <- N_at_age[y,1,s, sim] * wt_at_age[y,1,s,sim] 
      
      } # end sex loop
      
        # Now, calculate our SSB in the first year (only females matter in this case for calculating SSB)
        SSB[y,sim] <- sum(mat_at_age[y,,1,sim] * Biom_at_age[y,,1,sim], na.rm = TRUE)
        
        if(rec_type == "BH") { # do beverton holt recruitment
          # Now, calculate the number of recruits we get - this returns abundance - N at age 2
          rec_total[y,sim] <- beverton_holt_recruit_new(ssb = SSB[y,sim], h = h, r0 = r0, ssb0 = ssb0) * exp(rec_devs[y,sim] - ((sigma_rec^2)/2)) # Add lognormal correction and recdevs
        } 
        if(rec_type == "mean_rec") {
          rec_total[y,sim] <- exp(rnorm(1, mean = mu_rec, sd = sigma_rec))
        } # do mean recruitment
        
      }  # end if statement for if we are in the first year of the simulation

    if(y != 1) { # exiting the first year of the simulation
      
  # Ages Loop ---------------------------------------------------------------

      for(a in 1:length(ages)) {
        
      ### Sexes loop --------------------------------------------------------------
          
          for(s in 1:n_sex) {
            
            if(a != length(ages)) { # if we are not in the plus group nor are we in the recruit age

            # Calculate age, fleet, and sex specific mortality (returns vector fo fleet specific mortalities)
            fleet_mort <- sum(fish_mort[y-1,,sim] * Fish_selex_at_age[y-1,a,,s,sim], na.rm = TRUE)
            
            # Decrement population with Z = M + F
            N_at_age[y,a+1,s,sim] <- N_at_age[y-1,a,s,sim] * exp(-(Mort_at_age[y-1,a,sim] + fleet_mort))
          
            if(a == 1) {
              # Now, add in the recruits from previous year, assigned with the sex ratio
              N_at_age[y,1,s,sim] <- rec_total[y-1,sim] * sex_ratio[y-1,s]
            } # add recruits in at age-1 
          
        } # if we are not in the plus group
        
      ### Decrement population for our + group ---------------------------------------
        
        if(a == length(ages) & !is.na(N_at_age[y-1,length(ages),s,sim])) {
          
          # Calculate fishing mortality for the plus group
          fleet_mort_plus <- sum(fish_mort[y-1,,sim] * Fish_selex_at_age[y-1,a,,s,sim], na.rm = TRUE)
          
          # Applying mortality to plus group individuals last year, and add in recently recruited indviduals into the plus group
          N_at_age[y,a,s,sim] <-  N_at_age[y-1,a,s,sim] * exp(-(Mort_at_age[y-1,a,sim] + fleet_mort_plus)) +  N_at_age[y,a,s,sim] 
          
        } # if we are in the plus group 
          
      } # sexes loop
        
    } # ages loop
      
      ### Update Biomass values and Numbers + Generate Recruits -------------------
        
        # Update Biomass at age 
        Biom_at_age[y,,,sim] <- N_at_age[y,,,sim] * wt_at_age[y,,,sim]
        
        # Now, update SSB (only females matter so indexing 1 for the sex dimension)
        SSB[y,sim] <- sum(mat_at_age[y,,1,sim] * Biom_at_age[y,,1,sim], na.rm = TRUE)
        
        if(rec_type == "BH") { # do beverton holt recruitment
          # Now generate new recruits with the updated SSB
          rec_total[y,sim] <- beverton_holt_recruit_new(ssb = SSB[y,sim], h = h, r0 = r0, ssb0 = ssb0) * 
            exp(rec_devs[y,sim] - ((sigma_rec^2)/2))
        }
        if(rec_type == "mean_rec") {
          rec_total[y,sim] <- exp(rnorm(1, mean = mu_rec, sd = sigma_rec))
        } # do mean recruitment
        
  # Generate observations  ---------------------------------------------------
        
        if(check_equil == FALSE) { # end sampling if we want to check equilibrium
          
          # Get total mortality to calculate catch at age 
          if(n_fish_fleets > 1) { # need to row sum if > 1 fleet
            Z_s <- (Mort_at_age[y-1,,sim] + rowSums(fish_mort[y-1,,sim] * Fish_selex_at_age[y-1,,,s,sim]))
          } else{
            Z_s <- (Mort_at_age[y-1,,sim] + (fish_mort[y-1,,sim] * Fish_selex_at_age[y-1,,,s,sim]))
          } 
          
      ### Fishery fleet loop ------------------------------------------------------

          for(f in 1:n_fish_fleets) { # Loop for fishery fleets
            
            for(s in 1:n_sex) {
              
      ###  Get Catch at Age (Only F to C for now) -----------------------------------
              
              # Calculate instantaneous fishing mortality for a given fleet, sex, and age
              Fish_Fleet_Mort <- (fish_mort[y-1,f,sim] * Fish_selex_at_age[y-1,,f,s,sim])
              
              # Calculate our proportion of mortality via fishing
              Prop_fish_mort <- Fish_Fleet_Mort / Z_s
              
              # Now, get catch at age
              Catch_at_age[y-1,,f,s,sim] <- Prop_fish_mort * N_at_age[y-1,,s,sim] * (1-exp(-Z_s)) * wt_at_age[y,,s,sim]

      ### Sample Fishery Index and Comps ------------------------------------------
              
              # Only start sampling if y > Fish start year. 
              if(y > Fish_Start_yr[f]) { # Observation Model for Fishery
                
                # Generate a fishery index structured by fleet and sex (numbers based)
                Fishery_Index[y-1,f,s,sim] <- sample_index(Idx_Fleet = "Fishery")
                
                # Probability for fishery age comps
                Prob_Fish_Comps <- N_at_age[y-1,,s,sim] * Fish_selex_at_age[y-1,,f,s,sim]

                # Generate comps based 
                Fish_Age_Comps[y-1,,f,s,sim] <- sample_comps(Comp_Fleet = "Fishery",
                                                             error = "multinomial",
                                                             N_eff = fish_Neff[y,f], 
                                                             prob = Prob_Fish_Comps / sum(Prob_Fish_Comps))

              }  # Only start sampling if we are the start of the fish start year
              
            } # end sex index
            
            # Summarize this fishery index aggregated by sex and applying some error
            Fishery_Index_Agg[y-1,f,sim] <- sum(melt(Fishery_Index[y-1,f,,sim]), na.rm = TRUE) # Aggregate
            
            # Apply error here, index fish_CV vector
            Fishery_Index_Agg[y-1,f,sim] <- idx_obs_error(error = "log_normal", 
                                                          true_index = Fishery_Index_Agg[y-1,f,sim],
                                                          CV = fish_CV[f])
          } # end fishery fleet index and loop
          
  # Survey Index and Comps --------------------------------------------------

      ### Survey fleet loop -------------------------------------------------------

          for(sf in 1:n_srv_fleets) { # Loop for survey fleets
            
            for(s in 1:n_sex) {
              
              # Only start sampling if y > Survey Start Year.
              if(y > Surv_Start_yr[sf]) { 
                
                # Get survey index here (numbers based)
                Survey_Index[y-1,sf,s,sim] <- sample_index(Idx_Fleet = "Survey")
                
                # Get probability of sampling a given age class for use in multinomial
                Prob_Surv_at_age <- (N_at_age[y-1,,s,sim] * Surv_selex_at_age[y-1,,sf,s,sim])
                
                # Generate comps based on the expected CPUE at age
                Survey_Age_Comps[y-1,,sf,s,sim] <- sample_comps(Comp_Fleet = "Survey", 
                                                                error = "multinomial",
                                                                N_eff = srv_Neff[y,sf], 
                                                                prob = Prob_Surv_at_age)
                
              } # Only start sampling if we are at the start of the survey start year
              
              # Summarize this fishery index aggregated by sex and applying some error
              Survey_Index_Agg[y-1,sf,sim] <- sum(melt(Survey_Index[y-1,sf,,sim]), na.rm = TRUE) # Aggregate
              
              # Apply error here, index srv_CV vector
              Survey_Index_Agg[y-1,sf,sim] <- idx_obs_error(error = "log_normal", 
                                                            true_index = Survey_Index_Agg[y-1,sf,sim],
                                                            CV = srv_CV[sf])
              
            } # end sex loop for survey here
            
          } # end sf loop
          
        } # end check equilibrium loop (not sampling when checking equilibrium)
        
    } # if we are no longer in the first year

  } # end year loop
  
} # end simulation loop


plot_OM(path = here("figs", "Base_OM_Figs"), file_name = "OM_Check.pdf")


# WHAM checks ------------------------------------------------------------- 

# Set modelling structure here
selectivity <- list(model = c("logistic", "logistic"), 
                    re = c("none", "none"),
                    n_selblocks = 2)

# Fix M for now
M <- list(model = "constant", 
          initial_means = Mean_M) # need to specify est_ages for M to be estimated!

catchability <- list(re = c("none", "none"),
                     q_upper = c(0.5, 0.5)) # put bounds on q

NAA_model <- list(sigma = "rec", 
                  cor = "iid")

# Create objects to store stuff in
conv_vec <- vector(length = n_sims)
results_list <- list() # to store results
sdrep_list <- list()
wham_mod <- list()

# trace(fit_wham, edit = TRUE)

for(sim in 1:n_sims) {
  
  # set.seed(123)
  # Force our inputs into a list - so that it reads into wham
  basic <- make_input(n_fleets = 1, n_indices = 2, 
                      Catch_CV_Val = c(0.05), catch_error = FALSE,
                      n_sims = sim, bias_obs = TRUE, bias_process = TRUE,
                      units_indices = c(1,2), units_index_paa = c(2,2),
                      single_fleet = T, time_block = F, 
                      block_period_sel = list(c( 26), NA), 
                      block_period_idx = list(c( 26), NA))
  
  # Getting time blocks to work here
  # basic$index_cv[,] <- 0.2
  
  # Make WHAM inputs here 
  test_wham <- wham::prepare_wham_input(basic_info = basic, selectivity = selectivity, 
                                        recruit_model = 2, M = M, catchability = catchability,
                                        NAA_re = NULL)

  # Fit WHAM here
  tryCatch( {
    em_fit <- wham::fit_wham(input = test_wham, do.fit = T, do.osa = F, do.retro = F,
                             save.sdrep = TRUE, do.check = F, MakeADFun.silent = T)
    
    wham_mod[[sim]] <- em_fit
    
    convergence_check <- check_convergence(em_fit, ret = T) 
    
    results_list[[sim]] <- get_results(em_fit)
    sdrep_list[[sim]] <- em_fit$sdrep
    
    threp <- em_fit$report()
    sapply(grep("nll",names(threp),value=T), function(x) sum(threp[[x]]))
    
    # Get convergence status
    if(convergence_check$convergence == 0 &  convergence_check$is_sdrep == TRUE & convergence_check$na_sdrep == FALSE) {
      conv_vec[sim] <- "Converged"
    } else {
      conv_vec[sim] <- "Not Converged"
    }
} , error = function(error) {cat("ERROR :",conditionMessage(error), "\n")})
 
    cat(yellow("### Done with Simulation", sim, "###"))
      
} # end loop of number of simulations we want to run

sum(conv_vec == "Converged", na.rm = TRUE)
  
### SSB ---------------------------------------------------------------------

# Do some cleaning up and get RE
ssb_df <- om_em_results(results_list = results_list,
                        n_sims = n_sims, 
                        EM_variable = "SSB", OM_df = SSB,
                        conv_vec = conv_vec)

# Plot SSB with the trend across simulations
ssb_df[[1]] %>% 
  # filter(Sim  %in% c(43:49)) %>%
  ggplot(aes(x = Year, y = SSB, group = Sim, ymin = lwr, ymax = upr)) +
  geom_line() + 
  geom_ribbon(alpha = 0.5) + 
  geom_line(aes(x = Year, y = Truth), color = "red") +
  facet_wrap(~Sim, scales = "free") +
  theme_bw() +
  labs(x = "Year", y  = "SSB") 

# Plot relative error
ssb_df[[2]] %>% 
ggplot(aes(x = Year, y = mean)) +
  geom_line() +
  geom_vline(xintercept = 86, color = "blue", lty = 2, size = 1.2) +
  geom_ribbon(aes(ymin = lwr_95, ymax = up_95), alpha = 0.5) +
  geom_ribbon(aes(ymin = lwr_50, ymax = up_50), alpha = 0.5) +
  geom_hline(aes(yintercept = 0), col = "red", lty = 2, size = 1.2) +
  theme_bw() + 
  labs(x = "Year", y  ="Relative Error in SSB")  +
  ylim(-0.5, 0.5)

# 
# # Selex Comparison --------------------------------------------------------
# # Store values
# fish_sel <- data.frame()
# srv_sel <- data.frame()
# 
# for(i in 1:n_sims) {
#   
#   # Get truth sel
#   true_fish1 <- Fish_selex_at_age[1,,1,,i]
#   true_fish2 <- Fish_selex_at_age[1,,2,,i]
#   true_srv1 <- Surv_selex_at_age[1,,1,,i]
#   # true_srv2 <- Surv_selex_at_age[1,,2,,i]
#   # true_srv3 <- Surv_selex_at_age[1,,3,,i]
#   
#   # Get fishery sel
#   fish_sel_sub <- melt(wham_mod[[i]]$rep$selAA[[1]]) %>%
#     group_by(Var2) %>%
#     summarize(sel = mean(value)) %>%
#     mutate(Sim = i,
#            Conv = conv_vec[i],
#            F1 = true_fish1,
#            F2 = true_fish2,
#            fsh = "F1") %>%
#     rename(Age = Var2)
#   
#   # fish_sel_sub <- melt(wham_mod[[i]]$rep$selAA[[1]]) %>% 
#   #   # group_by(Var2) %>% 
#   #   # summarize(sel = mean(value)) %>% 
#   #   mutate(Sim = i,
#   #          Conv = conv_vec[i]) %>% 
#   #   rename(Age = Var2)
#   
#   # Get fishery sel
#   # fish_sel2_sub <- melt(wham_mod[[i]]$rep$selAA[[2]]) %>%
#   #   group_by(Var2) %>%
#   #   summarize(sel = mean(value)) %>%
#   #   mutate(Sim = i,
#   #          Conv = conv_vec[i],
#   #          F1 = true_fish2,
#   #          fsh = "F2") %>%
#   #   rename(Age = Var2)
#   # 
#   # Get survey sel 1
#   srv_sel_1sub <- melt(wham_mod[[i]]$rep$selAA[[2]]) %>%
#     group_by(Var2) %>%
#     summarize(sel = mean(value)) %>%
#     mutate(Sim = i,
#            Conv = conv_vec[i],
#            Truth = true_srv1,
#            Srv_Sel = "Srv_1") %>%
#     rename(Age = Var2)
#   # 
#   # # Get survey sel 2
#   # srv_sel_2sub <- melt(wham_mod[[i]]$rep$selAA[[3]]) %>% 
#   #   group_by(Var2) %>% 
#   #   summarize(sel = mean(value)) %>% 
#   #   mutate(Sim = i,
#   #          Conv = conv_vec[i],
#   #          Truth = true_srv2,
#   #          Srv_Sel = "Srv_2") %>% 
#   #   rename(Age = Var2)
#   # 
#   # # Get survey sel 3
#   # srv_sel_3sub <- melt(wham_mod[[i]]$rep$selAA[[4]]) %>% 
#   #   group_by(Var2) %>% 
#   #   summarize(sel = mean(value)) %>% 
#   #   mutate(Sim = i,
#   #          Conv = conv_vec[i],
#   #          Truth = true_srv3,
#   #          Srv_Sel = "Srv_3") %>% 
#   #   rename(Age = Var2)
#   
#   fish_sel <- rbind(fish_sel, fish_sel_sub)
#   srv_sel <- rbind(srv_sel, srv_sel_1sub)
# }
# 
# # Get Selex weighted by fishing mortality
# # Fmort rates
# c_prop <- melt(Catch_at_age[(Fish_Start_yr[1]:(n_years-1)),,,1,1]) %>% 
#   group_by(Var1, Var3) %>% 
#   summarize(value = sum(value)) %>% 
#   group_by(Var1) %>% 
#   mutate(sum_all = sum(value),
#          prop = value / sum_all )%>% 
#   select(Var1, Var3, prop)
# 
# 
# # fish_prop <- melt(fish_mort[(Fish_Start_yr[1]:(n_years-1)),,1]) %>% 
# #   group_by(Var1) %>% 
# #   mutate(sum = sum(value),
# #          prop = value / sum) %>% 
# #   select(Var1, prop, Var2)
# 
# # Gte selexe prop
# fish_sel_prop <- melt(Fish_selex_at_age[(Fish_Start_yr[1]:(n_years-1)),,,1,1]) %>% 
#   left_join(c_prop, by = c("Var1", "Var3" ))
# 
# sel_pr <- fish_sel_prop %>% 
#   mutate(value = value * prop,
#          Var2 = parse_number(paste(Var2)),
#          Var1 = parse_number(paste(Var1))) %>% 
#   group_by(Var1, Var2) %>% 
#   summarize(mean = mean(value)) %>% 
#   group_by(Var1) %>% 
#   mutate(max = max(mean),
#          mean = mean/max)
# 
# sel_pr %>% 
#   ggplot(aes(x = Var2, y = mean, color = Var1, group = Var1)) +
#   geom_line() +
#   scale_color_viridis_c()
# 
# # fish sel
# ggplot() +
#   # geom_line(sel_pr, 
#   #           mapping = aes(x = Var2, y = mean, group = Var1), color = "blue") +
#   geom_line(fish_sel %>% 
#               filter(Sim == 1), 
#             mapping = aes(x = Age, y = sel, group = Sim)) 
#   # scale_color_viridis_c() +
#   # facet_wrap(~Var1)
#   # geom_line(aes(x = Age, y = F1), col = "blue", lty = 2) +
#   # geom_line(aes(x = Age, y = F2), col = "orange", lty = 2) 
# 
# ggplot(fish_sel_sub, aes(x = Age, y = value, group = Var1,
#                          color = Var1)) +
#   geom_line() +
#   scale_color_viridis_c()
# 
# # srv sel
# ggplot(srv_sel, aes(x = Age, y = sel, group = Sim, color = Conv)) +
#   geom_line() +
#   geom_line(aes(x = Age, y = Truth), col = "black", lty = 2) +
#   facet_wrap(~Srv_Sel) +
#   scale_color_manual(values = c("grey", "red"))
# 
# # 
# # # Parameter checks ------------------------------------------------------
# # m_vec <- vector()
# # q1_vec <- vector()
# # q2_vec <- vector()
# # q3_vec <- vector()
# # q4_vec <- vector()
# # terminal_ssb <- vector()
# # terminal_F <- vector()
# # sigma_rec <- vector()
# # 
# # for(i in 1:n_sims) {
# #   
# #    # Get M
# #   # m_vec[i] <- (exp(sdrep_list[[i]]$par.fixed[names(sdrep_list[[i]]$par.fixed) == "M_a"]) - Mean_M) / Mean_M
# #   # Get q estiamtes
# #   q1 <- sdrep_list[[i]]$par.fixed[names(sdrep_list[[i]]$par.fixed) == "logit_q"][1]
# #   q2 <- sdrep_list[[i]]$par.fixed[names(sdrep_list[[i]]$par.fixed) == "logit_q"][2]
# #   # q3 <- sdrep_list[[i]]$par.fixed[names(sdrep_list[[i]]$par.fixed) == "logit_q"][3]
# #   # q4 <- sdrep_list[[i]]$par.fixed[names(sdrep_list[[i]]$par.fixed) == "logit_q"][4]
# #   
# #   # Get terminal ssb
# #   terminal_ssb[i] <- (ssb_df[[1]]$SSB[ssb_df[[1]]$Year == 90 & ssb_df[[1]]$Sim == i] - 
# #                         ssb_df[[1]]$Truth[ssb_df[[1]]$Year == 90 & ssb_df[[1]]$Sim == i]) / 
# #                         ssb_df[[1]]$Truth[ssb_df[[1]]$Year == 90 & ssb_df[[1]]$Sim == i]
# #   
# #   # Get terminal F
# #   terminal_F[i] <- (exp(summary(sdrep_list[[i]])[rownames(summary(sdrep_list[[i]])) == "log_F"][31]) - fish_mort[90,,i])/ ( fish_mort[90,,i])
# #   
# #   # Get sigma rec
# #   # sigma_rec[i] <- (exp(summary(sdrep_list[[i]])[rownames(summary(sdrep_list[[i]])) == "log_NAA_sigma"][2]) - 1.2) / 1.2
# #   
# #   # Get catchability estimates
# #   q1_vec[i] <- ((0 + (1-0)) / (1+exp(-1*q1)) - 0.035) / 0.035
# #   q2_vec[i] <- ((0 + (1-0)) / (1+exp(-1*q2)) - 0.05) / 0.05
# #   # q3_vec[i] <- ((0 + (1000-0)) / (1+exp(-1*q3)) - 0.01) / 0.01
# #   # q4_vec[i] <- ((0 + (1000-0)) / (1+exp(-1*q4)) - 0.035) / 0.035
# #   
# # } # end nsims loop
# # 
# # data_summary <- function(x) {
# #   m <- median(x)
# #   ymin <- m-sd(x)
# #   ymax <- m+sd(x)
# #   return(c(y=m,ymin=ymin,ymax=ymax))
# # }
# # 
# # # Put these into a df
# # par_df <- data.frame( q1 = q1_vec, q2 = q2_vec,
# #                      Term_SSB = terminal_ssb, Term_F = terminal_F) %>% 
# #   pivot_longer(names_to = c("par"), values_to = "val", cols = everything()) %>% 
# #   drop_na() %>% 
# #   group_by(par) %>% 
# #   mutate(median_col = abs(median(val)))
# # 
# # ggplot(par_df %>% 
# #          filter(par != "sigma_rec"), aes(x = par, y = val, fill = par)) +
# #   geom_violin(alpha = 0.5) +
# #   geom_hline(aes(yintercept = 0), col = "red", lty = 2, size = 1.5) +
# #   geom_boxplot(width = 0.1, alpha = 0.5) +
# #   # stat_summary(fun.data = data_summary, size = 0.8, alpha = 0.8) +
# #   scale_fill_brewer(palette = "Blues")+
# #   labs(x = "Parameter", y = "Relative Error") +
# #   theme_minimal() +
# #   theme(legend.position = "none") +
# #   ylim(-0.5, 0.5)

# Other plots -------------------------------------------------------------

plot_OM(path = here("figs", "Base_OM_Figs"), file_name = "OM_Check.pdf")
 
# Output plots
# Create directory to ouptut plots to
wham_out <- here("figs", "wham_checks")
# dir.create(wham_out)
plot_wham_output(wham_mod[[50]], dir.main = wham_out)



