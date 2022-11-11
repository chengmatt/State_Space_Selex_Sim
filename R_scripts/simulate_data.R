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

# Set up OM --------------------------------------------------------

# Read in excel sheet parameters and create OM objects to store values in
read_params_create_OM_objects(spreadsheet_path = spreadsheet_path)

# Specify fishing mortality pattern
get_Fs(Start_F = c(0.001, 0.001), Fish_Start_yr = c(150, 170), F_type = c("Contrast", "Increase_Plat"),
       n_years = n_years, max_rel_F_M = c(0.6, 0.75), desc_rel_F_M = c(0.5, 0), mean_nat_mort = Mean_M)

# Specify data scenarios here
fish_surv_data_scenarios(Fish_Start_yr = c(150, 175), Surv_Start_yr = c(160, 180),
                         fish_Neff = c(150, 150), srv_Neff = c(100, 75), 
                         fish_CV = c(0.2, 0), srv_CV = c(0.1, 0.1), 
                         Neff_Fish_Time = "F_Vary", fish_mort = fish_mort)

# Specify Natural Mortality
specify_nat_mort(Mort_Time = "Constant", Mean_M = Mean_M)

# Specify q for the fishery and survey
specify_q(q_Mean_Fish = c(0.001, 0.005), q_Mean_Surv = c(0.0001, 0.0005))

# Specify recruitment deviates here - loops though each simulation
specify_rec_devs(Rec_Dev_Type = "iid", rho_rec = NA) 

# Specify selectivity parameterizations here
specify_selex(fish_selex = c("logistic", "double_normal"), srv_selex = c("double_logistic", "uniform"), 
# Fishery parameters
fish_pars = list(Fleet_1_L <- matrix(data = c(2,5, 3,1), nrow = n_sex, byrow = TRUE),
                 Fleet_2_DN <- matrix(data = c(5,0.01,3,10,0.2,0.1,
                                               10,0.01,3,15,0.2,0.05), nrow = n_sex, byrow = TRUE)),
# Survey parameters
srv_pars = list(Fleet_1_DL <- matrix(data = c(3, 0.2, 5, 20, 
                                              4, 0.5, 8, 15), nrow = n_sex, byrow = TRUE),
                Fleet_2_U <- matrix()), bins = ages)

# Specify sex ratios
specify_sex(f_ratio = 0.8, m_ratio = 0.2)

check_equil <- TRUE

# Start Simulation to get to equilibrium conditions --------------------------------------------------------

# Start simulation loop here
for(sim in 1:n_sims) {
  
# Check equilibrium -------------------------------------------------------

  if(check_equil == TRUE) {
    
    # Turn fishing off
    rec_devs[,] <- 0 # Turn rec devs off
    Fish_selex_at_age[,,,,] <- 0 # Selectivity
    fish_mort[,,] <- 0
    
    print("### Checking whether equilibrium conditions have been met ###")
    
  } # checking equilibrium conditions

  # Print simulation iteration
  print(paste("### Simulation",sim,"out of", n_sims, "###"))
  
# Years loop  -------------------------------------------------------------

  for(y in 1:n_years) {
    
    #### Initialize Population here ----------------------------------------------

    if(y == 1) { 
      
      for(s in 1:n_sex) {
        
        # Initialize the population here first. We are going to seed the population with a starting number and the sex ratio
        N_at_age[y,1,s,sim] <- N_1 * sex_ratio[y,s] # Put that into our N_at_age array

        # Update Biomass at age after sex ratios have been assigned
        Biom_at_age[y,1,s,sim] <- N_at_age[y,1,s, sim] * wt_at_age[y,1,s,sim] 
      
      } # end sex loop
      
        # Now, calculate our SSB in the first year (only females matter in this case for calculating SSB)
        # Make sure we are not indexing any of the NAs in maturity at age and incorrectly multiplying vectors.
        SSB[y,sim] <- sum(mat_at_age[y,,1,sim][!is.na(Biom_at_age[y,,1,sim])] * Biom_at_age[y,,1,sim], na.rm = TRUE)
        
        # Now, calculate the number of recruits we get - this returns abundance - N at age 2
        rec_total[y,sim] <- beverton_holt_recruit_new(ssb = SSB[y,sim], h = h, r0 = r0, ssb0 = ssb0) * exp(rec_devs[y,sim] - ((sigma_rec^2)/2)) # Add lognormal correction and recdevs
      
      
        }  # if we are in the first year of the simulation

    if(y != 1) { # exiting the first year of the simulation
      
# Ages Loop ---------------------------------------------------------------

      for(a in 1:length(ages)) {
        
# Sexes loop --------------------------------------------------------------
          
          for(s in 1:n_sex) {
            
            if(a != length(ages)) { # if we are not in the plus group

    ### Decrement population with F and Z (!= + group) ---------------------------------------
            
          # Calculate age, fleet, and sex specific mortality (returns vector fo fleet specific mortalities)
          fleet_mort <- fish_mort[y-1,,sim] * Fish_selex_at_age[y-1,a,,s,sim]
          
          # Decrement population with Z = M + F
          N_at_age[y,a+1,s,sim] <- N_at_age[y-1,a,s,sim] * exp(-(Mort_at_age[y-1,a,sim] + sum(fleet_mort, na.rm = TRUE)))

          # Now, add in the recruits from previous year, assigned with the sex ratio
          N_at_age[y,1,s,sim] <- rec_total[y-1,sim] * sex_ratio[y-1,s]
          
        } # if we are not in the plus group
        
          
     ### Decrement population for our PLUS GROUP ---------------------------------------
        
        if(a == length(ages) & !is.na(N_at_age[y-1,length(ages),s,sim])) {
          
          # Calculate fishing mortality for the plus group
          fleet_mort_plus <- fish_mort[y-1,,sim] * Fish_selex_at_age[y-1,a,,s,sim]
          
          # Applying natural and fishing mortality to indivduals from the previous year  that were in the plus group
          N_at_age[y,a,s,sim] <-  (N_at_age[y-1,a,s,sim] * # Z = M + F
                                  exp(-(Mort_at_age[y-1,a,sim] + sum(fleet_mort_plus, na.rm = TRUE)))) +  
                                  N_at_age[y,a,s,sim] # Adding back the individuals that recently recruited into + group
          
        } # if we are in the plus group and we had a plus group in the previous year (so we are not accidentally
        # adding additional plus groups together)
          
        } # sexes loop
        
      } # ages loop
      
        ### Update Biomass values and Numbers + Generate Recruits -------------------
        
        # Update Biomass at age 
        Biom_at_age[y,,,sim] <- N_at_age[y,,,sim] * wt_at_age[y,,,sim]
        
        # Now, update SSB for year y and generate recruits for the new year with the updated SSB
        # Only females matter so indexing 1 for the sex dimension
        SSB[y,sim] <- sum(mat_at_age[y,,1,sim] * Biom_at_age[y,,1,sim], na.rm = TRUE)
          
        # Now generate new recruits with the updated SSB
        rec_total[y,sim] <- beverton_holt_recruit_new(ssb = SSB[y,sim], h = h, r0 = r0, ssb0 = ssb0) * exp(rec_devs[y,sim] - ((sigma_rec^2)/2))
        
  # Generate observations  ---------------------------------------------------
        
        if(check_equil == FALSE) { # end sampling if we want to check equilibrium
          
          for(f in 1:n_fish_fleets) { # Loop for fishery fleets
            
            for(s in 1:n_sex) {
              
              ###  Get Catch at Age (Only F to C for now) -----------------------------------
              
              # Calculate total mortality for a given sex (sum row-wise to get age specific total fishing mortality across fleets)
              Z_s <- (Mort_at_age[y-1,,sim] + rowSums(fish_mort[y-1,,sim] * Fish_selex_at_age[y-1,,,s,sim]))
              
              # Calculate instantaneous fishing mortality for a given fleet, sex, and age
              Fish_Fleet_Mort <- (fish_mort[y-1,f,sim] * Fish_selex_at_age[y-1,,f,s,sim])
              
              # Calculate our proportion of mortality via fishing
              Prop_fish_mort <- Fish_Fleet_Mort / Z_s
              
              # Now, get catch at age
              Catch_at_age[y-1,,f,s,sim] <- Prop_fish_mort * N_at_age[y-1,,s,sim] * (1-exp(-Z_s)) * wt_at_age[y,,s,sim]
              
              
              ### Sample Fishery Index and Comps ------------------------------------------
              
              # Only start sampling if y > Fish start year. If we want to cut out certain years, we can do this afterwards.
              if(y > min(Fish_Start_yr)) { # Observation Model for Fishery
                
                # Generate a fishery index structured by fleet and sex (numbers based)
                Fishery_Index[y-1,f,s,sim] <- sample_index(Idx_Fleet = "Fishery")
                
                # Generate comps based on the expected catch at age
                Fish_Age_Comps[y-1,,f,s,sim] <- sample_comps(Comp_Fleet = "Fishery", error = "multinomial",
                                                                 N_eff = fish_Neff[y,f])
                
              }  # Only start sampling if we are the start of the fish start year
              
            } # end sex index
            
            # Summarize this fishery index aggregated by sex and applying some error
            Fishery_Index_Agg[y-1,f,sim] <- sum(melt(Fishery_Index[y-1,f,,sim]), na.rm = TRUE) # Aggregate
            
            # Apply error here, index fish_CV vector
            Fishery_Index_Agg[y-1,f,sim] <- idx_obs_error(error = "log_normal", 
                                                          true_index = Fishery_Index_Agg[y-1,f,sim],
                                                          CV = fish_CV[f])
            
          } # end fishery fleet index and loop
          
          ### Survey Index and Comps --------------------------------------------------
          
          for(sf in 1:n_srv_fleets) { # Loop for survey fleets
            
            for(s in 1:n_sex) {
              
              # Only start sampling if y > Survey Start Year. We can cut out certain years if we 
              # want irregular sampling, etc. after the fact.
              if(y > min(Surv_Start_yr)) { 
                
                # Get survey index here (numbers based)
                Survey_Index[y-1,sf,s,sim] <- sample_index(Idx_Fleet = "Survey")
                
                # Generate comps based on the expected CPUE at age
                Survey_Age_Comps[y-1,,sf,s,sim] <- sample_comps(Comp_Fleet = "Survey", error = "multinomial",
                                                                    N_eff = srv_Neff[sf])
                
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
selectivity <- list(model = rep("logistic", 2))
M <- list(model = "constant", initial_means = 0.1, est_ages = 1) # need to specify est_ages for M to be estimated!
catchability <- list(re = c("none", "none")) # put bounds on q
NAA_model <- list(sigma = "rec", cor = "iid")

# Create objects to store stuff in
conv_vec <- vector()
results_list <- list() # to store results
sigma_rec_val <- vector()
m_val <- vector()

for(sim in 1:n_sims) {
  
  # Force our inputs into a list - so that it reads into wham -
  basic <- make_input(n_fleets = 1, n_indices = 2, Catch_CV_Val = 0.05, catch_error = FALSE,
                      n_sims = sim, bias_obs = TRUE, bias_process = TRUE, 
                      units_indices = c(2,2), units_index_paa = c(2,2)) 
  
  # Make WHAM inputs here 
  test_wham <- wham::prepare_wham_input(basic_info = basic, selectivity = selectivity, 
                                        recruit_model = 2, M = M, catchability = catchability,
                                        NAA_re = NAA_model)
  
  # Fit WHAM here
  tryCatch( {
    em_fit <- wham::fit_wham(input = test_wham, do.fit = T, do.osa = F, do.retro = F,
                             save.sdrep = TRUE, do.check = F, MakeADFun.silent = TRUE)
    
    convergence_check <- check_convergence(em_fit, ret = T) 
    
    results_list[[sim]] <- get_results(em_fit)
  } , error = function(error) {cat("ERROR :",conditionMessage(error), "\n")})
  
  # Get convergence status
  if(convergence_check$convergence == 0 & convergence_check$maxgr <= 0.001 &
     convergence_check$is_sdrep == TRUE & convergence_check$na_sdrep == FALSE) {
      conv_vec[sim] <- "Converged"
  } else conv_vec[sim] <- "Not Converged"
  
  # Get sigma rec re value
  sigma_rec_val[sim] <- em_fit$sdrep$value[which(names(em_fit$sdrep$value) == "NAA_sigma")]
  m_val[sim] <- exp(em_fit$sdrep$par.fixed[which(names(em_fit$sdrep$par.fixed) == "M_a")]) # Get M
  print(paste("### Done with Simulation", sim, "###")) 

} # end loop of number of simulations we want to run


# Plot checks -------------------------------------------------------------

# Recdevs RE
dir.rec.re <- here("output", "rec_dev_RE")
dir.create(dir.rec.re)
save(results_list, file = here(dir.rec.re, "rec_re.RData"))


### SSB ---------------------------------------------------------------------

# Do some cleaning up and get RE
ssb_df <- om_em_results(n_sims = n_sims, EM_variable = "SSB", OM_df = SSB)

# Plot SSB with the trend across simulations
ssb_df[[1]] %>% 
  ggplot(aes(x = Year, y = SSB, group = Sim, ymin = lwr, ymax = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) + 
  geom_line(aes(x = Year, y = Truth), color = "red") +
  facet_wrap(~Sim) +
  theme_bw() +
  labs(x = "Year", y  = "SSB") 

# Plot relative error
ggplot(ssb_df[[2]], aes(x = Year, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr_95, ymax = up_95), alpha = 0.5) +
  geom_ribbon(aes(ymin = lwr_50, ymax = up_50), alpha = 0.5) +
  geom_hline(aes(yintercept = 0), col = "red", lty = 2, size = 1.5) +
  theme_bw() + 
  labs(x = "Year", y  ="Relative Error in SSB") +
  ylim(-0.3,0.3)


### F -----------------------------------------------------------------------

# Do some cleaning up and get RE
f_df <- om_em_results(n_sims = n_sims, EM_variable = "F", OM_df = fish_mort)

# Plot F with the trend across simulations
f_df[[1]] %>% 
  ggplot(aes(x = Year, y = F_val, group = Sim, ymin = lwr, ymax = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) + 
  geom_line(aes(x = Year, y = Truth), color = "red") +
  facet_wrap(~Sim) +
  theme_bw() +
  labs(x = "Year", y  = "Fishing Mortality Rate") 

# Plot relative error
ggplot(f_df[[2]], aes(x = Year, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr_95, ymax = up_95), alpha = 0.5) +
  geom_ribbon(aes(ymin = lwr_50, ymax = up_50), alpha = 0.5) +
  geom_hline(aes(yintercept = 0), col = "red", lty = 2, size = 1.5) +
  theme_bw() +
  labs(x = "Year", y  ="Relative Error in Fishing Mortality") +
  ylim(-0.3,0.3)


### Catchability (Fishery) ------------------------------------------------------------

q_Fish_df <- om_em_results(n_sims = n_sims, EM_variable = "q_Fish", OM_df = q_Fish)

# Plot F with the trend across simulations
q_Fish_df[[1]] %>% 
  ggplot(aes(x = Year, y = q_Fish, group = Sim, ymin = q_Fish_lwr, ymax = q_Fish_upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  geom_line(aes(x = Year, y = Truth), color = "red") +
  facet_wrap(~Sim) +
  theme_bw() +
  labs(x = "Year", y  = "Fishing Mortality Rate") 

# Plot relative error
ggplot(q_Fish_df[[2]], aes(x = Year, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr_95, ymax = up_95), alpha = 0.5) +
  geom_ribbon(aes(ymin = lwr_50, ymax = up_50), alpha = 0.5) +
  geom_hline(aes(yintercept = 0), col = "red", lty = 2, size = 1.5) +
  theme_bw() +
  labs(x = "Year", y  ="Relative Error in Fishery Catchability") +
  ylim(-0.3,0.3)


# Sigma rec + Mortality ---------------------------------------------------------------
save(sigma_rec_val, file = here(dir.rec.re, "sigma_rec_values.RData"))
save(m_val, file = here(dir.rec.re, "m_values.RData"))

sig_rec_df <- data.frame(sigma_rec = sigma_rec_val, truth = sigma_rec) %>% 
  mutate(RE = (sigma_rec - truth )/ truth,
         par = "Sigma_Rec (0.5)") %>% 
  select(par, RE)

mort_df <- data.frame(m = m_val, truth = Mean_M) %>% 
  mutate(RE = (m_val - truth )/ truth,
         par = "Natural Mortality (0.1)") %>% 
  select(par, RE)

all_par <- rbind(sig_rec_df, mort_df)

ggplot(all_par, aes(x = par, y = RE)) +
  geom_violin(width = 0.3, alpha = 0.75) +
  geom_hline(aes(yintercept = 0), col = "red", lty = 2, size = 1.5) +
  theme_bw() + stat_summary(
    fun.min = function(z) { quantile(z,0.95) },
    fun.max = function(z) { quantile(z,0.05) },
    fun = median, size = 1, color = "black") +
  theme(axis.text = element_text(size = 15))
  
  


# Other plots -------------------------------------------------------------

plot_OM(path = here("figs", "Base_OM_Figs"), file_name = "OM_Check.pdf")
 
# Output plots
# Create directory to ouptut plots to
wham_out <- here("figs", "wham_checks")
# dir.create(wham_out)
plot_wham_output(em_fit, dir.main = wham_out)


