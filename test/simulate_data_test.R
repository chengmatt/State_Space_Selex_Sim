# Purpose: A customizable operating model that generates an age-structured population
# Options are to develop a fishery and survey index, as well as fishery and survey comps
# Creator: Matthew LH. Cheng
# Date: 10/30/22

# Set up ------------------------------------------------------------------

library(here)
library(reshape2)
library(tidyverse)

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
source(here("R_scripts", "functions", "make_wham_input.R")) # Put stuff into WHAM format
source(here("R_scripts", "functions", "get_results.R")) # Get WHAM data from model

# Path to general input biological parameters
spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")

# Set up OM --------------------------------------------------------

# Read in excel sheet parameters and create OM objects to store values in
read_params_create_OM_objects(spreadsheet_path = spreadsheet_path)

# Specify fishing mortality pattern
get_Fs(Start_F = c(0.01), 
     Fish_Start_yr = c(70), 
     F_type = c("Contrast"),
     n_years = n_years, 
     max_rel_F_M = c(1.5), 
     desc_rel_F_M = c(0.15), 
     mean_nat_mort = Mean_M,
     yr_chng = 86)

# Specify data scenarios here
fish_surv_data_scenarios(Fish_Start_yr = c(70), 
                       Surv_Start_yr = c(70),
                       fish_Neff_max = c(150),
                       srv_Neff_max = c(150), 
                       fish_CV = c(0.1), 
                       srv_CV = c(0.1), 
                       Neff_Fish_Time = "F_Vary", 
                       fish_mort = fish_mort,
                       fixed_Neff = 30)

# Specify Natural Mortality
specify_nat_mort(Mort_Time = "Constant", 
               Mean_M = Mean_M)

# Specify q for the fishery and survey
specify_q(q_Mean_Fish = c(0.08),
        q_Mean_Surv = c(0.01))

# Specify recruitment deviates here - loops though each simulation
specify_rec_devs(Rec_Dev_Type = "iid", 
               rho_rec = NA) 

# Specify selectivity parameterizations here
specify_selex(fish_selex = c("gamma"), srv_selex = c("logistic"), 
# Fishery parameters
fish_pars = list(Fleet_1_L = matrix(data = c(6,3), nrow = n_sex, byrow = TRUE)),
# Survey parameters
srv_pars = list(Fleet_3_SL = matrix(data = c(4,0.8), nrow = n_sex, byrow = TRUE)), 
bins = ages)

# Specify sex ratios
specify_sex(f_ratio = 1, m_ratio = 0) 

check_equil <- FALSE
rec_type <- "mean_rec"

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
      
      # Recruitment at first year = 0
      rec_total[y, sim] <- 0
      
    }  # end if statement for if we are in the first year of the simulation

  if(y != 1) { # exiting the first year of the simulation
    
    if(rec_type == "BH") { # do beverton holt recruitment
      # Now generate new recruits with the updated SSB
      rec_total[y,sim] <- beverton_holt_recruit_new(ssb = SSB[y,sim], h = h, r0 = r0, ssb0 = ssb0) * 
        exp(rec_devs[y,sim] - ((sigma_rec^2)/2))
    }
    if(rec_type == "mean_rec") {
      rec_total[y,sim] <- exp(mu_rec + rec_devs[y,sim] - ((sigma_rec^2)/2))
    } # do mean recruitment
    
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
            N_at_age[y,1,s,sim] <- rec_total[y,sim] * sex_ratio[y-1,s]
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

            # Now, get catch at age
            Catch_at_age[y-1,,f,s,sim] <- Fish_Fleet_Mort * N_at_age[y-1,,s,sim] * (1-exp(-Z_s)) * wt_at_age[y,,s,sim] / Z_s

    ### Sample Fishery Index and Comps ------------------------------------------
            
            # Only start sampling if y > Fish start year. 
            if(y > Fish_Start_yr[f]) { # Observation Model for Fishery
              
              # Generate a fishery index structured by fleet and sex (numbers based)
              Fishery_Index[y-1,f,s,sim] <- sample_index(Idx_Fleet = "Fishery")
              
              # Probability for fishery age comps
              Prob_Fish_Comps <- Fish_Fleet_Mort * N_at_age[y-1,,s,sim] * (1-exp(-Z_s)) / Z_s

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
                                                              prob = Prob_Surv_at_age / sum(Prob_Surv_at_age))
              
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
