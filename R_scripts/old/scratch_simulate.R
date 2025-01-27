# Purpose: Function to simulate data, which is meant to be called at the end
# of specifying the OM
# Date 1/2/23
# Creator: Matthew LH. Cheng

#' @param DM_Fish_Param Weight parameter for Dirichlet Multinomial in the fishery
#' @param DM_Srv_Param Weight parameter for Dirichlet Mutlinomial in the survey
#' @param fish_likelihood Fishery compositional likelihood
#' @param srv_likelihood Survey compositional likelihood

simulate_data <- function(fxn_path, 
                          spreadsheet_path, 
                          check_equil = FALSE,
                          rec_type = "mean_rec",
                          Start_F,
                          Fish_Start_yr,
                          Surv_Start_yr,
                          max_rel_F_M,
                          desc_rel_F_M,
                          F_type,
                          yr_chng, 
                          yr_chng_end,
                          fish_likelihood,
                          srv_likelihood = "multinomial",
                          Input_Fish_N_Max, 
                          Input_Srv_N_Max, 
                          fish_CV, 
                          srv_CV, 
                          DM_Fish_Param = NULL,
                          DM_Srv_Param = NULL,
                          catch_CV,
                          Input_N_Fish_Time,
                          Input_N_Fish_Fixed,
                          Mort_Time = "Constant",
                          q_Mean_Fish,
                          q_Mean_Surv,
                          Rec_Dev_Type = "iid",
                          rho_rec = NA,
                          fish_selex, 
                          srv_selex,
                          fish_pars,
                          srv_pars,
                          f_ratio, 
                          m_ratio) {
  
  require(tidyverse)
  require(reshape2)
  
  # Load in all functions from the functions folder
  files <- list.files(fxn_path)
  for(i in 1:length(files)) source(here(fxn_path, files[i]))
  
  # Read in parameters from our spreadsheet
  read_params_create_OM_objects(spreadsheet_path = spreadsheet_path)
  
  # Get fishing mortality pattern
  get_Fs(Start_F = Start_F, 
         Fish_Start_yr = Fish_Start_yr, 
         F_type = F_type,
         n_years = n_years, 
         max_rel_F_M = max_rel_F_M, 
         desc_rel_F_M = desc_rel_F_M, 
         mean_nat_mort = Mean_M,
         yr_chng = yr_chng,
         yr_chng_end = yr_chng_end)
  
  # Specify data scenarios here
  fish_surv_data_scenarios(Fish_Start_yr = Fish_Start_yr, 
                           Surv_Start_yr = Surv_Start_yr,
                           Input_Fish_N_Max = Input_Fish_N_Max,
                           Input_Srv_N_Max = Input_Srv_N_Max, 
                           fish_CV = fish_CV, 
                           srv_CV = srv_CV, 
                           Input_N_Fish_Time = Input_N_Fish_Time, 
                           fish_mort = fish_mort,
                           Input_N_Fish_Fixed = Input_N_Fish_Fixed)
  
  # Specify Natural Mortality
  specify_nat_mort(Mort_Time = Mort_Time, 
                   Mean_M = Mean_M)
  
  # Specify q for the fishery and survey
  specify_q(q_Mean_Fish = q_Mean_Fish,
            q_Mean_Surv = q_Mean_Surv)
  
  # Specify recruitment deviates here
  specify_rec_devs(Rec_Dev_Type = Rec_Dev_Type, rho_rec = rho_rec) 
  
  # Specify selectivity parameterizations here
  specify_selex(fish_selex = fish_selex, srv_selex = srv_selex, 
                fish_pars = fish_pars, srv_pars = srv_pars, bins = ages)
  
  # Specify sex ratios
  specify_sex(f_ratio = f_ratio, m_ratio = m_ratio) 
  
  
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
    
    for(y in 1:(n_years-1)) {
      
      if(y == 1) { 
        
        for(s in 1:n_sex) {
          
          # Initialize the population here first. We are going to seed the population with a starting number and the sex ratio
          N_at_age[y,1,s,sim] <- N_1 * sex_ratio[y,s] # Put that into our N_at_age array
          
          # Update Biomass at age after sex ratios have been assigned
          Biom_at_age[y,1,s,sim] <- N_at_age[y,1,s, sim] * wt_at_age[y,1,s,sim] 
          
        } # end sex loop
        
        # Recruitment at first year = 0
        rec_total[y, sim] <- 0
        
        # Now, calculate our SSB in the first year (only females matter in this case for calculating SSB)
        SSB[y,sim] <- sum(mat_at_age[y,,1,sim] * Biom_at_age[y,,1,sim], na.rm = TRUE)
        
      }  # end if statement for if we are in the first year of the simulation
      
      if(y != 1) { # exiting the first year of the simulation
        
        if(rec_type == "BH") { # do beverton holt recruitment
          # Now generate new recruits with the updated SSB
          rec_total[y,sim] <- beverton_holt_recruit(ssb = SSB[y-1,sim], h = h, r0 = r0) * 
            exp(rec_devs[y,sim] - ((sigma_rec^2)/2))
        }
        if(rec_type == "mean_rec") {
          rec_total[y,sim] <- exp(mu_rec + rec_devs[y,sim] - ((sigma_rec^2)/2))
        } # do mean recruitment
        
        # Ages Loop ---------------------------------------------------------------
        
        for(a in 1:length(ages)) {
          
          ### Sexes loop --------------------------------------------------------------
          
          for(s in 1:n_sex) {
            
            if(a == 1) {
              # Now, add in the recruits, assigned with the sex ratio
              N_at_age[y,1,s,sim] <- rec_total[y,sim] * sex_ratio[y-1,s]
            } # add recruits in at age-1 
            
            if(a < length(ages)) { # if we are not in the plus group nor are we in the recruit age
              
              # Calculate age, fleet, and sex specific mortality (returns vector fo fleet specific mortalities)
              fleet_mort <- sum(fish_mort[y,,sim] * Fish_selex_at_age[y,a,,s,sim], na.rm = TRUE)
              
              # Decrement population with Z = M + F
              N_at_age[y+1,a+1,s,sim] <- N_at_age[y,a,s,sim] * exp(-(Mort_at_age[y,a,sim] + fleet_mort))
              
            } # if we are not in the plus group
            
            ### Decrement population for our + group ---------------------------------------
            
            if(a == length(ages) & !is.na(N_at_age[y-1,length(ages),s,sim])) {
              
              # Calculate fishing mortality for the plus group
              fleet_mort_plus <- sum(fish_mort[y,,sim] * Fish_selex_at_age[y,a,,s,sim], na.rm = TRUE)
              
              # Applying mortality to plus group individuals last year, and add in recently recruited indviduals into the plus group
              N_at_age[y+1,a,s,sim] <-  N_at_age[y,a,s,sim] * exp(-(Mort_at_age[y,a,sim] + fleet_mort_plus)) +  N_at_age[y+1,a,s,sim] 
              
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
          
          ### Fishery fleet loop ------------------------------------------------------
          
          for(f in 1:n_fish_fleets) { # Loop for fishery fleets
            
            for(s in 1:n_sex) {
              # Calculate total mortality here
              if(n_fish_fleets > 1) { # If > 1 fishery fleet
                Z_s <- 0 # Reset Z here
                for(fl in 1:n_fish_fleets) {
                  Z_s_tmp <- fish_mort[y,fl,sim] * Fish_selex_at_age[y,,fl,s,sim] # get FAA here
                  Z_s <- Z_s_tmp + Z_s # increment
                  if(fl == n_fish_fleets) Z_s <- Z_s + Mort_at_age[y,,sim]
                } # end fl loop
                
                # Z_s = fish_mort[y,1,sim] * Fish_selex_at_age[y,,1,s,sim] +
                # fish_mort[y,2,sim] * Fish_selex_at_age[y,,2,s,sim] + Mort_at_age[y,,sim]
                
              } else{
                Z_s <- (fish_mort[y,f,sim] * Fish_selex_at_age[y,,f,s,sim]) +  Mort_at_age[y,,sim]
              } # if only 1 fishery fleet
              
              ###  Get Catch at Age (Only F to C for now) -----------------------------------
              
              # Calculate instantaneous fishing mortality for a given fleet, sex, and age
              Fish_Fleet_Mort <- fish_mort[y,f,sim] * Fish_selex_at_age[y,,f,s,sim]
              
              # Now, get catch at age in numbers
              Catch_at_age[y,,f,s,sim] <- N_at_age[y,,s,sim] * (1-exp(-Z_s)) * (Fish_Fleet_Mort / Z_s)
              
              if(s == n_sex) { # Get catch aggregated across ages and sexes, and add lognormal errors
                Catch_agg[y, f, sim] <- sum(Catch_at_age[y,,f,,sim] * wt_at_age[y,,,sim]) * 
                  exp( rnorm(1, 0, sqrt(log(catch_CV^2 + 1))) )
              } # if we are done w/ looping through sexes
              
              ### Sample Fishery Index and Comps ------------------------------------------
              
              # Only start sampling if y > Fish start year. 
              if(y >= Fish_Start_yr[f]) { # Observation Model for Fishery
                
                # Generate a fishery index structured by fleet and sex (numbers based)
                Fishery_Index[y,f,s,sim] <- q_Fish[y,f,sim]  *  sum(N_at_age[y,,s,sim] * wt_at_age[y,,s,sim] *  
                                                                      Fish_selex_at_age[y,,f,s,sim])
                
                # Probability for fishery age comps, using CAA as probability
                Prob_Fish_Comps <- Catch_at_age[y,,f,s,sim]
                
                # Generate comps based on catch at age
                Fish_Age_Comps[y,,f,s,sim] <- sample_comps(error = fish_likelihood,
                                                           Input_N = Input_N_Fish[y,f], 
                                                           DM_Param = DM_Fish_Param, 
                                                           prob = Prob_Fish_Comps / sum(Prob_Fish_Comps))
                
              }  # Only start sampling if we are the start of the fish start year
              
            } # end sex index
            
            # Summarize this fishery index aggregated by sex and applying some error
            Fishery_Index_Agg[y,f,sim] <- sum(melt(Fishery_Index[y,f,,sim]), na.rm = TRUE) # Aggregate
            
            # Apply error here, index fish_CV vector
            Fishery_Index_Agg[y,f,sim] <- idx_obs_error(error = "log_normal", 
                                                        true_index = Fishery_Index_Agg[y,f,sim],
                                                        CV = fish_CV[f])
          } # end fishery fleet index and loop
          
          # Survey Index and Comps --------------------------------------------------
          
          ### Survey fleet loop -------------------------------------------------------
          
          for(sf in 1:n_srv_fleets) { # Loop for survey fleets
            
            for(s in 1:n_sex) {
              
              # Only start sampling if y > Survey Start Year.
              if(y >= Surv_Start_yr[sf]) { 
                
                # Get survey index here (numbers based)
                Survey_Index[y,sf,s,sim] <- q_Surv[y,sf,sim]  * sum(N_at_age[y,,s,sim] *
                                                                      Surv_selex_at_age[y,,sf,s,sim])
                
                # Get probability of sampling a given age class for use in multinomial
                Prob_Surv_at_age <- (N_at_age[y,,s,sim] * Surv_selex_at_age[y,,sf,s,sim])
                
                # Generate comps based on the expected CPUE at age
                Survey_Age_Comps[y,,sf,s,sim] <- sample_comps(error = srv_likelihood,
                                                              Input_N = Input_N_Srv[y,sf], 
                                                              DM_Param = DM_Srv_Param, 
                                                              prob = Prob_Surv_at_age / sum(Prob_Surv_at_age))
                
              } # Only start sampling if we are at the start of the survey start year
              
            } # end sex loop for survey here
            
            # Summarize this fishery index aggregated by sex and applying some error
            Survey_Index_Agg[y,sf,sim] <- sum(melt(Survey_Index[y,sf,,sim]), na.rm = TRUE) # Aggregate
            
            # Apply error here, index srv_CV vector
            Survey_Index_Agg[y,sf,sim] <- idx_obs_error(error = "log_normal", 
                                                        true_index = Survey_Index_Agg[y,sf,sim],
                                                        CV = srv_CV[sf])
            
          } # end sf loop
          
        } # end check equilibrium loop (not sampling when checking equilibrium)
        
      } # if we are no longer in the first year
      
    } # end year loop
    
  } # end simulation loop
  
  # Return objects into environment
  N_at_age <<- N_at_age
  Biom_at_age <<- Biom_at_age
  SSB <<- SSB
  rec_total <<- rec_total
  Catch_at_age <<- Catch_at_age
  Catch_agg <<- Catch_agg
  Fishery_Index <<- Fishery_Index
  Fish_Age_Comps <<- Fish_Age_Comps
  Fishery_Index_Agg <<- Fishery_Index_Agg
  Survey_Index <<- Survey_Index
  Survey_Age_Comps <<- Survey_Age_Comps
  Survey_Index_Agg <<- Survey_Index_Agg
  rec_devs <<- rec_devs
  
  om_list <- list(N_at_age = N_at_age, Biom_at_age = Biom_at_age,
                  SSB = SSB, rec_total = rec_total, Catch_at_age = Catch_at_age, Catch_agg = Catch_agg,
                  Fishery_Index = Fishery_Index, Fish_Age_Comps = Fish_Age_Comps,
                  Fishery_Index_Agg = Fishery_Index_Agg, Survey_Index = Survey_Index,
                  Survey_Age_Comps = Survey_Age_Comps, Survey_Index_Agg = Survey_Index_Agg,
                  rec_devs = rec_devs)
  
} # end function