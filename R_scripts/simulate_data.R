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

# Path to general input biological parameters
spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")

# Set up OM --------------------------------------------------------

# Read in excel sheet parameters and create OM objects to store values in
read_params_create_OM_objects(spreadsheet_path = spreadsheet_path)

# Specify data scenarios here
fish_surv_data_scenarios(Fish_Start_yr = 125, Surv_Start_yr = 150, n_years = n_years,
                         Fish_freq = 1, Surv_freq = 1, Scenario = "Low")

# Specify selectivity parameterizations here
specify_selex(Fish_Funct_Form = "Asymp", Surv_Funct_Form = "Asymp", ages = ages,
              Fish_Time = "Constant", a50_fish = 4, k_fish = 3, a50_surv = 4, k_surv = 2)

# Specify Natural Mortality
specify_nat_mort(Mort_Time = "Constant", Mean_M = 0.1)

# Specify catchability for the fishery and survey
specify_q(q_Fish_Time = "Constant", q_Surv_Time = "Constant", q_Mean_Fish = 0.2, q_Mean_Surv = 0.03)

# Set 0 to fmort init
fish_mort[1:(Fish_Start_yr-1),] <- 0

# Increase
fish_mort[(Fish_Start_yr:(n_years-25)),] <- seq(0.005, 0.3, length.out = length(Fish_Start_yr:(n_years-25)))

# Ramp down
fish_mort[(n_years-25):n_years,] <- seq(0.15, 0.1, length.out = length((n_years-25):n_years))

fish_mort_df <- melt(fish_mort)
ggplot(fish_mort_df, aes(x = as.numeric(Var1), y = value, color = Var2)) +
  geom_line()

# True F to C
F_to_C <- TRUE

check_equil <- TRUE


# Start Simulation to get to equilibrium conditions --------------------------------------------------------

# Start simulation loop here
for(sim in 1:n_sims) {
  
# Check equilibrium -------------------------------------------------------

  if(check_equil == TRUE) {
    
    # Turn recruitment variation off
    sigma_rec <- 0
    # Turn fishing off
    Fish_selex_at_age[,,] <- 0 # Selectivity
    fish_mort[,] <- 0
    Fish_yrs <- NA # no fishing occuring
    Fish_Start_yr <- NA # no fishing occuring
    
  } # checking equilibrium conditions
  
  # Generate new rec devs for every simulation
  rec_devs[,] <- rnorm(n_years, 0, sigma_rec)
  
  # Print simulation iteration
  print(paste("Simulation",sim,"out of", n_sims))

# Years loop  -------------------------------------------------------------

  for(y in 1:n_years) {
    
    #### Initialize Population here ----------------------------------------------

    if(y == 1) { 
      # Initialize the population here first. We are going to seed the population with age-2s at 1600
      # Put that into our N_at_age array
      N_at_age[y,1,sim] <- N_1
      
      # Update Biomass at age 
      Biom_at_age[y,1,sim] <- N_at_age[y,1,sim] * wt_at_age[y,1,sim]

      # Now, calculate our SSB in the first year
      SSB[y,sim] <- sum(mat_at_age[y,,sim] * wt_at_age[y,,sim] * N_at_age[y,,sim], na.rm = TRUE)
      
      # Now, calculate the number of recruits we get - this returns abundance - N at age 2
      rec_total[y,sim] <- beverton_holt_recruit_new(ssb = SSB[y,sim], h = h, r0 = r0, ssb0 = ssb0) *
                                                                 exp(rec_devs[y,sim] - ((sigma_rec^2)/2)) # Add lognormal correction and recdevs

    }  # if we are in the first year of the simulation
    
    if(y != 1) {
      
# Ages Loop ---------------------------------------------------------------

      for(a in 1:length(ages)) {

    ### Decrement population with F and Z (!= + group) ---------------------------------------
        
        if(a != length(ages)) {
          
          # Apply a natural mortality rate and fishing mortality rate to decrement the population (Z = M + F)
          N_at_age[y,a+1,sim] <- N_at_age[y-1,a,sim] * # Indexing our numbers at age
                                 exp(-(Mort_at_age[y-1,a,sim] + (fish_mort[y-1,sim] * Fish_selex_at_age[y-1,a,sim]))) # Z = M + F
          
          # Now, add in the recruits generated from the previous year to the numbers at age matrix to recruit age
          N_at_age[y,1,sim] <- rec_total[y-1,sim]
          
        } # if we are not in the plus group
        
     ### Decrement population with F and Z (== + group) ---------------------------------------
        
        if(a == length(ages) & !is.na(N_at_age[y-1,length(ages),sim])) {
          
          # Decrement the plus group here by natural mortality and fishing mortality from year y - 1,
          # and add in back into the plus group at year y (the age bin previous to the plus group 
          # has already been decremeneted and added into the plus group here )
          
          # Applying natural and fishing mortality to indivduals from the previous year  that were in the plus group
          N_at_age[y,length(ages),sim] <- (N_at_age[y-1,length(ages),sim] * exp(-(Mort_at_age[y-1,a,sim] + (fish_mort[y-1,sim] * Fish_selex_at_age[y-1,a,sim])))) + # Z = M + F
                                             N_at_age[y,length(ages),sim] # Adding back the individuals that recently recruited into + group
                                           
        } # if we are in the plus group and we had a plus group in the previous year (so we are not accidentally
        # adding additional plus groups together)
        
      } # ages loop
      

    ### Update Values + Generate New Recruits + Get Catch at Age -----------------------------------
      
      if(F_to_C == TRUE) { #  catch equation 
        
        # Calculate Z at age - total mortality
        Z_at_age <- (Mort_at_age[y-1,,sim] + (fish_mort[y-1,sim] * Fish_selex_at_age[y-1,,sim]))
        
        # Calculate our proportion of mortality via fishing
        Prop_fish_mort <- (fish_mort[y-1,sim] * Fish_selex_at_age[y-1,a,sim]) / Z_at_age
        
        # Now, get catch at age
        Catch_at_age[y-1,,sim] <- Prop_fish_mort * N_at_age[y-1,,sim] * (1-exp(-Z_at_age)) * wt_at_age[y,,sim]
        
      } # Going from F to C
      
      # Update Biomass at age (Note that during the last year, N_at_age, Biomass_at_age,
      # and SSB has not had impacts from natural mortality or fishing mortality - it represents
      # what's left in the terminal year)
      Biom_at_age[y,,sim] <- N_at_age[y,,sim] * wt_at_age[y,,sim]
      
      # Now, update SSB for year y and generate recruits for the new year with the updated SSB
      SSB[y,sim] <- sum(mat_at_age[y,,sim] * wt_at_age[y,,sim] * N_at_age[y,,sim], na.rm = TRUE)
      
      # Now generate new recruits with the updated SSB
      rec_total[y,sim] <- beverton_holt_recruit_new(ssb = SSB[y,sim], h = h, r0 = r0, ssb0 = ssb0) *
                                                          exp(rec_devs[y,sim] - ((sigma_rec^2)/2)) # Add lognormal correction and recdev
      
# Observation Model (Survey and Fishery) ----------------------------------

    ### Fishery -----------------------------------------------------------------
      
        # Only start sampling if y > Fish start year and is in the specified fish years 
        # We only sample when y > Fish_Start_yr because we need to wait for the lag in Fmort to catch up
        if(y > Fish_Start_yr & y %in% c(Fish_yrs)) { # Observation Model for Fishery
          
          # Get proportions at a given age from catch
          Prob_Catch_at_age <- Catch_at_age[y-1,,sim]/sum(Catch_at_age[y-1,,sim], na.rm = TRUE)
          
          # Generate a fishery index based on I = q * sum(N * selex * weight) * log normal error
          Fishery_Index_at_age[y-1,,sim] <- (q_Fish[y-1,sim] * N_at_age[y-1,,sim] * Fish_selex_at_age[y-1,,sim] * wt_at_age[y-1,,sim]) * 
            exp(rnorm(1, 0, Fishery_CV))

          # Generate comps based on the expected catch at age
          Fish_Age_Comps[y-1,,sim] <- rmultinom(n = 1, size = N_eff_fish, prob = Prob_Catch_at_age)
          
        }  # Only start sampling if we are the start of the fish start year
      

    ### Survey ------------------------------------------------------------------
      
      # Only start sampling if y >= Survey start year + 1 (so it starts at the
      # correct yr ) and is in the specified survey years
      # Note that the survey fleet removals are not decremeneted in the population (we assume that they all
      # survive after being indexed.
      
      if(y >= (Surv_Start_yr+1) & y %in% c(Surv_yrs)) {

        # Get survey index here I = q * sum(N * selex * weight) * log normal error
        Survey_Index_at_age[y-1,,sim] <- (q_Surv[y-1,sim] * N_at_age[y-1,,sim] * Surv_selex_at_age[y-1,,sim] * wt_at_age[y-1,,sim]) *
                                          exp(rnorm(1, 0, Survey_CV))
        
        # Sample from the expected proportions at age in the index
        Prob_Surv_at_age <- Survey_Index_at_age[y-1,,sim]/sum(Survey_Index_at_age[y-1,,sim])
        
        # Generate comps based on the expected CPUE at age
        Survey_Age_Comps[y-1,,sim] <- rmultinom(n = 1, size = N_eff_surv, prob = Prob_Surv_at_age)
        
      } # Only start sampling if we are at the start of the survey start year
      
    } # if we are no longer in the first year

  } # end year loop
  
} # end simulation loop


# Plot checks Conditions ---------------------------------------------

# Create output folder to visualize OM figures
dir.figs <- here("figs", "Base_OM_Figs")
dir.create(dir.figs)

plot_OM(path = dir.figs)



