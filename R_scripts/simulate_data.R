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
source(here("R_scripts", "functions", "sample_index.R")) # Generate index data
source(here("R_scripts", "functions", "sample_age_comps.R")) # Generate age comp data
source(here("R_scripts", "functions", "make_input.R")) # Put stuff into WHAM format
source(here("R_scripts", "functions", "get_results.R")) # Get WHAM data from model

# Path to general input biological parameters
spreadsheet_path <- here("input", "Sablefish_Inputs.xlsx")

# Set up OM --------------------------------------------------------

# Read in excel sheet parameters and create OM objects to store values in
read_params_create_OM_objects(spreadsheet_path = spreadsheet_path)

# Specify data scenarios here
fish_surv_data_scenarios(Fish_Start_yr = 190, Surv_Start_yr = 190, n_years = n_years,
                           Fish_freq = 1, Surv_freq = 1, Scenario = "High")

# Specify selectivity parameterizations here
specify_selex(Fish_Funct_Form = "Asymp", Surv_Funct_Form = "Asymp", ages = ages,
              Fish_Time = "Constant", a50_fish = 7, k_fish = 1, a50_surv = 3, k_surv = 1)

# Specify Natural Mortality
specify_nat_mort(Mort_Time = "Constant", Mean_M = Mean_M)

# Specify catchability for the fishery and survey
specify_q(q_Fish_Time = "Constant", q_Surv_Time = "Constant", 
          q_Mean_Fish = 0.035, q_Mean_Surv = 0.01)

# Specify fishing mortality pattern
specify_F_pattern(Fish_Start_yr = Fish_Start_yr, F_type = "Contrast", Start_F = 0.001, F_sigma_dev = 0)

# Specify recruitment deviates here - loops though each simulation
specify_rec_devs(Rec_Dev_Type = "iid", rho_rec = NA) 

check_equil <- FALSE

# Start Simulation to get to equilibrium conditions --------------------------------------------------------

# Start simulation loop here
for(sim in 1:n_sims) {
  
# Check equilibrium -------------------------------------------------------

  if(check_equil == TRUE) {
    
    # Turn fishing off
    Fish_selex_at_age[,,] <- 0 # Selectivity
    fish_mort[,] <- 0
    Fish_yrs <- NA # no fishing occuring
    Fish_Start_yr <- NA # no fishing occuring
    Surv_yrs <- NA # no survey occuring
    Surv_Start_yr <- NA # no survey occuring
    
    print("### Checking whether equilibrium conditions have been met ###")
    
  } # checking equilibrium conditions

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
      SSB[y,sim] <- sum(mat_at_age[y,,sim] * Biom_at_age[y,1,sim], na.rm = TRUE)

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
                                 exp(-(Mort_at_age[y-1,a,sim] + (fish_mort[y-1,sim] * Fish_selex_at_age[y-1,a,sim])))
                                     # Z = M + F
          
          # Now, add in the recruits generated from the previous year to the numbers at age matrix to recruit age
          N_at_age[y,1,sim] <- rec_total[y-1,sim]
          
        } # if we are not in the plus group
        
     ### Decrement population with F and Z (== + group) ---------------------------------------
        
        if(a == length(ages) & !is.na(N_at_age[y-1,length(ages),sim])) {
          
          # Decrement the plus group here by natural mortality and fishing mortality from year y - 1,
          # and add in back into the plus group at year y (the age bin previous to the plus group 
          # has already been decremeneted and added into the plus group here )
          
          # Applying natural and fishing mortality to indivduals from the previous year  that were in the plus group
          N_at_age[y,length(ages),sim] <-  (N_at_age[y-1,length(ages),sim] * exp(-(Mort_at_age[y-1,a,sim] + (fish_mort[y-1,sim] * Fish_selex_at_age[y-1,a,sim])))) + # Z = M + F
                                            N_at_age[y,length(ages),sim] # Adding back the individuals that recently recruited into + group
          
                                           
        } # if we are in the plus group and we had a plus group in the previous year (so we are not accidentally
        # adding additional plus groups together)
        
      } # ages loop
      
    ###  Get Catch at Age (Only F to C for now) -----------------------------------
      
        # Calculate Z at age - total mortality
        Z_at_age <- (Mort_at_age[y-1,,sim] + (fish_mort[y-1,sim] * Fish_selex_at_age[y-1,,sim]))
        
        # Calculate our proportion of mortality via fishing
        Prop_fish_mort <- (fish_mort[y-1,sim] * Fish_selex_at_age[y-1,,sim]) / Z_at_age
        
        # Now, get catch at age
        Catch_at_age[y-1,,sim] <- Prop_fish_mort * N_at_age[y-1,,sim] * (1-exp(-Z_at_age)) * wt_at_age[y,,sim]
        
# Observation Model (Survey and Fishery) ----------------------------------

    ### Fishery -----------------------------------------------------------------
      
        # Only start sampling if y > Fish start year and is in the specified fish years 
        # We only sample when y > Fish_Start_yr because we need to wait for the lag in Fmort to catch up
        if(y > Fish_Start_yr & y %in% c(Fish_yrs)) { # Observation Model for Fishery
          
          # Generate a fishery index using our function w/ normal error (numbers based)
          Fishery_Index[y-1,sim] <- sample_index(Idx_Fleet = "Fishery", error = "log_normal")

          # Generate comps based on the expected catch at age
          Fish_Age_Comps[y-1,,sim] <- sample_age_comps(Comp_Fleet = "Fishery", error = "multinomial")
          
        }  # Only start sampling if we are the start of the fish start year
      

    ### Survey ------------------------------------------------------------------
      
      # Only start sampling if y >= Survey start year + 1 (so it starts at the
      # correct yr ) and is in the specified survey years
      # Note that the survey fleet removals are not decremeneted in the population (we assume that they all
      # survive after being indexed.
      
      if(y >= (Surv_Start_yr+1) & y %in% c(Surv_yrs)) {

        # Get survey index here (numbers based) 
        Survey_Index[y-1,sim] <- sample_index(Idx_Fleet = "Survey", error = "log_normal")
        
        # Generate comps based on the expected CPUE at age
        Survey_Age_Comps[y-1,,sim] <- sample_age_comps(Comp_Fleet = "Survey", error = "multinomial")
        
      } # Only start sampling if we are at the start of the survey start year
        
        
        ### Update Biomass values and Numbers + Generate Recruits -------------------
        
        # Update Biomass at age 
        Biom_at_age[y,,sim] <- N_at_age[y,,sim] * wt_at_age[y,,sim]
        
        # Now, update SSB for year y and generate recruits for the new year with the updated SSB
        SSB[y,sim] <- sum(mat_at_age[y,,sim] * Biom_at_age[y,,sim], na.rm = TRUE)
          
          # Now generate new recruits with the updated SSB
          rec_total[y,sim] <- beverton_holt_recruit_new(ssb = SSB[y,sim], h = h, r0 = r0, ssb0 = ssb0) *
                                                        exp(rec_devs[y,sim] - ((sigma_rec^2)/2))
          
    } # if we are no longer in the first year

  } # end year loop
  
} # end simulation loop


# WHAM checks ------------------------------------------------------------- (changed recruit var to be lower)

# Set modelling structure here
selectivity <- list(model = rep("logistic", 2))

M <- list(model = "constant", initial_means = 0.1, est_ages = 1) # need to specify est_ages for M to be estimated!

catchability <- list(re = "none") # put bounds on q

for(sim in 1:n_sims) {
  
  # Force our inputs into a list - so that it reads into wham -
  basic <- make_input(n_fleets = 1, n_indices = 2, Catch_CV_Val = 0.05, catch_error = FALSE,
                      n_sims = sim, bias_obs = TRUE, bias_process = TRUE) 
  
  # Make WHAM inputs here 
  test_wham <- wham::prepare_wham_input(basic_info = basic, selectivity = selectivity, 
                                        recruit_model = 2, M = M)
  
  # Fit WHAM here
  em_fit <- wham::fit_wham(input = test_wham, do.fit = T, do.osa = F, do.retro = F,
                           save.sdrep = TRUE, do.check = T)
  
} # end loop of number of simulations we want to run

# Debugging
em_fit$sdrep
rep <- em_fit$report()
sapply(grep("nll",names(rep),value=T), function(x) sum(rep[[x]]))


# Fishery selectivity
fish_sel_df <- melt(Fish_selex_at_age)
names(fish_sel_df) <- c("Year", "Age", "Sim", "Selex")
fish_sel_df <- fish_sel_df %>% 
  group_by(Age) %>% 
  summarize(Selex = mean(Selex)) %>% 
  mutate(Type = "Fishery",
         Age = parse_number(as.character(Age))-1)

surv_sel_df <- melt(Surv_selex_at_age)
names(surv_sel_df) <- c("Year", "Age", "Sim", "Selex")
surv_sel_df <- surv_sel_df %>% 
  group_by(Age) %>% 
  summarize(Selex = mean(Selex)) %>% 
  mutate(Type = "Survey",
         Age = parse_number(as.character(Age))-1)
  
# Output plots
# Create directory to ouptut plots to
wham_out <- here("figs", "wham_checks")
dir.create(wham_out)
plot_wham_output(em_fit, dir.main = wham_out)

# Check other stuff
results <- get_results(em_fit)

# Plot ssb
ssb_df <- results[[1]]

# SSB_TRUTH
ssb_truth <- c(SSB[190:200,sim])
ssb_df <- ssb_df %>% 
  mutate(SSB_Truth = ssb_truth)

ggplot(ssb_df, mapping = aes(x = Year, y = SSB, ymin = lwr, ymax = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  geom_line(mapping = aes(x = Year, y = SSB_Truth), col = "red") 

# Plot F
f_df <- results[[2]]

# SSB_TRUTH
f_truth <- c(fish_mort[190:200,20])
f_df <- f_df %>% 
  mutate(f_truth = f_truth)

ggplot(f_df, mapping = aes(x = Year, y = F_val, ymin = lwr, ymax = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  geom_line(mapping = aes(x = Year, y = f_truth), col = "red")



truth_selex <- rbind(fish_sel_df, surv_sel_df) %>% 
  rename(Selex_Truth = Selex)
# Plot selex
selex_df <- results[[3]] %>% 
  group_by( Type, Age) %>% 
  summarize(selex = mean(Selex)) %>% 
  left_join(truth_selex, by  = c("Type", "Age"))

ggplot(selex_df, aes(x = Age, y = selex)) +
  geom_line() +
  geom_line(selex_df, mapping = aes(x = Age, y = Selex_Truth), col = "red") +
  facet_wrap(~Type, scales = "free")


# plot_OM(path = here("figs", "Base_OM_Figs"), file_name = "OM_Check.pdf")
