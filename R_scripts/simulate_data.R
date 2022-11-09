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
fish_surv_data_scenarios(Fish_Start_yr = 180, Surv_Start_yr = 180, n_years = n_years,
                           Fish_freq = 1, Surv_freq = 1, Scenario = "High")

# Specify selectivity parameterizations here
specify_selex(Fish_Funct_Form = "Asymp", Surv_Funct_Form = "Asymp", ages = ages,
              Fish_Time = "Constant", a50_fish = 7, k_fish = 1, a50_surv = 3, k_surv = 1)

# Specify Natural Mortality
specify_nat_mort(Mort_Time = "Constant", Mean_M = Mean_M)

# Specify catchability for the fishery and survey
specify_q(q_Fish_Time = "Constant", q_Surv_Time = "Constant", 
          q_Mean_Fish = 0.01, q_Mean_Surv = 0.01, 
          q_Fish_rho = 0.5, q_Fish_sigma = 0.005)

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
  print(paste("### Simulation",sim,"out of", n_sims, "###"))
  
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


