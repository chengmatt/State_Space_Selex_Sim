# Purpose: To collate results from simulation study - meant to be used
# after OM-EM models have been run! 
# Within this script, we have the following functions:
# 1) get_sdrep_values (to get sdrep values), 2) get_results (to get results from models - uses
# the get_sdrep_values function)
# Creator: Matthew LH. Cheng
# Date 11/4/22


# Get sdrep values --------------------------------------------------------

#' @param parameter Parameter name from sdreport
get_sdrep_values <- function(sd_rep = sd_rep, parameter) {
  
  idx <- which(rownames(sd_rep) == parameter) # find out where in the sdrep the 
  # parameter is
  
  if(parameter == "log_FAA_tot") { # Fishing mortality by age
    
    # The code here is adapted from WHAM - from plot.SSB.F.trend()
    
    # Stick these values into a matrix (nrow = years, length = ages)
    log_faa <- matrix(sd_rep[idx,1], length(Fish_Start_yr:(n_years-1)), length(ages))
    
    # Get CIs for the above
    faa_ci <- matrix(sd_rep[idx,2], length(Fish_Start_yr:(n_years-1)), length(ages))
    
    # Get max F - fully selected F 
    full_selex_F_index <- apply(log_faa,1, function(x) max(which(x == max(x)))) 
    
    # Now, get fully selected F values and associated SEs
    # Create empty vectors to store vlaues in
    value <- vector()
    se <- vector()
    
    # Loop through this because we might get years where the fully selected F differs by yrs
    for(y in 1:length(full_selex_F_index)) {
      
      # Get year index for fully selected F
      f_index <- full_selex_F_index[y]
      
      # search for the yth index and put into our vec
      value[y] <- log_faa[y,f_index]
      se[y] <- faa_ci[y,f_index]
      
    } # end y loop

  } else{ 
    value <- sd_rep[idx,1] # Get estimated values from sdrep
    se <- sd_rep[idx,2] # get standard errors
  } # else - not log_FAA_tot
  
  return(list(value, se)) # return values
  
} # function to get values from sdrep


# Get general results from sdrep -----------------------------------------------------

#' @param model Model object estimated from WHAM

get_results <- function(model) {
  
  require(tidyverse)
  
  # Get report file here from model
  report <- model$report()
  
  # Get sdreport file here from model
  sd_rep <- summary(model$sdrep)
  
# SSB, F, and q -----------------------------------------------------------

  # Get SSB values from sdrep
  ssb_values <- get_sdrep_values(sd_rep = sd_rep, parameter = "log_SSB")
  
  # Stick these into a dataframe and calculate CIs
  SSB_rep <- data.frame(Year = Fish_Start_yr:(n_years-1), log_SSB = ssb_values[[1]], SE = ssb_values[[2]]) %>% 
    mutate(SSB = exp(log_SSB),
           lwr = exp(log_SSB - (1.96*SE)),
           upr = exp(log_SSB + (1.96*SE)))
  
  # Get fishing mortality rates
  f_values <- get_sdrep_values(sd_rep = sd_rep, parameter = "log_FAA_tot")
  F_rep <- data.frame(Year = Fish_Start_yr:(n_years-1), log_F = f_values[[1]], SE = f_values[[2]]) %>% 
    mutate(F_val = exp(log_F),
           lwr = exp(log_F - (1.96*SE)),
           upr = exp(log_F + (1.96*SE)))

# Selectivity -------------------------------------------------------------
  
  # Get fishery selectivity
  fish_sel <- melt(report$selAA[[1]]) %>% 
    rename(Year = Var1, Age = Var2, Selex = value) %>% 
    mutate(Type = "Fishery")
  
  # Get survey selectivity
  surv_sel <- melt(report$selAA[[2]]) %>% 
    rename(Year = Var1, Age = Var2, Selex = value) %>% 
    mutate(Type = "Survey")
  
  selex_rep <- rbind(fish_sel, surv_sel) # bind these together!
  
  return(list(SSB = SSB_rep, F = F_rep, Selex = selex_rep))
}
