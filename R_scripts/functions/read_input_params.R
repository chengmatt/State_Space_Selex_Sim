# Purpose: To read in input parameters for the OM from an excel spreadsheet. The function calls another
# function to which creates OM objects to store values within it.
# Creator: Matthew LH. Cheng
# Date: 10/30/22

#' @param spreadsheet_path Path to input spreadsheet

read_params_create_OM_objects <- function(spreadsheet_path) {
  
  require(readxl)
  require(tidyverse)
  require(here)
  
  source(here("R_scripts", "functions", "create_OM_objects.R")) # Create OM objects 

# Controls ----------------------------------------------------------------

  ctl <- read_xlsx(spreadsheet_path, sheet = "Controls")
  
  # Read in and save objects in our environment
  n_sims <<- ctl$Value[ctl$Par == "n_sims"] # Number of simulations
  n_years <<- ctl$Value[ctl$Par == "n_years"] # Number of years
  N_1 <<- ctl$Value[ctl$Par == "N_1"] # Numbers at year 1 in first age of recruitment
  n_sex <<- ctl$Value[ctl$Par == "n_sex"] # Numbers of sexes
  n_fish_fleets <<- ctl$Value[ctl$Par == "n_fish_fleets"] # Numbers of fishery fleets
  n_srv_fleets <<- ctl$Value[ctl$Par == "n_srv_fleets"] # Numbers of survey fleets

# Age Bins ----------------------------------------------------------------

  age_pars <- read_xlsx(spreadsheet_path, sheet = "Age_Bins")
  
  # Read in ages
  ages <<- age_pars$ages
  
  # Create OM objects here
  create_OM_objects(n_years = n_years, n_sims = n_sims, ages = ages, n_sex = n_sex,
                    n_fish_fleets = n_fish_fleets)

# Maturity at age ---------------------------------------------------------

  maturity <- read_xlsx(spreadsheet_path, sheet = "Maturity_At_Age")
  
  # Munge maturity at age - pivot longer
  maturity_long <- maturity %>% 
    pivot_longer(!c(Year, Time, Sex), names_to = "Ages", values_to = "Maturity")
  
  # If this is time-invariant 
  if(unique(maturity_long$Time) == "Time_Inv") {
    
    # Fill this in
    for(s in 1:length(unique(maturity_long$Sex))) {
      
      for(y in 1:nrow(mat_at_age)) {
        
        mat_at_age[y,,s,] <-  maturity_long$Maturity[maturity_long$Sex == s]
        
      } # end s sex loop

    } # end i
    
  } # end time invariant
  
  # If maturity is time-varying
  if(unique(maturity_long$Time) == "Time_Var") {
    
    if(length(unique(maturity_long$Year)) != nrow(mat_at_age)) {
      stop("Years of Maturity At Age in Excel Sheet != the number of rows in the array")
    } 
    
    # Fill this in
    for(s in 1:length(unique(maturity_long$Sex))) {
      
      for(y in 1:nrow(mat_at_age)) {
        
        mat_at_age[y,,s,] <-  maturity_long$Maturity[maturity_long$Sex == s & maturity_long$Year == y]
        
      } # end s sex loop
      
    } # end i
    
  } # end time varying
  
  mat_at_age <<- mat_at_age # Output this to environment
  
# Weight at Age -----------------------------------------------------------

  weight <- read_xlsx(spreadsheet_path, sheet = "Weight_At_Age")
  
  # Munge maturity at age - pivot longer
  weight_long <- weight %>% 
    pivot_longer(!c(Year, Time, Sex), names_to = "Ages", values_to = "Weight")
  
  # If this is time-invariant 
  if(unique(weight_long$Time) == "Time_Inv") {
    
    # Fill this in
    for(s in 1:length(unique(weight_long$Sex))) {
      
      for(y in 1:nrow(wt_at_age)) {
        
        wt_at_age[y,,s,] <-  weight_long$Weight[weight_long$Sex == s]
        
      } # end s sex loop
      
    } # end i
    
  } # end time invariant
  
  # If maturity is time-varying
  if(unique(weight_long$Time) == "Time_Var") {
    
    if(length(unique(weight_long$Year)) != nrow(wt_at_age)) {
      stop("Years of Weight At Age in Excel Sheet != the number of rows in the array")
    } 
    
    # Fill this in
    for(s in 1:length(unique(weight_long$Sex))) {
      
      for(y in 1:nrow(wt_at_age)) {
        
        wt_at_age[y,,s,] <-  weight_long$Weight[weight_long$Sex == s & weight_long$Year == y]
        
      } # end s sex loop
      
    } # end i
    
  } # end time varying
  
  wt_at_age <<- wt_at_age # Output to environment
  
  
# Recruitment + Mortality-------------------------------------------------------------
  
  recruitment_pars <- read_xlsx(spreadsheet_path, sheet = "Recruitment_Mortality")
  h <<- recruitment_pars$Value[recruitment_pars$Par == "h"] # Steepness (Recruitment at 20% of SSB0)
  r0 <<- recruitment_pars$Value[recruitment_pars$Par == "r0"] # Virgin Recruitment
  ssb0 <<- recruitment_pars$Value[recruitment_pars$Par == "ssb0"] # Virgin SSB
  sigma_rec <<- recruitment_pars$Value[recruitment_pars$Par == "sigma_rec"] # Recruitment variability
  mu_rec <<- recruitment_pars$Value[recruitment_pars$Par == "mu_rec"] # Mean recruitment
  Mean_M <<- recruitment_pars$Value[recruitment_pars$Par == "M"] # Mortality
  

  print("### Input parameters have been read in and OM objects have been created ###")

} # end function
