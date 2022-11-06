# Purpose: To specify F pattern scenarios for the fishery
# Creator: Matthew LH. Cheng
# Date: 10/31/22

#' @param Fish_Start_yr Start fishing year
#' @param Start_F Start value for F mort
#' @param F_type Type of fishing mortality pattern (Contrast, Constant, Increase, and Increase_Plat)
#' @param F_sigma_dev If we want to add some noise around the fishing mortality pattern, default
#' is 0. 

specify_F_pattern <- function(Fish_Start_yr, Start_F, F_type, F_sigma_dev = 0) {
  
  # Get mean natural mortality rate
  mean_nat_mort <- mean(Mort_at_age[,,])
  
  # Set fishing mortality to be 0 before the fishery starts
  fish_mort[1:(Fish_Start_yr-1),] <- 0
  
  if(F_type == "Contrast") {
    
    # Determine the mid point of when to ramp back down
    midpoint <- round(mean(Fish_Start_yr:n_years))
    
    for(sim in 1:n_sims) {
      
      # Increase fishing mortality (ramp up); here we are specifying fishing mortality relative to natural mortality.
      # We will allow it go over natural mortality first by 1.5x
      fish_mort[Fish_Start_yr:midpoint,sim] <- abs(
        seq(Start_F, (1.5 * mean_nat_mort), length.out = length(Fish_Start_yr:midpoint)) +
          # Add some normal random error around this
          rnorm(length(Fish_Start_yr:Fish_Start_yr:midpoint), 0, F_sigma_dev)
      ) # absolute values so they're always positive
      
      # Now, decrease fihsing mortality (ramp down) - decrease to 0.75 of the mean natural mortality
      fish_mort[midpoint:n_years,sim] <- abs(
        seq((1.5 * mean_nat_mort), (0.75 * mean_nat_mort), length.out = length(midpoint:n_years)) +
          # Add some noraml random error around this
          rnorm(length(midpoint:n_years), 0, F_sigma_dev)
      ) # absolute values so they're always positive
      
    } # end sim loop

  } # if statement for F pattern that is contrasting
  
  if(F_type == "Constant") {
    
    # does not use the Start_F argument - automatically scales it up to 0.75 of the natural mortality 
    
    for(sim in 1:n_sims) { 
      
      # Constant fishing mortality - 0.75 of natural mortality
      fish_mort[Fish_Start_yr:n_years,] <- abs(
        rep((0.75 * mean_nat_mort), times = length(Fish_Start_yr:n_years)) +
          # Add some normal random error around this
          rnorm(length(Fish_Start_yr:n_years), 0, F_sigma_dev)
      ) # Absolute values so they're always positive
      
    } # end sim loop
    
  } # if statement for F pattern that is Constant
  
  if(F_type == "Increase") {
    
    for(sim in 1:n_sims) {
      
      fish_mort[Fish_Start_yr:n_years,sim] <- abs(
        seq(Start_F, (1.5 * mean_nat_mort), length.out = length(Fish_Start_yr:n_years)) +
          # Add some normal random error around this
          rnorm(length(Fish_Start_yr:n_years), 0, F_sigma_dev)
      ) # Making absolute so they're always positive numbers
      
    } # end sims loop
    
  } # if statement for F pattern that is increasing
  
  
  if(F_type == "Increase_Plat") { # if increas and plateaus out
    
    # Calculate plateau year
    plat_year <- round(0.8 * n_years)  # make it 0.5 of the number of years in the simulation
    
    for(sim in 1:n_sims) {
      
      # Calculate increasing f mort
      increasing <- abs(seq(Start_F, (0.75 * mean_nat_mort), length.out = length(Fish_Start_yr:plat_year)) +
          rnorm(length(Fish_Start_yr:plat_year), 0, F_sigma_dev) # Add some normal random error around this
      ) # Making absolute so they're always positive numbers
      
      # Get vector of plateau fmorts
      plat <- rep((0.75 * mean_nat_mort), length.out = length((plat_year+1):n_years)) + 
        rnorm(length((plat_year+1):n_years), 0, F_sigma_dev)
      
      fish_mort[Fish_Start_yr:n_years,sim] <- c(increasing, plat)
      
    } # end sims loop
    
  } # if statement for F pattern that is increasing
  
  # Output into environment
  fish_mort <<- fish_mort
  
  print(paste("Fishing mortality specified as:" , F_type))
  
} # end function



