# Purpose: To specify F pattern scenarios for the fishery
# Creator: Matthew LH. Cheng
# Date: 10/31/22


# Specify Fs --------------------------------------------------------------

#' @param Start_F Starting F vector
#' @param Fish_Start_yr Starting Fishing years for fleets
#' @param n_years Number of years in simulation
#' @param max_rel_F_M Maximum F value relative to the M
#' @param desc_rel_F_M Descending F value relative to the M
#' @param F_type Specify type of fishing mortality pattern
#' @param mean_nat_mort Mean natural mortality
#' @param yr_chng If we are simulating a fleet structure change, when do we want this change to occur at 
#' Note that this is only applicable for (Const_Inc and Contrast)

specify_F_pattern <- function(Start_F, Fish_Start_yr, F_type, n_years, max_rel_F_M, desc_rel_F_M,
                              mean_nat_mort, yr_chng = NULL) {
  
  if(F_type == "Contrast") { # Contrast F type
    
    if(desc_rel_F_M > max_rel_F_M) stop("Descending limb of fishing mortality contrast is larger than ascending limb! Specify another pattern or adjust the relative mortality rates to reflect a contrast pattern!")
    
    if(!is.null(yr_chng)) { # if this is a fleet sturcture change,
      # calculate when this change should occur
      chngpoint <- yr_chng
      print(paste("F mort contrast decreasing at user specified value:", chngpoint))
    } else{
      # Determine the mid point of when to ramp back down
      chngpoint <- round(mean(Fish_Start_yr:(n_years)))
      print(paste("F mort contrast decreasing at fxn calculated value:", chngpoint))
    } # else statement for change point
    
    # Increase fishing mortality (ramp up); here we are specifying fishing mortality relative to natural mortality.
    # We will allow it go over natural mortality first by 1.5x
    fish_mort_1 <- seq(Start_F, 
                       (max_rel_F_M * mean_nat_mort), 
                       length.out = length(Fish_Start_yr:chngpoint)) 
    
    # Now, decrease fihsing mortality (ramp down) - decrease to 0.75 of the mean natural mortality
    fish_mort_2 <- seq(((max_rel_F_M * mean_nat_mort)), 
                       (desc_rel_F_M * mean_nat_mort), 
                       length.out = length(chngpoint:n_years))
    
    # Now bind all of these fish morts into a vector
    F_vec <- c(fish_mort_1, fish_mort_2[-1]) # Remove first element so it connects smoothly
    
  } # if statement for F pattern that is contrasting
  
  if(F_type == "Constant") { # Constant F type
    
    # Constant fishing mortality - 0.75 of natural mortality
    F_vec <- rep(Start_F, times = length(Fish_Start_yr:n_years)) 
    
  } # if statement for F pattern that is Constant
  
  if(F_type == "Increase") { # Increasing F type
    
    F_vec <- seq(Start_F, (max_rel_F_M * mean_nat_mort), length.out = length(Fish_Start_yr:n_years))
    
  } # if statement for F pattern that is increasing

  if(F_type == "Increase_Plat") { # if increases and plateaus out
    
    # Calculate plateau year
    plat_year <- round(0.95 * n_years)  # make it 0.5 of the number of years in the simulation
    
    # Calculate increasing f mort
    increasing <- seq(Start_F, (max_rel_F_M * mean_nat_mort), length.out = length(Fish_Start_yr:(plat_year-1))) 
    
    # Get vector of plateau fmorts
    plat <- rep((max_rel_F_M * mean_nat_mort), length.out = length((plat_year):n_years))
    
    F_vec <- c(increasing, plat)
    
  } # if statement for F pattern that is increasing
  
  if(F_type == "Const_Inc") { # Constant and Increasing F type
    
    if(!is.null(yr_chng)) { # if this is a fleet sturcture change,
      # calculate when this change should occur
      chngpoint <- yr_chng
      print(paste("F mort Const_Inc increasing at user specified value:", chngpoint))
    } else{
      # Determine the mid point of when to ramp back down
      chngpoint <- round(mean(Fish_Start_yr:(n_years)))
      print(paste("F mort Const_Inc increasing at fxn calculated value:", chngpoint))
    } # else statement for change point
    
    # F stays constant initially
    constant_F <- seq(Start_F, Start_F, length.out = length(Fish_Start_yr:chngpoint))
    # Change point for F
    change_F <- seq(Start_F, (max_rel_F_M * mean_nat_mort), length.out = length(chngpoint:n_years))
    # F vector 
    F_vec <- c(constant_F, change_F[-1]) # Removing the first year from this sequence
    
  } # if fleet structure is constant and increases after
  
  return(F_vec)
  
} # end function


# Get fishing mortality rates ---------------------------------------------

#' @param Start_F Starting F vector
#' @param Fish_Start_yr Starting Fishing years for fleets
#' @param n_years Number of years in simulation
#' @param max_rel_F_M Maximum F value relative to the M
#' @param desc_rel_F_M Descending F value relative to the M
#' @param F_type Specify type of fishing mortality pattern
#' @param mean_nat_mort Mean natural mortality

get_Fs <- function(Start_F, Fish_Start_yr, F_type, n_years, max_rel_F_M, desc_rel_F_M,
                   mean_nat_mort, yr_chng = NULL) {
  
  if(n_fish_fleets != length(Start_F) & n_fish_fleets != length(Fish_Start_yr) &
     n_fish_fleets != length(F_type) & n_fish_fleets != length(max_rel_F_M) &
     n_fish_fleets != length(desc_rel_F_M)
     ) stop("Vector length of Start_F, Fish_Start_yr, F_type, max_rel_F_M, desc_rel_F_M does not equal the number of fleets specified.")
  
  # Loop through fleets to specify Fs
  for(f in 1:n_fish_fleets) {
    
    # Set Fs to 0 before fishery starts
    fish_mort[1:(Fish_Start_yr-1)[f],f,] <- 0
    
    # Set F pattern to fishery - feed in inputs to arguments above
    fish_mort[Fish_Start_yr[f]:n_years,f,] <- specify_F_pattern(Start_F = Start_F[f], # Starting F
                                              Fish_Start_yr = Fish_Start_yr[f], # When fishery starts
                                              F_type = F_type[f], # What type of F pattern
                                              n_years = n_years, # Number of model years
                                              max_rel_F_M = max_rel_F_M[f], # Apical F relative to M
                                              desc_rel_F_M = desc_rel_F_M[f], # Descending F after apical relative to M
                                              mean_nat_mort = mean_nat_mort, # Mean Natural Mortality
                                              yr_chng = yr_chng) # what year fleet structure begins to change
  } # end f loop
  
  # Output into environment
  fish_mort <<- fish_mort
  Fish_Start_yr <<- Fish_Start_yr 
  
} # end function



